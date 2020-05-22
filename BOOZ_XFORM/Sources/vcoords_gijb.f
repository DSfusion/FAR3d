      SUBROUTINE vcoords_gijb(rmnc, zmns, pmns, gmnc, rmns, zmnc, pmnc,
     1    gmns, xmb, xnb, cosm_boz, sinm_boz, cosn_boz, sinn_boz, mboz,
     2    nboz, mnmax, js, jr, nsd, gi, grr, grt, gtt, grroj, grtoj,
     3    gttoj, jbgrr, jbgrt, jbgtt, nznt, nfp, lasym_b)
      USE stel_kinds
      USE booz_params, ONLY: ns, ohs
      USE booz_persistent, ONLY: sfull
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: js, jr, nsd, mnmax, nznt, nfp, mboz, nboz
      REAL(rprec), DIMENSION(mnmax,nsd), INTENT(in) :: 
     1  rmnc, zmns, pmns, gmnc    !BOOZER COORDINATES
      REAL(rprec), DIMENSION(mnmax,nsd), INTENT(in) ::
     1  rmns, zmnc, pmnc, gmns
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: xmb, xnb
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: gi, grr, grt, gtt,
     1  grroj, grtoj, gttoj, jbgrr, jbgrt, jbgtt
      REAL(rprec), DIMENSION(nznt,0:mboz) :: cosm_boz, sinm_boz
      REAL(rprec), DIMENSION(nznt,0:nboz) :: cosn_boz, sinn_boz
      LOGICAL, INTENT(in) :: lasym_b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, m, n
      REAL(rprec) :: gc, rc, drc, dzs, dps, sgn
      REAL(rprec), DIMENSION(nznt) :: tsin, tcos, g, r, drr, dzr, dpr,
     1   drt, dzt, dpt
C-----------------------------------------------

      g   = 0
      r   = 0
      drr = 0
      dzr = 0
      dpr = 0
      drt = 0
      dzt = 0
      dpt = 0

!
!       Compute Reven, Rodd and Zeven, Zodd in Real BOOZER Space
!       on half radial grid (rmncb, zmnsb, etc ARE ALREADY on half mesh
!

      DO mn = 1, mnmax
         m = NINT(xmb(mn))
         n = NINT(ABS(xnb(mn)/nfp))
         sgn = SIGN(one,xnb(mn))

         tcos = cosm_boz(:,m)*cosn_boz(:,n)
     1        + sinm_boz(:,m)*sinn_boz(:,n)*sgn
         tsin = sinm_boz(:,m)*cosn_boz(:,n)
     1        - cosm_boz(:,m)*sinn_boz(:,n)*sgn
         IF (jr < ns) THEN
          gc  = 0.5*(gmnc(mn,js+1)+gmnc(mn,js))
          rc  = 0.5*(rmnc(mn,js+1)+rmnc(mn,js))
          drc = (rmnc(mn,js+1)-rmnc(mn,js))*ohs  
          dzs = (zmns(mn,js+1)-zmns(mn,js))*ohs
          dps = (pmns(mn,js+1)-pmns(mn,js))*ohs
          g   = g + tcos*gc
          r   = r + tcos*rc
          drr = drr + tcos*drc
          dzr = dzr + tsin*dzs
          dpr = dpr + tsin*dps
          drc = -0.5*xmb(mn)*(rmnc(mn,js+1)+rmnc(mn,js))  
          dzs =  0.5*xmb(mn)*(zmns(mn,js+1)+zmns(mn,js))
          dps =  0.5*xmb(mn)*(pmns(mn,js+1)+pmns(mn,js))
          drt = drt + tsin*drc
          dzt = dzt + tcos*dzs
          dpt = dpt + tcos*dps
          IF (lasym_b) THEN
            gc  = 0.5*(gmns(mn,js+1)+gmns(mn,js))
            rc  = 0.5*(rmns(mn,js+1)+rmns(mn,js))
            drc = (rmns(mn,js+1)-rmns(mn,js))*ohs  
            dzs = (zmnc(mn,js+1)-zmnc(mn,js))*ohs
            dps = (pmnc(mn,js+1)-pmnc(mn,js))*ohs
            g   = g + tsin*gc
            r   = r + tsin*rc
            drr = drr + tsin*drc
            dzr = dzr + tcos*dzs
            dpr = dpr + tcos*dps
            drc =  0.5*xmb(mn)*(rmns(mn,js+1)+rmns(mn,js))  
            dzs = -0.5*xmb(mn)*(zmnc(mn,js+1)+zmnc(mn,js))
            dps = -0.5*xmb(mn)*(pmnc(mn,js+1)+pmnc(mn,js))
            drt = drt + tcos*drc
            dzt = dzt + tsin*dzs
            dpt = dpt + tsin*dps
          END IF
         ELSE
          gc  = 0.5*(3*gmnc(mn,js)-gmnc(mn,js-1))
          rc  = 0.5*(3*rmnc(mn,js)-rmnc(mn,js-1))
          drc = (2*rmnc(mn,js)-3*rmnc(mn,js-1)+rmnc(mn,js-2))*ohs  
          dzs = (2*zmns(mn,js)-3*zmns(mn,js-1)+zmns(mn,js-2))*ohs
          dps = (2*pmns(mn,js)-3*pmns(mn,js-1)+pmns(mn,js-2))*ohs
          g   = g + tcos*gc
          r   = r + tcos*rc
          drr = drr + tcos*drc
          dzr = dzr + tsin*dzs
          dpr = dpr + tsin*dps
          drc = -0.5*xmb(mn)*(3*rmnc(mn,js)-rmnc(mn,js-1))  
          dzs =  0.5*xmb(mn)*(3*zmns(mn,js)-zmns(mn,js-1))
          dps =  0.5*xmb(mn)*(3*pmns(mn,js)-pmns(mn,js-1))
          drt = drt + tsin*drc
          dzt = dzt + tcos*dzs
          dpt = dpt + tcos*dps
          IF (lasym_b) THEN
            gc  = 0.5*(3*gmns(mn,js)-gmns(mn,js-1))
            rc  = 0.5*(3*rmns(mn,js)-rmns(mn,js-1))
            drc = (2*rmns(mn,js)-3*rmns(mn,js-1)+rmns(mn,js-2))*ohs  
            dzs = (2*zmnc(mn,js)-3*zmnc(mn,js-1)+zmnc(mn,js-2))*ohs
            dps = (2*pmnc(mn,js)-3*pmnc(mn,js-1)+pmnc(mn,js-2))*ohs
            g   = g + tsin*gc
            r   = r + tsin*rc
            drr = drr + tsin*drc
            dzr = dzr + tcos*dzs
            dpr = dpr + tcos*dps
            drc =  0.5*xmb(mn)*(3*rmns(mn,js)-rmns(mn,js-1))  
            dzs = -0.5*xmb(mn)*(3*zmnc(mn,js)-zmnc(mn,js-1))
            dps = -0.5*xmb(mn)*(3*pmnc(mn,js)-pmnc(mn,js-1))
            drt = drt + tcos*drc
            dzt = dzt + tsin*dzs
            dpt = dpt + tsin*dps
          END IF
         END IF
      END DO

      gi  = 1.0_dp/g
      grr = 4*sfull(jr)*sfull(jr)*(drr*drr+dzr*dzr+r*r*dpr*dpr)
      grt = 2*(drr*drt+dzr*dzt+r*r*dpr*dpt)
      gtt = (drt*drt+dzt*dzt+r*r*dpt*dpt)/(sfull(jr)*sfull(jr))
      grroj = -grr/g
      grtoj = -grt/g
      gttoj = -gtt/g
      jbgrr = -g*grr
      jbgrt = -g*grt
      jbgtt = -g*gtt

      END SUBROUTINE vcoords_gijb
