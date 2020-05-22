      SUBROUTINE boozer_gij (thgrd, ztgrd, gib, grrb, grtb, gttb,
     1     grrojb, grtojb, gttojb, jbgrrb, jbgrtb, jbgttb, xmb, xnb,
     1     gimncb, grrmncb, grtmnsb, gttmncb, grrojmncb, grtojmnsb,
     2     gttojmncb, jbgrrmncb, jbgrtmnsb, jbgttmncb,
     1     gimnsb, grrmnsb, grtmncb, gttmnsb, grrojmnsb, grtojmncb,
     2     gttojmnsb, jbgrrmnsb, jbgrtmncb, jbgttmnsb,
     2     scl, cosm_boz, sinm_boz, cosn_boz, sinn_boz,
     3     mboz, nboz, mnmax, nfp, nznt, nu2, nv, lasym_b)
C...MODIFIED 6/98 by A. WARE to speed up by factor of 8
C
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mnmax, nznt, mboz, nboz, nfp, nu2, nv
      REAL(rprec), DIMENSION(nznt), INTENT(in) ::
     1   thgrd, ztgrd, gib, grrb, grtb, gttb, grrojb, grtojb, gttojb,
     2   jbgrrb, jbgrtb, jbgttb
      REAL(rprec), DIMENSION(mnmax), INTENT(in) ::
     1   xmb, xnb, scl
      REAL(rprec), DIMENSION(nznt,0:mboz) :: cosm_boz, sinm_boz
      REAL(rprec), DIMENSION(nznt,0:nboz) :: cosn_boz, sinn_boz
      REAL(rprec), DIMENSION(mnmax), INTENT(out) ::
     1   gimncb, grrmncb, grtmnsb, gttmncb, grrojmncb, grtojmnsb,
     2   gttojmncb, jbgrrmncb, jbgrtmnsb, jbgttmncb, gimnsb, grrmnsb,
     3   grtmncb, gttmnsb, grrojmnsb, grtojmncb, gttojmnsb, jbgrrmnsb,
     4   jbgrtmncb, jbgttmnsb
      LOGICAL, INTENT(in) :: lasym_b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one=1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, m, n, imax, i
      REAL(rprec) :: sgn
      REAL(rprec), DIMENSION(nznt) :: tsin, tcos
C-----------------------------------------------
      IF (.not.lasym_b) THEN
!     ONLY INTEGRATE IN U HALF WAY AROUND (FOR LASYM=F)
         i = nv*(nu2-1)+1                                      !u=pi interval: i:imax
         imax = i-1+nv
         DO m = 0,mboz
            cosm_boz(1:nv,m)    = 0.5*cosm_boz(1:nv,m)         !u=0
            cosm_boz(i:imax,m)  = 0.5*cosm_boz(i:imax,m)       !u=pi
            sinm_boz(1:nv,m)    = 0.5*sinm_boz(1:nv,m)         !should be zeroes
            sinm_boz(i:imax,m)  = 0.5*sinm_boz(i:imax,m)       !should be zeroes
         END DO
      END IF

      DO mn = 1, mnmax
         m = NINT(xmb(mn))
         n = NINT(ABS(xnb(mn)/nfp))
         sgn = SIGN(one,xnb(mn))

         tcos = cosm_boz(:,m)*cosn_boz(:,n)
     1        + sinm_boz(:,m)*sinn_boz(:,n)*sgn
         tsin = sinm_boz(:,m)*cosn_boz(:,n)
     1        - cosm_boz(:,m)*sinn_boz(:,n)*sgn

         gimncb(mn) = DOT_PRODUCT(gib,tcos)
         grrmncb(mn) = DOT_PRODUCT(grrb,tcos)
         grtmnsb(mn) = DOT_PRODUCT(grtb,tsin)
         gttmncb(mn) = DOT_PRODUCT(gttb,tcos)
         grrojmncb(mn) = DOT_PRODUCT(grrojb,tcos)
         grtojmnsb(mn) = DOT_PRODUCT(grtojb,tsin)
         gttojmncb(mn) = DOT_PRODUCT(gttojb,tcos)
         jbgrrmncb(mn) = DOT_PRODUCT(jbgrrb,tcos)
         jbgrtmnsb(mn) = DOT_PRODUCT(jbgrtb,tsin)
         jbgttmncb(mn) = DOT_PRODUCT(jbgttb,tcos)

         IF (.not.lasym_b) CYCLE

         gimnsb(mn) = DOT_PRODUCT(gib,tsin)
         grrmnsb(mn) = DOT_PRODUCT(grrb,tsin)
         grtmncb(mn) = DOT_PRODUCT(grtb,tcos)
         gttmnsb(mn) = DOT_PRODUCT(gttb,tsin)
         grrojmnsb(mn) = DOT_PRODUCT(grrojb,tsin)
         grtojmncb(mn) = DOT_PRODUCT(grtojb,tcos)
         gttojmnsb(mn) = DOT_PRODUCT(gttojb,tsin)
         jbgrrmnsb(mn) = DOT_PRODUCT(jbgrrb,tsin)
         jbgrtmncb(mn) = DOT_PRODUCT(jbgrtb,tcos)
         jbgttmnsb(mn) = DOT_PRODUCT(jbgttb,tsin)

      END DO

      gimncb = scl*gimncb
      grrmncb = scl*grrmncb
      grtmnsb = scl*grtmnsb
      gttmncb = scl*gttmncb
      grrojmncb = scl*grrojmncb
      grtojmnsb = scl*grtojmnsb
      gttojmncb = scl*gttojmncb
      jbgrrmncb = scl*jbgrrmncb
      jbgrtmnsb = scl*jbgrtmnsb
      jbgttmncb = scl*jbgttmncb

      IF (lasym_b) THEN

        gimnsb = scl*gimnsb
        grrmnsb = scl*grrmnsb
        grtmncb = scl*grtmncb
        gttmnsb = scl*gttmnsb
        grrojmnsb = scl*grrojmnsb
        grtojmncb = scl*grtojmncb
        gttojmnsb = scl*gttojmnsb
        jbgrrmnsb = scl*jbgrrmnsb
        jbgrtmncb = scl*jbgrtmncb
        jbgttmnsb = scl*jbgttmnsb

      ELSE

!     RECOVER cosm AND sinm FOR u=0 AND u=pi
         i = nv*(nu2-1)+1                                      !u=pi interval: i:imax
         imax = i-1+nv
         DO m = 0,mboz
            cosm_boz(1:nv,m)    = 2.0*cosm_boz(1:nv,m)         !u=0
            cosm_boz(i:imax,m)  = 2.0*cosm_boz(i:imax,m)       !u=pi
            sinm_boz(1:nv,m)    = 2.0*sinm_boz(1:nv,m)         !should be zeroes
            sinm_boz(i:imax,m)  = 2.0*sinm_boz(i:imax,m)       !should be zeroes
         END DO

      END IF

      END SUBROUTINE boozer_gij
