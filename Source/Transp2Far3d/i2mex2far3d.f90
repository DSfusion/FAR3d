SUBROUTINE i2mex2far3d(feq, rzplot, rhovar)
  USE i2mex_mod
  USE far3d_wout
  IMPLICIT NONE

  TYPE(fardata), INTENT(INOUT) :: feq
  LOGICAL, INTENT(IN)          :: rzplot
  INTEGER, INTENT(IN)          :: rhovar

  REAL(rp), DIMENSION(:), ALLOCATABLE :: psi, rho, tht, buf, dpdr
  REAL(rp), DIMENSION(:,:), ALLOCATABLE :: buf2d, jac, gtt, gtp, gpp
  REAL(rp), PARAMETER :: twopi = 6.283185307179586476925286766559D0
  REAL(rp), PARAMETER :: bigno = 1.0d+12
  REAL(rp), PARAMETER :: mu0inv = 7.95774715459476678813d+5
  REAL(rp) :: phip, psi_min, psi_max, rerror, rho_b
  INTEGER :: nsi, nt1, ier, js, nth, nsurf, ntheta, jth, rhoexp
  LOGICAL, PARAMETER :: compute=.FALSE., cleanup=.TRUE.

  ! Query dimensions
  CALL i2mex_getOriNs(nsi, ier)
  CALL i2mex_getOriNt1(nt1, ier)
  PRINT *,'i2mex ns x nt1 = ',nsi,nt1
  IF ((nsi.LT.1).OR.(nt1.LT.2)) RETURN

  ! Temporary local storage
  nth = nt1 - 1
  ALLOCATE(psi(nsi), rho(nsi), buf(nsi), dpdr(nsi), tht(nth))
  ALLOCATE(buf2d(nth,nsi))
  ALLOCATE(jac(nth,nsi), gtt(nth,nsi), gtp(nth,nsi), gpp(nth,nsi))

  ! Allocate space for FAR3D arrays
  CALL far3d_init(feq, nsi, nt1)

  ! Uniform theta spacing
  tht = twopi*(/ (REAL(js-1,rp)/REAL(nth,rp), js=1,nth) /)

  ! Original poloidal flux profile psi
  CALL i2mex_getOriPsi(nsi, buf, ier);  CALL i2mex_error(ier)
  psi_min = MINVAL(buf);  psi_max = MAXVAL(buf)
  PRINT *,'i2mex psi range: ',psi_min,' to ',psi_max

  ! Adjust psi spacing?
  SELECT CASE (rhovar)
  CASE (ROOTFLUX) ! Uniform spacing in sqrt(phi)
     PRINT *,'Radial variable is sqrt(normalized toroidal flux).'
     rhoexp = 2
     CALL fluxmap(nsi, rhoexp, psi, phip, ier)
     feq%phip = phip  ! Phip = (1/twopi) * dPhi/ds (constant)
     !PRINT *,'phip = ',phip
     CALL i2mex_getPhi(nsi, psi, rho, ier);  CALL i2mex_error(ier)
     rho_b = rho(nsi)
     rho = SQRT(rho/rho_b)
     dpdr = 2.0_rp*rho_b*rho
  CASE(FLUX)
     PRINT *,'Radial variable is normalized toroidal flux.'
     rhoexp = 1
     CALL fluxmap(nsi, rhoexp, psi, phip, ier)
     feq%phip = phip  ! Phip = (1/twopi) * dPhi/ds (constant)
     CALL i2mex_getPhi(nsi, psi, rho, ier);  CALL i2mex_error(ier)
     rho_b = rho(nsi)
     rho = rho/rho_b
     dpdr = rho_b
  CASE DEFAULT ! Uniform spacing in psi
     PRINT *,'Radial variable is poloidal flux.'
     rhoexp = 0
     psi = psi_min + (psi_max - psi_min)* &
          (/ (REAL(js-1,rp)/REAL(nsi-1,rp), js=1,nsi) /)
     CALL i2mex_getPhiP(nsi, psi, feq%phip, ier);  CALL i2mex_error(ier)
     rho = psi
     dpdr = 1.0d0
  END SELECT

  !PRINT *,'psi spacing =',psi(2:nsi)-psi(1:nsi-1)
  !PRINT *,'delta-rho = ',rho(2:nsi)-rho(1:nsi-1)

  ! Rotational transform
  CALL i2mex_getQ(nsi, psi, buf, ier);  CALL i2mex_error(ier)
  PRINT *,'i2mex Q from ',MINVAL(buf),' to ',MAXVAL(buf)
  WHERE (buf.NE.0.0)
     feq%iota =  1.0d0/buf !twopi/buf
  ELSEWHERE
     feq%iota = bigno
  END WHERE
  !PRINT *,'iota =',feq%iota

  ! dpsi/drho
  IF ((rhovar.EQ.ROOTFLUX) .OR. (rhovar.EQ.FLUX)) &
       dpdr = dpdr*feq%iota
  !PRINT *,'Analytic dpsi/drho = ',dpdr

  ! Current profile
  CALL i2mex_getPlasmaCurProf(nth, nsi, tht, psi, feq%cur, ier)
  CALL i2mex_error(ier)
  feq%cur = feq%cur/twopi
  !PRINT *,'Ip =',feq%cur

  ! FF' = R B_toroidal
  CALL i2mex_getG(nsi, psi, feq%f, ier)
  CALL i2mex_error(ier)
  feq%f = -feq%f
  !PRINT *,'f = ',feq%f

  ! Pressure profile
  CALL i2mex_getP(nsi, psi, feq%pres, ier);  CALL i2mex_error(ier)
  PRINT *,'i2mex p range: ',feq%pres(1),' to ',feq%pres(nsi)
  feq%pres = mu0inv * feq%pres   ! Convert to Pa

  ! Cylindrical coordinates
  CALL i2mex_getX(nth, nsi, tht, psi, buf2d, ier);  CALL i2mex_error(ier)
  PRINT *,'i2mex R range: ',MINVAL(buf2d),' to ',MAXVAL(buf2d)
  IF (rzplot) THEN
     nsurf = 15;  ntheta = 32

     OPEN(17,FILE='fs.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
     IF (ier /= 0) STOP
     write(17,'(2I10)') nsurf,nth
     DO js=1,nsurf
        WRITE(17,'(4E22.12)') buf2d(:,js*nsi/nsurf)
     ENDDO

     OPEN(18,FILE='ct.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
     IF (ier /= 0) STOP
     write(18,'(2I10)') nsi,ntheta
     DO jth=1,ntheta
        WRITE(18,'(4E22.12)') buf2d(jth*nth/ntheta,:)
     ENDDO
  ENDIF !rzplot
  CALL scxform_many(buf2d, nth, nsi, feq%RMNC, feq%RMNS, compute)
  PRINT *,'R = ',feq%RMNC(0,nsi),' + ',feq%RMNC(1,nsi),' cos theta +...'

  CALL i2mex_getZ(nth, nsi, tht, psi, buf2d, ier);  CALL i2mex_error(ier)
  PRINT *,'i2mex Z range: ',MINVAL(buf2d),' to ',MAXVAL(buf2d)
  IF (rzplot) THEN
     DO js=1,nsurf
        WRITE(17,'(4E22.12)') buf2d(:,js*nsi/nsurf)
     ENDDO
     CLOSE(17)

     DO jth=1,ntheta
        WRITE(18,'(4E22.12)') buf2d(jth*nth/ntheta,:)
     ENDDO
     CLOSE(18)
  ENDIF !rzplot
  CALL scxform_many(buf2d, nth, nsi, feq%ZMNC, feq%ZMNS, compute)
  PRINT *,'Z = ',feq%ZMNC(0,nsi),' + ',feq%ZMNS(1,nsi),' sin theta +...'


  ! Mod B
  CALL get_i2mex_modB(nth, nsi, tht, psi, buf2d, ier);  CALL i2mex_error(ier)
  PRINT *,'i2mex mod B range: ',MINVAL(buf2d),' to ',MAXVAL(buf2d)
  CALL scxform_many(buf2d, nth, nsi, feq%BMNC, feq%BMNS, compute)

  ! Jacobian, metric terms
  CALL far3d_getJmet(nth, nsi, tht, psi, SPREAD(dpdr,1,nth), &
       jac, gtt, gtp, gpp, ier)
  CALL i2mex_error(ier)
  psi_min = MINVAL(jac);  psi_max = MAXVAL(jac)
  PRINT *,'Jacobian range: ',psi_min,' to ',psi_max
  IF (ANY(jac.eq.0.)) THEN
     PRINT *,'ERROR: Zero Jacobian detected in i2mex2far3d!'
     STOP
  ENDIF
  CALL scxform_many(jac, nth, nsi, feq%GMNC, feq%GMNS, compute)

  ! Metric terms
  CALL invert_matrix(nth, nsi, gtt, gtp, gpp) !Upper to lower
  gtp = 2.0d0*gtp
  DO js=2,nsi
     gpp(:,js) = 4.0d0*(rho(js)**rhoexp)*gpp(:,js)
     gtt(:,js) = gtt(:,js)/(rho(js)**rhoexp)
  ENDDO
  PRINT *,'gtt range: ',MINVAL(gtt),MAXVAL(gtt)
  CALL scxform_many(gtt, nth, nsi, feq%GTTMNC, feq%GTTMNS, compute)
  buf2d = -jac * gtt
  CALL scxform_many(buf2d, nth, nsi, feq%JBGTTC, feq%JBGTTS, compute)
  CALL scxform_many(gtp, nth, nsi, feq%GRTMNC, feq%GRTMNS, compute)
  buf2d = -jac * gtp
  CALL scxform_many(buf2d, nth, nsi, feq%JBGRTC, feq%JBGRTS, compute)
  CALL scxform_many(gpp, nth, nsi, feq%GRRMNC, feq%GRRMNS, compute)
  buf2d = -jac * gpp
  CALL scxform_many(buf2d, nth, nsi, feq%JBGRRC, feq%JBGRRS, compute)

  ! Inverse Jacobian
  jac = 1.0_rp/jac
  CALL scxform_many(jac, nth, nsi, feq%GIMNC, feq%GIMNS, compute)

  ! Metric terms
  buf2d = -jac * gtt
  PRINT *,'gtt/J range: ',MINVAL(buf2d),MAXVAL(buf2d)
  CALL scxform_many(buf2d, nth, nsi, feq%GTTOJC, feq%GTTOJS, compute)
  buf2d = -jac * gtp
  CALL scxform_many(buf2d, nth, nsi, feq%GRTOJC, feq%GRTOJS, compute)
  buf2d = -jac * gpp
  CALL scxform_many(buf2d, nth, nsi, feq%GRROJC, feq%GRROJS, compute)

  ! Check Grad-Shafranov error
  CALL i2mex_getGsError(nth, nsi, tht, psi, buf2d, ier)
  CALL i2mex_error(ier)
  PRINT *, 'Min gs err = ',MINVAL(buf2d)
  PRINT *,'Max gs err = ',MAXVAL(buf2d)
  PRINT *,'Mean gs err = ',SUM(buf2d)/(REAL(nth,rp)*REAL(nsi,rp))
  rerror = SQRT(SUM(buf2d**2)/(REAL(nth,rp)*REAL(nsi,rp)))
  PRINT *,'RMS gs err = ',rerror
  IF (ABS(rerror).GT.1.0d0) THEN
     PRINT *,' WARNING: HUGE GS error > 100%!!!'
  ELSE IF (ABS(rerror).GT.1.0d-2) THEN
     PRINT *,'Warning: LARGE rel GS error > 1%!'
  ENDIF

  ! Free up temporary storage
  CALL scxform_many(buf2d, nth, nsi, buf2d, buf2d, cleanup)
  DEALLOCATE(psi, rho, buf, dpdr, buf2d, tht, jac, gtt, gtp, gpp)
END SUBROUTINE i2mex2far3d

!-----------------------------------------------------------------------
SUBROUTINE get_i2mex_modB(nt1, ns, the, psi, modb, ier)
  USE i2mex_mod
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt1, ns
  REAL(i2mex_r8), INTENT(IN) :: the(nt1)
  REAL(i2mex_r8), INTENT(IN) :: psi(ns)
  REAL(i2mex_r8), INTENT(OUT) :: modb(nt1,ns)
  INTEGER, INTENT(OUT) :: ier

  REAL(i2mex_r8), DIMENSION(:,:), ALLOCATABLE :: dum1, dum2
  REAL(i2mex_r8) :: g(ns)
  INTEGER iok

  ALLOCATE(dum1(nt1,ns), dum2(nt1,ns))

  CALL i2mex_getMetric(nt1, ns, the, psi, dum1, dum2, modb, ier)
  CALL i2mex_getG(ns, psi, g, iok)
  modb = modb + SPREAD(g*g, dim=1, ncopies=nt1)
  CALL i2mex_getX(nt1, ns, the, psi, dum1, ier)
  modb = SQRT(modb)/dum1

  DEALLOCATE(dum1, dum2)
  IF (iok /= 0) ier = 109
END SUBROUTINE get_i2mex_modB

!-----------------------------------------------------------------------
! Determine a set of ns monotonically increasing values of poloidal flux
! psi that are uniformly spaced in square-root of toroidal flux phi.
SUBROUTINE fluxmap(ns, exp_s, psiarr, phip, ier)
  USE i2mex_mod
  USE ISO_C_BINDING, ONLY : rp=>C_DOUBLE
  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)                  :: ns, exp_s
  REAL(rp), DIMENSION(ns), INTENT(OUT) :: psiarr
  REAL(rp), INTENT(OUT)                :: phip
  INTEGER, INTENT(OUT)                 :: ier

  REAL(rp), EXTERNAL :: phi2psi_i2mex

  ! Local variables
  REAL(rp), ALLOCATABLE :: oripsi(:), oriphi(:)
  REAL(rp), PARAMETER :: rtol = 1.0d-12
  REAL(rp) :: phimin, phimax, phitarg, lb, ub, phib, res
  INTEGER, PARAMETER :: maxrec = 64
  INTEGER :: nsi, js, ks, irec

  ier = 1; phip=0.0d0

  IF (ns.LT.1) RETURN

  PRINT *
  PRINT *,'In fluxmap.'

  ! Find original data dimensions, allocate
  CALL i2mex_getOriNs(nsi, ier); CALL i2mex_error(ier)
  IF (nsi.LT.1) THEN
     PRINT *,'Fluxmap: Original data has too few surfaces.'
     ier = 2
     RETURN
  ENDIF
  PRINT *,'nsi,ns=',nsi,ns
  ALLOCATE(oripsi(nsi), oriphi(nsi))

  ! Get psi and phi on original mesh
  CALL i2mex_getOriPsi(nsi, oripsi, ier);  CALL i2mex_error(ier)
  PRINT *,'first to last psi:',oripsi(1),oripsi(nsi)
  PRINT *,'min,max:',MINVAL(oripsi),MAXVAL(oripsi)
  PRINT *,'Getting phi...'
  CALL i2mex_getPhi(nsi, oripsi, oriphi, ier)
  CALL i2mex_error(ier)
  PRINT *,'first to last phi:',oriphi(1),oriphi(nsi)
  phimin = MINVAL(oriphi);  phimax = MAXVAL(oriphi)
  PRINT *,'min,max:',phimin,phimax
  phip = oriphi(nsi) - oriphi(1)
  IF (ANY(phip*(oriphi(2:nsi) - oriphi(1:nsi-1)).LE.0.0)) THEN
     PRINT *,'Fluxmap: non-monotonic toroidal flux profile detected.'
     ier = 3
     RETURN
  ENDIF

  ! Set target toroidal flux values
  psiarr = phimin + (phimax - phimin)* &
       (/ ((REAL(js-1,rp)/REAL(ns-1,rp))**exp_s, js=1,ns) /)

  ! Search for corresponding psi values
  psiarr(1) = oripsi(1);  psiarr(ns) = oripsi(nsi)
  ks = 1
  DO js=2,ns-1
     phitarg = psiarr(js)
     DO
        IF (phitarg.LE.oriphi(ks+1)) EXIT
        ks = ks + 1
        IF (ks.GT.nsi-1) THEN
           PRINT *,'fluxmap: target phi value out of range!'
           ier = 5
           RETURN
        ENDIF
     ENDDO

     ! Recursively bracket the phi value
     lb = oripsi(ks);  ub = oripsi(ks+1);  phib = oriphi(ks)
     DO irec=1,maxrec
        psiarr(js) = phi2psi_i2mex(phitarg, lb, ub, phib, 64, res)
        IF (ABS(res/phip).LE.rtol) THEN
           EXIT
        ENDIF
     ENDDO
     IF (irec.GE.maxrec) PRINT *,'Warning: convergence failure in fluxmap.'
  ENDDO

  DEALLOCATE(oripsi, oriphi)
  ier = 0
  PRINT *,'Mapping complete.'
END SUBROUTINE fluxmap

!-----------------------------------------------------------------------
REAL*8 FUNCTION phi2psi_i2mex(phitarg, lb, ub, phib, nsb, resid)
  USE ISO_C_BINDING, ONLY : rp=>C_DOUBLE
  USE i2mex_mod
  IMPLICIT NONE
  REAL(rp), INTENT(IN)    :: phitarg
  REAL(rp), INTENT(INOUT) :: lb, ub, phib
  INTEGER, INTENT(IN)     :: nsb
  REAL(rp), INTENT(OUT)   :: resid

  REAL(rp),DIMENSION(:),ALLOCATABLE :: psif,phif
  INTEGER ii, ier

  ALLOCATE(psif(nsb), phif(nsb))
  psif = lb + (ub - lb) * &
       (/ ( REAL(ii-1,rp)/REAL(nsb-1,rp), ii=1,nsb) /)
  CALL i2mex_getPhi(nsb, psif, phif, ier); CALL i2mex_error(ier)
  ii = MINLOC(ABS((phif + phib) - phitarg),1)
  phi2psi_i2mex = psif(ii)
  resid = phif(ii) + phib - phitarg
  lb = psif(MAX(ii-1,1))
  ub = psif(MIN(ii+1,nsb))
  phib = phib + phif(MAX(ii-1,1))
  DEALLOCATE(psif,phif)
END FUNCTION phi2psi_i2mex

!-----------------------------------------------------------------------
! Invert 2x2 metric matrix at each (rho,theta) point to convert between
! covariant and contravariant coordinate systems.
SUBROUTINE invert_matrix(nt, ns, gtt, gtp, gpp)
  USE ISO_C_BINDING, ONLY : rp=>C_DOUBLE
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt, ns
  REAL(rp), DIMENSION(nt,ns), INTENT(INOUT) :: gtt, gtp, gpp

  REAL(rp) :: di,tmp
  INTEGER  :: js,kt

  DO js=1,ns
     DO kt=1,nt
        di = 1.0d0/(gpp(kt,js)*gtt(kt,js) - gtp(kt,js)**2)
        tmp = di*gtt(kt,js)
        gtt(kt,js) = di*gpp(kt,js)
        gpp(kt,js) = tmp
        gtp(kt,js) = -di*gtp(kt,js)
     ENDDO
  ENDDO
END SUBROUTINE invert_matrix

!-----------------------------------------------------------------------
! Use i2mex routines to get Jacobian, metric terms with respect to theta and
!  Far3D normalized minor radius rho ~ sqrt(Phi), where Phi is toroidal flux.
!
! gtt = |grad the|^2
! gtp = grad the . grad psi
! gpp = |grad psi|^2
! 
SUBROUTINE far3d_getJmet(nt, ns, the, psi, dpsidrho, &
                         jac, gtt, gtp, gpp, ier)
  USE ISO_C_BINDING, ONLY : rp=>C_DOUBLE
  USE i2mex_mod
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: nt, ns
  REAL(rp), INTENT(IN)  :: the(nt), psi(ns), dpsidrho(nt,ns)
  REAL(rp), INTENT(OUT) :: jac(nt,ns), gtt(nt,ns), gtp(nt,ns), gpp(nt,ns)
  INTEGER, INTENT(OUT)  :: ier

  REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xt, xp, zt, zp, xa
  REAL(rp) :: signj
  INTEGER iok

  ier = 0

  ALLOCATE(xt(nt,ns), xp(nt,ns), zt(nt,ns), zp(nt,ns), xa(nt,ns), STAT=iok)
  IF (iok/=0) THEN
     ier=44
     RETURN
  ENDIF

  CALL i2mex_getX(nt, ns, the, psi, xa, iok)
  CALL i2mex_getGradX(nt, ns, the, psi, xt, xp, iok)
  xp = xp * dpsidrho
  CALL i2mex_getGradZ(nt, ns, the, psi, zt, zp, iok)
  zp = zp * dpsidrho

  signj = 2.0_rp*(i2mex_o%isThetaClockwise-0.5_rp)
  jac = (xp*zt - xt*zp)*xa * signj

  call i2mex_toAxis(nt, ns, the, psi, 0.0_rp, jac, iok)    
  call i2mex_error(iok)
  jac(:,1) = abs(jac(:,1)*jac(:,2))/jac(:,2) ! enforce sign

  gpp =  xa*xa*(xt**2 + zt**2)/jac**2
  gtp = -xa*xa*(zt*zp + xt*xp)/jac**2
  gtt =  xa*xa*(zp*zp + xp*xp)/jac**2

  call i2mex_toAxis(nt, ns, the, psi, -1.0_rp, gtt, iok)
  call i2mex_error(iok)
  gtt(:,1) = abs(gtt(:,1)) ! enforce positive definiteness
  call i2mex_toAxis(nt, ns, the, psi, 0.0_rp, gtp, iok)
  call i2mex_error(iok)
  call i2mex_toAxis(nt, ns, the, psi, +1.0_rp, gpp, iok)
  call i2mex_error(iok)
  gpp(:,1) = abs(gpp(:,1)) ! enforce positive definiteness

  DEALLOCATE(xt, xp, zt, zp, xa)

  IF (iok /= 0) ier = 41
END SUBROUTINE far3d_getJmet
