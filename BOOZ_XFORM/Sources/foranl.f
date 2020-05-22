      SUBROUTINE foranl(nu, nv, nfp, nunv, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: cosm_b, cosn_b, sinm_b, sinn_b,
     1    cosm_nyq, cosn_nyq, sinm_nyq, sinn_nyq, thgrd, ztgrd
      USE booz_params, ONLY: mpol1, ntor, mpol_nyq, ntor_nyq
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nu, nv, nfp, nunv
      LOGICAL, INTENT(in) :: lasym
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lk, lt, lz, istat=0
      REAL(rprec) :: dth, dzt, twopi
C-----------------------------------------------
      IF (.not.ALLOCATED(cosm_b)) ALLOCATE (
     1    cosm_b(nunv,0:mpol1), sinm_b(nunv,0:mpol1),
     2    cosn_b(nunv,0:ntor),  sinn_b(nunv,0:ntor),
     3    cosm_nyq(nunv,0:mpol_nyq), sinm_nyq(nunv,0:mpol_nyq),
     4    cosn_nyq(nunv,0:ntor_nyq),  sinn_nyq(nunv,0:ntor_nyq),
     6    thgrd(nunv), ztgrd(nunv), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in foranl'

      twopi = 8*ATAN(1.0_dp)
!
!     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
!
      IF (lasym) THEN
         dth = twopi/nu                  !USE THIS FOR FULL 2-pi
      ELSE
         dth = twopi/(2*(nu-1))          !Half-around in theta
      END IF

      dzt = twopi/(nv*nfp)
      lk = 0

      DO lt = 1, nu
         DO lz=1, nv
           lk = lk + 1
           thgrd(lk) = (lt-1)*dth
           ztgrd(lk) = (lz-1)*dzt
          END DO
      END DO

      CALL trigfunc (thgrd, ztgrd, cosm_b, sinm_b, cosn_b, sinn_b, 
     1               mpol1, ntor, nunv)

      CALL trigfunc (thgrd, ztgrd, cosm_nyq, sinm_nyq, cosn_nyq, 
     1               sinn_nyq, mpol_nyq, ntor_nyq, nunv)

      END SUBROUTINE foranl
