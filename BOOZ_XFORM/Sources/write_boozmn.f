      SUBROUTINE write_boozmn (extension)
      USE booz_params
      USE booz_persistent, ONLY: xmb, xnb 
      USE read_wout_mod, ONLY: rmax_surf, rmin_surf, betaxis, aspect
      USE safe_open_mod
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) :: extension
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: version =
     1   "Boozer Transformation Code Version 2.0"
!-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, iunit, istat, js
C-----------------------------------------------
      iunit = unit_booz
      CALL safe_open (iunit, istat, 'boozmn.' // extension, 'replace',
     1     'unformatted')
      IF (istat .ne. 0) THEN
         PRINT *,' istat = ', istat
         STOP 'Error opening boozmn file in XBOOZ_XFORM!'
      END IF

!
!     Write out surface quantities needed by ballooning code
!
      WRITE(iunit, iostat=istat, err=100)
     1   nfp, ns, aspect, rmax_surf, rmin_surf, betaxis

      DO js = 2, ns
         WRITE(iunit, iostat=istat, err=100) hiota(js), pres(js),
     1   beta_vol(js), phip(js), phi(js), bvco(js), buco(js)
      END DO

!     SPH 070909: ADDED lasym to dump
!     CRCook 10/8/12: ADDED lrfp to dump
      WRITE (iunit, iostat=istat, err=100) mboz, nboz, mnboz, jsize
      WRITE (iunit, iostat=istat, err=100) version, lasym_b, lrfp_b

      WRITE (iunit, iostat=istat, err=100) NINT(xnb(:mnboz)), 
     1                                     NINT(xmb(:mnboz))

!
!     Write packed (in radial coordinate) 2D arrays
!
      DO i = 1, jsize
        js = jlist(i)
        IF (js.le.0 .or. js.gt.ns) CYCLE
        WRITE (iunit, iostat=istat, err=100) js
        WRITE (iunit, iostat=istat, err=100) bmncb(:mnboz,i),
     1       rmncb(:mnboz,i), zmnsb(:mnboz,i), pmnsb(:mnboz,i),
     2       gmncb(:mnboz,i)
!SPH070909: WRITE OUT ASYMMETRIC PARTS
        IF (lasym_b) THEN
        WRITE (iunit, iostat=istat, err=100) bmncb(:mnboz,i),
     1       rmnsb(:mnboz,i), zmncb(:mnboz,i), pmncb(:mnboz,i),
     2       gmnsb(:mnboz,i)
        ENDIF	
      END DO


 100  CONTINUE
      IF (istat .gt. 0)
     1    PRINT *,' Error writing in subroutine write_boozmn:',
     2            ' istat = ', istat

      CLOSE(iunit)

      END SUBROUTINE write_boozmn
