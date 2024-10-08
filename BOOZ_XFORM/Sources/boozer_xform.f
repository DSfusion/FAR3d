      PROGRAM boozer_xform
      USE booz_params
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, jrad, iread, numargs
      REAL(rprec) :: t1, t2
      CHARACTER(LEN=50) :: arg1, arg2
      CHARACTER(LEN=120) :: extension
      LOGICAL :: lwritecdf
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the boozer coordinates are required
!     writes the boozer coordinates to a file, boozmn.extension
!
!     call this as follows:
!
!          xbooz_xform input.boz [T or F]
!
!     WHERE input.boz CONTAINS the mboz, nboz, wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     mboz   nboz
!     FILE_EXTENSION (does NOT have to include the .nc or .txt extension for netcdf,text file)
!     1  3   5   10  12
!
!     The OPTIONAL command line argument, (T) or (F), Allows the user to turn off screen
!     output IF set to F.
!
!     CALL xbooz_xform -h brings up a help screen
!
      lscreen = .true.          !!Default, write to screen
      far = .false.             !!Default, don't calculate metric elements
      lwritecdf = .false.       !!Default, write binary, not NetCDF output

!
!     Read command line argument to get input file name
!
      CALL getcarg(1, arg1, numargs)
      IF (numargs .gt. 1) CALL getcarg(2, arg2, numargs)

      IF (numargs .lt. 1) THEN
         PRINT *,'Invalid command line in calling xbooz_xform'
         PRINT *,'Type xbooz_xform -h to get more information'
         STOP
      ELSE IF (arg1 .eq. '-h' .or. arg1 .eq. '/h') THEN
         PRINT *,' ENTER INPUT FILE NAME ON COMMAND LINE'
         PRINT *,' For example: xbooz_xform in_booz.ext'
         PRINT *
         PRINT *,' WHERE in_booz.ext is the input file'
         PRINT *
         PRINT *,' Optional command line argument'
         PRINT *,' xbooz_xform <infile> (T or F)'
         PRINT *
         PRINT *,' where F suppresses output to the screen'
         STOP
      ELSE IF (numargs .gt. 1) THEN
         IF (TRIM(arg2).eq.'f' .or. TRIM(arg2).eq.'F') lscreen = .false.
         IF (TRIM(arg2).eq.'far' .or. TRIM(arg2).eq.'FAR') far = .true.
      ENDIF

      iread = unit_booz-1
      CALL safe_open (iread, istat, TRIM(arg1), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'Error opening input file in booz_xform'

      READ (iread, *, iostat = istat) mboz, nboz
      READ (iread, *, iostat = istat) extension
      IF (istat .ne. 0) STOP 'Error reading input file in booz_xform'

!
!     READ IN PARAMETERS, DATA FROM WOUT FILE
!
      CALL read_wout_booz(extension, iread, istat)
      CALL second0(t1)
      IF (istat .ne. 0) THEN
         PRINT *,' ierr_vmec !=0 in booz_xform read_wout_booz'
         GOTO 1010
      END IF

      CLOSE (unit=iread)

!
!     COMPUTE BOOZER TRANSFORM, SURFACE BY SURFACE
!
      DO jrad = 1, ns
        IF (lsurf_boz(jrad)) CALL boozer_coords(jrad)
      END DO

! CRCook for debugging purposes
      OPEN(UNIT = 16, FILE = 'sum_gmncb.txt')
      DO jrad = 1, jsize
        WRITE(16,*) jrad, SUM(gmncb(:,jrad))
      ENDDO
      CLOSE(16)

      IF (far) call boozer_metric

! CRCook if LRFP = T, compute chip and chi
      IF (lrfp_b) THEN
        chip = hiota*phip
        DO jrad = 1, ns
          chi(jrad) = chip(jrad)*(jrad-1)/(ns-1)
        ENDDO
      ENDIF

 !
 !    WRITE OUT CONVERTED RESULTS
 !
      IF (lwritecdf) THEN
         CALL write_boozmn_cdf(extension)
      ELSE
         CALL write_boozmn(extension)
      ENDIF

 1010 CONTINUE
!
!     FREE MEMORY : USER MUST CALL BOOZER_COORDS AT END WITH LDEALLOC = TRUE
!     OTHERWISE MEMORY ALLOCATED WILL NOT BE AVAILABLE
!

      CALL free_mem_boozer

      CALL second0(t2)

      IF (lscreen) PRINT 120, t2-t1
 120  FORMAT(/,' TIME IN BOOZER TRANSFORM CODE:',1pe12.2,' SEC')

      END PROGRAM boozer_xform
