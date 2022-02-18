!-----------------------------------------------------------------------
! Use FFTW-3 C Interface to compute sin, cos transform of multiple arrays
!  of real periodic input data.  Retain buffers, plan between calls to
!  save time.
SUBROUTINE scxform_many_(inarr, nreal, narr, outcos, outsin, lcleanup)
  USE, INTRINSIC :: iso_c_binding !For fftw3 C interface
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  !Arguments
  INTEGER, INTENT(IN)                               :: nreal, narr
  REAL(C_DOUBLE), DIMENSION(nreal,narr), INTENT(IN) :: inarr
  REAL(C_DOUBLE), DIMENSION(0:nreal/2,*), INTENT(OUT) :: outcos, outsin
  LOGICAL, INTENT(IN)                               :: lcleanup

  !Local variables
  TYPE(C_PTR), SAVE                        :: fftwplan, pin, pout
  REAL(C_DOUBLE), POINTER, SAVE            :: rin(:,:)
  COMPLEX(C_DOUBLE_COMPLEX), POINTER, SAVE :: cout(:,:)
  REAL(C_DOUBLE), SAVE                     :: scale=0.0
  INTEGER, SAVE                            :: rsize=0, csize=0, nrows=0
  LOGICAL, SAVE                            :: first=.TRUE.

  IF (lcleanup.OR.&
       ((.NOT.first).AND.((nreal.NE.rsize).OR.(narr.NE.nrows)))) THEN
     CALL fftw_destroy_plan(fftwplan)
     CALL fftw_free(pin);  CALL fftw_free(pout)
     first = .TRUE.
     IF (lcleanup) THEN
        CALL fftw_cleanup
        RETURN
     ENDIF
  ENDIF

  IF ((nreal.LT.1).OR.(narr.LT.1)) RETURN

  IF (first) THEN
     rsize = nreal;  csize = rsize/2 + 1;  nrows = narr
     scale = 2.0d0/REAL(rsize,C_DOUBLE)

     pin = fftw_alloc_real(INT(rsize*nrows, C_SIZE_T))
     CALL C_F_POINTER(pin, rin, [rsize, nrows])

     pout = fftw_alloc_complex(INT(csize*nrows, C_SIZE_T))
     CALL C_F_POINTER(pout, cout, [csize, nrows])

     fftwplan = fftw_plan_many_dft_r2c(1, (/ rsize /), nrows, &
          rin, (/ rsize /), 1, rsize, &
          cout, (/ csize /), 1, csize, FFTW_MEASURE)

     first = .FALSE.
  ENDIF

  rin = inarr
  CALL fftw_execute_dft_r2c(fftwplan, rin, cout)

  ! Assume zeta = pi + theta, thus signs of odd coefs are flipped.
  outcos(0,1:nrows)   = (5.0d-1*scale)* REAL(cout(1,1:nrows))
  outcos(1:csize-1:2,1:nrows) = -scale* REAL(cout(2:csize:2,1:nrows))
  outcos(2:csize-1:2,1:nrows) =  scale* REAL(cout(3:csize:2,1:nrows))
  outsin(0:csize:2,1:nrows)   = -scale*AIMAG(cout(1:csize+1:2,1:nrows))
  outsin(1:csize:2,1:nrows)   =  scale*AIMAG(cout(2:csize+1:2,1:nrows))

END SUBROUTINE scxform_many_
