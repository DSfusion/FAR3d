MODULE transp_eq
  use iso_c_binding, only: c_double
  implicit none

  INTEGER, PARAMETER, private :: r8 = c_double

  type transpeq
     integer :: version                                         !Version number for tracking changes
     integer :: nsm1, ntm1                                      !R_minor, poloidal angle dimensions respectively
     integer :: nr, nz                                          !Cylindrical dims for inverse equilibrium
     real(r8), dimension(:), allocatable   :: psis, qs, fs, ps  !Flux, safety factor, R*B_T, pressure (0:nsm1)
     real(r8), dimension(:,:), allocatable :: rg, zg            !Cylindrical (R,Z) coords (0:nsm1,0:ntm1)
     real(r8), dimension(:), allocatable   :: theta             !Poloidal angle values (0:ntm1)
     real(r8), dimension(:,:), allocatable :: psirz             !Inverse equilibrium psi(R,z) (1:nr,1:nz)
     real(r8), dimension(:,:), allocatable :: BR,Bphi,Bz        !Cylindrical components of magnetic field (1:nr,1:nz) (T)
     real(r8), dimension(:), allocatable :: rr, zz              !R coords on R grid (1:nr) and z coords on z grid (1:nz)
  end type transpeq

CONTAINS
  SUBROUTINE readtranspeq(eqdat, ta, ier)
    IMPLICIT NONE

    TYPE(transpeq), INTENT(OUT) :: eqdat
    REAL(r8), INTENT(OUT)       :: ta
    INTEGER, INTENT(OUT)        :: ier

    CHARACTER*128 :: inname
    INTEGER iunit, ia, nt, ns

    IF (ALLOCATED(eqdat%psis)) DEALLOCATE(eqdat%psis)
    IF (ALLOCATED(eqdat%ps)) DEALLOCATE(eqdat%ps)
    IF (ALLOCATED(eqdat%qs)) DEALLOCATE(eqdat%qs)
    IF (ALLOCATED(eqdat%fs)) DEALLOCATE(eqdat%fs)
    IF (ALLOCATED(eqdat%rg)) DEALLOCATE(eqdat%rg)
    IF (ALLOCATED(eqdat%zg)) DEALLOCATE(eqdat%zg)
    IF (ALLOCATED(eqdat%psirz)) DEALLOCATE(eqdat%psirz)
    IF (ALLOCATED(eqdat%BR)) DEALLOCATE(eqdat%BR)
    IF (ALLOCATED(eqdat%Bphi)) DEALLOCATE(eqdat%Bphi)
    IF (ALLOCATED(eqdat%Bz)) DEALLOCATE(eqdat%Bz)
    IF (ALLOCATED(eqdat%rr)) DEALLOCATE(eqdat%rr)
    IF (ALLOCATED(eqdat%zz)) DEALLOCATE(eqdat%zz)

    !Read data from file
10  FORMAT(1p,3e24.16)
    iunit = 10;  inname = 'dconeq.dat'
    OPEN(iunit, FILE=inname, STATUS='OLD', IOSTAT=ier)
    IF (ier.eq.0) THEN
       READ(iunit,'(1e14.6)')ta
       READ(iunit,'(2I7)')nt,ns
       eqdat%ntm1 = nt - 1;  eqdat%nsm1 = ns - 1
       READ(iunit,'(2I7)')eqdat%nr, eqdat%nz
       IF ((eqdat%nsm1.gt.0).and.(eqdat%ntm1.gt.0)) THEN
          ALLOCATE(eqdat%psis(0:eqdat%nsm1), eqdat%ps(0:eqdat%nsm1), eqdat%qs(0:eqdat%nsm1), &
               eqdat%fs(0:eqdat%nsm1), eqdat%rg(0:eqdat%nsm1,0:eqdat%ntm1), &
               eqdat%zg(0:eqdat%nsm1,0:eqdat%ntm1))
          READ(iunit,10)eqdat%psis
          READ(iunit,10)eqdat%ps
          READ(iunit,10)eqdat%qs
          READ(iunit,10)eqdat%fs
          READ(iunit,10)(eqdat%rg(:,ia),ia=0,eqdat%ntm1)
          READ(iunit,10)(eqdat%zg(:,ia),ia=0,eqdat%ntm1)
       ENDIF
       IF ((eqdat%nr.gt.0).and.(eqdat%nz.gt.0)) THEN
          ALLOCATE(eqdat%psirz(eqdat%nr,eqdat%nz))
          READ(iunit,10)(eqdat%psirz(:,ia),ia=1,eqdat%nz)
          ALLOCATE(eqdat%BR(eqdat%nr,eqdat%nz), &
                   eqdat%Bphi(eqdat%nr,eqdat%nz), &
                   eqdat%Bz(eqdat%nr,eqdat%nz))
          ALLOCATE(eqdat%rr(eqdat%nr), eqdat%zz(eqdat%nz))
          READ(iunit,10,END=100)(eqdat%BR(:,ia),ia=1,eqdat%nz)
          READ(iunit,10)(eqdat%Bphi(:,ia),ia=1,eqdat%nz)
          READ(iunit,10)(eqdat%Bz(:,ia),ia=1,eqdat%nz)
          READ(iunit,10)eqdat%rr
          READ(iunit,10)eqdat%zz
          CLOSE(iunit)
          RETURN
100       eqdat%BR = 0.0;  eqdat%Bphi = 0.0;  eqdat%Bz = 0.0
          eqdat%rr = 0.0;  eqdat%zz = 0.0
       ENDIF
       CLOSE(iunit)
    ELSE
       PRINT *,'Error: could not open file '//TRIM(inname)//' for reading.'
    ENDIF
  END SUBROUTINE readtranspeq
END MODULE transp_eq
