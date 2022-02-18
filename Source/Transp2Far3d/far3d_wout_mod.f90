MODULE far3d_wout
  USE ISO_C_BINDING, ONLY : rp => C_DOUBLE
  IMPLICIT NONE

  INTEGER, PARAMETER :: ROOTFLUX=0, FLUX=1
  INTEGER, PARAMETER :: BMSRC_PS=1, BMSRC_TRCOM=2

  ! Datatype for VMEC equilibrium
  TYPE fardata
     INTEGER :: ns, nt1  ! Dimensions
     REAL(rp), DIMENSION(:),   ALLOCATABLE :: phip, iota, cur, f, pres
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: RMNC, RMNS, ZMNC, ZMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: BMNC, BMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GMNC, GMNS, GIMNC, GIMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GRRMNC, GRRMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GRROJC, GRROJS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: JBGRRC, JBGRRS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GRTMNC, GRTMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GRTOJC, GRTOJS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: JBGRTC, JBGRTS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GTTMNC, GTTMNS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: GTTOJC, GTTOJS
     REAL(rp), DIMENSION(:,:), ALLOCATABLE :: JBGTTC, JBGTTS
  END type fardata

  ! Datatype for beam, minority ion profiles
  TYPE bdstruc
     INTEGER  :: nr      ! Number of radial zones for profile data
     REAL(rp) :: rgcen   ! Major radius of geometric center (cm)
     REAL(rp) :: minrad  ! Average minor radius (cm)
     REAL(rp) :: rm      ! R_max
     REAL(rp) :: triavg  ! Average plasma triangularity
     REAL(rp) :: bt      ! Vacuum toroidal B on axis in Tesla
     REAL(rp) :: elong   ! Average plasma elongation
     REAL(rp) :: beta0   ! Target beta_0
     REAL(rp) :: rcenrm  ! Radial position of magnetic axis in cm
     REAL(rp) :: ionmass ! Atomic mass number of main ion species
     REAL(rp) :: omcy    ! Bulk ion cyclotron frequency normalized to 1/t_A
     REAL(rp) :: iflr    ! Bulk ion gyroradius normalized to minor radius
     REAL(rp) :: epflr   ! Beam ion gyroradius normalized to minor radius
     REAL(rp) :: bet0_f  ! Beam ion beta
     REAL(rp) :: omcyalp ! Alpha cyclotron frequency normalized to 1/t_A
     REAL(rp) :: bet0_alp ! Alpha particle beta
     REAL(rp) :: dpres   ! Electron pressure as a fraction of eqbm pressure
     REAL(rp), DIMENSION(:), ALLOCATABLE :: rho, q, bdens, idens
     REAL(rp), DIMENSION(:), ALLOCATABLE :: edens, zdens, adens
     REAL(rp), DIMENSION(:), ALLOCATABLE :: tbeam, ti, te, ta
     REAL(rp), DIMENSION(:), ALLOCATABLE :: pbeam, pt, pe, trot, prot
     CHARACTER(LEN=32) :: mcontsp ! Main contaminant species
     LOGICAL  :: lalpha  ! =TRUE if alpha species present
  END type bdstruc

CONTAINS
  SUBROUTINE far3d_init(feq, nsi, nti)
    IMPLICIT NONE
    TYPE(fardata), INTENT(INOUT) :: feq
    INTEGER, INTENT(IN) :: nsi, nti
    INTEGER mmax

    CALL far3d_free(feq)
    IF ((nsi.LT.1).OR.(nti.LT.1)) RETURN

    feq%ns = nsi;  feq%nt1 = nti
    mmax = (nti - 1)/2

    ALLOCATE(feq%phip(nsi), feq%iota(nsi), feq%cur(nsi), &
         feq%f(nsi), feq%pres(nsi))
    ALLOCATE(feq%RMNC(0:mmax,nsi), feq%RMNS(0:mmax,nsi))
    ALLOCATE(feq%ZMNC(0:mmax,nsi), feq%ZMNS(0:mmax,nsi))
    ALLOCATE(feq%BMNC(0:mmax,nsi), feq%BMNS(0:mmax,nsi))
    ALLOCATE(feq%GMNC(0:mmax,nsi), feq%GMNS(0:mmax,nsi))
    ALLOCATE(feq%GIMNC(0:mmax,nsi), feq%GIMNS(0:mmax,nsi))
    ALLOCATE(feq%GRRMNC(0:mmax,nsi), feq%GRRMNS(0:mmax,nsi))
    ALLOCATE(feq%GRTMNC(0:mmax,nsi), feq%GRTMNS(0:mmax,nsi))
    ALLOCATE(feq%GTTMNC(0:mmax,nsi), feq%GTTMNS(0:mmax,nsi))
    ALLOCATE(feq%GRROJC(0:mmax,nsi), feq%GRROJS(0:mmax,nsi))
    ALLOCATE(feq%GRTOJC(0:mmax,nsi), feq%GRTOJS(0:mmax,nsi))
    ALLOCATE(feq%GTTOJC(0:mmax,nsi), feq%GTTOJS(0:mmax,nsi))
    ALLOCATE(feq%JBGRRC(0:mmax,nsi), feq%JBGRRS(0:mmax,nsi))
    ALLOCATE(feq%JBGRTC(0:mmax,nsi), feq%JBGRTS(0:mmax,nsi))
    ALLOCATE(feq%JBGTTC(0:mmax,nsi), feq%JBGTTS(0:mmax,nsi))
  END SUBROUTINE far3d_init

!-----------------------------------------------------------------------
  SUBROUTINE far3d_free(feq)
    IMPLICIT NONE
    TYPE(fardata), INTENT(INOUT) :: feq

    IF (ALLOCATED(feq%phip)) DEALLOCATE(feq%phip)
    IF (ALLOCATED(feq%iota)) DEALLOCATE(feq%iota)
    IF (ALLOCATED(feq%cur)) DEALLOCATE(feq%cur)
    IF (ALLOCATED(feq%f)) DEALLOCATE(feq%f)
    IF (ALLOCATED(feq%pres)) DEALLOCATE(feq%pres)
    IF (ALLOCATED(feq%RMNC)) DEALLOCATE(feq%RMNC)
    IF (ALLOCATED(feq%RMNS)) DEALLOCATE(feq%RMNS)
    IF (ALLOCATED(feq%ZMNC)) DEALLOCATE(feq%ZMNC)
    IF (ALLOCATED(feq%ZMNS)) DEALLOCATE(feq%ZMNS)
    IF (ALLOCATED(feq%BMNC)) DEALLOCATE(feq%BMNC)
    IF (ALLOCATED(feq%BMNS)) DEALLOCATE(feq%BMNS)
    IF (ALLOCATED(feq%GMNC)) DEALLOCATE(feq%GMNC)
    IF (ALLOCATED(feq%GMNS)) DEALLOCATE(feq%GMNS)
    IF (ALLOCATED(feq%GIMNC)) DEALLOCATE(feq%GIMNC)
    IF (ALLOCATED(feq%GIMNS)) DEALLOCATE(feq%GIMNS)
    IF (ALLOCATED(feq%GRRMNC)) DEALLOCATE(feq%GRRMNC)
    IF (ALLOCATED(feq%GRRMNS)) DEALLOCATE(feq%GRRMNS)
    IF (ALLOCATED(feq%GRROJC)) DEALLOCATE(feq%GRROJC)
    IF (ALLOCATED(feq%GRROJS)) DEALLOCATE(feq%GRROJS)
    IF (ALLOCATED(feq%JBGRRC)) DEALLOCATE(feq%JBGRRC)
    IF (ALLOCATED(feq%JBGRRS)) DEALLOCATE(feq%JBGRRS)
    IF (ALLOCATED(feq%GRTMNC)) DEALLOCATE(feq%GRTMNC)
    IF (ALLOCATED(feq%GRTMNS)) DEALLOCATE(feq%GRTMNS)
    IF (ALLOCATED(feq%GRTOJC)) DEALLOCATE(feq%GRTOJC)
    IF (ALLOCATED(feq%GRTOJS)) DEALLOCATE(feq%GRTOJS)
    IF (ALLOCATED(feq%JBGRTC)) DEALLOCATE(feq%JBGRTC)
    IF (ALLOCATED(feq%JBGRTS)) DEALLOCATE(feq%JBGRTS)
    IF (ALLOCATED(feq%GTTMNC)) DEALLOCATE(feq%GTTMNC)
    IF (ALLOCATED(feq%GTTMNS)) DEALLOCATE(feq%GTTMNS)
    IF (ALLOCATED(feq%GTTOJC)) DEALLOCATE(feq%GTTOJC)
    IF (ALLOCATED(feq%GTTOJS)) DEALLOCATE(feq%GTTOJS)
    IF (ALLOCATED(feq%JBGTTC)) DEALLOCATE(feq%JBGTTC)
    IF (ALLOCATED(feq%JBGTTS)) DEALLOCATE(feq%JBGTTS)

    feq%ns=0;  feq%nt1=0
  END SUBROUTINE far3d_free

!-----------------------------------------------------------------------
  SUBROUTINE beamdata_init(bd, nr)
    IMPLICIT NONE

    TYPE(bdstruc), INTENT(INOUT) :: bd
    INTEGER, INTENT(IN) :: nr

    INTEGER j

    CALL beamdata_free(bd)

    ! Initialize scalars to zero
    bd%rcenrm = 0.0;  bd%bt = 0.0;  bd%rgcen = 0.0;  bd%minrad = 0.0
    bd%elong = 0.0;  bd%triavg = 0.0;  bd%mcontsp = '';  bd%ionmass = 0.0
    bd%beta0 = 0.0;  bd%rm = 0.0;  bd%omcy = 0.0;  bd%iflr = 0.0
    bd%epflr = 0.0;  bd%bet0_f = 0.0;  bd%omcyalp = 0.0;  bd%bet0_alp = 0.0
    bd%dpres = 0.0;  bd%lalpha = .FALSE.

    bd%nr = nr
    IF (nr.GT.0) THEN
       ALLOCATE(bd%rho(nr), bd%q(nr), bd%bdens(nr), bd%idens(nr), &
            bd%edens(nr), bd%adens(nr), bd%zdens(nr), bd%tbeam(nr), &
            bd%ti(nr), bd%te(nr), bd%ta(nr), bd%pbeam(nr), bd%pt(nr), &
            bd%pe(nr), bd%trot(nr), bd%prot(nr))

       ! Initialize arrays to zero
       bd%q = 0.0;  bd%bdens = 0.0
       bd%idens = 0.0;  bd%edens = 0.0;  bd%adens = 0.0;  bd%zdens = 0.0
       bd%tbeam = 0.0;  bd%ti = 0.0;  bd%te = 0.0;  bd%ta = 0.0
       bd%pbeam = 0.0;  bd%pt = 0.0;  bd%pe = 0.0
       bd%trot  = 0.0;  bd%prot = 0.0

       ! Initialize radial variable
       bd%rho = (/ ( j/REAL(nr-1,rp), j=0,nr-1) /)
    ENDIF
  END SUBROUTINE beamdata_init

!-----------------------------------------------------------------------
  SUBROUTINE beamdata_free(bd)
    IMPLICIT NONE
    TYPE(bdstruc), INTENT(INOUT) :: bd

    IF (ALLOCATED(bd%rho)) DEALLOCATE(bd%rho)
    IF (ALLOCATED(bd%q)) DEALLOCATE(bd%q)
    IF (ALLOCATED(bd%bdens)) DEALLOCATE(bd%bdens)
    IF (ALLOCATED(bd%idens)) DEALLOCATE(bd%idens)
    IF (ALLOCATED(bd%edens)) DEALLOCATE(bd%edens)
    IF (ALLOCATED(bd%zdens)) DEALLOCATE(bd%zdens)
    IF (ALLOCATED(bd%adens)) DEALLOCATE(bd%adens)
    IF (ALLOCATED(bd%tbeam)) DEALLOCATE(bd%tbeam)
    IF (ALLOCATED(bd%ti)) DEALLOCATE(bd%ti)
    IF (ALLOCATED(bd%te)) DEALLOCATE(bd%te)
    IF (ALLOCATED(bd%ta)) DEALLOCATE(bd%ta)
    IF (ALLOCATED(bd%pbeam)) DEALLOCATE(bd%pbeam)
    IF (ALLOCATED(bd%pt)) DEALLOCATE(bd%pt)
    IF (ALLOCATED(bd%pe)) DEALLOCATE(bd%pe)
    IF (ALLOCATED(bd%trot)) DEALLOCATE(bd%trot)
    IF (ALLOCATED(bd%prot)) DEALLOCATE(bd%prot)

    bd%nr = 0
  END SUBROUTINE beamdata_free
END MODULE far3d_wout
