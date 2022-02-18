!\
! Collection of data and methods to access IDSs for FAR3D interface
!/
module dcon_ids_mod
#ifdef IMAS
  
  USE ids_routines, only : MDSPLUS_BACKEND, HDF5_BACKEND, ASCII_BACKEND  
  USE ids_schemas, only : ids_equilibrium, ids_core_profiles
  use iso_c_binding, only : c_double
  implicit none

  integer, parameter, private :: R8 = c_double

  type (ids_equilibrium), private :: eq ! IDS equilibrium
  type (ids_core_profiles), private :: cp ! core profiles

  character(len=5), parameter, private     :: treename = 'imas'
  character(len=12), parameter, private    :: imas_version='3.34.0'
  character(len=12), parameter, private    :: al_version='4.9.3'
  integer, parameter, private :: ids_backend = HDF5_BACKEND ! MDSPLUS_BACKEND / ASCII_BACKEND
     
  logical, private :: first_slice = .TRUE.
#endif    

  contains
  !\
  ! Save beam data in IMAS
  !/ 
  subroutine ids_beamPut(tr_beam, ztime, ierr)
    use far3d_wout 
    
    type(bdstruc), intent(in)  :: tr_beam
    real(R8), intent(in)  :: ztime
    integer, intent(out)  :: ierr

    ierr = 0
#ifdef IMAS
    !
    ! add information to equilibrium IDS
    if (.not.associated(eq%time_slice)) then 
      allocate(eq%time_slice(1))
      eq%time_slice(1)%time = ztime
      if (.not.associated(eq%time)) allocate(eq%time(1))
      eq%time(1) = ztime
    endif
    eq%time_slice(1)%boundary%triangularity = tr_beam%triavg
    eq%time_slice(1)%boundary%elongation    = tr_beam%elong
    eq%time_slice(1)%boundary%minor_radius  = tr_beam%minrad    ! should be in [m]
    eq%time_slice(1)%boundary%geometric_axis%r  = tr_beam%rgcen ! rm is the same
    eq%time_slice(1)%global_quantities%magnetic_axis%r = tr_beam%rcenrm
    eq%time_slice(1)%global_quantities%magnetic_axis%b_field_tor = tr_beam%bt ! [T]
    eq%time_slice(1)%global_quantities%beta_normal = tr_beam%beta0 ! normalized beta tor time_s 
    !
    ! fill core_profiles IDS
    if (.not.associated(cp%profiles_1d)) allocate(cp%profiles_1d(1))
    cp%profiles_1d(1)%time = ztime
    if (.not.associated(cp%time)) allocate(cp%time(1))
    cp%time(1) = ztime
    if (.not.associated(cp%profiles_1d(1)%grid%rho_tor_norm))  &
       allocate(cp%profiles_1d(1)%grid%rho_tor_norm(tr_beam%nr))
    cp%profiles_1d(1)%grid%rho_tor_norm(1:tr_beam%nr) = tr_beam%rho(1:tr_beam%nr)
    if (.not.associated(cp%profiles_1d(1)%q)) allocate(cp%profiles_1d(1)%q(tr_beam%nr))
    cp%profiles_1d(1)%q(1:tr_beam%nr) = tr_beam%q(1:tr_beam%nr)
    ! main ion species
    if (.not.associated(cp%profiles_1d(1)%ion)) then
       if (tr_beam%lalpha) then
         allocate(cp%profiles_1d(1)%ion(3))
         ! alphas 
         if (.not.associated(cp%profiles_1d(1)%ion(3)%density)) &
           allocate(cp%profiles_1d(1)%ion(3)%density(tr_beam%nr))
         cp%profiles_1d(1)%ion(3)%density(1:tr_beam%nr) = tr_beam%zdens(1:tr_beam%nr)*1d20 ! in [m^-3]
         if (.not.associated(cp%profiles_1d(1)%ion(3)%temperature)) &
           allocate(cp%profiles_1d(1)%ion(3)%temperature(tr_beam%nr))
         cp%profiles_1d(1)%ion(3)%temperature(1:tr_beam%nr) = tr_beam%ta(1:tr_beam%nr)*1d3 ! in [eV]
       else
         allocate(cp%profiles_1d(1)%ion(2))
       endif
    endif
    if (.not.associated(cp%profiles_1d(1)%ion(1)%element)) &
       allocate(cp%profiles_1d(1)%ion(1)%element(1))
    cp%profiles_1d(1)%ion(1)%element(1)%a = tr_beam%ionmass ! in [amu]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%density)) &
       allocate(cp%profiles_1d(1)%ion(1)%density(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%density(1:tr_beam%nr) = tr_beam%idens(1:tr_beam%nr)*1d20 ! in [m^-3]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%density_fast)) &
       allocate(cp%profiles_1d(1)%ion(1)%density_fast(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%density_fast(1:tr_beam%nr) = tr_beam%bdens(1:tr_beam%nr)*1d20 ! in [m^-3]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%temperature)) &
       allocate(cp%profiles_1d(1)%ion(1)%temperature(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%temperature(1:tr_beam%nr) = tr_beam%ti(1:tr_beam%nr)*1d3 ! in [eV]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%pressure)) &
       allocate(cp%profiles_1d(1)%ion(1)%pressure(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%pressure(1:tr_beam%nr) = tr_beam%pt(1:tr_beam%nr)*1d3 ! in [Pa]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%pressure_fast_perpendicular)) &
       allocate(cp%profiles_1d(1)%ion(1)%pressure_fast_perpendicular(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%pressure_fast_perpendicular(1:tr_beam%nr) = tr_beam%pbeam(1:tr_beam%nr)*1d3 ! in [Pa]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%pressure_fast_parallel)) &
       allocate(cp%profiles_1d(1)%ion(1)%pressure_fast_parallel(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%pressure_fast_parallel(1:tr_beam%nr) = tr_beam%pbeam(1:tr_beam%nr)*1d3 ! in [Pa]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%velocity%toroidal)) &
       allocate(cp%profiles_1d(1)%ion(1)%velocity%toroidal(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%velocity%toroidal(1:tr_beam%nr) = tr_beam%trot(1:tr_beam%nr)*1d3 ! in [m/s]
    if (.not.associated(cp%profiles_1d(1)%ion(1)%velocity%poloidal)) &
       allocate(cp%profiles_1d(1)%ion(1)%velocity%poloidal(tr_beam%nr))
    cp%profiles_1d(1)%ion(1)%velocity%poloidal(1:tr_beam%nr) = tr_beam%prot(1:tr_beam%nr)*1d3 ! in [m/s]
    ! electrons
    if (.not.associated(cp%profiles_1d(1)%electrons%density)) &
       allocate(cp%profiles_1d(1)%electrons%density(tr_beam%nr))
    cp%profiles_1d(1)%electrons%density(1:tr_beam%nr) = tr_beam%edens(1:tr_beam%nr)*1d20 ! in [m^-3]
    if (.not.associated(cp%profiles_1d(1)%electrons%temperature)) &
       allocate(cp%profiles_1d(1)%electrons%temperature(tr_beam%nr))
    cp%profiles_1d(1)%electrons%temperature(1:tr_beam%nr) = tr_beam%te(1:tr_beam%nr)*1d3 ! in [eV]
    if (.not.associated(cp%profiles_1d(1)%electrons%pressure)) &
       allocate(cp%profiles_1d(1)%electrons%pressure(tr_beam%nr))
    cp%profiles_1d(1)%electrons%pressure(1:tr_beam%nr) = tr_beam%pe(1:tr_beam%nr)*1d3 ! in [Pa]
    ! impurities
    if (.not.associated(cp%profiles_1d(1)%ion(2)%label)) &
       allocate(cp%profiles_1d(1)%ion(2)%label(1))
    cp%profiles_1d(1)%ion(2)%label(1) = tr_beam%mcontsp
    if (.not.associated(cp%profiles_1d(1)%ion(2)%density)) &
       allocate(cp%profiles_1d(1)%ion(2)%density(tr_beam%nr))
    cp%profiles_1d(1)%ion(2)%density(1:tr_beam%nr) = tr_beam%zdens(1:tr_beam%nr)*1d20 ! in [m^-3]
#else
    write(*,*) ' ids_beamPut: imas libraries are not available.'     
    ierr = -1
#endif
  end subroutine ids_beamPut

  !\
  ! Get beam data from IMAS
  ! This routine assumes that only one time slice exist
  !/ 
  subroutine ids_beamGet(tr_beam, znonlin, ierr)
    use far3d_wout,  only : bdstruc
    
    type(bdstruc), intent(inout)  :: tr_beam

    integer, intent(in)   :: znonlin
    integer, intent(out)  :: ierr
    !
    ierr = 0
    !
#ifdef IMAS
    !
    ! get information from equilibrium IDS
    if (size(eq%time_slice) < 1) then
       ierr = 1
       return
    endif
    tr_beam%triavg = eq%time_slice(1)%boundary%triangularity
    tr_beam%elong = eq%time_slice(1)%boundary%elongation
    tr_beam%minrad = eq%time_slice(1)%boundary%minor_radius    ! should be in [m]
    tr_beam%rgcen = eq%time_slice(1)%boundary%geometric_axis%r  ! [m]
    tr_beam%rm = eq%time_slice(1)%boundary%geometric_axis%r     ! [m]
    tr_beam%rcenrm = eq%time_slice(1)%global_quantities%magnetic_axis%r ! [m]
    tr_beam%bt = eq%time_slice(1)%global_quantities%magnetic_axis%b_field_tor ! [T]
    tr_beam%beta0 = eq%time_slice(1)%global_quantities%beta_normal ! normalized beta tor time_s 
    !
    ! get info from core_profiles IDS
    if (size(cp%profiles_1d) < 1) then 
      ierr = 1
      return 
    endif
    tr_beam%nr = size(cp%profiles_1d(1)%grid%rho_tor_norm)
    if (.not.allocated(tr_beam%rho)) allocate(tr_beam%rho(tr_beam%nr))
    tr_beam%rho(1:tr_beam%nr) = cp%profiles_1d(1)%grid%rho_tor_norm(1:tr_beam%nr)
    if (.not.allocated(tr_beam%q)) allocate(tr_beam%q(tr_beam%nr))
    tr_beam%q(1:tr_beam%nr) = cp%profiles_1d(1)%q(1:tr_beam%nr) 
    if (size(cp%profiles_1d(1)%ion) == 3) then
      tr_beam%lalpha = .TRUE.
      if (.not.allocated(tr_beam%adens)) allocate(tr_beam%adens(tr_beam%nr))
      tr_beam%zdens(1:tr_beam%nr) = cp%profiles_1d(1)%ion(3)%density(1:tr_beam%nr)*1d-20 ! in [10^20 m^-3]
      if (.not.allocated(tr_beam%ta)) allocate(tr_beam%ta(tr_beam%nr))
      tr_beam%ta(1:tr_beam%nr) = cp%profiles_1d(1)%ion(3)%temperature(1:tr_beam%nr)*1d-3 ! in [keV]
    else
      tr_beam%lalpha = .FALSE.
    endif
    ! main ion species
    tr_beam%ionmass = cp%profiles_1d(1)%ion(1)%element(1)%a ! in [amu]
    if (.not.allocated(tr_beam%idens)) allocate(tr_beam%idens(tr_beam%nr))
    tr_beam%idens(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%density(1:tr_beam%nr)*1d-20 ! in [10^20 m^-3]
    if (.not.allocated(tr_beam%bdens)) allocate(tr_beam%bdens(tr_beam%nr))
    tr_beam%bdens(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%density_fast(1:tr_beam%nr)*1d-20 ! in [10^20 m^-3]
    if (.not.allocated(tr_beam%ti)) allocate(tr_beam%ti(tr_beam%nr))
    tr_beam%ti(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%temperature(1:tr_beam%nr)*1d-3 ! in [keV]
    if (.not.allocated(tr_beam%pt)) allocate(tr_beam%pt(tr_beam%nr))
    tr_beam%pt(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%pressure(1:tr_beam%nr)*1d-3 ! in [kPa]
    if (.not.allocated(tr_beam%pbeam)) allocate(tr_beam%pbeam(tr_beam%nr))
    tr_beam%pbeam(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%pressure_fast_perpendicular(1:tr_beam%nr)*1d-3 ! in [kPa]
    if (.not.allocated(tr_beam%tbeam)) allocate(tr_beam%tbeam(tr_beam%nr))
    tr_beam%tbeam = 0_r8
    where (tr_beam%bdens.NE.0.0)
      tr_beam%tbeam = tr_beam%pbeam/(1.602d+01*tr_beam%bdens)
    end where
    if (.not.allocated(tr_beam%trot)) allocate(tr_beam%trot(tr_beam%nr))
    tr_beam%trot(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%velocity%toroidal(1:tr_beam%nr)*1d-3 ! in [km/s]
    if (.not.allocated(tr_beam%prot)) allocate(tr_beam%prot(tr_beam%nr))
    tr_beam%prot(1:tr_beam%nr) = cp%profiles_1d(1)%ion(1)%velocity%poloidal(1:tr_beam%nr)*1d-3 ! in [km/s]
    ! electrons
    if (.not.allocated(tr_beam%edens)) allocate(tr_beam%edens(tr_beam%nr))
    tr_beam%edens(1:tr_beam%nr) = cp%profiles_1d(1)%electrons%density(1:tr_beam%nr)*1d-20 ! in [10^20 m^-3]
    if (.not.allocated(tr_beam%te)) allocate(tr_beam%te(tr_beam%nr))
    tr_beam%te(1:tr_beam%nr) = cp%profiles_1d(1)%electrons%temperature(1:tr_beam%nr)*1d-3 ! in [keV]
    if (.not.allocated(tr_beam%pe)) allocate(tr_beam%pe(tr_beam%nr))
    tr_beam%pe(1:tr_beam%nr) = cp%profiles_1d(1)%electrons%pressure(1:tr_beam%nr)*1d-3 ! in [kPa]
    ! impurities
    tr_beam%mcontsp = cp%profiles_1d(1)%ion(2)%label(1)
    if (.not.allocated(tr_beam%zdens)) allocate(tr_beam%zdens(tr_beam%nr))
    tr_beam%zdens(1:tr_beam%nr) = cp%profiles_1d(1)%ion(2)%density(1:tr_beam%nr)*1d-20 ! in [10^20 m^-3]
#else
    write(*,*) ' ids_beamGet: imas libraries are not available.'     
    ierr = -1
#endif
  end subroutine ids_beamGet

  
  !\
  ! Save eq data in IMAS
  !/ 
  subroutine ids_eqPut(tr_eq, ztime, ierr)
    use transp_eq 
    
    type(transpeq), intent(in)  :: tr_eq
    real(R8), intent(in)  :: ztime
    integer, intent(out)  :: ierr

    ierr = 0
#ifdef IMAS
    if (.not.associated(eq%time_slice))                         &
             allocate(eq%time_slice(1))
    if (.not.associated(eq%time_slice(1)%profiles_1d%psi))      &
             allocate(eq%time_slice(1)%profiles_1d%psi(tr_eq%nsm1+1))
    if (.not.associated(eq%time_slice(1)%profiles_1d%q))        &
             allocate(eq%time_slice(1)%profiles_1d%q(tr_eq%nsm1+1))
    if (.not.associated(eq%time_slice(1)%profiles_1d%f))        &
             allocate(eq%time_slice(1)%profiles_1d%f(tr_eq%nsm1+1))
    if (.not.associated(eq%time_slice(1)%profiles_1d%pressure)) &
             allocate(eq%time_slice(1)%profiles_1d%pressure(tr_eq%nsm1+1))
    eq%time_slice(1)%profiles_1d%psi(1:tr_eq%nsm1+1)      = tr_eq%psis(0:tr_eq%nsm1) ! in Wb
    eq%time_slice(1)%profiles_1d%q(1:tr_eq%nsm1+1)        = tr_eq%qs(0:tr_eq%nsm1)   
    eq%time_slice(1)%profiles_1d%f(1:tr_eq%nsm1+1)        = tr_eq%fs(0:tr_eq%nsm1)   ! in T.m
    eq%time_slice(1)%profiles_1d%pressure(1:tr_eq%nsm1+1) = tr_eq%ps(0:tr_eq%nsm1)   ! in Pa
    eq%time_slice(1)%time = ztime 
    if (.not.associated(eq%time_slice(1)%profiles_2d))           &
             allocate(eq%time_slice(1)%profiles_2d(1))
    ! (R,Z) coords (1:nr) (1:nz)
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%grid%dim1))   &
             allocate(eq%time_slice(1)%profiles_2d(1)%grid%dim1(tr_eq%nr))
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%grid%dim2))   &
             allocate(eq%time_slice(1)%profiles_2d(1)%grid%dim2(tr_eq%nz))
    eq%time_slice(1)%profiles_2d(1)%grid%dim1(1:tr_eq%nr)   = tr_eq%rr(1:tr_eq%nr)
    eq%time_slice(1)%profiles_2d(1)%grid%dim2(1:tr_eq%nz)   = tr_eq%zz(1:tr_eq%nz)
    ! Magnetic field in T
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%b_field_r))   &
             allocate(eq%time_slice(1)%profiles_2d(1)%b_field_r(1:tr_eq%nr,1:tr_eq%nz))
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%b_field_z))   &
             allocate(eq%time_slice(1)%profiles_2d(1)%b_field_z(1:tr_eq%nr,1:tr_eq%nz))
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%b_field_tor)) &
             allocate(eq%time_slice(1)%profiles_2d(1)%b_field_tor(1:tr_eq%nr,1:tr_eq%nz))
    eq%time_slice(1)%profiles_2d(1)%b_field_r(1:tr_eq%nr,1:tr_eq%nz)   = tr_eq%Br(1:tr_eq%nr,1:tr_eq%nz)
    eq%time_slice(1)%profiles_2d(1)%b_field_z(1:tr_eq%nr,1:tr_eq%nz)   = tr_eq%Bz(1:tr_eq%nr,1:tr_eq%nz)
    eq%time_slice(1)%profiles_2d(1)%b_field_tor(1:tr_eq%nr,1:tr_eq%nz) = tr_eq%Bphi(1:tr_eq%nr,1:tr_eq%nz)
    ! psi in Wb
    if (.not.associated(eq%time_slice(1)%profiles_2d(1)%psi))   &
             allocate(eq%time_slice(1)%profiles_2d(1)%psi(1:tr_eq%nr,1:tr_eq%nz))
    eq%time_slice(1)%profiles_2d(1)%psi(1:tr_eq%nr,1:tr_eq%nz) = tr_eq%psirz(1:tr_eq%nr,1:tr_eq%nz)
    ! Cylindrical (R,Z) coords (0:nsm1,0:ntm1) in m
    if (.not.associated(eq%time_slice(1)%coordinate_system%r))    &
             allocate(eq%time_slice(1)%coordinate_system%r(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1))
    if (.not.associated(eq%time_slice(1)%coordinate_system%z))    &
             allocate(eq%time_slice(1)%coordinate_system%z(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1))
    eq%time_slice(1)%coordinate_system%r(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1) = tr_eq%rg(0:tr_eq%nsm1,0:tr_eq%ntm1)
    eq%time_slice(1)%coordinate_system%z(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1) = tr_eq%zg(0:tr_eq%nsm1,0:tr_eq%ntm1)
    
    if (.not.associated(eq%time)) allocate(eq%time(1))
    eq%time(1) = ztime

#else
    write(*,*) ' ids_eqPut: imas libraries are not available.'     
    ierr = -1
#endif

  end subroutine ids_eqPut
 
  !\
  ! Get eq data from IMAS
  ! This routine assumes that only one time slice exist
  !/ 
  subroutine ids_eqGet(tr_eq, znonlin, ierr)
    use transp_eq 
    
    type(transpeq), intent(out)  :: tr_eq

    integer, intent(in)   :: znonlin
    integer, intent(out)  :: ierr
    !
    ierr = 0
    !
#ifdef IMAS
    tr_eq%nsm1 = size(eq%time_slice(1)%profiles_1d%psi) - 1
    if (.not.allocated(tr_eq%psis)) allocate(tr_eq%psis(0:tr_eq%nsm1))
    tr_eq%psis(0:tr_eq%nsm1) = eq%time_slice(1)%profiles_1d%psi(1:tr_eq%nsm1+1)      ! in Wb
    if (.not.allocated(tr_eq%qs)) allocate(tr_eq%qs(0:tr_eq%nsm1))
    tr_eq%qs(0:tr_eq%nsm1) = eq%time_slice(1)%profiles_1d%q(1:tr_eq%nsm1+1)
    if (.not.allocated(tr_eq%fs)) allocate(tr_eq%fs(0:tr_eq%nsm1))
    tr_eq%fs(0:tr_eq%nsm1) = eq%time_slice(1)%profiles_1d%f(1:tr_eq%nsm1+1)          ! in T.m
    if (.not.allocated(tr_eq%ps)) allocate(tr_eq%ps(0:tr_eq%nsm1))
    tr_eq%ps(0:tr_eq%nsm1) = eq%time_slice(1)%profiles_1d%pressure(1:tr_eq%nsm1+1)   ! in Pa
    if (size(eq%time_slice(1)%profiles_2d)<1) then
      write(znonlin,*) ' ids_eqGet: equilibrium IDS does not include 2d profiles'
      ierr = -2
      return
    endif
    ! (R,Z) coords (1:nr) (1:nz)
    tr_eq%nr = size(eq%time_slice(1)%profiles_2d(1)%grid%dim1)
    tr_eq%nz = size(eq%time_slice(1)%profiles_2d(1)%grid%dim2)
    if (.not.allocated(tr_eq%rr)) allocate(tr_eq%rr(1:tr_eq%nr))
    tr_eq%rr(1:tr_eq%nr) = eq%time_slice(1)%profiles_2d(1)%grid%dim1(1:tr_eq%nr)
    if (.not.allocated(tr_eq%zz)) allocate(tr_eq%zz(1:tr_eq%nz))
    tr_eq%zz(1:tr_eq%nz) = eq%time_slice(1)%profiles_2d(1)%grid%dim2(1:tr_eq%nz)
    ! Magnetic field in T
    if (.not.allocated(tr_eq%Br)) allocate(tr_eq%Br(1:tr_eq%nr,1:tr_eq%nz))
    tr_eq%Br(1:tr_eq%nr,1:tr_eq%nz) = eq%time_slice(1)%profiles_2d(1)%b_field_r(1:tr_eq%nr,1:tr_eq%nz)
    if (.not.allocated(tr_eq%Bz)) allocate(tr_eq%Bz(1:tr_eq%nr,1:tr_eq%nz))
    tr_eq%Bz(1:tr_eq%nr,1:tr_eq%nz) = eq%time_slice(1)%profiles_2d(1)%b_field_z(1:tr_eq%nr,1:tr_eq%nz)
    if (.not.allocated(tr_eq%Bphi)) allocate(tr_eq%Bphi(1:tr_eq%nr,1:tr_eq%nz))
    tr_eq%Bphi(1:tr_eq%nr,1:tr_eq%nz) = eq%time_slice(1)%profiles_2d(1)%b_field_tor(1:tr_eq%nr,1:tr_eq%nz)
    ! psi in Wb
    if (.not.allocated(tr_eq%psirz)) allocate(tr_eq%psirz(1:tr_eq%nr,1:tr_eq%nz))
    tr_eq%psirz(1:tr_eq%nr,1:tr_eq%nz) = eq%time_slice(1)%profiles_2d(1)%psi(1:tr_eq%nr,1:tr_eq%nz)
    ! Cylindrical (R,Z) coords (0:nsm1,0:ntm1) in m
    tr_eq%ntm1 = size(eq%time_slice(1)%coordinate_system%r(tr_eq%nsm1,:))-1
    if (.not.allocated(tr_eq%rg)) allocate(tr_eq%rg(0:tr_eq%nsm1,0:tr_eq%ntm1))
    tr_eq%rg(0:tr_eq%nsm1,0:tr_eq%ntm1) = eq%time_slice(1)%coordinate_system%r(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1)
    if (.not.allocated(tr_eq%zg)) allocate(tr_eq%zg(0:tr_eq%nsm1,0:tr_eq%ntm1))
    tr_eq%zg(0:tr_eq%nsm1,0:tr_eq%ntm1) = eq%time_slice(1)%coordinate_system%z(1:tr_eq%nsm1+1,1:tr_eq%ntm1+1)
    eq%time_slice(1)%profiles_2d(1)%psi(1:tr_eq%nr,1:tr_eq%nz) = tr_eq%psirz(1:tr_eq%nr,1:tr_eq%nz)

#else
    write(znonlin,*) ' ids_eqGet: imas libraries are not available.'     

#endif

  end subroutine ids_eqGet

  !\
  ! write dcon data to ids
  !/
  subroutine dcon2ids(saveLastOnly, runid, zctok, ierr)
#ifdef IMAS   
    use ids_routines             ! Access Layer routines + management of IDS structures
    use ids_schemas              ! Fortran type definitions for the Physics Data Model
#endif
    logical, intent(in)   :: saveLastOnly   ! True if only the last time slice is to be saved
    integer, intent(out)  :: ierr
    character(len=*), intent(in) :: runid, zctok
#ifdef IMAS   
    !
    !  Variables that define the record in the IMAS database 
    ! real(ids_real)  :: time_ids
    integer :: idx = -1

    integer              :: run, shot
    character(len=25)    :: usr
    character(len=25) :: zrunid
    integer   :: zl

    zrunid = adjustl(runid)
    zl = len_trim(zrunid) 
    read(zrunid(1:zl-3), *) shot
    read(zrunid(zl-1:zl), *) run
    run = ichar(zrunid(zl-2:zl-2))*100 + run

    IF (TRIM(usr) .EQ. '') call get_environment_variable("USER",usr)

    call ual_begin_pulse_action(ids_backend, shot, run, TRIM(usr), TRIM(zctok), imas_version, idx, ierr)
    if (first_slice.OR.saveLastOnly) then
      if (ierr == 0) then
        call ual_open_pulse(idx, FORCE_CREATE_PULSE, "", ierr)
      end if
      if (ierr /= 0) then
        write(*,*) ' dcon2ids: can not create new IMAS environment'
        ierr = -1
        return
      endif
      ! basic information
      eq%ids_properties%homogeneous_time = 0 ! the time values are stored in the various time fields at lower levels in the tree
      cp%ids_properties%homogeneous_time = 0 ! the time values are stored in the various time fields at lower levels in the tree
      allocate(eq%ids_properties%version_put%data_dictionary(1))
      allocate(eq%ids_properties%version_put%access_layer(1))
      allocate(eq%ids_properties%version_put%access_layer_language(1))
      allocate(eq%ids_properties%comment(1))
      allocate(eq%ids_properties%provider(1))
      eq%ids_properties%comment(1) = 'Interface equilibirum data for coupling with FAR3D'
      eq%ids_properties%provider(1) = 'TRANSP'
      eq%ids_properties%version_put%data_dictionary(1) = imas_version
      eq%ids_properties%version_put%access_layer(1)    = al_version
      eq%ids_properties%version_put%access_layer_language(1) = 'Fortran'
      allocate(cp%ids_properties%version_put%data_dictionary(1))
      allocate(cp%ids_properties%version_put%access_layer(1))
      allocate(cp%ids_properties%version_put%access_layer_language(1))
      allocate(cp%ids_properties%comment(1))
      allocate(cp%ids_properties%provider(1))
      cp%ids_properties%comment(1) = 'Interface equilibirum data for coupling with FAR3D'
      cp%ids_properties%provider(1) = 'TRANSP'
      cp%ids_properties%version_put%data_dictionary(1) = imas_version
      cp%ids_properties%version_put%access_layer(1)    = al_version
      cp%ids_properties%version_put%access_layer_language(1) = 'Fortran'
    else 
      if (ierr == 0) then
        call ual_open_pulse(idx, OPEN_PULSE, "", ierr)
      end if
      if (ierr /= 0) then
        write(*,*) ' dcon2ids: can not create new IMAS environment'
        ierr = -1
        return
      endif
    endif
 
    if (first_slice.OR.saveLastOnly) then
      CALL ids_put(idx,"equilibrium",eq,ierr)
      if (ierr==0) CALL ids_put(idx,"core_profiles",cp,ierr)
      first_slice = .False. 
    else
      CALL ids_put_slice(idx,"equilibrium",eq,ierr)
      if (ierr==0) CALL ids_put_slice(idx,"core_profiles",cp,ierr)
    endif
    if (ierr /= 0) then
      write(*,*) ' dcon2ids: can not save equilibrium and core_profiles IDSs'
      ierr = -1
      return
    endif
    CALL imas_close(idx, ierr)
    if (ierr /= 0) then
      write(*,*) ' dcon2ids: can not close IMAS files'
      ierr = -1
      return
    endif
#else
    write(*,*) ' dcon2ids: imas libraries are not available.'     
    ierr = -1
#endif
  end subroutine dcon2ids

  !\
  ! Open and read IDSs for the FAR3D interface
  !/
  subroutine ids2dcon(znonlin, zrunid, zctok, ierr)
#ifdef IMAS   
    use ids_routines             ! Access Layer routines + management of IDS structures
    use ids_schemas              ! Fortran type definitions for the Physics Data Model
#endif
    integer, intent(in)   :: znonlin
    integer, intent(out)  :: ierr
    character(len=*), intent(in) :: zrunid, zctok
#ifdef IMAS    
    !  Variables that define the record in the IMAS database 
    ! real(ids_real)  :: time_ids
    integer :: idx = -1
    !
    integer              :: run, shot
    character(len=25)    :: usr
    !
    character(len=25) :: ztmp
    integer   :: zl

    ztmp = adjustl(zrunid)
    zl = len_trim(ztmp) 
    read(ztmp(1:zl-3), *) shot
    read(ztmp(zl-1:zl), *) run
    run = ichar(ztmp(zl-2:zl-2))*100 + run

    IF (TRIM(usr) .EQ. '') call get_environment_variable("USER",usr)

    call ual_begin_pulse_action(ids_backend, shot, run, TRIM(usr), zctok, imas_version, idx, ierr)

    call ual_open_pulse(idx, OPEN_PULSE, "", ierr)
    if (ierr /= 0) then
      write(znonlin,*) ' ids2dcon: can not open IMAS environment'
      return
    endif
    ! read the equilibrium IDS
    CALL ids_get(idx,"equilibrium",eq,ierr)
    if (ierr /= 0) then
      write(znonlin,*) ' ids2dcon: can not get the equilibrium IDS'
      return
    endif
    if (size(eq%time_slice)<1) then
      write(znonlin,*) ' ids2dcon: equilibrium IDS does not include time slices'
      ierr = -1
      return
    endif

    ! read the core_profiles IDS
    CALL ids_get(idx,"core_profiles",cp,ierr)
    if (ierr /= 0) then
      write(znonlin,*) ' ids2dcon: can not get the core_profiles IDS'
      return
    endif
    if (size(cp%profiles_1d)<1) then
      write(znonlin,*) ' ids2dcon: core_profiles IDS does not include any profiles'
      ierr = -1
      return
    endif

    CALL imas_close(idx, ierr)
    if (ierr /= 0) then
      write(znonlin,*) ' ids2dcon: can not close the IDSs'
      return
    endif
#else
    write(*,*) ' ids2dcon: imas libraries are not available.'     
    ierr = -1
#endif
  end subroutine ids2dcon

end module dcon_ids_mod
