program driver
  use iso_c_binding, only : rp => c_double
  implicit none

  integer, parameter :: maxmodes=11
  integer nmodes
  integer, dimension(maxmodes) :: ntor
  real(rp), dimension(maxmodes) :: gamma,sig_gam,omega,sig_om

#ifdef NAMELIST_INPUT
  call far3d_main(.TRUE.,nmodes,ntor,gamma,omega,sig_gam,sig_om)
#else
  call far3d_main(.FALSE.,nmodes,ntor,gamma,omega,sig_gam,sig_om)
#endif

end program driver
