program driver
  implicit none

  integer, parameter :: r8 = kind(1.0d0)
  character(len=80) :: arg
  real(r8) :: center_freq, growth_rate
  integer :: n_eigenvals

  CALL getarg(1, arg)
  read(arg,*) center_freq
  CALL getarg(2, arg)
  read(arg,*) growth_rate
  call getarg(3, arg)
  if (len_trim(arg).gt.0) then
     read(arg,*) n_eigenvals
  else
     n_eigenvals = 7
  endif

  if (n_eigenvals.gt.0) &
       call far3d_eigen(center_freq, growth_rate, n_eigenvals)
end program driver
