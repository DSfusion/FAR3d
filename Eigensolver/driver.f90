program driver
  implicit none

  integer, parameter :: r8 = kind(1.0d0)
  character(len=80) :: arg1, arg2
  real(r8) :: center_freq, growth_rate

  CALL getarg(1, arg1)
  read(arg1,*) center_freq
  CALL getarg(2, arg2)
  read(arg2,*) growth_rate

  call far3d_eigen(center_freq, growth_rate)
  stop
end program driver
