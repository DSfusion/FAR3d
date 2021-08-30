program driver
  implicit none

#ifdef NAMELIST_INPUT
  call far3d_main(.TRUE.)
#else
  call far3d_main(.FALSE.)
#endif

end program driver
