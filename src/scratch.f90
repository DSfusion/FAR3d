module scratch
  use param
  implicit none
  save
  real(IDP), dimension(:,:), allocatable :: sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10,sc11,sc12
  real(IDP), dimension(:,:), allocatable :: sceq1,sceq2,sceq3,sceq4,sceq5,sceq6,sceq7
  real(IDP), dimension(:), allocatable :: sd1,sd2,sd3,sd4,sd5,sd6,sd7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     !  sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8 dummy dynamic variables with radial and angular dependency
     !  sceq1,sceq2,sceq3,sceq4,sceq5,sceq6,sceq7 dummy equilibrium variables with radial and angular dependency
     !  sd1,sd2,sd3,sd4,sd5,sd6 dummy variables with only radial dependency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module scratch
