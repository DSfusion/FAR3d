subroutine grid

! this sub sets up the r grid
! and calculates the time independent geometric scale factors

  use param
  use processor
  use cotrol
  use domain

  implicit none

  integer :: mjp, j, jm, jp

  if (Auto_grid_on == 0) then
     mjp=ne+nis+ni
     if(mj /= mjp) then
        if (myPE == 0) write(6,'("  mj=",i5,"  but ni+nis+ne=",i5)') mj,mjp
        stop
     end if

! compute the r values

     call findr
  else
     do j=0,mj
        r(j)=1.0_IDP*j/mj
     end do
  end if
  mjm1=mj-1
  mjm2=mj-2     
  
  do j=1,mjm1
     jm=j-1
     jp=j+1
     rinv(j)=1./r(j)
     dc1m(j)=-(r(jp)-r(j))/((r(j)-r(jm))*(r(jp)-r(jm)))
     dc1p(j)=(r(j)-r(jm))/((r(jp)-r(j))*(r(jp)-r(jm)))
     dc2m(j)=2./((r(j)-r(jm))*(r(jp)-r(jm)))
     dc2p(j)=2./((r(jp)-r(j))*(r(jp)-r(jm)))
     del2cm(j)=(3.*r(j)-r(jp))/(r(j)*(r(j)-r(jm))*(r(jp)-r(jm)))
     del2cp(j)=(3.*r(j)-r(jm))/(r(j)*(r(jp)-r(j))*(r(jp)-r(jm)))
  end do
  rinv(mj)=1./r(mj)
  rinv(0)=0.

end subroutine grid
