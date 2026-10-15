subroutine bigf(f,g,itypeg,h,itypeh,sx1,sx2,c1,c2)

  use param
  use processor
  use var_para
  use domain
  use dbyd
  use mult_mod

  implicit none

  integer :: itypeg,itypeh,l
  real(IDP) :: c1,c2
  real(IDP), dimension(mj_start:,0:) :: f,g,h,sx1,sx2

  call dbydr_par(sx2,g,0.0_IDP,-1.0_IDP,0)
  call mult(sx1,h,itypeh,sx2,itypeg,0.0_IDP,1.0_IDP)
  call dbydth_par(sx2,sx1,itypeg*itypeh,0.0_IDP,1.0_IDP,1)
  if (myPE == 0) then
     do l=1,lmax
        sx2(0,l)=0.
     end do
  end if
!$OMP PARALLEL DO
  do l=1,lmax
     f(mj_start:mj_end,l)=c1*f(mj_start:mj_end,l)+c2*sx2(mj_start:mj_end,l)
  end do
!$OMP END PARALLEL DO
  call dbydth_par(sx1,g,itypeg,0.0_IDP,1.0_IDP,0)
  call mult(sx2,h,itypeh,sx1,-itypeg,0.0_IDP,1.0_IDP)
  if (myPE == 0) then
     do l=1,lmax
        if (mm(l) == 0) f(0,l)=f(0,l)+c2*2.*rinv(1)*sx2(1,l)
     end do
  end if
!$OMP PARALLEL DO
  do l=1,lmax
     sx2(mj_start:mj_end,l)=r(mj_start:mj_end)*sx2(mj_start:mj_end,l)
  end do
!$OMP END PARALLEL DO
  call dbydr_par(sx1,sx2,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
  do l=1,lmax
     f(mj_start:mj_end,l)=f(mj_start:mj_end,l)+c2*rinv(mj_start:mj_end)*sx1(mj_start:mj_end,l)
  end do
!$OMP END PARALLEL DO

end subroutine bigf
