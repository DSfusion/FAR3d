subroutine divv(f,c1,c2)

  use param
  use var_para
  use domain
  use equil
  use dynamo
  use dbyd
  use mult_mod
  use scratch

  implicit none

  integer :: l
  real(IDP) :: c1,c2
  real(IDP), dimension(mj_start:,0:) :: f

  sd2=feq/(feq-qqinv*cureq)
  sd3=cureq/(feq-qqinv*cureq)
  do l=1,leqmax
     sceq1(:,l)=-(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
  end do
  call dbydr_par(sc1,phi,0.0_IDP,1.0_IDP,0)
  call multed(sc3,sceq1,-1,sc1,-1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-(r*dbsjzoj(:,l)-sd2*djroj(:,l))
  end do
  call dbydth_par(sc1,phi,-1,0.0_IDP,1.0_IDP,0)
  call multed(sc3,sceq1,1,sc1,1,1.0_IDP,1.0_IDP)
  call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
  do l=1,lmax
     sc3(:,l)=sc3(:,l)+sd4(mj_start:mj_end)*sc1(:,l)
  end do
!$OMP END PARALLEL DO
  do l=1,leqmax
     sceq1(:,l)=-(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
  end do
  call dbydzt_par(sc1,phi,-1,0.0_IDP,1.0_IDP)
  call multed(sc3,sceq1,1,sc1,1,1.0_IDP,1.0_IDP)
  call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
!$OMP PARALLEL DO
  do l=1,lmax
     sc3(:,l)=sc3(:,l)-rinv(mj_start:mj_end)*sd4(mj_start:mj_end)*sc1(:,l)
  end do
!$OMP END PARALLEL DO

!  parallel thermal velocity term  

  do l=1,leqmax
     sceq1(:,l)=bmod(:,l)/(feq-qqinv*cureq)
  end do
  call grdpar(sc1,vthprlf,-1,0.0_IDP,-1.0_IDP)
  call multed(sc3,sceq1,1,sc1,1,1.0_IDP,1.0_IDP)
  call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=sceq2(:,l)/(feq-qqinv*cureq)
  end do
  call multed(sc3,sceq1,-1,vthprlf,-1,1.0_IDP,1.0_IDP)

  call mult(f,pr,1,sc3,1,c1,c2)

end subroutine divv
