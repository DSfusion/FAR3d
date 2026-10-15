subroutine lincheck

  use mpi
  use param
  use processor
  use var_para
  use cotrol
  use domain
  use equil
  use dynamo
  use om_mod
  use block_mod
  use gam_mod
  use dbyd
  use scratch

  implicit none

  integer :: l,i,mnum3,l1,lp,l1t,l2t,imt1,imt2,ivar,j,lnum3,nbxdim,ierr,iPE,tag
  real(IDP) :: epsq,oneos,betfc,betfc_f,betfc_alp,omcyd,omcydalp,beteom,betiom,coef
  character(len=80) :: formato
  integer, dimension(MPI_STATUS_SIZE) :: status
  real(IDP), dimension(:), allocatable :: ytrhs
  real(IDP), dimension(:), allocatable :: slhs,srhs,crhs,grwth,omega

  epsq=eps*eps
  oneos=0.0_IDP
  if (s > 0.0_IDP) oneos=1.0_IDP/s
  betfc=bet0/(2.*epsq)
  betfc_f=LcA2*bet0_f/(2.*epsq)
  betfc_alp=LcA2alp*bet0_alp/(2.*epsq)
  omcyd = omcy
  omcydalp = omcyalp
  beteom=dpres*bet0/(2*epsq*omcyd)
  betiom=(1-dpres)*bet0/(2*epsq*omcyd)

  yt=0.0_IDP

! put in r.h.s. of equations

! psi equation

  sd1=1.0_IDP
  call omc0(sd1,1,2,0,0,0,1.0_IDP)
  call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*sceq2(:,l)
  end do
  call b2lx(sceq1,-1,1,1,1,0,0,oneos)
  call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*sceq2(:,l)
  end do
  call b2lx(sceq1,1,1,1,0,1,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*grroj(:,l)
  end do
  call b2lx(sceq1,1,1,1,2,0,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=-2.*eta*feq*grtoj(:,l)
  end do
  call b2lx(sceq1,-1,1,1,1,1,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*gttoj(:,l)
  end do
  call b2lx(sceq1,1,1,1,0,2,0,oneos)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call b2lx(sceq2,1,1,1,0,1,0,oneos)

! Ion FLR effects

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then

     if (ext_prof == 1) then
        if (iflr_on == 1) then
           sd1=iflr*iflr
        else
           sd1=vthi*vthi*tieq/(epsq*omcyd*omcyd)
        end if
        do l=1,leqmax
           sceq5(:,l)=1.2533*sd1*vAlfven*vAlfven*bmod(:,l)/(vtherm_elecP*(feq-qqinv*cureq))
        end do
     else
        do l=1,leqmax
           sceq5(:,l)=1.2533*iflr*iflr*bmod(:,l)/(denseq*vtherm_elc*(feq-qqinv*cureq))
        end do
     end if
     call b2lx_landau_grad_parallel(sceq5,1,1,iq,0,0,0,1.0_IDP)

  end if

! Two fluid terms

  if (twofl_on == 1) then
     sd1=beteom/denseq
     call omc0(sd1,1,3,0,0,0,-1.0_IDP)
  end if

! u-zeta equation

  call dbydreq(sceq1,sqg,0.0_IDP,-1.0_IDP,0)
  call b2lx(sceq1,1,2,3,1,0,0,betfc)
  call dbydtheq(sceq2,sqg,1,0.0_IDP,1.0_IDP,0)
  call b2lx(sceq2,-1,2,3,0,1,0,betfc)

! fast ion coupling

  call b2lx(sceq1,1,2,5,1,0,0,betfc_f)
  call b2lx(sceq2,-1,2,5,0,1,0,betfc_f)

  if (alpha_on == 1) then
     call b2lx(sceq1,1,2,8,1,0,0,betfc_alp)
     call b2lx(sceq2,-1,2,8,0,1,0,betfc_alp)
  end if

! Shared equilibrium toroidal flow velocity for u-zeta equation

  call b2lx0(vzt_eq,2,4,0,0,1,-1.0_IDP)

  call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
  call omc(sceq2,-1,2,1,1,0,0,1.0_IDP)
  call grpareq(sceq1,sceq2,-1,0.0_IDP,1.0_IDP)
  call b2lx(sceq1,1,2,1,1,0,0,1.0_IDP)
  call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
  call omc(sceq2,1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,sceq2,1,0.0_IDP,1.0_IDP)
  call b2lx(sceq1,-1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,grtoj,-1,0.0_IDP,-2.0_IDP)
  call b2lx(sceq1,1,2,1,1,1,0,1.0_IDP)
  call grpareq(sceq1,gttoj,1,0.0_IDP,1.0_IDP)
  call b2lx(sceq1,-1,2,1,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call b2lx(sceq2,-1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,grroj,1,0.0_IDP,1.0_IDP)
  call b2lx(sceq1,-1,2,1,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-grtoj(:,l)
  end do
  call omc(sceq1,-1,2,1,1,1,0,2.0_IDP*1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=grroj(:,l)
  end do
  call omc(sceq1,1,2,1,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=gttoj(:,l)
  end do
  call omc(sceq1,1,2,1,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call omc(sceq2,1,2,1,0,1,0,1.0_IDP)
  call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)
  sd2=rinv*sd1
  call d2bydr20(sd2,cureq,1.0_IDP,-1.0_IDP,0)
  sd1=rinv*sd2/epsq
  call b2lx0(sd1,2,1,1,0,0,1.0_IDP)
  call dbydtheq(sceq2,bst,-1,0.0_IDP,-1.0_IDP,0)
  do l=1,leqmax
     sceq1(:,l)=r*sceq2(:,l)/epsq
  end do
  call dbydreq(sceq2,sceq1,0.0_IDP,-1.0_IDP,0)
  call b2lx(sceq2,1,2,1,1,0,0,1.0_IDP)
  call dbydtheq(sceq2,sceq1,1,0.0_IDP,1.0_IDP,0)
  call b2lx(sceq2,-1,2,1,0,1,0,1.0_IDP)

! Ion FLR effects

  ! if (iflr_on == 1) then
  !    coef=omegar*iflr*iflr
  !    call b2lx_dlsq(2,4,coef)
  ! end if

! Electron-ion Landau damping

  if (ieldamp_on == 1) then   
     call b2lxjl(eildrr,1,2,2,0,2,0,eilnd)
     call b2lxjl(eildrt,-1,2,2,1,1,0,eilnd)
     call b2lxjl(eildrz,-1,2,2,0,1,1,eilnd)
     call b2lxjl(eildtt,1,2,2,2,0,0,eilnd)
     call b2lxjl(eildtz,1,2,2,1,0,1,eilnd)
     call b2lxjl(eildzz,1,2,2,0,0,2,eilnd)
     call b2lxjl(eildr,1,2,2,0,1,0,eilnd)
     call b2lxjl(eildt,-1,2,2,1,0,0,eilnd)
     call b2lxjl(eildz,-1,2,2,0,0,1,eilnd)
  end if

! Two fluid terms

  if (twofl_on == 1) then
     sd2=feq-qqinv*cureq
     do l=1,leqmax
        sceq1(:,l)=feq*dpreqdr*jbgrr(:,l)/sd2
        sceq2(:,l)=feq*dpreqdr*jbgrt(:,l)/sd2
        sceq3(:,l)=feq*dpreqdr*jbgtt(:,l)/sd2
     end do
     call b2lx(sceq1,1,2,2,3,0,0,-betiom)
     call b2lx(sceq2,-1,2,2,2,1,0,2*betiom)
     call b2lx(sceq3,1,2,2,1,2,0,-betiom)
     call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,3)
     call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,3)
     do l=1,leqmax
        sceq4(:,l)=sceq4(:,l)+rinv*sceq2(:,l)
     end do
     call b2lx(sceq4,-1,2,2,2,0,0,-betiom)
     call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,3)
     call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,3)
     call b2lx(sceq4,1,2,2,1,1,0,betiom)
     do l=1,leqmax
        sceq1(:,l)=rinv*cureq*dpreqdr*jbgrr(:,l)/sd2
        sceq2(:,l)=rinv*cureq*dpreqdr*jbgrt(:,l)/sd2
        sceq3(:,l)=rinv*cureq*dpreqdr*jbgtt(:,l)/sd2
     end do
     call b2lx(sceq1,1,2,2,2,0,1,betiom)
     call b2lx(sceq2,-1,2,2,1,1,1,-2*betiom)
     call b2lx(sceq3,1,2,2,0,2,1,betiom)
     call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,4)
     call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,4)
     call b2lx(sceq4,-1,2,2,1,0,1,betiom)
     call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,4)
     call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,4)
     do l=1,leqmax
        sceq4(:,l)=sceq4(:,l)-rinv*sceq3(:,l)
     end do
     call b2lx(sceq4,1,2,2,0,1,1,-betiom)
     do l=1,leqmax
        sceq1(:,l)=rinv*cureq*dpreqdr*dgrrz(:,l)
        sceq2(:,l)=rinv*cureq*dpreqdr*dgrtz(:,l)
        sceq3(:,l)=rinv*cureq*dpreqdr*dgttz(:,l)
     end do
     call b2lx(sceq1,-1,2,2,2,0,0,0.5*betiom)
     call b2lx(sceq2,1,2,2,1,1,0,-betiom)
     call b2lx(sceq3,-1,2,2,0,2,0,0.5*betiom)
     sd3=rinv*feq*dpreqdr
     sd3(0)=2*feq(0)*(preq(1)-preq(0))/(r(1)*r(1))
     do l=1,leqmax
        sceq1(:,l)=sd3*dgrrt(:,l)-rinv*cureq*dpreqdr*dgrrz(:,l)
        sceq2(:,l)=feq*dpreqdr*dgttr(:,l)+2*sd3*jbgtt(:,l)/sd2-rinv*cureq*dpreqdr*dgrtz(:,l)
     end do
     call dbydtheq(sceq3,sceq1,-1,0.0_IDP,1.0_IDP,2)
     call dbydreq(sceq3,sceq2,1.0_IDP,-1.0_IDP,2)
     call dbydr0(sd4,qqinv,0.0_IDP,1.0_IDP,0)
     call dbydr0(sd5,cureq,0.0_IDP,1.0_IDP,0)
     call dbydreq(sceq2,jsq,0.0_IDP,0.5_IDP,0)
     do l=1,leqmax
        sceq1(:,l)=cureq*dpreqdr*(qqinv*(rinv*dgrtt(:,l)-dgttr(:,l))-(sd4+2*rinv*qqinv)*jbgtt(:,l)/sd2+ &
                   rinv*(dbsjtbj(:,l)-rinv*(cureq*sceq2(:,l)+sd5*jsq(:,l))/sd2)/(eps*eps))
     end do
     call dbydreq(sceq3,sceq1,1.0_IDP,-1.0_IDP,4)
     call b2lx(sceq3,1,2,2,1,0,0,-0.5*betiom)
     call dbydtheq(sceq3,sceq1,1,0.0_IDP,1.0_IDP,4)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*(feq*dgttt(:,l)-cureq*dgttz(:,l))
        sceq2(:,l)=2*sd3*(dgrtt(:,l)-jbgtt(:,l)/sd2)-dpreqdr*(feq*dgttr(:,l)+rinv*cureq*dgrtz(:,l))
     end do
     call dbydreq(sceq5,sceq1,0.0_IDP,1.0_IDP,3)
     do l=1,leqmax
        sceq3(:,l)=sceq3(:,l)+rinv*sceq5(:,l)
     end do
     call dbydtheq(sceq3,sceq2,1,1.0_IDP,-1.0_IDP,2)
     call b2lx(sceq3,-1,2,2,0,1,0,-0.5*betiom)
  end if

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call b2lx0_dlsq(sd1,2,4,stdifu)
  else
     call b2lx0_dlsqnr(fctr_dif,2,4,stdifun,dfctr_difdr)
  end if

! p equation

  call b2lx0(dpreqdr,3,2,1,0,0,1.0_IDP)

!  do l=1,leqmax
!     sceq1(:,l)=preq*djroj(:,l)
!  end do
!  call b2lx(sceq1,1,3,2,1,0,0,gamma)
!  do l=1,leqmax
!     sceq1(:,l)=-preq*djtoj(:,l)
!  end do
!  call b2lx(sceq1,-1,3,2,0,1,0,gamma)

  sd2=feq/(feq-qqinv*cureq)
  sd3=cureq/(feq-qqinv*cureq)
  do l=1,leqmax
     sceq1(:,l)=-preq*(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
  end do
  call b2lx(sceq1,-1,3,2,0,1,0,gamma)
  do l=1,leqmax
     sceq1(:,l)=-preq*(r*dbsjzoj(:,l)-sd2*djroj(:,l))
  end do
  call b2lx(sceq1,1,3,2,1,0,0,gamma)
  do l=1,leqmax
     sceq1(:,l)=-preq*(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
  end do
  call b2lx(sceq1,1,3,2,0,0,1,gamma)
  call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
  sd2=preq*sd4
  call b2lx0(sd2,3,2,1,0,0,gamma)
  call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
  sd2=-rinv*preq*sd4
  call b2lx0(sd2,3,2,0,0,1,gamma)

! parallel thermal velocity term  

  do l=1,leqmax
     sceq1(:,l)=preq*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call omc(sceq1,1,3,7,0,0,0,-gamma)
  call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=preq*sceq2(:,l)/(feq-qqinv*cureq)
  end do
  call b2lx(sceq1,-1,3,7,0,0,0,gamma)

! Shared equilibrium toroidal flow velocity for pressure equation

  call b2lx0(vzt_eq,3,3,0,0,1,-1.0_IDP)

! Two fluid terms

  if (twofl_on == 1) then
     sd2=gamma*betiom*preq/denseq
     sd3=sd2/(feq-qqinv*cureq)
     call dbydr0(sd4,feq,0.0_IDP,1.0_IDP,0)
     sd4=sd3*sd4
     call b2lx0(sd4,3,3,1,0,0,1.0_IDP)
     call dbydr0(sd4,cureq,0.0_IDP,1.0_IDP,0)
     sd4=rinv*sd3*sd4
     call b2lx0(sd4,3,3,0,0,1,-1.0_IDP)
     call dbydzteq(sceq1,bst,-1,0.0_IDP,1.0_IDP)
     call dbydtheq(sceq2,bst,-1,0.0_IDP,1.0_IDP,0)
     do l=1,leqmax
        sceq1(:,l)=r*sd3*sceq1(:,l)
        sceq2(:,l)=r*sd3*sceq2(:,l)
     end do
     call b2lx(sceq1,1,3,3,1,0,0,-1.0_IDP)
     call b2lx(sceq2,1,3,3,0,0,1,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=sd2*omdr(:,l)
        sceq2(:,l)=sd2*omdt(:,l)
        sceq3(:,l)=sd2*omdz(:,l)
     end do
     call b2lx(sceq1,-1,3,3,0,1,0,-2.0_IDP)
     call b2lx(sceq2,1,3,3,1,0,0,-2.0_IDP)
     call b2lx(sceq3,1,3,3,0,0,1,-2.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*sd3*grtoj(:,l)
        sceq2(:,l)=dpreqdr*sd3*gttoj(:,l)
     end do
     call omc(sceq1,-1,3,1,1,0,0,-epsq)
     call omc(sceq2,1,3,1,0,1,0,epsq)
     sd4=rinv*dpreqdr*sd3*cureq
     call b2lx0(sd4,3,1,1,1,0,-1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=r*dpreqdr*sd3*bst(:,l)
     end do
     call b2lx(sceq1,-1,3,1,2,0,0,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*sd2*dgrtp(:,l)
        sceq2(:,l)=dpreqdr*sd2*dgttp(:,l)
     end do
     call b2lx(sceq1,1,3,1,1,0,0,-epsq)
     call b2lx(sceq2,-1,3,1,0,1,0,epsq)
     do l=1,leqmax
        sceq1(:,l)=rinv*dpreqdr*sd3*cureq*djtoj(:,l)
        sceq2(:,l)=dpreqdr*sd2*dbsjtoj(:,l)
     end do
     call b2lx(sceq1,-1,3,1,0,1,0,-1.0_IDP)
     call b2lx(sceq2,1,3,1,1,0,0,1.0_IDP)
  end if

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call b2lx0_dlsq(sd1,3,3,stdifp)
  else
     call b2lx0_dlsqnr(fctr_dif,3,3,stdifpn,dfctr_difdr)
  end if

! u-zeta expression

  do l=1,leqmax
     sceq2(:,l)=denseq*jbgrt(:,l)
  end do
  call b2lx(sceq2,-1,4,2,1,1,0,-2.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=denseq*jbgrr(:,l)
  end do
  call b2lx(sceq2,1,4,2,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=denseq*jbgtt(:,l)
  end do
  call b2lx(sceq2,1,4,2,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=rinv*sceq2(:,l)
  end do
  call b2lx(sceq1,1,4,2,0,1,0,1.0_IDP)
  call clgam(sceq2,1,1,2,0.0_IDP,-1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=denseq*sceq2(:,l)-denseqr*jbgrt(:,l)
  end do
  call b2lx(sceq1,-1,4,2,1,0,0,1.0_IDP)
  call clgam(sceq2,1,2,2,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=denseq*sceq2(:,l)+denseqr*jbgtt(:,l)
  end do
  call b2lx(sceq1,1,4,2,0,1,0,1.0_IDP)
  sd1=-1.0_IDP
  call b2lx0(sd1,4,4,0,0,0,1.0_IDP)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NBI particle effects   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! Load 1st Omega-d terms in fast ion density equation:

  if (trapped_on .eq. 0) then
     do l=1,leqmax
        sceq1(:,l)=vfova2*omdr(:,l)/(epsq*omcyd)
        sceq2(:,l)=vfova2*omdt(:,l)/(epsq*omcyd)
        sceq3(:,l)=vfova2*omdz(:,l)/(epsq*omcyd)
     end do
  else
     do l=1,leqmax
        sceq1(:,l)=vfova2*omdrprp(:,l)/(epsq*omcyd)
        sceq2(:,l)=vfova2*omdtprp(:,l)/(epsq*omcyd)
        sceq3(:,l)=vfova2*omdzprp(:,l)/(epsq*omcyd)
     end do
  end if
  
  call b2lx(sceq1,-1,5,5,0,1,0,-1.0_IDP)
  call b2lx(sceq2,1,5,5,1,0,0,-1.0_IDP)
  call b2lx(sceq3,1,5,5,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion density equation

  call b2lx0(vzt_eq,5,5,0,0,1,-1.0_IDP)

! Load 1st Omega-d terms in fast ion parallel velocity equation:

  call b2lx(sceq1,-1,6,6,0,1,0,-1.0_IDP)
  call b2lx(sceq2,1,6,6,1,0,0,-1.0_IDP)
  call b2lx(sceq3,1,6,6,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

  call b2lx0(vzt_eq,6,6,0,0,1,-1.0_IDP)

! Load 2nd Omega-d terms in fast ion density equation:

  if (trapped_on .eq. 0) then
     do l=1,leqmax
        sceq1(:,l)=nfeq(:)*omdr(:,l)
        sceq2(:,l)=nfeq(:)*omdt(:,l)
        sceq3(:,l)=nfeq(:)*omdz(:,l)
     end do
  else
     do l=1,leqmax
        sceq1(:,l)=nfeq(:)*omdrprp(:,l)/vfova
        sceq2(:,l)=nfeq(:)*omdtprp(:,l)/vfova
        sceq3(:,l)=nfeq(:)*omdzprp(:,l)/vfova
     end do
  end if

  call b2lx(sceq1,-1,5,2,0,1,0,-1.0_IDP)
  call b2lx(sceq2,1,5,2,1,0,0,-1.0_IDP)
  call b2lx(sceq3,1,5,2,0,0,1,-1.0_IDP)

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call b2lx0_dlsq(sd1,5,5,stdifnf)
  else
     call b2lx0_dlsqnr(fctr_dif,5,5,stdifnfn,dfctr_difdr)
  end if

! Load remaining terms in fast ion density equation

! Parallel gradient term

  do l=1,leqmax
     sceq1(:,l)=nfeq*bmod(:,l)/(feq-qqinv*cureq)
  end do
  
  call omc(sceq1,1,5,6,0,0,0,-1.0_IDP)

! Omega* term

  sd1=dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
  sd2=dnfeqdr*feq/(feq-qqinv*cureq)
  
  call b2lx0(sd1,5,2,0,0,1,1.0_IDP)
  call b2lx0(sd2,5,2,1,0,0,-1.0_IDP)

! EP FLR effects

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then
     sd1=-dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
     sd2=-dnfeqdr*feq/(feq-qqinv*cureq)
     call b2lx0(sd1,5,iw,0,0,1,1.0_IDP)
     call b2lx0(sd2,5,iw,1,0,0,-1.0_IDP)
     ! sd1=epsq*omcyd*omegar*nfeq/vfova2
     ! call b2lx0(sd1,5,iw,0,0,0,1.0_IDP)
  end if

! Load remaining terms in fast ion parallel velocity equation

! Landau closure term

  do l=1,leqmax
     sceq1(:,l)=1.414213*LcA1*vfova*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call b2lx_landau_grad_parallel(sceq1,1,6,6,0,0,0,1.0_IDP)

! Parallel gradient terms

  do l=1,leqmax
     sceq1(:,l)=2*LcA0*vfova2*bmod(:,l)/(nfeq*(feq-qqinv*cureq))
  end do
  call omc(sceq1,1,6,5,0,0,0,-1.0_IDP)

  do l=1,leqmax
     sceq1(:,l)=2*LcA0*epsq*omcyd*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call omc(sceq1,1,6,2,0,0,0,-1.0_IDP)

! Omega* term

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then

!  EP FLR effects

     sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
     sd2=vfova2*dnfeqdr*rinv*feq/(nfeq*(feq-qqinv*cureq))
     call b2lx0(sd1,6,ix1,0,0,0,1.0_IDP)
     call b2lx0(sd2,6,ix2,0,0,0,-1.0_IDP)
   
  else 

     sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
     sd2=vfova2*dnfeqdr*feq/(nfeq*(feq-qqinv*cureq))
     call b2lx0(sd1,6,1,0,0,1,1.0_IDP)
     call b2lx0(sd2,6,1,1,0,0,-1.0_IDP)
   
  end if 

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call b2lx0_dlsq(sd1,6,6,stdifvf)
  else
     call b2lx0_dlsqnr(fctr_dif,6,6,stdifvfn,dfctr_difdr)
  end if

! End of fast ion moment equations

! Load terms of the thermal moment of energetic particles

! Pressure gradient term

  do l=1,leqmax
     sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
  end do
  call omc(sceq1,1,7,3,0,0,0,-1.0_IDP)
  
! Magnetic field perturbation term  
  
  do l=1,leqmax
     sceq2(:,l)=sceq1(:,l)*dpreqdr
  end do  
  call b2lx(sceq2,1,7,1,1,0,0,1.0_IDP)  
  
! Shared equilibrium toroidal flow velocity for the thermal moment of energetic particles equation

  call b2lx0(vzt_eq,7,7,0,0,1,-1.0_IDP)

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call b2lx0_dlsq(sd1,7,7,stdifv)
  else
     call b2lx0_dlsqnr(fctr_dif,7,7,stdifvn,dfctr_difdr)
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Alpha particle effects   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (alpha_on == 1) then

! Load 1st Omega-d terms in fast ion density equation:

     do l=1,leqmax
        sceq1(:,l)=valphaova2*omdr(:,l)/(epsq*omcydalp)
        sceq2(:,l)=valphaova2*omdt(:,l)/(epsq*omcydalp)
        sceq3(:,l)=valphaova2*omdz(:,l)/(epsq*omcydalp)
     end do
  
     call b2lx(sceq1,-1,8,8,0,1,0,-1.0_IDP)
     call b2lx(sceq2,1,8,8,1,0,0,-1.0_IDP)
     call b2lx(sceq3,1,8,8,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion density equation

     call b2lx0(vzt_eq,8,8,0,0,1,-1.0_IDP)

! Load 1st Omega-d terms in fast ion parallel velocity equation:

     call b2lx(sceq1,-1,9,9,0,1,0,-1.0_IDP)
     call b2lx(sceq2,1,9,9,1,0,0,-1.0_IDP)
     call b2lx(sceq3,1,9,9,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

     call b2lx0(vzt_eq,9,9,0,0,1,-1.0_IDP)

! Load 2nd Omega-d terms in fast ion density equation:

     do l=1,leqmax
        sceq1(:,l)=nalpeq(:)*omdr(:,l)
        sceq2(:,l)=nalpeq(:)*omdt(:,l)
        sceq3(:,l)=nalpeq(:)*omdz(:,l)
     end do

     call b2lx(sceq1,-1,8,2,0,1,0,-1.0_IDP)
     call b2lx(sceq2,1,8,2,1,0,0,-1.0_IDP)
     call b2lx(sceq3,1,8,2,0,0,1,-1.0_IDP)

! diffusion term added
     sd1=1.0_IDP
     if (difnr_on == 0) then
        call b2lx0_dlsq(sd1,8,8,stdifnalp)
     else
        call b2lx0_dlsqnr(fctr_dif,8,8,stdifnalpn,dfctr_difdr)
     end if

! Load remaining terms in fast ion density equation

! Parallel gradient term

     do l=1,leqmax
        sceq1(:,l)=nalpeq*bmod(:,l)/(feq-qqinv*cureq)
     end do
  
     call omc(sceq1,1,8,9,0,0,0,-1.0_IDP)

! Omega* term

     sd1=dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
     sd2=dnalpeqdr*feq/(feq-qqinv*cureq)
  
     call b2lx0(sd1,8,2,0,0,1,1.0_IDP)
     call b2lx0(sd2,8,2,1,0,0,-1.0_IDP)
   

! EP FLR effects

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then
        sd1=-dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
        sd2=-dnalpeqdr*feq/(feq-qqinv*cureq)
        call b2lx0(sd1,8,iwa,0,0,1,1.0_IDP)
        call b2lx0(sd2,8,iwa,1,0,0,-1.0_IDP)
        ! sd1=epsq*omcydalp*omegar*nalpeq/valphaova2
        ! call b2lx0(sd1,8,iwa,0,0,0,1.0_IDP)
     end if

! Load remaining terms in fast ion parallel velocity equation

! Landau closure term

     do l=1,leqmax
        sceq1(:,l)=1.414213*LcA1alp*valphaova*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call b2lx_landau_grad_parallel(sceq1,1,9,9,0,0,0,1.0_IDP)

! Parallel gradient term

     do l=1,leqmax
        sceq1(:,l)=2*LcA0alp*valphaova2*bmod(:,l)/(nalpeq*(feq-qqinv*cureq))
     end do
     call omc(sceq1,1,9,8,0,0,0,-1.0_IDP)

     do l=1,leqmax
        sceq1(:,l)=2*LcA0alp*epsq*omcydalp*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call omc(sceq1,1,9,2,0,0,0,-1.0_IDP)

! Omega* term

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then

!  EP FLR effects

        sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
        sd2=valphaova2*dnalpeqdr*rinv*feq/(nalpeq*(feq-qqinv*cureq))
        call b2lx0(sd1,9,ix1a,0,0,0,1.0_IDP)
        call b2lx0(sd2,9,ix2a,0,0,0,-1.0_IDP)

     else

        sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
        sd2=valphaova2*dnalpeqdr*feq/(nalpeq*(feq-qqinv*cureq))
        call b2lx0(sd1,9,1,0,0,1,1.0_IDP)
        call b2lx0(sd2,9,1,1,0,0,-1.0_IDP)

     end if 

! diffusion term added
     sd1=1.0_IDP
     if (difnr_on == 0) then
        call b2lx0_dlsq(sd1,9,9,stdifvalp)
     else
        call b2lx0_dlsqnr(fctr_dif,9,9,stdifvalpn,dfctr_difdr)
     end if

! End of alpha particles moment equations

  end if

  nbxdim=lmx*mj
  if (myPE < numPElm1) then
     allocate (ytrhs(nskpxn(n_start)+1:nskpxn(n_end+1)))
  else if (myPE== numPElm1) then
     allocate (ytrhs(nskpxn(n_start)+1:nbxdim))
  else
     allocate (ytrhs(1))
  end if

  ytrhs=yt
  yt=0.0_IDP

! put in l.h.s. of equation

! psi equation

  sd1=1.0_IDP
  call b2lx0(sd1,1,1,0,0,0,1.0_IDP)

! u-zeta equation

  call b2lx0(sd1,2,4,0,0,0,1.0_IDP)

! Ion FLR effects

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then
     if (iflr_on == 1) then
        sd2=-iflr*iflr
     else
        sd2=-vthi*vthi*tieq/(epsq*omcyd*omcyd)
     end if
     call b2lx_dlsq(2,4,sd2)
  end if

! p equation

  call b2lx0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density moment equation

  call b2lx0(sd1,5,5,0,0,0,1.0_IDP)

! EP FLR effects

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then
     if (epflr_on == 1) then
        sd2=nfeq/(r_epflr*r_epflr*omcyd)
     else if (epflr_on == 2) then
        sd2=epsq*omcyd*nfeq/vfova2
     end if
     call b2lx0(sd2,5,iw,0,0,0,1.0_IDP)
  end if
  
! fast ion parallel velocity moment equation

  call b2lx0(sd1,6,6,0,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-epsq*omcyd*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call b2lx(sceq1,1,6,1,0,0,0,1.0_IDP)

! thermal moment equation

  call b2lx0(sd1,7,7,0,0,0,1.0_IDP)  

  if (alpha_on == 1) then

! fast ion density moment equation

     call b2lx0(sd1,8,8,0,0,0,1.0_IDP)

! EP FLR effects

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then
        if (epflr_on == 1) then
           sd2=nalpeq/(r_epflralp*r_epflralp*omcydalp)
        else if (epflr_on == 2) then
           sd2=epsq*omcydalp*nalpeq/valphaova2
        end if
        call b2lx0(sd2,8,iwa,0,0,0,1.0_IDP)
     end if

! fast ion parallel velocity moment equation

     call b2lx0(sd1,9,9,0,0,0,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=-epsq*omcydalp*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call b2lx(sceq1,1,9,1,0,0,0,1.0_IDP)
  
  end if

  allocate (slhs(nvar*lmaxn),srhs(nvar*lmaxn),crhs(nvar*lmaxn),grwth(nvar*lmaxn),omega(nvar*lmaxn))

  do i=n_start,n_end
     mnum3=noeqn*mnumn(i)
     do l1t=1,nvar*mnumn(i)
        l=l1t+nvar*lnumn(i-1)
        slhs(l)=0.
        srhs(l)=0.
        crhs(l)=0.
        do j=1,mjm1
           imt1=l1t+mnum3*(j-1)+nskpxn(i)
           srhs(l)=srhs(l)+ytrhs(imt1)*yt(imt1)
           slhs(l)=slhs(l)+yt(imt1)**2
        end do
        ivar=l1t/mnumn(i)+1
        lp=mod(l1t,mnumn(i))
        if (lp == 0) then
           ivar=ivar-1
           lp=mnumn(i)
        end if
        l1=lp+lnumn(i-1)
        if (lo(l1) == 0) cycle
        l2t=(ivar-1)*mnumn(i)+lo(l1)-lnumn(i-1)
        do j=1,mjm1
           imt1=l1t+mnum3*(j-1)+nskpxn(i)
           imt2=l2t+mnum3*(j-1)+nskpxn(i)
           crhs(l)=crhs(l)+ytrhs(imt2)*yt(imt1)
        end do
        if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == 9)  crhs(l)=-crhs(l)
     end do
  end do
  do i=n_start,n_end
     mnum3=nvar*mnumn(i)
     do l1t=1,mnum3
        l=l1t+nvar*lnumn(i-1)
        grwth(l)=0.
        omega(l)=0.
        ivar=l1t/mnumn(i)+1
        l2t=mod(l1t,mnumn(i))
        if (l2t == 0) then
           ivar=ivar-1
           l2t=mnumn(i)
        end if
        l1=l2t+lnumn(i-1)
        if (lo(l1) == 0) then
           if (slhs(l) /= 0) grwth(l)=srhs(l)/slhs(l)
        else if (signl(lln(l1)) > 0) then
           lp=(ivar-1)*mnumn(i)+lo(l1)+(nvar-1)*lnumn(i-1)
           if (slhs(l) /= 0.0 .or. slhs(lp) /= 0.0) then
              grwth(l)=(srhs(l)+srhs(lp))/(slhs(l)+slhs(lp))
              omega(l)=(crhs(l)-crhs(lp))/(slhs(l)+slhs(lp))
           end if
        end if
     end do
     if (myPE /= 0) then
        tag=i
        call MPI_SEND(grwth(nvar*lnumn(i-1)+1),nvar*mnumn(i),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
        tag=nnum+i
        call MPI_SEND(omega(nvar*lnumn(i-1)+1),nvar*mnumn(i),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
     end if
  end do
  if (myPE == 0) then
     do iPE=1,numPElm1
        do i=n_st(iPE),n_nd(iPE)
           tag=i
           call MPI_RECV(grwth(nvar*lnumn(i-1)+1),nvar*mnumn(i),MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,status,ierr)
           tag=nnum+i
           call MPI_RECV(omega(nvar*lnumn(i-1)+1),nvar*mnumn(i),MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,status,ierr)
        end do
     end do
  end if

  if (myPE == 0) then
     write(6,'(/)')
     do i=1,nnum
        do l1=1,mnumn(i)
           lnum3=nvar*lnumn(i-1)
           lp=l1+lnumn(i-1)
           if (signl(lln(lp)) < 0) cycle
           do ivar=1,nvar
              l=lnum3+(ivar-1)*mnumn(i)+l1
              select case (ivar)
                 case (1)
                    formato='(" psi     : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (2)
                    formato='(" phi     : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (3)
                    formato='(" pr      : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (4)
                    formato='(" uzt     : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (5)
                    formato='(" nfast   : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (6)
                    formato='(" vfast   : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (7)
                    formato='(" vthfast : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (8)
                    formato='(" nalpha  : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
                 case (9)
                    formato='(" valpha  : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)'
              end select
              if (ivar == 4 .or. grwth(l) == 0.0_IDP) cycle
              if (ivar == 7 .and. s == 0.0_IDP .and. gamma == 0.0_IDP .and. twofl_on == 0 .and. iflr_on == 0) cycle
              write(6,formato) mm(lln(lp)),nn(lln(lp)),grwth(l),omega(l)
           end do
        end do
     end do
  end if
  
  deallocate (omega,grwth,crhs,srhs,slhs)
  deallocate (ytrhs)

end subroutine lincheck
