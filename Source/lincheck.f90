	subroutine lincheck

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none
		integer :: j,l,l1,ivar,lp,nvar
		real(IDP) :: epsq,oneos,betfc,betfc_f,omcyd,grwth_avg,omega_r_avg,sum_iter,devation_growth,devation_freq, &
                     betfc_alp,omcydalp,beteom,betiom,coef,deviation_freq,deviation_growth		
		real(IDP), dimension(:), allocatable :: slhs,srhs,crhs,grwth,omega_r
		real(IDP), dimension(:,:), allocatable :: ytrhs			

		interface
			subroutine clgam(rgam,na,nb,nc,c1,c2)
				use param
				implicit none
				integer :: na,nb,nc
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: rgam
			end subroutine clgam
			subroutine omc(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine omc
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
			subroutine b2lx_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx_landau_grad_parallel
			subroutine omc0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine omc0
			subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine b2lx0
			subroutine b2lx_dlsq(ieqn,ivar,coef)
				use param
				implicit none
				integer :: ieqn,ivar
				real(IDP) :: coef
			end subroutine b2lx_dlsq
			subroutine dbydth(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydth
			subroutine dbydzt(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzt			
			subroutine dbydr(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydr			
			subroutine dbydtheq(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydtheq
			subroutine grpareq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine grpareq
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydzteq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzteq
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine d2bydr20(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine d2bydr20		
		end interface

!		In this subroutine we calculate the instability growth rate and frequency
!		using the subroutines b2lx and b2lx0		
		
		epsq=eps*eps
		oneos=1.0_IDP/s
		betfc=bet0/(2.*epsq)
		betfc_f=LcA2*bet0_f/(2.*epsq)
		omcyd = omcy
		beteom=dpres*bet0/(2*epsq*omcyd)
		betiom=(1-dpres)*bet0/(2*epsq*omcyd)	

                if(alpha_on .eq. 1) then
		  betfc_alp=LcA2alp*bet0_alp/(2.*epsq)
		  omcydalp = omcyalp
		end if	

		yt=0.0_IDP

!  put in r.h.s. of equations


!  psi equation

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

!  Ion FLR effects

		if (iflr_on == 1) then

			if (ext_prof == 1) then
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*vAlfven*vAlfven*bmod(:,l)/(vtherm_ionP*(feq-qqinv*cureq))
				end do
			else
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*bmod(:,l)/(denseq*vtherm_ion*(feq-qqinv*cureq))
				end do
			end if
			call b2lx_landau_grad_parallel(sceq5,1,1,iq,0,0,0,1.0_IDP)

		end if
		 
! Two fluid terms

        if(twofl_on .eq. 1) then   

			sd1=beteom/denseq
			call omc0(sd1,1,3,0,0,0,-1.0_IDP)					
			
		end if	
		
!  u-zeta equation

		call dbydreq(sceq1,sqg,0.0_IDP,-1.0_IDP,0)
		call b2lx(sceq1,1,2,3,1,0,0,betfc)
		call dbydtheq(sceq2,sqg,1,0.0_IDP,1.0_IDP,0)
		call b2lx(sceq2,-1,2,3,0,1,0,betfc)

!  fast ion coupling

		call b2lx(sceq1,1,2,5,1,0,0,betfc_f)
		call b2lx(sceq2,-1,2,5,0,1,0,betfc_f)

		if (alpha_on == 1) then
			call b2lx(sceq1,1,2,8,1,0,0,betfc_alp)
			call b2lx(sceq2,-1,2,8,0,1,0,betfc_alp)
		end if	
		
!  Sheared equilibrium toroidal flow velocity for u-zeta equation

        call b2lx0(vzt_eq,2,4,0,0,1,-1.0_IDP) 

!  Sheared equilibrium poloidal flow velocity for u-zeta equation

                vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,2,4,1,0,0,-1.0_IDP) 			
		
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
		
!  Ion FLR effects

		if (iflr_on == 1) then

			coef=omegar*iflr*iflr
			call b2lx_dlsq(2,4,coef)

		end if

!   Electron-ion Landau damping

        if(ieldamp_on .eq. 1 .and. ext_prof .eq. 1) then   
		  
		    norm_eildump=-(1-dpres)*bet0/(2.*epsq*omegar)		
			call b2lx(eildrr,1,2,2,0,2,0,norm_eildump)
			call b2lx(eildrt,-1,2,2,1,1,0,norm_eildump)
			call b2lx(eildrz,-1,2,2,0,1,1,norm_eildump)
			call b2lx(eildtt,1,2,2,2,0,0,norm_eildump)
			call b2lx(eildtz,1,2,2,1,0,1,norm_eildump)
			call b2lx(eildzz,1,2,2,0,0,2,norm_eildump)
			call b2lx(eildr,1,2,2,0,1,0,norm_eildump)
			call b2lx(eildt,-1,2,2,1,0,0,norm_eildump)
			call b2lx(eildz,-1,2,2,0,0,1,norm_eildump)			
			
		end if

! Two fluid terms

        if(twofl_on .eq. 1) then   

			call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
			sd2=feq-qqinv*cureq
			do l=1,leqmax
				sceq1(:,l)=feq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=feq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=feq*sd1*jbgtt(:,l)/sd2
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
				sceq1(:,l)=rinv*cureq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=rinv*cureq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=rinv*cureq*sd1*jbgtt(:,l)/sd2
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
				sceq1(:,l)=rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=rinv*cureq*sd1*dgrtz(:,l)
				sceq3(:,l)=rinv*cureq*sd1*dgttz(:,l)
			end do
			call b2lx(sceq1,-1,2,2,2,0,0,0.5*betiom)
			call b2lx(sceq2,1,2,2,1,1,0,-betiom)
			call b2lx(sceq3,-1,2,2,0,2,0,0.5*betiom)
			sd3=rinv*feq*sd1
			sd3(0)=2*feq(0)*(preq(1)-preq(0))/(r(1)*r(1))
			do l=1,leqmax
				sceq1(:,l)=sd3*dgrrt(:,l)-rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=feq*sd1*dgttr(:,l)+2*sd3*jbgtt(:,l)/sd2-rinv*cureq*sd1*dgrtz(:,l)
			end do
			call dbydtheq(sceq3,sceq1,-1,0.0_IDP,1.0_IDP,2)
			call dbydreq(sceq3,sceq2,1.0_IDP,-1.0_IDP,2)
			call dbydr0(sd4,qqinv,0.0_IDP,1.0_IDP,0)
			call dbydr0(sd5,cureq,0.0_IDP,1.0_IDP,0)
			call dbydreq(sceq2,jsq,0.0_IDP,0.5_IDP,0)
			do l=1,leqmax
				sceq1(:,l)=cureq*sd1*(qqinv*(rinv*dgrtt(:,l)-dgttr(:,l))-(sd4+2*rinv*qqinv)*jbgtt(:,l)/sd2+ &
						      rinv*(dbsjtbj(:,l)-rinv*(cureq*sceq2(:,l)+sd5*jsq(:,l))/sd2)/(eps*eps))
			end do
			call dbydreq(sceq3,sceq1,1.0_IDP,-1.0_IDP,4)
			call b2lx(sceq3,1,2,2,1,0,0,-0.5*betiom)
			call dbydtheq(sceq3,sceq1,1,0.0_IDP,1.0_IDP,4)
			do l=1,leqmax
				sceq1(:,l)=sd1*(feq*dgttt(:,l)-cureq*dgttz(:,l))
				sceq2(:,l)=2*sd3*(dgrtt(:,l)-jbgtt(:,l)/sd2)-sd1*(feq*dgttr(:,l)+rinv*cureq*dgrtz(:,l))
			end do
			call dbydreq(sceq5,sceq1,0.0_IDP,1.0_IDP,3)
			do l=1,leqmax
				sceq3(:,l)=sceq3(:,l)+rinv*sceq5(:,l)
			end do
			call dbydtheq(sceq3,sceq2,1,1.0_IDP,-1.0_IDP,2)
			call b2lx(sceq3,-1,2,2,0,1,0,-0.5*betiom)
			
		end if
		
! diffusion term added

		call b2lx_dlsq(2,4,stdifu)

! end diffusion term

! p equation

		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
		sd2=-rinv*sd1*cureq(:)/(feq(:)-qqinv(:)*cureq(:))	
		call b2lx0(sd2,3,2,0,0,1,1.0_IDP)		
		sd2=sd1*feq(:)/(feq(:)-qqinv(:)*cureq(:))			
		call b2lx0(sd2,3,2,1,0,0,1.0_IDP)			
		
		call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)		
                sd2=-rinv*sd1*preq/(feq-qqinv*cureq)
                call b2lx0(sd2,3,2,0,0,1,gamma)		
		call dbydr0(sd1,feq,0.0_IDP,1.0_IDP,0)		
                sd2=sd1*preq/(feq-qqinv*cureq)
                call b2lx0(sd2,3,2,1,0,0,gamma)		
                call dbydtheq(sceq1,bst,-1,0.0_IDP,1.0_IDP,0)					
		do l=1,leqmax
			sceq1(:,l)=sceq1(:,l)*r*preq/(feq(:)-qqinv(:)*cureq(:))
		end do	
		call b2lx(sceq1,1,3,2,0,0,1,gamma)		
		call dbydzteq(sceq1,bst,-1,0.0_IDP,1.0_IDP)			
		do l=1,leqmax
			sceq1(:,l)=-r*sceq1(:,l)*preq/(feq(:)-qqinv(:)*cureq(:))
		end do	
		call b2lx(sceq1,1,3,2,1,0,0,gamma)		

		sd1=1/(feq(:)-qqinv(:)*cureq(:))		
		call dbydr0(sd2,sd1,0.0_IDP,1.0_IDP,0)			
		sd3=-cureq*rinv*sd2*preq
		call b2lx0(sd3,3,2,0,0,1,gamma)	
		sd3=feq*sd2*preq
		call b2lx0(sd3,3,2,1,0,0,gamma)	
		do l=1,leqmax
			sceq1(:,l)=-sqgdroj(:,l)*cureq*rinv*preq/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,1,3,2,0,0,1,gamma)	
		do l=1,leqmax
			sceq1(:,l)=sqgdroj(:,l)*feq*preq/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,1,3,2,1,0,0,gamma)			
		do l=1,leqmax
			sceq1(:,l)=sqgdthojbst(:,l)*r*preq/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,1,3,2,0,0,1,gamma)	
		do l=1,leqmax
			sceq1(:,l)=-sqgdthoj(:,l)*feq*preq/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,-1,3,2,0,1,0,gamma)	
		do l=1,leqmax
			sceq1(:,l)=-sqgdztojbst(:,l)*r*preq/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,1,3,2,1,0,0,gamma)	
		do l=1,leqmax
			sceq1(:,l)=sqgdztoj(:,l)*cureq*preq*rinv/(feq(:)-qqinv(:)*cureq(:))
		end do		
		call b2lx(sceq1,-1,3,2,0,1,0,gamma)
		
!  parallel thermal velocity terms		

                 do l=1,leqmax
			sceq1(:,l)=preq*r*qqinv*bmod(:,l)/(feq(:)-qqinv(:)*cureq(:))
		end do
		call b2lx(sceq1,1,3,7,1,0,0,-gamma)	
                         do l=1,leqmax
			sceq1(:,l)=preq*bmod(:,l)/(feq(:)-qqinv(:)*cureq(:))
		end do
		call b2lx(sceq1,1,3,7,0,0,1,-gamma)	

		do l=1,leqmax
			sceq1(:,l)=preq*r*qqinv*sqgibmodith(:,l)/(feq-qqinv*cureq)
		end do
		call b2lx(sceq1,-1,3,7,0,0,0,-gamma)
		do l=1,leqmax
			sceq1(:,l)=preq*sqgibmodizt(:,l)/(feq-qqinv*cureq)
		end do
		call b2lx(sceq1,-1,3,7,0,0,0,-gamma)
		
!  Sheared equilibrium toroidal flow velocity for pressure equation

        call b2lx0(vzt_eq,3,3,0,0,1,-1.0_IDP) 	

!  Sheared equilibrium poloidal flow velocity for pressure equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,3,3,1,0,0,-1.0_IDP) 		

! Two fluid terms		
		
        if(twofl_on .eq. 1) then   

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
				sceq1(:,l)=sd1*sd3*grtoj(:,l)
				sceq2(:,l)=sd1*sd3*gttoj(:,l)
			end do
			call omc(sceq1,-1,3,1,1,0,0,-epsq)
			call omc(sceq2,1,3,1,0,1,0,epsq)
			sd4=rinv*sd1*sd3*cureq
			call b2lx0(sd4,3,1,1,1,0,-1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=r*sd1*sd3*bst(:,l)
			end do
			call b2lx(sceq1,-1,3,1,2,0,0,1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd1*sd2*dgrtp(:,l)
				sceq2(:,l)=sd1*sd2*dgttp(:,l)
			end do
			call b2lx(sceq1,1,3,1,1,0,0,-epsq)
			call b2lx(sceq2,-1,3,1,0,1,0,epsq)
			do l=1,leqmax
				sceq1(:,l)=rinv*sd1*sd3*cureq*djtoj(:,l)
				sceq2(:,l)=sd1*sd2*dbsjtoj(:,l)
			end do
			call b2lx(sceq1,-1,3,1,0,1,0,-1.0_IDP)
			call b2lx(sceq2,1,3,1,1,0,0,1.0_IDP)
		
	end if
		
! diffusion term added

		call b2lx_dlsq(3,3,stdifp)

! end diffusion term

!  u-zeta expression

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

!  Load 1st Omega-d terms in fast ion density equation:

        if(Trapped_on .eq. 0) then 
		do l=1,leqmax
                        write(87,'("ep-den&ep-vprl omega-d sym-break coeffs: omdr omdt omdz")')
			sceq1(:,l)=vfova2(:)*omdr(:,l)/(epsq*omcyd)
			sceq2(:,l)=vfova2(:)*omdt(:,l)/(epsq*omcyd)
			sceq3(:,l)=vfova2(:)*omdz(:,l)/(epsq*omcyd)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,3(2x,e15.7))') r(j), omdr(j,l), omdt(j,l), omdz(j,l)
			end do
		end do
	else
		do l=1,leqmax
                                                write(87,'("ep-den&ep-vprl omega-d sym-break coeffs: omdrprp omdtprp omdzprp")')
			sceq1(:,l)=vfova2(:)*omdrprp(:,l)/(epsq*omcyd)
			sceq2(:,l)=vfova2(:)*omdtprp(:,l)/(epsq*omcyd)
			sceq3(:,l)=vfova2(:)*omdzprp(:,l)/(epsq*omcyd)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,3(2x,e15.7))') r(j), omdrprp(j,l), omdtprp(j,l), omdzprp(j,l)
			end do
		end do
	end if                 

	       call b2lx(sceq1,-1,5,5,0,1,0,-1.0_IDP)
	       call b2lx(sceq2,1,5,5,1,0,0,-1.0_IDP)
	       call b2lx(sceq3,1,5,5,0,0,1,-1.0_IDP)
		   
!  Sheared equilibrium toroidal flow velocity for fast ion density equation

        call b2lx0(vzt_eq,5,5,0,0,1,-1.0_IDP) 		

!  Sheared equilibrium poloidal flow velocity for fast ion density equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,5,5,1,0,0,-1.0_IDP) 		   

!  Load 1st Omega-d terms in fast ion parallel velocity equation:

	       call b2lx(sceq1,-1,6,6,0,1,0,-1.0_IDP)
	       call b2lx(sceq2,1,6,6,1,0,0,-1.0_IDP)
	       call b2lx(sceq3,1,6,6,0,0,1,-1.0_IDP)
		   
!  Sheared equilibrium toroidal flow velocity for fast ion parallel velocit equation

        call b2lx0(vzt_eq,6,6,0,0,1,-1.0_IDP) 		

!  Sheared equilibrium poloidal flow velocity for fast ion parallel velocity equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,6,6,1,0,0,-1.0_IDP) 			   

!  Load 2nd Omega-d terms in fast ion density equation:

        if(Trapped_on .eq. 0) then 
		write(87,'("ep-den 2nd omega-d non-sym-break coeffs: omdr omdt omdz")')
		do l=1,leqmax
			sceq1(:,l)=nfeq(:)*omdr(:,l)
			sceq2(:,l)=nfeq(:)*omdt(:,l)
			sceq3(:,l)=nfeq(:)*omdz(:,l)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,3(2x,e15.7))') r(j), sceq1(j,l), sceq2(j,l), sceq3(j,l)
			end do
		end do
	else
		write(87,'("ep-den 2nd omega-d non-sym-break coeffs: omdrprp omdtprp omdzprp")')
		do l=1,leqmax
			sceq1(:,l)=nfeq(:)*omdrprp(:,l)
			sceq2(:,l)=nfeq(:)*omdtprp(:,l)
			sceq3(:,l)=nfeq(:)*omdzprp(:,l)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,3(2x,e15.7))') r(j), sceq1(j,l), sceq2(j,l), sceq3(j,l)
			end do
		end do
	end if 

		call b2lx(sceq1,-1,5,2,0,1,0,-1.0_IDP)
		call b2lx(sceq2,1,5,2,1,0,0,-1.0_IDP)
		call b2lx(sceq3,1,5,2,0,0,1,-1.0_IDP)

! diffusion term added

		call b2lx_dlsq(5,5,stdifnf)

! end diffusion term

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

		write(87,'("ep-den par-grad non-sym-break coeff")')
		do l=1,leqmax
			sceq1(:,l)=nfeq*bmod(:,l)/(feq-qqinv*cureq)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,2x,e15.7)') r(j), sceq1(j,l)
			end do
		end do
		call omc(sceq1,1,5,6,0,0,0,-1.0_IDP)

!    Omega* term

		write(87,'("ep-den omega-star non-sym-break 1d coeffs")')
		
		sd1=dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
		sd2=dnfeqdr*feq/(feq-qqinv*cureq)
			do j=0,mj
			 write(87,'(e15.7,2x,e15.7)') r(j), sd1(j), sd2(j)
			end do
		  call b2lx0(sd1,5,2,0,0,1,1.0_IDP)
		  call b2lx0(sd2,5,2,1,0,0,-1.0_IDP)
		
!!  EP FLR effects		

		if (epflr_on == 1) then
			sd1=-dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=-dnfeqdr*feq/(feq-qqinv*cureq)
			call b2lx0(sd1,5,iw,0,0,1,1.0_IDP)
			call b2lx0(sd2,5,iw,1,0,0,-1.0_IDP)
			sd1=epsq*omcyd*omegar*nfeq/vfova2
			call b2lx0(sd1,5,iw,0,0,0,1.0_IDP)
		end if	
		 
!  Load remaining terms in fast ion parallel velocity equation

!    Landau closure term

		write(87,'("ep-vprl Landau closure sym-break coeff")')
		do l=1,leqmax
			sceq1(:,l)=1.414213*LcA1*vfova*bmod(:,l)/(feq-qqinv*cureq)
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,2x,e15.7)') r(j), sceq1(j,l)
			end do
		end do
		call b2lx_landau_grad_parallel(sceq1,1,6,6,0,0,0,1.0_IDP)

!    Parallel gradient term

		write(87,'("ep-vprl prl-grad non-sym-break coeff")')
		do l=1,leqmax
			sceq1(:,l)=2*LcA0*vfova2*bmod(:,l)/(nfeq*(feq-qqinv*cureq))
			write(87,'(i4,2x,i4)') mmeq(l), nneq(l)
			do j=0,mj
			 write(87,'(e15.7,2x,e15.7)') r(j), sceq1(j,l)
			end do
		end do
		call omc(sceq1,1,6,5,0,0,0,-1.0_IDP)

!    Omega* term

        if(epflr_on .eq. 0) then 

		  write(87,'("ep-vprl omega-star non-sym-break coeff")')
			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*feq/(nfeq*(feq-qqinv*cureq))
			do j=0,mj
			 write(87,'(e15.7,2x,e15.7)') r(j), sd1(j), sd2(j)
			end do
			do j=0,mj
			 write(87,'(e15.7,3(2x,e15.7))') r(j), rinv(j), cureq(j), feq(j)
			end do
		  call b2lx0(sd1,6,1,0,0,1,1.0_IDP)
		  call b2lx0(sd2,6,1,1,0,0,-1.0_IDP)
		  
		else 	

!!   EP FLR effects	

			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*rinv*feq/(nfeq*(feq-qqinv*cureq))
			call b2lx0(sd1,6,ix1,0,0,0,1.0_IDP)
			call b2lx0(sd2,6,ix2,0,0,0,-1.0_IDP)
		 
		end if 	
		 
! diffusion term added

		call b2lx_dlsq(6,6,stdifvf)

! end diffusion term
		close(unit=87)    !close file for fast ion coefficient checks
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

!  Load terms of the thermal moment of energetic particles
		
!    Perturbed pressure gradient term

		do l=1,leqmax
			sceq1(:,l)=bet0*r*qqinv*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
		end do
		call b2lx(sceq1,1,7,3,1,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=-bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
		end do
		call b2lx(sceq1,1,7,3,0,0,1,1.0_IDP)
		
!    Perturbed magnetic field term		
		
		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)		
		do l=1,leqmax
			sceq2(:,l)=-sceq1(:,l)*sd1
		end do		
		call b2lx(sceq2,1,7,1,1,0,0,1.0_IDP)		

!  Sheared equilibrium toroidal flow velocity for the thermal moment of energetic particles equation

        call b2lx0(vzt_eq,7,7,0,0,1,-1.0_IDP) 

!  Sheared equilibrium poloidal flow velocity for the thermal moment of energetic particles equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,7,7,1,0,0,-1.0_IDP) 	

! diffusion term added
		call b2lx_dlsq(7,7,stdifv)			
		
!  End of fast ion moment equations		 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Alpha particle effects	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
        if(alpha_on .eq. 1) then  
		
!  Load 1st Omega-d terms in fast alpha density equation:

		  do l=1,leqmax
				sceq1(:,l)=valphaova2*omdr(:,l)/(epsq*omcydalp)
				sceq2(:,l)=valphaova2*omdt(:,l)/(epsq*omcydalp)
				sceq3(:,l)=valphaova2*omdz(:,l)/(epsq*omcydalp)
		  end do
		
	       call b2lx(sceq1,-1,8,8,0,1,0,-1.0_IDP)
	       call b2lx(sceq2,1,8,8,1,0,0,-1.0_IDP)
	       call b2lx(sceq3,1,8,8,0,0,1,-1.0_IDP)
		   
!  Sheared equilibrium toroidal flow velocity for fast alpha density equation

          call b2lx0(vzt_eq,8,8,0,0,1,-1.0_IDP) 

!  Sheared equilibrium poloidal flow velocity for fast alpha density equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,8,8,1,0,0,-1.0_IDP) 			   

!  Load 1st Omega-d terms in fast alpha parallel velocity equation:

	       call b2lx(sceq1,-1,9,9,0,1,0,-1.0_IDP)
	       call b2lx(sceq2,1,9,9,1,0,0,-1.0_IDP)
	       call b2lx(sceq3,1,9,9,0,0,1,-1.0_IDP)
		   
!  Sheared equilibrium toroidal flow velocity for fast alpha parallel velocity equation

          call b2lx0(vzt_eq,9,9,0,0,1,-1.0_IDP) 

!  Sheared equilibrium poloidal flow velocity for fast alpha parallel velocity equation

        vth_eq1=rinv*vth_eq/eps
        call b2lx0(vth_eq1,9,9,1,0,0,-1.0_IDP) 					   

!  Load 2nd Omega-d terms in fast alpha density equation:

		  do l=1,leqmax
		       
				sceq1(:,l)=nalpeq(:)*omdr(:,l)
				sceq2(:,l)=nalpeq(:)*omdt(:,l)
				sceq3(:,l)=nalpeq(:)*omdz(:,l)
		  end do
		
		  call b2lx(sceq1,-1,8,2,0,1,0,-1.0_IDP)
		  call b2lx(sceq2,1,8,2,1,0,0,-1.0_IDP)
		  call b2lx(sceq3,1,8,2,0,0,1,-1.0_IDP)

! diffusion term added

			call b2lx_dlsq(8,8,stdifnalp)

! end diffusion term

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

		  do l=1,leqmax
				sceq1(:,l)=nalpeq*bmod(:,l)/(feq-qqinv*cureq)
		  end do
		
		  call omc(sceq1,1,8,9,0,0,0,-1.0_IDP)

!    Omega* term

			sd1=dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=dnalpeqdr*feq/(feq-qqinv*cureq)		
		    call b2lx0(sd1,8,2,0,0,1,1.0_IDP)
		    call b2lx0(sd2,8,2,1,0,0,-1.0_IDP)

!   EP FLR effects

			if (epflr_on == 1) then
				sd1=-dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
				sd2=-dnalpeqdr*feq/(feq-qqinv*cureq)
				call b2lx0(sd1,8,iwa,0,0,1,1.0_IDP)
				call b2lx0(sd2,8,iwa,1,0,0,-1.0_IDP)
				sd1=epsq*omcydalp*omegar*nalpeq/valphaova2
				call b2lx0(sd1,8,iwa,0,0,0,1.0_IDP)
			end if
			
!  Load remaining terms in fast alpha parallel velocity equation

!    Landau closure term

		  do l=1,leqmax
				sceq1(:,l)=1.414213*LcA1alp*valphaova*bmod(:,l)/(feq-qqinv*cureq)
		  end do
		  call b2lx_landau_grad_parallel(sceq1,1,9,9,0,0,0,1.0_IDP)

!    Parallel gradient term

		  do l=1,leqmax
				sceq1(:,l)=2*LcA0alp*valphaova2*bmod(:,l)/(nalpeq*(feq-qqinv*cureq))
		  end do
		  call omc(sceq1,1,9,8,0,0,0,-1.0_IDP)

!    Omega* term

          if(epflr_on .eq. 0) then 
		  
				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*feq/(nalpeq*(feq-qqinv*cureq))
		    call b2lx0(sd1,9,1,0,0,1,1.0_IDP)
		    call b2lx0(sd2,9,1,1,0,0,-1.0_IDP)
		  
			else
		  
!!   EP FLR effects	

				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*rinv*feq/(nalpeq*(feq-qqinv*cureq))
				call b2lx0(sd1,9,ix1a,0,0,0,1.0_IDP)
				call b2lx0(sd2,9,ix2a,0,0,0,-1.0_IDP)		  
		  end if			  
		
! diffusion term added

			call b2lx_dlsq(9,9,stdifvalp)

! end diffusion term		
				
		end if
		
		allocate (ytrhs(lmx,mjm1))
		ytrhs=yt
		yt=0.0_IDP

!  put in l.h.s. of equation

!  psi equation

		sd1=1.0_IDP
		call b2lx0(sd1,1,1,0,0,0,1.0_IDP)

!  u-zeta equation

		call b2lx0(sd1,2,4,0,0,0,1.0_IDP)

! p equation

		call b2lx0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density moment equation

		call b2lx0(sd1,5,5,0,0,0,1.0_IDP)
		
! fast ion parallel velocity moment equation

		call b2lx0(sd1,6,6,0,0,0,1.0_IDP)
		
! thermal moment equation

		call b2lx0(sd1,7,7,0,0,0,1.0_IDP)	

        if(alpha_on .eq. 1) then 	

! fast alpha density moment equation

		  call b2lx0(sd1,8,8,0,0,0,1.0_IDP)
		
! fast alpha parallel velocity moment equation

		  call b2lx0(sd1,9,9,0,0,0,1.0_IDP)

		end if

		nvar=iq
		if (iflr_on == 1) nvar=nvar-1

		allocate (slhs(nvar*lmaxn))
		allocate (srhs(nvar*lmaxn))
		allocate (crhs(nvar*lmaxn))
		allocate (grwth(nvar*lmaxn))
		allocate (omega_r(nvar*lmaxn))

		slhs=0.0_IDP
		srhs=0.0_IDP
		crhs=0.0_IDP
		grwth=0.0_IDP
		omega_r=0.0_IDP
	
!		Here we calculate the instability growth rate (grwth) and frequency (omega_r)
	
		open(unit=98,file="egn_eval_1.out",status="unknown")
                write(98,'("lmx, lmaxn, noeqn =",3(2x,i5))') nvar*lmaxn, lmaxn, nvar
		write(98,'(//,"l ivar l1 lo(l1) lp lln(l1) mm nn slhs(l) srhs(l) crhs(l)")')
		do l=1,nvar*lmaxn
		        srhs(l) = 0.0_IDP; slhs(l) = 0.0_IDP; crhs(l) = 0.0_IDP
			do j=1,mjm1
				srhs(l)=srhs(l)+ytrhs(l,j)*yt(l,j)
				slhs(l)=slhs(l)+yt(l,j)**2
			end do
			ivar=l/lmaxn+1
			l1=mod(l,lmaxn)
			if (l1 == 0) then
				ivar=ivar-1
				l1=lmaxn
			end if
			if (ivar == 4) cycle
			if (lo(l1) == 0) cycle
			lp=(ivar-1)*lmaxn+lo(l1)
			do j=1,mjm1
				crhs(l)=crhs(l)+ytrhs(lp,j)*yt(l,j)
			end do
			if (ivar == 2 .or. ivar == 6 .or. ivar == 7 .or. ivar == 9) crhs(l)=-crhs(l)
			write(98,'(8(i4,2x),3(e16.8,2x))') l, ivar, l1, lo(l1), lp, lln(l1), mm(lln(l1)), nn(lln(l1)), slhs(l), srhs(l), crhs(l)
		end do

		do l=1,nvar*lmaxn
		        grwth(l) = 0.0_IDP; omega_r(l) = 0.0_IDP
			ivar=l/lmaxn+1
			l1=mod(l,lmaxn)
			if (l1 == 0) then
				ivar=ivar-1
				l1=lmaxn
			end if
			if (ivar == 4) cycle
			if (lo(l1) == 0) then
				if (slhs(l) /= 0.0_IDP) grwth(l)=srhs(l)/slhs(l)
			else if (signl(lln(l1)) > 0) then
				lp=(ivar-1)*lmaxn+lo(l1)
				if (slhs(l) == 0.0_IDP .and. slhs(lp) == 0.0_IDP) cycle
					grwth(l)=(srhs(l)+srhs(lp))/(slhs(l)+slhs(lp))
					omega_r(l)=(crhs(l)-crhs(lp))/(slhs(l)+slhs(lp))
			end if
			write(98,'(i4,3(3x,i4))') l1, lp, (mm(lln(l1))+mm(lln(lo(l1)))), (nn(lln(l1))+nn(lln(lo(l1))))
		end do
                close(unit=98)

		grwth_avg = 0.; omega_r_avg = 0.
		grwth_avg = 0.; omega_r_avg = 0.
		sum_iter = 0.
		open (unit=6,file="farprt",status="old",POSITION="APPEND")
		do l1=1,lmaxn
			if (signl(lln(l1)) < 0) cycle
			do ivar=1,nvar
				l=(ivar-1)*lmaxn+l1
				if (ivar == 1) write(6,'(" psi   : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 2) write(6,'(" phi   : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 3) write(6,'(" pr    : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 4) cycle
				if (ivar == 5) write(6,'(" nfast : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 6) write(6,'(" vfast : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 7) write(6,'(" vth   : m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)						
				if (ivar == 8) write(6,'(" nalpha: m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)
				if (ivar == 9) write(6,'(" valpha: m=",i4," n=",i4," gam=",1pe13.5," om_r=",1pe13.5)') &
						mm(lln(l1)),nn(lln(l1)),grwth(l),omega_r(l)						
				sum_iter = sum_iter + 1.
				grwth_avg = grwth_avg + grwth(l)
				omega_r_avg = omega_r_avg + omega_r(l)
			end do
		end do
		close(6)
                grwth_avg = grwth_avg/sum_iter
		omega_r_avg = omega_r_avg/sum_iter
		
		deviation_growth = 0.; deviation_freq = 0.
		do l1=1,lmaxn
			if (signl(lln(l1)) < 0) cycle
			do ivar=1,nvar
				l=(ivar-1)*lmaxn+l1
				if (ivar == 4) cycle
				deviation_growth = deviation_growth + (grwth(l) - grwth_avg)**2/sum_iter
				deviation_freq = deviation_freq + (omega_r(l) - omega_r_avg)**2/sum_iter		
			end do
		end do
		open (unit=6,file="farprt",status="old",POSITION="APPEND")
		write(6, '(/,30x,"Avg. gam:",10x,"Avg. om_r:")')
		write(6, '(27x,1pe13.5,6x,1pe13.5,/)') grwth_avg, omega_r_avg
		write(6, '(27x,1pe13.5,6x,1pe13.5,/)') deviation_growth, deviation_freq
		close(6)
		open(unit=88,file="temp_grwth_omega",status="unknown")
		write(88,'(1pe13.5,2x,1pe13.5)') grwth_avg, omega_r_avg
		write(88,'(1pe13.5,2x,1pe13.5)') deviation_growth, deviation_freq
		close(unit=88)
		deallocate (omega_r)
		deallocate (grwth)
		deallocate (crhs)
		deallocate (srhs)
		deallocate (slhs)
		
	end subroutine lincheck
