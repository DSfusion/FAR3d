	subroutine inital

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none
		
		integer :: j
		real(IDP) :: qlambw,coeff

		interface
			subroutine setmod
			end subroutine setmod
			subroutine mmlims
			end subroutine mmlims
			subroutine grid
			end subroutine grid
			subroutine seteq
			end subroutine seteq
			subroutine etachi
			end subroutine etachi
			subroutine pert
			end subroutine pert
		end interface

!		This subroutine is called if the simulation is new, not a continuation		
		
		numruns=numrun
		
		rewind(5)
		if (rc == 0.) rc=rs(1)
		if (delta >= rc .or. delta >= (1.-rc)) delta=min(rc,1.-rc)
		
!		Subroutine setmod set up the modes of the model  		
		call setmod  
		
!		n value in the complex exponential representation.	 	
		call mmlims  
		
!		Subroutine grid set up the time independent geometric scale factors	
		call grid	
		
!		Subroutine seteq set up	the equilibria	
		call seteq
		
!		Subroutine etachi set up the magnetic diffusivity of the model	
		call etachi		
		
!		Subroutine pert set up the perturbation in the equilibria			
		call pert     

	end subroutine inital

	subroutine resume

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

!		namelist/nam/maxstp,ndump,nprint,lplots,itime,dt0,nonlin,mj,lmax,leqmax,mm,nn,mmeq,nneq,Adens,Bdens,ni,nis,ne, &
!			     delta,rc,fti,fte,stdifp,stdifu,ngeneq,eps,bet0,rs,leq,LcA0,LcA1,LcA2,LcA3, &
!                 ext_prof,omegar,iflr,r_epflr,iflr_on,ieldamp_on,epflr_on,twofl_on,dpres,tfl_Un,tfl_Psin,tfl_pn, &
!			     ndat,xr,xdat,etascl,reta,eta0,etalmb,ietaeq, &
!			     s,gamma,ipert,widthi,gammai,pertscl,m0dy,nocpl,cnep,ctep,cvep,omcy,bet0_f,cnfp,cvfp,stdifnf,stdifvf,LcA0alp,LcA1alp,LcA2alp, &
!				 LcA3alp,alpha_on,stdifnfalp,stdifvfalp,Adensalp,Bdensalp,bet0_falp,cnfpalp,cvfpalp,omcyalp,DIIID_u,r_epflralp, &
!				 EP_dens_on,Alpha_dens_on,ext_prof_len,betath_factor,EP_vel_on,Alpha_vel_on,deltaq,deltaiota,matrix_out, &
!				leqdim,ldim,jdim,nstres,etascl,reta,eta0,etalmb,q_prof_on,Trapped_on,omcyb,rbound

		integer :: i

		interface
			subroutine rddump
			end subroutine rddump
			subroutine setmod
			end subroutine setmod
			subroutine mmlims
			end subroutine mmlims
			subroutine etachi
			end subroutine etachi
			subroutine pert
			end subroutine pert
		end interface

!		This subroutine is called if the simulation is a continuation run		

!		Subroutine setmod set up the modes of the model  		
		call setmod
		
!		n value in the complex exponential representation.	 			
		call mmlims

!		Subroutine rddump reads the data required to continue the run		
		call rddump

!		Change the thermal beta
		bet0=bet0*betath_factor

!		Displacment included to the safety factor / iota
		qq=(1.0_IDP/qqinv) + deltaq
		qqinv = (1.0_IDP/qq) + deltaiota	

!		Subroutine ae_profiles adds the external profiles in the model			
                if(ext_prof .eq. 1) then	
                         call ae_profiles
		end if
		
!		Subroutine etachi set up the magnetic diffusivity of the model			
		call etachi
		
!		Subroutine pert set up the perturbation in the equilibria			
		call pert	  

	end subroutine resume