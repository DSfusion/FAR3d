!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GNU LINCENSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    FAR3d gyrofluid code ver 1.0                                              !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Copyright (C) 2018  D. Spong, L. Garcia and J. Varela                     !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  											  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is free software: you can redistribute it and/or modify      !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    it under the terms of the GNU General Public License as published by      !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    the Free Software Foundation, either version 3 of the License, or         !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    (at your option) any later version.										  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!																			  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is distributed in the hope that it will be useful,           !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of            !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the			  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    GNU General Public License for more details.							  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!																			  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!   You should have received a copy of the GNU General Public License          !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    along with this program.  If not, see <https://www.gnu.org/licenses/>     !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
	subroutine far3d_main(input_namelist)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

                logical, intent(in) :: input_namelist

!		namelist/nam/maxstp,ndump,nprint,lplots,itime,dt0,nonlin,mj,lmax,leqmax,mm,nn,mmeq,nneq,Adens,Bdens,ni,nis,ne, &
!			     delta,rc,fti,fte,stdifp,stdifu,ngeneq,eps,bet0,rs,leq,tfl_Un,tfl_Psin,tfl_pn, &
!			     ndat,xr,xdat,etascl,reta,eta0,etalmb,ietaeq,LcA0,LcA1,LcA2,LcA3,ext_prof,omegar, &
!				 iflr,r_epflr,iflr_on,s,gamma,ipert,widthi,gammai,pertscl,m0dy,nocpl,cnep,ctep,cvep,omcy,bet0_f,cnfp,cvfp,stdifnf,stdifvf,dpres, &
!			     ieldamp_on,epflr_on,twofl_on,LcA0alp,LcA1alp,LcA2alp,LcA3alp,alpha_on,stdifnfalp,stdifvfalp, &
!				 Adensalp,Bdensalp,bet0_falp,cnfpalp,cvfpalp,omcyalp,DIIID_u,r_epflralp,EP_dens_on,Alpha_dens_on,ext_prof_len,betath_factor, &
!				 EP_vel_on,Alpha_vel_on,deltaq,deltaiota,matrix_out,leqdim,ldim,jdim,nstres,etascl,reta,eta0,etalmb, &
!                                                                                      q_prof_on,Trapped_on,omcyb,rbound
		character(len=10) :: timew,datew,zone_h
		integer, dimension(8) :: values_s,values_e
		character(len=8) :: confil
		integer :: i,l,j,iend
		real(IDP) :: d,dk
		character*1 :: cdum0
		character(len=132) :: char5			

		interface
			subroutine inputlist
			end subroutine inputlist
                        subroutine inputlist_namelist
                        end subroutine inputlist_namelist
			subroutine dfault
			end subroutine dfault
			subroutine far3d_inital
			end subroutine far3d_inital
			subroutine far3d_resume
			end subroutine far3d_resume
			subroutine linstart
			end subroutine linstart
			subroutine output
			end subroutine output
			subroutine wrdump
			end subroutine wrdump
			subroutine far3d_endrun
			end subroutine far3d_endrun
			subroutine numinc
			end subroutine numinc
			subroutine energy(i)
				implicit none
				integer :: i
			end subroutine energy
			subroutine elapsed_time(values_s,values_e)
				implicit none
				integer, dimension(8) :: values_s,values_e
			end subroutine elapsed_time
		end interface

!		Input read	 
                                                   
		print *, "====================/ WELCOME TO \===================="	
		print *, "======================================================"			
		print *, "===============================  ______========_======"		
		print *, "=====_____====______====_____===|_____ |=======||====="	
		print *, "====| ____|==| ____ |==| ___ |========||=======||====="	
		print *, "====||____===||====||==||===||===_____||=======||====="	
		print *, "====| ____|==||____||==||___||==|_____ |== ____||====="			
		print *, "====||=======| ____ |==||=\\==========||==| ___ |====="	
		print *, "====||=======||====||==||==\\====_____||==||___||====="	
		print *, "====||=======||====||==||== \\==|______|==|_____|====="	
		print *, "======================================================"	
		print *, "======================\ ver1.0 /======================"		
			                                     
		write(0,'(" ====> Checking input list ... ")')		
	
                if (input_namelist) then
                   call inputlist_namelist
                else
                   call inputlist
                endif

		write(0,'(" ====> Input list check DONE !! ")')		
		
		lmax=ldim
		leqmax=leqdim	
		mj=jdim	
		
!		Date and time		
		call date_and_time(datew,timew,zone_h,values_s)

!		open (unit=6,file="farprt",status="old",POSITION="APPEND")		
		write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
											     datew(7:8),datew(5:6),datew(1:4)	

!		close(6)
												
!		Identify the number of equations in the model												 
	
		noeqn=7	
		inalp=0
		ivalp=0
		iq=0
		iw=0
		ix1=0
		ix2=0
		iwa=0
		ix1a=0
		ix2a=0
		if (alpha_on == 1) then
			inalp=noeqn+1
			noeqn=noeqn+2
			ivalp=noeqn
		end if
		if (iflr_on == 1) noeqn=noeqn+1
		iq=noeqn
		if (epflr_on == 1) then
			iw=noeqn+1
			ix1=noeqn+2
			noeqn=noeqn+3
			ix2=noeqn
			if (alpha_on == 1) then
				iwa=noeqn+1
				ix1a=noeqn+2
				noeqn=noeqn+3
				ix2a=noeqn
			end if
		end if

		if (ext_prof == 0) ieldamp_on = 0	
	
!		Initialize the perturbed variables	

		psi=0.0_IDP
		phi=0.0_IDP
		pr=0.0_IDP
		nf=0.0_IDP
		vprlf=0.0_IDP
		vthprlf=0.0_IDP
		if (alpha_on == 1) then
			nalp=0.0_IDP
			vprlalp=0.0_IDP
		end if
		ihist=ihist+1

		allocate (numhist(ihist))			  	

!		nstres indicates if the run is a new run or a continuation
!		If it is a new run the far3d_inital subroutine is called
!		If it is a continuation the far3d_resume subroutine is called
		if (nstres == 0) call far3d_inital
		if (nstres == 0) write(0,'(" ====> Preparing new run ...")')	
		if (nstres /= 0) call far3d_resume
		if (nstres /= 0) write(0,'(" ====> Preparing run continuation ...")')		
	
!		Modification of the equlibria thermal beta		
		if (betath_factor /= 1) write(0,'(" WARNING THERMAL BETA FACTOR ACTIVE = ",1pe13.6)') betath_factor 
		
!		Modification of the safety factor / iota profiles
		if (deltaq /= 0) write(0,'(" WARNING SAFETY FACTOR PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaq
		if (deltaiota /= 0) write(0,'(" WARNING IOTA PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaiota		

!		Time step normalization		
		if (itime == 1) then
			dk=0.
			do l=1,lmax
				do j=0,mj
					d=nn(l)-mm(l)*qqinv(j)
					dk=max(dk,abs(d))
				end do
			end do
			dt=dt0*2./(dk)
			dtd2=dt*.5
		else
			dt=dt0
			dtd2=dt*.5
		end if

		confil="fs"//numrun(1)//numrun(2)
		numhist(ihist)=confil

!		The model parameter are included in farprt output file		
		call output

!		The model parameter are saved in fs####z file	
		if (nstres == 0) call wrdump
		call numinc
		
!		From here the model is advanced in time	

		write(0,'(" ====> Time stepping begins ...")')

		do i = 1,10
		  nstep_count(i)=maxstp*i/10
		end do
		if (nstres /= 0) nstep_count(:)=nstep_count(:)+nstep		
		i=1
		if (maxstp > 0) then
			
			nstep1=nstep

!		Subroutine linstart creates the tridiagonal matrix where the right and left side of the
!		model equations are added  	
		
			call linstart

			iend = 0
			do while (nstep < nstep1+maxstp)

				nstep=nstep+1
				
				if (nstep_count(i) == nstep) then		
					write(0,'(" ====> Time step = ",i8)') nstep
					i=i+1
				end if					
				time=time+dt
				if (nstep == nstep1+maxstp .and. mod(nstep,nprint) /= 0) iend=1

!		Subroutine listep calculates perturvation variable time advance
!		Subroutine energy calculates the radial and poloidal component of the magentic and velocity fields
!		as well as kinetic and magnetic energy of the system.					
				
				if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(1)
				call linstep
				if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(2)

				if (mod(nstep,ndump) == 0 .and. nstep /= nstep1+maxstp) then
				
!		Subroutine wrdump writes an output file to continue the run				
					call wrdump
					call numinc
				end if

			end do				
			
!  		Stepping done.
!		Subroutine lincheck calculates the growth rate and the frequency of the instability	
			if (nonlin == 0) call lincheck
			numrun(3)="z"	

!		Subroutine far3d_endrun writes the eigenfuctions output
			call wrdump
			call far3d_endrun

		end if

		write(0,'(" ====> Simulation DONE !! ")')		
		
		call date_and_time(datew,timew,zone_h,values_e)
		write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
											     datew(7:8),datew(5:6),datew(1:4)

		call elapsed_time(values_s,values_e)

		close(6)

       end subroutine far3d_main
