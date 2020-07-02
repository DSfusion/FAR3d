	subroutine seteq

		use param
		use cotrol
		use domain
		use equil
		use scratch
		implicit none

		character(len=1) :: t
		integer :: mjp,l2p,j,l,i
		real(IDP) :: betfac,beta,pbar,rsq
		real(IDP), dimension(0:jdim) :: ppeq,ffpeq,psieqr,preqr,feqr
		real(IDP), dimension(:), allocatable :: aa,bb,cc,dd,ee,xint,expo

		interface
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydrl(d,a,c1,c2,k,l)
				use param
				implicit none
				integer :: k,l
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydrl
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
			subroutine quadq(n,x,y,x0,b,c,d,result)
				use param
				implicit none
				integer :: n
				real(IDP) :: x0,result
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine quadq
			subroutine hbeq
			end subroutine hbeq
			subroutine vmec
			end subroutine vmec
		end interface

!		Equilibrium set up	
		
		t=char(9)

!		User defined thermal plasma density profile 		
		  if (cnep(0) > 0.0_IDP) then
			denseq=cnep(0)
			do j=1,mj
				rsq=r(j)*r(j)
				do i=1,10
					denseq(j)=denseq(j)+cnep(i)*rsq**i
				end do
			end do
			denseq=denseq/cnep(0)
		  end if
		  
!		User defined thermal electron plasma temperature profile
		  if (ctep(0) > 0.0_IDP) then
			teeq=ctep(0)
			do j=1,mj
				rsq=r(j)*r(j)
				do i=1,10
					teeq(j)=teeq(j)+ctep(i)*rsq**i
				end do
			end do
			teeq=teeq/ctep(0)
			do j=0,mj
				if (teeq(j) < 1.e-4_IDP) teeq(j)=1.e-4_IDP
			end do
		  end if

!		User defined energetic particles density profile		  
		  if (cnfp(0) > 0.0_IDP) then
			nfeq=cnfp(0)
			do j=0,mj
				rsq=r(j)*r(j)
				do i=1,10
					nfeq(j)=nfeq(j)+cnfp(i)*rsq**i
				end do				
			end do
			nfeq=nfeq/cnfp(0)
			do j=0,mj
				if (nfeq(j) < 1.e-4_IDP) nfeq(j)=1.e-4_IDP
			end do
		  end if

!		User defined energetic particles parallel velocity profile		  	  
		if (cvfp(0) > 0.0_IDP) then
			vfova=cvfp(0)
			do j=1,mj
				rsq=r(j)*r(j)
				do i=1,10
					vfova(j)=vfova(j)+cvfp(i)*rsq**i
				end do
			end do
			vfova=vfova/LcA3
			do j=0,mj
				if (vfova(j) < 1.e-4_IDP) vfova(j)=1.e-4_IDP
			end do
		end if	

!		User defined 2nd species energetic particles density profile		  
          if (alpha_on .eq. 1) then		  
		   if (cnfpalp(0) > 0.0_IDP) then
			nalpeq=cnfpalp(0)
			do j=0,mj
				rsq=r(j)*r(j)
				do i=1,10
					nalpeq(j)=nalpeq(j)+cnfpalp(i)*rsq**i
				end do
			end do
			nalpeq=nalpeq/cnfpalp(0)
			do j=0,mj
				if (nalpeq(j) < 1.e-4_IDP) nalpeq(j)=1.e-4_IDP
			end do
		   end if	

!		User defined 2nd species energetic particles parallel velocity profile		
	
		  if (cvfpalp(0) > 0.0_IDP) then
			valphaova=cvfpalp(0)
			do j=1,mj
				rsq=r(j)*r(j)
				do i=1,10
					valphaova(j)=valphaova(j)+cvfpalp(i)*rsq**i
				end do
			end do
			valphaova=valphaova/LcA3alp
			do j=0,mj
				if (valphaova(j) < 1.e-4_IDP) valphaova(j)=1.e-4_IDP
			end do
		  end if	
        end if	
		  
!!!!!!!  User defined profile for NBI density, parameter EP_dens_on, Adens and Bdens in input  !!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
		  
		  if(EP_dens_on .eq. 1) then
		    do j=0,mj
			    nfeq(j) = (.5*(1.+tanh(Adens*(Bdens-r(j)))) + 0.02)/(.5*(1.+tanh(Adens*Bdens)) + 0.02) 
				if (nfeq(j) < 1.e-4_IDP) nfeq(j)=1.e-4_IDP
			end do	
		  end if

!!!!!!!  User defined profile for 2nd EP species density, parameter Alpha_dens_on, Adensalp and Bdensalp in input  !!!!!!!			

               if (alpha_on .eq. 1) then	
		  if(Alpha_dens_on .eq. 1) then
		    do j=0,mj				
			    nalpeq(j) = (.5*(1.+tanh(Adensalp*(Bdensalp-r(j)))) + 0.02)/(.5*(1.+tanh(Adensalp*Bdensalp)) + 0.02)		  
				if (nalpeq(j) < 1.e-4_IDP) nalpeq(j)=1.e-4_IDP
			end do	
		  end if
                       end if
				
!       Ion flr components: normalized thermal ion velocity profile if no external profile

	   	    if (cvep(0) > 0.0_IDP) then
		      vtherm_ion=cvep(0)
			    do j=1,mj
				  rsq=r(j)*r(j)
				  do i=1,10
					  vtherm_ion(j)=vtherm_ion(j)+cvep(i)*rsq**i
				  end do
			    end do
			    do j=0,mj
				  if (vtherm_ion(j) < 1.e-4_IDP) vtherm_ion(j)=1.e-4_IDP
			    end do
		    end if	

!       Thermal ion toroidal velocity profile if no external profile

	   	    if (eqvt(0) > 0.0_IDP) then
		      vzt_eq=eqvt(0)
			    do j=1,mj
				  rsq=r(j)*r(j)
				  do i=1,10
					  vzt_eq(j)=vzt_eq(j)+eqvt(i)*rsq**i
				  end do
			    end do
			    do j=0,mj
				  if (vzt_eq(j) < 1.e-4_IDP) vzt_eq(j)=1.e-4_IDP
			    end do
		    end if	

!       Thermal ion poloidal velocity profile if no external profile

	   	    if (eqvp(0) > 0.0_IDP) then
		      vth_eq=eqvp(0)
			    do j=1,mj
				  rsq=r(j)*r(j)
				  do i=1,10
					  vth_eq(j)=vth_eq(j)+eqvp(i)*rsq**i
				  end do
			    end do
			    do j=0,mj
				  if (vth_eq(j) < 1.e-4_IDP) vth_eq(j)=1.e-4_IDP
			    end do
		    end if		  

!		The equilibria main parameters are calculated	
		
		if (ngeneq == 1) call vmec		  	

	end subroutine seteq