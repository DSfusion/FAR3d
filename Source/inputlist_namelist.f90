!
!     This version of inputlist.f90 reads all data in Input_Model using namelist statements. This is useful for
!      setting up automated parameter scans and integrated modeling applications where it is necessary to read and
!      write the main FAR3d input file from other codes. The namelist format is a standard way of doing this and
!      allows the reads and writes to be carried out using just a few lines of code.
!     The Input_Model file is separated into two separate namelists: far3d_init and far3d_params. far3d_init
!       is read first and contains some basic information needed to allocate arrays, such as nubmer of equilibrium
!       and dynamic Fourier modes and radial grid resolution; this is needed for some of the data in far3d_params.
!       Once this is done, a call is made to set default values. Then Input_Model is rewound and far3d_init and
!       far3d_params are read a second time so that uneeded defaults are overwritten.
!
	subroutine inputlist_namelist
	
		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		real(IDP) :: widthix
                namelist /far3d_init/ nstres, numrun, numruno, numvac, nonlin, ngeneq, eq_name, maxstp, dt0, ldim, leqdim, jdim    &		
		/far3d_params/ ext_prof,ext_prof_name,ext_prof_len,iflr_on,epflr_on,ieldamp_on,twofl_on,alpha_on,                  &
		Trapped_on,matrix_out,m0dy,mm,nn,mmeq,nneq,ipert,widthix,Auto_grid_on,ni,nis,ne,delta,rc,Edge_on,edge_p,gamma,     &
		s,betath_factor,ietaeq,bet0_f,bet0_alp,omcy,omcyb,rbound,omcyalp,itime,dpres,stdifp,stdifu,stdifv,                 &
		stdifnf,stdifvf,stdifnalp,stdifvalp,LcA0,LcA1,LcA2,LcA3,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omegar,                    &
		iflr,r_epflr,r_epflralp,lplots,nprint,ndump,DIIID_u,EP_dens_on,Adens,Bdens,Alpha_dens_on,Adensalp,                 &
		Bdensalp,EP_vel_on,Alpha_vel_on,q_prof_on,Eq_vel_on,Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,deltaq,                &
		deltaiota,etascl,eta0,reta,etalmb,cnep,ctep,cnfp,cvep,cvfp,cnfpalp,cvfpalp,eqvt,eqvp

		character(len=132) :: char5		
		integer :: i
		character(len=8) :: confil	

!		Open and read the input file and open the main output farprt
		
		open (unit=5,file="Input_Model",status="old",form="formatted")
		open (unit=6,file="farprt",status="unknown")
		
		write(6,'("===================================================================================================")')			
		write(6,'("========================================================______========_============================")')		
		write(6,'("============================_____====______====_____===|_____ |=======||===========================")')	
		write(6,'("===========================| ____|==| ____ |==| ___ |========||=======||===========================")')	
		write(6,'("===========================||____===||====||==||===||===_____||=======||===========================")')	
		write(6,'("===========================| ____|==||____||==||___||==|_____ |== ____||===========================")')			
		write(6,'("===========================||=======| ____ |==||=\\==========||==| ___ |===========================")')	
		write(6,'("===========================||=======||====||==||==\\====_____||==||___||===========================")')	
		write(6,'("===========================||=======||====||==||== \\==|______|==|_____|=======VER1.0==============")')	
		write(6,'("===================================================================================================")')			

		write(6,'(" Namelist output of Input_Model:")')	
		
!              The following lines can be uncommented if a verbatim
!              copy (including comments)of Input_Model is desired in farprt.	
!		do
!			read(5,'(a)',end=20) char5
!			write(6,*) char5
!		end do
!		
!   20           rewind(5)		
		
		read(5,nml=far3d_init)
		write(6,'(/,"FAR3d_INIT Namelist:",/)')
		write(6,nml=far3d_init)										 
		
		if (nstres /= 0) then

			confil="fs"//numruno(1)//numruno(2)//numruno(3)
			open(unit=8,file=confil,status='unknown',form='unformatted')

			read(8) ihist
			rewind(8)

		end if	
		allocate (sgnleq(leqdim))		!  Vector with the sign of the equilibrium modes
		allocate (mmeq(leqdim))			!  Vector with the equilibrium polidal mode numbers
		allocate (nneq(leqdim))			!  Vector with the equilibrium toroidal mode numbers
		allocate (mheq(leqdim))			!  Vector with the helicities equilibrium polidal mode numbers
		allocate (nheq(leqdim))			!  Vector with the helicities equilibrium toroidal mode numbers	
		allocate (mm(ldim))				!  Vector with the dynamic polidal mode numbers
		allocate (nn(ldim))				!  Vector with the dynamic toroidal mode numbers
		allocate (mh(ldim))				!  Vector with the helicities dynamic polidal mode numbers
		allocate (nh(ldim))				!  Vector with the helicities dynamic toroidal mode numbers	
		allocate (signl(ldim))			!  Vector with the sign of the dynamic modes
		allocate (rs(ldim))		
		allocate (widthi(ldim))			!  Perturbation width
		allocate (gammai(ldim))			!  Perturbation parameter
		allocate (dc1m(jdim))	        !  Derivate weigths to be used in grid subroutine
		allocate (dc1p(jdim))	
		allocate (dc2m(jdim))	
		allocate (dc2p(jdim))	
		allocate (del2cm(jdim))	
		allocate (del2cp(jdim))		
		allocate (r(0:jdim))			!	Normalized minor radius
		allocate (rinv(0:jdim))			!	Inverse of the normalized minor radius
		allocate (wt1m(jdim,2))			!	Derivative weights to be used in blockj, block0, b2lx and b2lx0 subroutines
		allocate (wt10(jdim,2))		
		allocate (wt1p(jdim,2))		
		allocate (wt2m(jdim,2))		
		allocate (wt20(jdim,2))		
		allocate (wt2p(jdim,2))	
		allocate (eta(0:jdim))			!	Normalized resistivity (to the value in the axis)
		allocate (etann(0:jdim))		!	resistivity (m^3 Kg / (s C^2))
		allocate (sd1(0:jdim))			!	Dummy vectors
		allocate (sd2(0:jdim))		
		allocate (sd3(0:jdim))		
		allocate (sd4(0:jdim))		
		allocate (sd5(0:jdim))		
		allocate (sd6(0:jdim))				
		
		allocate (qq(0:jdim))			!	Safety factor
		allocate (qqinv(0:jdim))		!	Inverse of the safety factor
		allocate (qqinvp(0:jdim))		!	Radial derivate of the safety factor inverse
		allocate (denseq(0:jdim))		!	Normalized thermal plasma density
		allocate (denseqr(0:jdim))		!	Radial derivate of the normalized thermal plasma density
		allocate (psieq(0:jdim))		!	Normalized equilibrium poloidal flux
		allocate (chieq(0:jdim))		!	Normalized equilibrium toroidal flux
		allocate (preq(0:jdim))			!	Normalized equilibrium thermal plasma pressure
		allocate (feq(0:jdim))			!	Normalized equilibrium current density
		allocate (cureq(0:jdim))		!	Normalized equilibrium electric currect
		allocate (teeq(0:jdim))			!	Normalized equilibrium plasma temperature
		allocate (nfeq(0:jdim))			!	Normalized equilibrium fast particle density
		allocate (dnfeqdr(0:jdim))		!	Radial derivate of the normalized fast particle density
		allocate (vfova(0:jdim))		!	Fast particle thermal velocity normalized to the Alfven velocity
		allocate (vfova2(0:jdim))		!	vfova^2
		allocate (vtherm_ion(0:jdim))	!	User defined normalized thermal ion velocity
		allocate (vzt_eq(0:jdim))	!	User defined equilibrium thermal toroidal velocity
		allocate (vth_eq(0:jdim))	!	User defined equilibrium thermal poloidal velocity
		allocate (vth_eq1(0:jdim))	!	User defined equilibrium thermal poloidal velocity
		allocate (nalpeq(0:jdim))		!	Normalized equilibrium density of the second fast particle species	
		allocate (valphaova(0:jdim))	!	Fast particle thermal velocity normalized to the Alfven velocity of the second fast particle species
		allocate (valphaova2(0:jdim))	!	vfalphaova^2
		allocate (dnalpeqdr(0:jdim))	!	Radial derivate of the normalized density of the second fast particle species
		
		allocate (psi(0:jdim,0:ldim))			!	Normalized fluctuating poloidal flux
		allocate (phi(0:jdim,0:ldim))			!	Normalized fluctuating electrostatic potential
		allocate (pr(0:jdim,0:ldim))			!	Normalized fluctuating pressure
		allocate (uzt(0:jdim,0:ldim))			!	Normalized fluctuating vorticity (toroidal component)
		allocate (nf(0:jdim,0:ldim))			!	Normalized fluctuating fast particle density
		allocate (vprlf(0:jdim,0:ldim))			!	Normalized fluctuating fast particle parallel velocity
		allocate (vthprlf(0:jdim,0:ldim))		!	Normalized fluctuating parallel thermal velocity
	               allocate (nalp(0:jdim,0:ldim))			!	Normalized fluctuating density of the second fast particle species
		allocate (vprlalp(0:jdim,0:ldim))		!	Normalized fluctuating parallel velocity of the second fast particle species
		
		allocate (sc1(0:jdim,0:ldim))			! dummy variable for dynamic modes
		allocate (sc2(0:jdim,0:ldim))	
		allocate (sc3(0:jdim,0:ldim))	
		allocate (sc4(0:jdim,0:ldim))	
		allocate (sc5(0:jdim,0:ldim))	
		allocate (sc6(0:jdim,0:ldim))	
		allocate (sc7(0:jdim,0:ldim))	
		allocate (sc8(0:jdim,0:ldim))		
		allocate (sc9(0:jdim,0:ldim))	
		allocate (sc10(0:jdim,0:ldim))	
		allocate (sc11(0:jdim,0:ldim))											
		allocate (sceq1(0:jdim,0:leqdim))		! dummy variable for equilibrium modes
		allocate (sceq2(0:jdim,0:leqdim))	
		allocate (sceq3(0:jdim,0:leqdim))	
		allocate (sceq4(0:jdim,0:leqdim))	
		allocate (sceq5(0:jdim,0:leqdim))	
		allocate (sceq6(0:jdim,0:leqdim))	
		allocate (sceq7(0:jdim,0:leqdim))		

		allocate (epsi(ldim,2),ephi(ldim,2),epr(ldim,2),eprnc(ldim,2),ekenc(ldim,2),eke(ldim,2), &
			  emenc(ldim,2),eme(ldim,2),ealp(ldim,2),ealpnc(ldim,2))
			  
		call dfault

!		close(unit=5)	
!		open (unit=5,file="Input_Model",status="old",form="formatted")
		rewind(5)
		read(5,nml=far3d_init)
		read(5,nml=far3d_params)
											 
		write(6,'(/,"FAR3d_PARAMS Namelist:",/)')
		write(6,nml=far3d_params)										 
		write(6,'(///)')
		
		widthi(:)=widthix

		if (Auto_grid_on == 1) then

		  ni =  (jdim/2) - 1
		  nis = (jdim/4) + 1
		  ne =  jdim/4
                  
		end if	

!             Uncomment to test numrun, numruno, numvac
!		write(6,'(2a2)') (numrun(i),i=1,2)	
!		write(6,'(2a2,a1)') (numruno(i),i=1,3)
!		write(6,'(a5)') numvac
		
		close(5)		
!		close(6)
	end subroutine inputlist_namelist
