!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!                                PROGRAM vtk_format (v1.0)
!!!                            FAR3d distribution graphics addon
!!!                        
!!!  Description: the module generates vtk files using the FAR3d output.
!!!  Required files from FAR3d distribution: 
!!!      1) Equilibria file
!!!      2) Model input file (Input_Model)
!!!      3) Dump file of the FAR3d run (fs####z)
!!!  Module output: vtk files with a dimension (mj_r,lth_r,lzt_r) for the variables selected
!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                  J. Varela, D. Spong, L. Garcia and Y. Ghai (01/07/2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!              MODULE LIST              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	module param
		implicit none
		integer, parameter :: IDP = kind(1.d0)
		integer :: mj_r = 100
		integer :: lth_r = 50
		integer :: lzt_r = 50
	end module param												
													
	module cotrol
		use param
		implicit none
		save
		character(len=8), dimension(:), allocatable :: numhist
		character(len=2), dimension(3) :: numrun,numruno,numruns
		character(len=5) :: numvac
		
		real(IDP) :: stdifp,stdifu,stdifnf,stdifvf,stdifv,stdifnalp,stdifvalp,dt0,Adens,Bdens,Adensalp,Bdensalp, &
			     LcA0,LcA1,LcA2,LcA3,LcA0alp,LcA1alp,LcA2alp,LcA3alp,ext_prof,omegar,iflr,r_epflr,r_epflralp, &
			     dpres,DIIID_u,betath_factor,etascl,reta,eta0,etalmb,deltaq,deltaiota
		integer :: ihist,nocpl,maxstp,nstep,ndump,nprint,lplots,itime,nstep1,nonlin,noeqn,edge_p, &
                           iflr_on,epflr_on,ieldamp_on,twofl_on,alpha_on,inalp,ivalp,iq,iw,ix1,ix2,iwa,ix1a,ix2a, &
		           EP_dens_on,Alpha_dens_on,EP_vel_on,Alpha_vel_on,Trapped_on,ext_prof_len,q_prof_on,Eq_vel_on, & 
                           Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,Auto_grid_on,Edge_on
                                                                         
		character(len=40) :: eq_name,ext_prof_name		
		logical :: matrix_out
		integer :: leqdim,ldim,jdim,nstres	
	
	end module cotrol

	module domain
		use param
		implicit none
		save
		integer :: mjeq,nfp
		real(IDP), dimension(:), allocatable :: rfar,qfar
		integer :: mj,lmax,leqmax,lbmax,mmin,mmax,nmin,nmax,mmineq,mmaxeq,nmineq,nmaxeq,mjm1,mjm2,lmax0,lmaxn,lmx,l0,leq0,lhmax, &
				   lheqmx,nnum,m0dy,mxmband,nnd,nst,mminb,mmaxb,nminb,nmaxb,mxmbandb,mmaxx,mmaxxb,mbandeq,mmaxxeq
		integer, dimension(:), allocatable :: mm,nn,mh,nh
		integer, dimension(:), allocatable :: mmeq,nneq,mheq,nheq
		integer, dimension(:), allocatable :: mmb,nnb
		integer, dimension(:), allocatable :: mmstart,mmend
		integer, dimension(:), allocatable :: mmstartb,mmendb
		real(IDP), dimension(:), allocatable :: r,rinv,r_red
		real(IDP), dimension(:), allocatable :: dc1m,dc1p,dc2m,dc2p,del2cm,del2cp
		real(IDP), dimension(:,:), allocatable :: wt1m,wt10,wt1p,wt2m,wt20,wt2p
		integer, dimension(:,:), allocatable :: ll
		integer, dimension(:,:), allocatable :: llb
		integer, dimension(:,:), allocatable :: lleq
		integer, dimension(:), allocatable :: m1n,mrang
		integer, dimension(:), allocatable :: ll0
		integer, dimension(:), allocatable :: lln,lo,llno
		integer, dimension(:), allocatable :: signl
		integer, dimension(:), allocatable :: sgnleq
		integer, dimension(:), allocatable :: lnumn
		integer, dimension(:), allocatable :: mnumn												
		integer :: ni,nis,ne
		real(IDP) :: delta,rc,fti,fte,bmodn,mu0,uion,vthi,vthe,xnuelc0,coul_log
		real(IDP), dimension(:), allocatable :: dnnbi,dne,dni,temp_epnn,ti,te,temp_ep,qprofile, &
		                                        dnnbinn,dnenn,dninn,pthermalnn,tinn,tenn,tbn,tbnnn,pepnn,ptotnn, &
                                                        pthermal,pep,ptot,vthermalep,vAlfven,vtherm_ionP,vtherm_elecP, &
							dnalpha,dnalphann,talpha,talphann,vzt_eqp,vth_eqp
												                                            												
 	end module domain

	module equil
		use param
		implicit none
		save
		real(IDP), dimension(:), allocatable :: qq,qqinv,qqinvp,denseq,denseqr,psieq,chieq,preq,feq,cureq,teeq,nfeq,dnfeqdr,vfova, &
							vfova2,vtherm_ion,vzt_eq,vth_eq,vth_eq1,nalpeq,valphaova,valphaova2,dnalpeqdr
		character(len=8), dimension(2) :: ndevice
		integer :: ngeneq,leq
		real(IDP) :: eps,bet0
		real(IDP) :: omcy,bet0_f,bet0_alp,omcyalp,omcyb,rbound,norm_eildump
		real(IDP), dimension(:), allocatable :: rs
		integer, dimension(1:10) :: nstep_count		
		real(IDP), dimension(0:10) :: cnep,ctep,cnfp,cvfp,cvep,cnfpalp,cvfpalp,eqvt,eqvp
		real(IDP), dimension(:,:), allocatable :: grr,grt,gtt,grz,gtz,gzz,sqgi,sqg,rmn,zmn,pmn, &
							  grroj,grtoj,gttoj,grzoj,gtzoj,gzzoj,grrup,grtup,grzup,gttup,gtzup,gzzup, &
							  bmod,bst,jbgrr,jbgrt,jbgtt,jbgrz,jbgtz,omdrprp,omdtprp,omdzprp, &
							  omdr,omdt,omdz,djroj,djtoj,djzoj,dbsjtoj,dbsjzoj,dbsjtbj,dgttr,dgrrt,dgrtt,dgttt,dgrrz, &
							  dgrtz,dgttz,dgrtp,dgttp,jsq,bsgrt,bsgtt,bsq,bsqgtt,lplrr,lplrt,lplrz,lplr,lpltt, &
							  lpltz,lplt,lplz,lplzz,eildr,eildt,eildz,eildrr,eildrt,eildrz,eildtt,eildtz,eildzz, &
							  sqgdroj,sqgdthoj,sqgdztoj,sqgdthojbst,sqgdztojbst,sqgibmodith,sqgibmodizt, &
                                                          test,testr,testt,testrr,testtt,testrt
		real(IDP), dimension(:,:), allocatable :: rmn_r,zmn_r,pmn_r
		LOGICAL :: lasym	
						   			
	end module equil	

	module dynamo
		use param
		implicit none
		save
		integer :: ietaeq,ipert,neta,ith,izt
		real(IDP) :: dt,dtd2,time,s,gamma,pertscl		
		real(IDP), dimension(:), allocatable :: eta,etann
		real(IDP), dimension(:), allocatable :: widthi,gammai		
		real(IDP), dimension(:,:), allocatable :: psi,phi,pr,uzt,nf,vprlf,vthprlf,nalp,vprlalp
		real(IDP), dimension(:,:,:), allocatable :: cmamm,cmamp,cmapm,cmapp
		real(IDP), dimension(:,:,:), allocatable :: amat,bmat,cmat,amatw,bmatw,cmatw,amatwalp,bmatwalp,cmatwalp
		real(IDP), dimension(:,:), allocatable :: xt,yt,xw
		integer, dimension(:,:), allocatable :: ipc,ipcw,ipcwalp	
		real(IDP), dimension(:), allocatable :: cs,sn
		real(IDP), dimension(:,:,:), allocatable :: psi_b,phi_b,pr_b,uzt_b,nf_b,vprlf_b, &
                                           	vthprlf_b,Vr_b,Vth_b,Br_b,Bth_b,Bzt_b,sqgi_b, &
                                                nalp_b,vprlalp_b,rmn_c,zmn_c,pmn_c,xx1,xx2,xx3		
	end module dynamo

	module scratch
		use param
		implicit none
		save
		real(IDP), dimension(:,:), allocatable :: sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10,sc11
 		real(IDP), dimension(:,:), allocatable :: sceq1,sceq2,sceq3,sceq4,sceq5,sceq6,sceq7
		real(IDP), dimension(:), allocatable :: sd1,sd2,sd3,sd4,sd5,sd6
		real(IDP), dimension(:), allocatable :: sd1_r,sd2_r
		
	end module scratch

	module findrf
		use param
		implicit none
		save
		integer :: n0,nt
		real(IDP) :: d,dx
				
	end module findrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	program vtk_format

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer, dimension(8) :: values_s,values_e
		character(len=8) :: confil
		integer :: i,l,j,iend,n,m
		real(IDP) :: d,dk
		character*1 :: cdum0
		character(len=132) :: char5	

		integer :: lth,lzt
		real(IDP) :: scmn,scmx,sclmn,sclmx,scnorm
		character(len=1) :: t
		real(IDP) :: pi = 3.1415926

		interface
			subroutine inputlist
			end subroutine inputlist		
			subroutine dfault
			end subroutine dfault
			subroutine setmod
			end subroutine setmod
			subroutine mmlims
			end subroutine mmlims
			subroutine vmec
			end subroutine vmec
			subroutine sum_var
			end subroutine sum_var
			subroutine vtk_file
			end subroutine vtk_file
		end interface


!		Input read	 
                                                   
		write(0,'("====================/ WELCOME TO \====================")')	
		write(0,'("======================================================")')			
		write(0,'("===============================  ______========_======")')		
		write(0,'("=====_____====______====_____===|_____ |=======||=====")')	
		write(0,'("====| ____|==| ____ |==| ___ |========||=======||=====")')	
		write(0,'("====||____===||====||==||===||===_____||=======||=====")')	
		write(0,'("====| ____|==||____||==||___||==|_____ |== ____||=====")')			
		write(0,'("====||=======| ____ |==||=\\==========||==| ___ |=====")')	
		write(0,'("====||=======||====||==||==\\====_____||==||___||=====")')	
		write(0,'("====||=======||====||==||== \\==|______|==|_____|=====")')
		write(0,'("======================================================")')	
		write(0,'("===============GRAPHICS ADDON MODULE==================")')
		write(0,'("======================================================")')	
		write(0,'("======================\ ver1.0 /======================")')		
			                                     
		write(0,'(" ====> Checking input list ... ")')		
	
		call inputlist	

		write(0,'(" ====> Input list check DONE !! ")')		
		
		lmax=ldim
		leqmax=leqdim	
		mj=jdim	

!		Subroutine setmod set up the modes of the model  		
		call setmod
		
!		n value in the complex exponential representation.	 			
		call mmlims

!		Read dump file
		write(0,'(" ====> Reading dump file ... ")')
                call rddump
		write(0,'(" ====> Reading dump file DONE ... ")')

!		Read equilibria file
		write(0,'(" ====> Reading equilibria file ... ")')
                call vmec
		write(0,'(" ====> Reading equilibria file DONE ... ")')	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		The variables are calculated for a given set of radial and angular points
!               in Boozer coordinates (rho,theta,zeta) and transformed to Cartesian coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(0,'(" ====> Adding Fourier terms for each variable ... ")')
                call sum_var
		write(0,'(" ====> Adding Fourier terms for each variable DONE ... ")')	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		The vtk files are created
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

		write(0,'(" ====> Creating vtk files ... ")')
                call vtk_file
		write(0,'(" ====> Creating vtk files DONE ")')	

	end program vtk_format

	subroutine inputlist
	
		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none
		 
		character*1 :: cdum0
		character(len=132) :: char5		
		integer :: i
		character(len=8) :: confil	
		real(IDP) :: widthix	

!		Read the input files
		
		open (unit=5,file="Input_Model",status="old",form="formatted")
		
		do
			read(5,'(a)',end=20) char5
		end do
20    rewind(5)		
		
		read(5,'(a1)') cdum0										 
		read(5,'(a1)') cdum0
		read(5,*) cdum0
		read(5,'(a1)') cdum0
		read(5,'(2a2)') (numrun(i),i=1,2)
		read(5,'(a1)') cdum0
		read(5,'(2a2,a1)') (numruno(i),i=1,3)
		read(5,'(a1)') cdum0
		read(5,'(a5)') numvac	
                read(5,'(a1)') cdum0
                read(5,*) nonlin	
                read(5,'(a1)') cdum0
                read(5,*) ngeneq
                read(5,'(a1)') cdum0
                read(5,'(a40)') eq_name			
                read(5,'(a1)') cdum0					
                read(5,*) maxstp	
                read(5,'(a1)') cdum0
                read(5,*) dt0	
                read(5,'(a1)') cdum0		
                read(5,*) ldim		
                read(5,'(a1)') cdum0
                read(5,*) leqdim	
                read(5,'(a1)') cdum0	
                read(5,*) jdim		

                nstres=1
		confil="fs"//numruno(1)//numruno(2)//numruno(3)
		open(unit=8,file=confil,status='unknown',form='unformatted')

		read(8) ihist
		rewind(8)	
		
        close (5)		
		
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

		allocate (psi_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (phi_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (pr_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (uzt_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (nf_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (vprlf_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (vthprlf_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (Vr_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (Vth_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (Br_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (Bth_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (Bzt_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (sqgi_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (nalp_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (vprlalp_b(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (rmn_c(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (zmn_c(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (pmn_c(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (xx1(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (xx2(0:mj_r,0:lth_r,0:lzt_r))	
		allocate (xx3(0:mj_r,0:lth_r,0:lzt_r))	

		allocate (r_red(0:mj_r))
                allocate (sd1_r(0:mj_r))	
                allocate (sd2_r(0:mj_r))

		call dfault	

		open (unit=5,file="Input_Model",status="old",form="formatted")		

		read(5,'(a1)') cdum0										 
		read(5,'(a1)') cdum0
		read(5,*) nstres
		read(5,'(a1)') cdum0
		read(5,'(2a2)') (numrun(i),i=1,2)
		read(5,'(a1)') cdum0
		read(5,'(2a2,a1)') (numruno(i),i=1,3)
		read(5,'(a1)') cdum0
		read(5,'(a5)') numvac	
		read(5,'(a1)') cdum0
		read(5,*) nonlin	
		read(5,'(a1)') cdum0
		read(5,*) ngeneq	
		read(5,'(a1)') cdum0
		read(5,'(a40)') eq_name			
		read(5,'(a1)') cdum0	
		read(5,*) maxstp	
		read(5,'(a1)') cdum0
		read(5,*) dt0	
		read(5,'(a1)') cdum0		
		read(5,*) ldim		
		read(5,'(a1)') cdum0
		read(5,*) leqdim	
		read(5,'(a1)') cdum0	
		read(5,*) jdim	
		read(5,'(a1)') cdum0	
		read(5,*) ext_prof	
		read(5,'(a1)') cdum0	
		read(5,'(a40)') ext_prof_name	
		read(5,'(a1)') cdum0	
		read(5,*) ext_prof_len		
		read(5,'(a1)') cdum0	
		read(5,*) iflr_on	
		read(5,'(a1)') cdum0		
		read(5,*) epflr_on	
		read(5,'(a1)') cdum0	
		read(5,*) ieldamp_on	
		read(5,'(a1)') cdum0	
		read(5,*) twofl_on	
		read(5,'(a1)') cdum0
		read(5,*) alpha_on	
		read(5,'(a1)') cdum0	
		read(5,*) Trapped_on	
		read(5,'(a1)') cdum0	
		read(5,*) matrix_out
		read(5,'(a1)') cdum0
		read(5,*) m0dy
		read(5,'(a1)') cdum0	
					
		read(5,'(a1)') cdum0										 
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0		
		read(5,*) (mm(i),i=1,ldim)
		read(5,'(a1)') cdum0		
		read(5,*) (nn(i),i=1,ldim)				
		read(5,'(a1)') cdum0		
		read(5,*) (mmeq(i),i=1,leqdim)				
		read(5,'(a1)') cdum0		
		read(5,*) (nneq(i),i=1,leqdim)	  
		read(5,'(a1)') cdum0	
		read(5,'(a1)') cdum0		
		read(5,*) ipert	
		read(5,'(a1)') cdum0		
		read(5,*) widthix
		read(5,'(a1)') cdum0		
		read(5,*) Auto_grid_on		
		read(5,'(a1)') cdum0		
		read(5,*) ni
		read(5,'(a1)') cdum0		
		read(5,*) nis
		read(5,'(a1)') cdum0		
		read(5,*) ne
		read(5,'(a1)') cdum0		
		read(5,*) delta
		read(5,'(a1)') cdum0		
		read(5,*) rc	
		read(5,'(a1)') cdum0
		read(5,*) Edge_on	
		read(5,'(a1)') cdum0
		read(5,*) edge_p	
		read(5,'(a1)') cdum0			
		read(5,'(a1)') cdum0		
		read(5,*) gamma			
		read(5,'(a1)') cdum0		
		read(5,*) s	
		read(5,'(a1)') cdum0		
		read(5,*) betath_factor
		read(5,'(a1)') cdum0		
		read(5,*) ietaeq		
		read(5,'(a1)') cdum0		
		read(5,*) bet0_f	
		read(5,'(a1)') cdum0		
		read(5,*) bet0_alp		
		read(5,'(a1)') cdum0		
		read(5,*) omcy			
		read(5,'(a1)') cdum0		
		read(5,*) omcyb			
		read(5,'(a1)') cdum0
		read(5,*) rbound			
		read(5,'(a1)') cdum0		
		read(5,*) omcyalp		
		read(5,'(a1)') cdum0		
		read(5,*) itime	
		read(5,'(a1)') cdum0		
		read(5,*) dpres		
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0		
		read(5,*) stdifp			
		read(5,'(a1)') cdum0		
		read(5,*) stdifu	
		read(5,'(a1)') cdum0		
		read(5,*) stdifv	
		read(5,'(a1)') cdum0	
		read(5,*) stdifnf		
		read(5,'(a1)') cdum0		
		read(5,*) stdifvf	
		read(5,'(a1)') cdum0		
		read(5,*) stdifnalp	
		read(5,'(a1)') cdum0		
		read(5,*) stdifvalp
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0		
		read(5,*) LcA0			
		read(5,'(a1)') cdum0		
		read(5,*) LcA1	
		read(5,'(a1)') cdum0		
		read(5,*) LcA2		
		read(5,'(a1)') cdum0		
		read(5,*) LcA3	
		read(5,'(a1)') cdum0		
		read(5,*) LcA0alp	
		read(5,'(a1)') cdum0		
		read(5,*) LcA1alp
		read(5,'(a1)') cdum0		
		read(5,*) LcA2alp	
		Read(5,'(a1)') cdum0		
		read(5,*) LcA3alp		
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0		
		read(5,*) omegar			
		read(5,'(a1)') cdum0		
		read(5,*) iflr	
		read(5,'(a1)') cdum0		
		read(5,*) r_epflr		
		read(5,'(a1)') cdum0		
		read(5,*) r_epflralp	
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0
		read(5,*) lplots			
		read(5,'(a1)') cdum0		
		read(5,*) nprint	
		read(5,'(a1)') cdum0		
		read(5,*) ndump	
		read(5,'(a1)') cdum0		
		read(5,'(a1)') cdum0
		read(5,*) DIIID_u	
		
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0
		read(5,'(a1)') cdum0		
		read(5,*) EP_dens_on			
		read(5,'(a1)') cdum0		
		read(5,*) Adens			
		read(5,'(a1)') cdum0		
		read(5,*) Bdens	
		read(5,'(a1)') cdum0		
		read(5,*) Alpha_dens_on		
		read(5,'(a1)') cdum0		
		read(5,*) Adensalp		
		read(5,'(a1)') cdum0		
		read(5,*) Bdensalp		
		read(5,'(a1)') cdum0		
		read(5,*) EP_vel_on	
		read(5,'(a1)') cdum0	
		read(5,*) Alpha_vel_on	
		read(5,'(a1)') cdum0		
		read(5,*) q_prof_on
		read(5,'(a1)') cdum0		
		read(5,*) Eq_vel_on
		read(5,'(a1)') cdum0	
		read(5,*) Eq_velp_on				
		read(5,'(a1)') cdum0	
		read(5,*) Eq_Presseq_on			
		read(5,'(a1)') cdum0
		read(5,*) Eq_Presstot_on			
		read(5,'(a1)') cdum0
		read(5,*) deltaq			
		read(5,'(a1)') cdum0	
		read(5,*) deltaiota	
		read(5,'(a1)') cdum0	
		read(5,*) etascl
		read(5,'(a1)') cdum0	
		read(5,*) eta0
		read(5,'(a1)') cdum0	
		read(5,*) reta
		read(5,'(a1)') cdum0	
		read(5,*) etalmb		
		read(5,'(a1)') cdum0			
		read(5,*) (cnep(i),i=0,10)		
		read(5,'(a1)') cdum0		
		read(5,*) (ctep(i),i=0,10)	
		read(5,'(a1)') cdum0	
		read(5,*) (cnfp(i),i=0,10)			
		read(5,'(a1)') cdum0		
		read(5,*) (cvep(i),i=0,10)
		read(5,'(a1)') cdum0		
		read(5,*) (cvfp(i),i=0,10)	
		read(5,'(a1)') cdum0		
		read(5,*) (cnfpalp(i),i=0,10)	
		read(5,'(a1)') cdum0		
		read(5,*) (cvfpalp(i),i=0,10)	
		read(5,'(a1)') cdum0		
		read(5,*) (eqvt(i),i=0,10)
		read(5,'(a1)') cdum0		
		read(5,*) (eqvp(i),i=0,10)	

		widthi(:)=widthix

		if (Auto_grid_on == 1) then

		  ni =  (jdim/2) - 1
		  nis = (jdim/4) + 1
		  ne =  jdim/4
                  
		end if	
		
		close (5)
		
	end subroutine inputlist

	subroutine dfault

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

!		Default values of the model variables and parameters						
		
		numrun(1)="zz"
		numrun(2)="zz"
		numrun(3)="0"
		numruno(1)="fs"
		numruno(2)="##"
		numruno(3)="#"
		numruns(1)="zz"
		numruns(2)="zz"
		numruns(3)="0"
		numvac="vvvvv"
		ihist=0
		maxstp=100000
		nstep=0
		ndump=100000
		nprint=100
		lplots=0
		itime=1
		dt0=.66			
		nstep1=0
		nonlin=1	
		mj=0
		lmax=0
		leqmax=0
		mmin=0
		mmax=0
		nmin=0
		nmax=0
		mmineq=0
		mmaxeq=0 				
		nmineq=0
		nmaxeq=0
		mjm1=mj-1		
		mjm2=mj-2
		r=0.
		mm=0 	
		nn=0
		mmeq=0
		nneq=0
		l0=0
		leq0=0
		rinv=0. 	
		dc1m=0.
		dc1p=0.
		dc2m=0.
		dc2p=0.
		del2cm=0.
		del2cp=0.
		signl=0
		sgnleq=0
		lhmax=0		
		mh=0
		nh=0
		lheqmx=0
		mheq=0
		nheq=0
		m0dy=0
		nocpl=0
		ni=mj/3
		ne=mj/3
		nis=mj-ni-ne
		delta=.333
		Adens=7.0
		Bdens=0.5	
		LcA0=0.5
		LcA1=-0.88623
		LcA2=1
		LcA3=1	
		ext_prof=0	
		iflr_on=0
		ieldamp_on=0		
		epflr_on=0
		twofl_on=0	
		dpres=0.5
		DIIID_u=0
		omegar=0	
		iflr=0	
		r_epflr=0	
		r_epflralp=0		
		rc=.5
		fti=.95
		fte=.95  				
		qq=0.
		qqinv=0.
		psieq=0.
		chieq=0.
		preq=0.
		feq=0.
		cureq=0.
		denseq=1.
		teeq=1.
                vzt_eq=0.
                vtherm_ion=0.
		ngeneq=0.
		ndevice(1)="        "
		ndevice(2)="        "
		eps=0.		
		bet0=1.
		rs=0.
		leq=1
		nfp=1
		eta=0.
		etann=0.
		dt=0.
		time=0.
		etascl=1.
		reta=0.
		eta0=0.
		etalmb=0.
		ietaeq=0
		s=1.e+5
		gamma=0.
		ipert=0
		widthi=0.
		gammai=0.
		pertscl=1.
		stdifp=0.
		stdifu=0.
		stdifv=0.
		stdifnf=0.
		stdifvf=0.
		cnep=0.
		ctep=0.
		cvep=0.
                LcA0alp=0.
		LcA1alp=0.
		LcA2alp=0.
		LcA3alp=0.
		alpha_on=0.
		stdifnalp=0.
		stdifvalp=0.
		Adensalp=0.
		Bdensalp=0.
		bet0_alp=0.		
                cnfpalp=0.
		cvfpalp=0.		
		omcyalp=0.
		betath_factor=1.
		EP_dens_on=0.
		Alpha_dens_on=0.
		EP_vel_on=0.
		Alpha_vel_on=0.
		ext_prof_len=100
                deltaq=0.
		deltaiota=0.	
                Eq_vel_on=0.
                q_prof_on=0.
                omcyb=0.
	        rbound=0.
                Trapped_on=0.

	end subroutine dfault

	subroutine setmod

		use param
		use cotrol
		use domain
		use dynamo
		implicit none

		interface
			subroutine inputlist
			end subroutine inputlist	
		end interface	
		
		integer :: l,n,m,lh,mnxxx,mt,nt,ngcd,ndiv,lll,nxeq,i,lm,lq,mcnt,mcnt1,n1,nl,lp,lpp,lppn,lop,l1,mband,neq,np,nm,icou,ncou,nlin
		integer, dimension(:), allocatable :: iflag,nfm,nfmo,mcnh		

!		Set up the modes distribution and the mode couplings of the model  
		
		mmin=0
		mmax=0
		nmin=0
		nmax=0
		mmineq=0
		mmaxeq=0
		nmineq=0
		nmaxeq=0
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			signl(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) signl(l)=-1
			if (n == 0 .and. m == 0) signl(l)=0
		end do
		mmax=maxval(mm(1:lmax))
		mmin=minval(mm(1:lmax))
		nmax=maxval(nn(1:lmax))
		nmin=minval(nn(1:lmax))
		mmaxx=max(mmax,abs(mmin))
		do l=1,leqmax
			m=mmeq(l)
			n=nneq(l)
			sgnleq(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) sgnleq(l)=-1
			if (n == 0 .and. m == 0) sgnleq(l)=0
		end do
		mmaxeq=maxval(mmeq(1:leqmax))
		mmineq=minval(mmeq(1:leqmax))
		nmaxeq=maxval(nneq(1:leqmax))
		nmineq=minval(nneq(1:leqmax))
		mmaxxeq=max(mmaxeq,abs(mmineq))

		allocate (ll(mmin:mmax,nmin:nmax))

		ll=0
		do l=1,lmax
			ll(mm(l),nn(l))=l
		end do
		l0=0
		if (0 >= mmin .and. 0 <= mmax .and. 0 >= nmin .and. 0 <= nmax) l0=ll(0,0)
		if (l0 == 0) then
			write (6,'("  setmod: l0=0")')
			stop
		end if

		allocate (lleq(mmineq:mmaxeq,nmineq:nmaxeq))

		lleq=0
		do l=1,leqmax
			lleq(mmeq(l),nneq(l))=l
		end do
		leq0=0
		if (0 >= mmineq .and. 0 <= mmaxeq .and. 0 >= nmineq .and. 0 <= nmaxeq) leq0=lleq(0,0)
		if (leq0 == 0) then
			write (6,'("  setmod: leq0=0")')
			stop
		end if

!  find prime harmonics.

		lh=0
		mnxxx=1
		do l=1,lmax

			mt=mm(l)*signl(l)
			nt=nn(l)*signl(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) then
					mnxxx=max(mnxxx,nt)
					nt=1
				end if
				if (mt /= 0) then
					mnxxx=max(mnxxx,mt)
					mt=1
				end if
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mnxxx=max(mnxxx,ngcd)
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mh(lll) .and. nt == nh(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mh(lh)=mt
				nh(lh)=nt
			end if
		end do
		lhmax=lh

		lh=0
		do l=1,leqmax

			mt=mmeq(l)*sgnleq(l)
			nt=nneq(l)*sgnleq(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) nt=1
				if (mt /= 0) mt=1
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mheq(lll) .and. nt == nheq(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mheq(lh)=mt
				nheq(lh)=nt
			end if
		end do
		lheqmx=lh

		nxeq=nmaxeq
		if (nmaxeq > 0) then
			do l=1,leqmax
				if (nneq(l) > 0) nxeq=min(nxeq,nneq(l))
			end do
			do l=1,leqmax
				if (nneq(l) /= nxeq*(nneq(l)/nxeq)) then
					write (6,'("  setmod: nneq(",i2,") =",i4," is not multiple of nxeq =",i3)') l,nneq(l),nxeq
					stop
				end if
			end do
		else
			nxeq=0
		end if

		lmaxn=lmax-leqmax
		if (m0dy < 0) m0dy=0
		lmaxn=lmaxn+m0dy
		lmx=noeqn*lmaxn

		allocate (m1n(nmin:nmax))
		allocate (mrang(nmin:nmax))

		do n=nmin,nmax
			m1n(n)=mmax+1
			mrang(n)=mmin-1
			do l=1,lmaxn
				if (nn(l) /= n) cycle
				m1n(n)=min(m1n(n),mm(l))
				mrang(n)=max(mrang(n),mm(l))
			end do
			if (mrang(n) == mmin-1) then
				mrang(n)=0
			else
				mrang(n)=mrang(n)-m1n(n)+1
			endif
		end do
		if (nmin < 0 .and. nmax > 0) then
			if (nmin /= -nmax) then
				write (6,'("  setmod: nmin =",i4," nmax =",i3)') nmin,nmax
				stop
			end if			
			do n=1,nmax
				if (mrang(n)+mrang(-n) > 0 .and. (mrang(-n) /= mrang(n) .or. m1n(n) /= -(m1n(-n)+mrang(-n)-1))) then
					write (6,'("  setmod: n =",i4," mrang =",i3,i4," m1n =",i3,i5)') n,mrang(n),mrang(-n),m1n(n),m1n(-n)+mrang(-n)-1
					stop
				end if
			end do
			if (mrang(0) /= 0 .and. (mod(mrang(0),2) == 0 .or. m1n(0)+(mrang(0)-1)/2 /= 0)) then
				write (6,'("  setmod: n = 0  mrang =",i3," m1n =",i3,i5)') mrang(0),m1n(0),m1n(0)+mrang(0)-1
				stop
			end if
		end if

		lmax0=0
		do i=1,mrang(0)
			m=m1n(0)+i-1
			l=ll(m,0)
			if (l == 0) cycle
			lmax0=lmax0+1
		end do

		if (lmax0 > 0) then

			allocate (ll0(lmax0))

			lm=0
			do i=1,mrang(0)
				m=m1n(0)+i-1
				l=ll(m,0)
				if (l == 0) cycle
				lm=lm+1
				ll0(lm)=l
			end do

		end if

		allocate (lln(lmaxn))
		allocate (iflag(nmin:nmax))

		iflag=0
		nnum=0
		n1=0

		if (nocpl == 0) then

			if (nxeq > 0) then

!	helical couplings.

				allocate (nfm(nmax+1))
				allocate (nfmo(0:nmax))
				allocate (mcnh(nmax+1))
				nfmo=0

				do n=1,nmax

					if (mrang(n) == 0 .or. iflag(n) == 1) cycle
					ncou=0
					do i=0,n-1
						ncou=ncou+iflag(i)
					end do
					if (ncou > 0) then
						icou=0
						do neq=nxeq,nmaxeq,nxeq
							nm=abs(neq-n)
							icou=icou+iflag(nm)
							if (iflag(nm) == 1) n1=nfmo(nm) 
						end do
					end if
					if (ncou == 0 .or. icou == 0) then
						nnum=nnum+1
						nfm(nnum)=n
						nfmo(n)=nnum
						mcnh(nnum)=0
						n1=nnum
					end if
					mcnt1=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnh(n1)=mcnh(n1)+1
						mcnt1=mcnt1+1
					end do
					if (mcnh(nnum) == 0) then
						nnum=nnum-1
					else
						iflag(n)=1
					end if

					if (mcnt1 == 0) cycle
					nfmo(n)=n1

					do neq=nxeq,nmaxeq,nxeq

						nm=abs(neq-n)
						if (mrang(nm) > 0 .and. iflag(nm) == 0) then
							do i=1,mrang(nm)
								m=m1n(nm)+i-1
								l=ll(m,nm)
								lq=0
								if (nm >= nmineq .and. nm <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,nm)
								if (l == 0 .or. lq > m0dy) cycle
								mcnh(n1)=mcnh(n1)+1
							end do
							iflag(nm)=1
							nfmo(nm)=n1
						end if

						np=abs(neq+n)
						if (np > nmax) cycle
						if (mrang(np) > 0 .and. iflag(np) == 0) then
							do i=1,mrang(np)
								m=m1n(np)+i-1
								l=ll(m,np)
								lq=0
								if (np >= nmineq .and. np <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,np)
								if (l == 0 .or. lq > m0dy) cycle
								mcnh(n1)=mcnh(n1)+1
							end do
							iflag(np)=1
							nfmo(np)=n1
						end if

					end do

				end do

				allocate (mnumn(nnum))
				allocate (lnumn(0:nnum))

				nnd=0
				lnumn(0)=0
				lp=0
				do nl=1,nnum

					lpp=0
					do n=1,nmax
						if (nfmo(n) /= nl) cycle
						do i=1,mrang(n)
							m=m1n(n)+i-1
							l=ll(m,n)
							lq=0
							if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
							if (l == 0 .or. lq > m0dy) cycle
							lp=lp+1
							lln(lp)=l
							lpp=lpp+1
						end do
					end do
					if (nfmo(0) == nl) then
						if (nmin < 0) then
							do i=(mrang(0)+3)/2,mrang(0)
								m=m1n(0)+i-1
								l=ll(m,0)
								lq=0
								if (0 >= nmineq .and. 0 <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,0)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
							l=ll(0,0)
							lq=0
							if (0 >= nmineq .and. 0 <= nmaxeq .and. 0 >= mmineq .and. 0 <= mmaxeq) lq=lleq(0,0)
							if (l /= 0 .and. lq <= m0dy) then
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end if
							do i=(mrang(0)-1)/2,1,-1
								m=m1n(0)+i-1
								l=ll(m,0)
								lq=0
								if (0 >= nmineq .and. 0 <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,0)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
						else
							do i=2,mrang(0)
								m=m1n(n)+i-1
								l=ll(m,n)
								lq=0
								if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
							l=ll(0,0)
							lq=0
							if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
							if (l /= 0 .and. lq <= m0dy) then
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end if
						end if
					end if
					if (nmin < 0) then
						do n=-1,nmin,-1
							if (nfmo(-n) /= nl) cycle
							do i=mrang(n),1,-1
								m=m1n(n)+i-1
								l=ll(m,n)
								lq=0
								if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
						end do
					end if
					mnumn(nl)=lpp
					lnumn(nl)=lp

				end do

			else

!	toroidal couplings

				do n=nmin,nmax
					if (mrang(n) == 0) cycle
					if (n < 0 .and. nmax > 0) cycle
					nnum=nnum+1
					mcnt=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnt=mcnt+1
					end do
					if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
						do i=1,mrang(-n)
							m=m1n(-n)+i-1
							l=ll(m,-n)
							lq=0
							if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
							if (l == 0 .or. lq > m0dy) cycle
							mcnt=mcnt+1
						end do
					end if
					if (mcnt == 0) nnum=nnum-1
				end do

				allocate (mnumn(nnum))
				allocate (lnumn(0:nnum))

				nnd=0
				lnumn(0)=0
				nl=0
				lp=0
				do n=nmin,nmax
					if (mrang(n) == 0) cycle
					if (n < 0 .and. nmax > 0) cycle
					nl=nl+1
					lpp=0
					lppn=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						lp=lp+1
						lln(lp)=l
						lpp=lpp+1
					end do
					if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
						do i=mrang(-n),1,-1
							m=m1n(-n)+i-1
							l=ll(m,-n)
							lq=0
							if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
							if (l == 0 .or. lq > m0dy) cycle
							lp=lp+1
							lln(lp)=l
							lppn=lppn+1
						end do
					end if
					if (lpp > 0) then
						mnumn(nl)=lpp
						lnumn(nl)=lp
						if (nmin == 0 .or. nmax == 0 .or. n == 0) then
							nnd=nnd+1
						else if (lppn /= lpp) then
							write (6,'("  setmod: n =",i4," mnumn =",i3,i4)') n,lpp,lppn
							stop
						end if
					else
						nl=nl-1
					end if
				end do
				nst=nnd+1
			
			end if

		else

!	no couplings, cylinder.

			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if(l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nnum=nnum+1
					mcnt=1
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) mcnt=mcnt+1
				end do
			end do

			allocate (mnumn(nnum))
			allocate (lnumn(0:nnum))

			nnd=0
			lnumn(0)=0
			nl=0
			lp=0
			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if (l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nl=nl+1
					lp=lp+1
					lpp=1
					lln(lp)=l
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) then
						lp=lp+1
						lpp=lpp+1
						lln(lp)=lop
					end if
					mnumn(nl)=lpp
					lnumn(nl)=lp
					if (lpp == 1) nnd=nnd+1
				end do
			end do
			nst=nnd+1
			
			write(6,'(/" number of n-values with one parity:   ",i3)') nnd
			write(6,'(" number of n-values with both parities:",i3/)') nnum-nnd

		end if

		allocate (lo(lmaxn))

		lo=0
		do l=1,lmaxn
			m=-mm(lln(l))
			n=-nn(lln(l))
			if (m < mmin .or. m > mmax .or. n < nmin .or. n > nmax) cycle
			do l1=1,lmaxn
				if (mm(lln(l1)) == m .and. nn(lln(l1)) == n) exit
			end do
			if (l1 > lmaxn) then
				lo(l)=0
			else
				lo(l)=l1
			end if
		end do

		mband=0
		do m=0,mmax
			if (ll(m,0) /= 0) mband=mband+1
		end do
		mxmband=2*mband-1
		do n=1,nmax
			mband=0
			do m=mmin,mmax
				if (ll(m,n) /= 0) mband=mband+1
			end do
			if (mband > mxmband) mxmband=mband
		end do
	end subroutine setmod

	subroutine mmlims

	!   Calculate the minimum and maximum m values for each n value
	!   in the complex exponential representation.

		use param
		use domain
		implicit none

		integer :: n,l

		allocate (mmstart(-nmax:nmax))
		allocate (mmend(-nmax:nmax))

		do n = 0,nmax
			mmstart(n) = 1000
			mmend(n) = -1000
		end do

		do l = 1,lmax

			n = nn(l)
			if (n >= 0) then
				if (mm(l) < mmstart(n)) mmstart(n) = mm(l)
				if (mm(l) > mmend(n)) mmend(n) = mm(l)
			else
				if (-mm(l) < mmstart(-n)) mmstart(-n) = -mm(l)
				if (-mm(l) > mmend(-n)) mmend(-n) = -mm(l)
			end if

		end do

		!   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
		!   so limits for negative n rows are just negatives of corresponding
		!   positive n limits.

		do n = 0,nmax

			mmstart(-n) = -mmend(n)
			mmend(-n) = -mmstart(n)

		end do

	end subroutine mmlims

	subroutine rddump

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

		integer :: i,l,j,lmaxo,m,n,ihisto,l1
		integer, dimension(ldim) :: mmo,nno
		integer, dimension(ldim) :: lmap
		real(IDP) :: pertsclo
		real(IDP), dimension(ldim) :: widthio,gammaio
		character(len=2), dimension(3) :: numrunp
		character(len=8) :: confil		

!		The data of the previous run is read to continue the simulation		

		confil="fs"//numruno(1)//numruno(2)//numruno(3)		
		open(unit=8,file=confil,status='unknown',form='unformatted',POSITION="REWIND")
		allocate (numhist(ihist))

		read(8) ihisto,numruno,numrunp,numruns,(numhist(i),i=1,ihisto),nstep,maxstp, &
				ndump,nprint,lplots,itime,dt0,nstep1,nonlin,m0dy,nocpl
		read(8) mj,lmaxo,leqmax,mjm1,mjm2,leqdim,ldim,jdim
		read(8) (r(j),j=0,mj),(mmo(l),l=1,lmaxo),(nno(l),l=1,lmaxo),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &
				(dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
				
		if (lmax == 0) lmax=lmaxo
		do l=1,lmaxo
			if (mm(l) /= 0 .or. nn(l) /= 0) exit
			mm(l)=mmo(l)
			nn(l)=nno(l)
		end do

		lmap=0
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			do l1=1,lmaxo
				if (mmo(l1) == m .and. nno(l1) == n) exit
			end do
			if (l1 <= lmaxo) lmap(l1)=l
		end do

		allocate (sqgi(0:mj,0:leqmax))
		allocate (sqg(0:mj,0:leqmax))
		allocate (bst(0:mj,0:leqmax))
		allocate (grr(0:mj,0:leqmax))
		allocate (grt(0:mj,0:leqmax))
		allocate (gtt(0:mj,0:leqmax))
		allocate (grz(0:mj,0:leqmax))
		allocate (gtz(0:mj,0:leqmax))
		allocate (gzz(0:mj,0:leqmax))				
		allocate (grroj(0:mj,0:leqmax))
		allocate (grtoj(0:mj,0:leqmax))
		allocate (gttoj(0:mj,0:leqmax))		
		allocate (grzoj(0:mj,0:leqmax))
		allocate (gtzoj(0:mj,0:leqmax))
		allocate (gzzoj(0:mj,0:leqmax))
		allocate (bmod(0:mj,0:leqmax))				
		allocate (grrup(0:mj,0:leqmax))	
		allocate (grtup(0:mj,0:leqmax))	
		allocate (grzup(0:mj,0:leqmax))	
		allocate (gttup(0:mj,0:leqmax))	
		allocate (gtzup(0:mj,0:leqmax))	
		allocate (gzzup(0:mj,0:leqmax))					
		allocate (jbgrr(0:mj,0:leqmax))
		allocate (jbgrt(0:mj,0:leqmax))
		allocate (jbgtt(0:mj,0:leqmax))
		allocate (lplr(0:mj,0:leqmax))
		allocate (lplt(0:mj,0:leqmax))		
		allocate (lplz(0:mj,0:leqmax))	
		allocate (djroj(0:mj,0:leqmax))	
		allocate (djtoj(0:mj,0:leqmax))	
		allocate (djzoj(0:mj,0:leqmax))	
		allocate (omdr(0:mj,0:leqmax))
		allocate (omdt(0:mj,0:leqmax))
		allocate (omdz(0:mj,0:leqmax))
		
		allocate (dbsjtoj(0:mj,0:leqmax))	
		allocate (dbsjzoj(0:mj,0:leqmax))
		allocate (dbsjtbj(0:mj,0:leqmax))
		allocate (dgttr(0:mj,0:leqmax))
		allocate (dgrrt(0:mj,0:leqmax))
		allocate (dgrtt(0:mj,0:leqmax))
		allocate (dgttt(0:mj,0:leqmax))
		allocate (dgrrz(0:mj,0:leqmax))
		allocate (dgrtz(0:mj,0:leqmax))
		allocate (dgttz(0:mj,0:leqmax))
		allocate (dgrtp(0:mj,0:leqmax))
		allocate (dgttp(0:mj,0:leqmax))
		allocate (jsq(0:mj,0:leqmax))
		allocate (bsgrt(0:mj,0:leqmax))
		allocate (bsgtt(0:mj,0:leqmax))
		allocate (bsq(0:mj,0:leqmax))
		allocate (bsqgtt(0:mj,0:leqmax))
		allocate (lplrr(0:mj,0:leqmax))
		allocate (lplrt(0:mj,0:leqmax))
		allocate (lplrz(0:mj,0:leqmax))
		allocate (lpltt(0:mj,0:leqmax))
		allocate (lpltz(0:mj,0:leqmax))
		allocate (lplzz(0:mj,0:leqmax))
					
		allocate (sqgdroj(0:mj,0:leqmax))	
		allocate (sqgdthoj(0:mj,0:leqmax))	
		allocate (sqgdztoj(0:mj,0:leqmax))							 			 
		allocate (sqgdthojbst(0:mj,0:leqmax))
		allocate (sqgdztojbst(0:mj,0:leqmax))	
		allocate (sqgibmodith(0:mj,0:leqmax))
		allocate (sqgibmodizt(0:mj,0:leqmax))		
							 		
        if(ieldamp_on .eq. 1 .and. ext_prof .eq. 1) then		
		  allocate (eildrr(0:mj,0:leqmax))		
		  allocate (eildrt(0:mj,0:leqmax))		
		  allocate (eildrz(0:mj,0:leqmax))			
		  allocate (eildtt(0:mj,0:leqmax))		
		  allocate (eildtz(0:mj,0:leqmax))		
		  allocate (eildzz(0:mj,0:leqmax))			
		  allocate (eildr(0:mj,0:leqmax))		
		  allocate (eildt(0:mj,0:leqmax))		
		  allocate (eildz(0:mj,0:leqmax))			  			  
		end if	

        if(Trapped_on .eq. 1) then	 
		allocate (omdrprp(0:mj,0:leqmax))
		allocate (omdtprp(0:mj,0:leqmax))
		allocate (omdzprp(0:mj,0:leqmax))	
		end if	
		
		sqgi=0.0_IDP
		sqg=0.0_IDP
		bst=0.0_IDP
		grr=0.0_IDP
		grt=0.0_IDP
		gtt=0.0_IDP
		grz=0.0_IDP
		gtz=0.0_IDP
		gzz=0.0_IDP		
		grroj=0.0_IDP
		grtoj=0.0_IDP
		gttoj=0.0_IDP
		grzoj=0.0_IDP
		gtzoj=0.0_IDP	
		gzzoj=0.0_IDP	
		bmod=0.0_IDP
		grrup=0.0_IDP
		grtup=0.0_IDP
		grzup=0.0_IDP
		gttup=0.0_IDP
		gtzup=0.0_IDP	
		gzzup=0.0_IDP		
		jbgrr=0.0_IDP
		jbgrt=0.0_IDP
		jbgtt=0.0_IDP
		lplr=0.0_IDP
		lplt=0.0_IDP
		omdr=0.0_IDP
		omdt=0.0_IDP
		omdz=0.0_IDP
		sqgdroj=0.0_IDP
		sqgdthoj=0.0_IDP
		sqgdztoj=0.0_IDP
		sqgdthojbst=0.0_IDP
		sqgdztojbst=0.0_IDP
		sqgibmodith=0.0_IDP
		sqgibmodizt=0.0_IDP

		dbsjtoj=0.0_IDP
		dbsjzoj=0.0_IDP
		dbsjtbj=0.0_IDP
		dgttr=0.0_IDP
		dgrrt=0.0_IDP
		dgrtt=0.0_IDP
		dgttt=0.0_IDP
		dgrrz=0.0_IDP
		dgrtz=0.0_IDP
		dgttz=0.0_IDP
		dgrtp=0.0_IDP
		dgttp=0.0_IDP
		jsq=0.0_IDP
		bsgrt=0.0_IDP
		bsgtt=0.0_IDP
		bsq=0.0_IDP
		bsqgtt=0.0_IDP
		lplrr=0.0_IDP
		lplrt=0.0_IDP
		lplrz=0.0_IDP
		lpltt=0.0_IDP
		lpltz=0.0_IDP
		lplzz=0.0_IDP

        if(Trapped_on .eq. 1) then
		omdrprp=0.0_IDP
		omdtprp=0.0_IDP
		omdzprp=0.0_IDP
		end if	

        if(ieldamp_on .eq. 1 .and. ext_prof .eq. 1) then		
		  eildrr=0.0_IDP
		  eildrt=0.0_IDP	
		  eildrz=0.0_IDP	
		  eildtt=0.0_IDP	
		  eildtz=0.0_IDP	
		  eildzz=0.0_IDP	
		  eildr=0.0_IDP	
		  eildt=0.0_IDP	
		  eildz=0.0_IDP		  
		end if	 				

		read(8) ni,nis,ne,delta,rc,fti,fte
		read(8) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(psieq(j),j=0,mj),(chieq(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
				(cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=1,mj),(vfova(j),j=1,mj), &
                                (vtherm_ion(j),j=1,mj),(vzt_eq(j),j=1,mj)

                if(alpha_on .eq. 1) then
		  read(8) (nalpeq(j),j=1,mj),(valphaova(j),j=1,mj)
		end if	

		read(8) ngeneq,ndevice(1),ndevice(2),eps,bet0,(rs(lmap(l)),l=1,lmaxo),leq,mjeq,nfp  

		read(8) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
				 ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
				 ((grz(j,l),j=0,mj),l=1,leqmax),((gtz(j,l),j=0,mj),l=1,leqmax),((gzz(j,l),j=0,mj),l=1,leqmax), &
				 ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
				 ((grzoj(j,l),j=0,mj),l=1,leqmax),((gtzoj(j,l),j=0,mj),l=1,leqmax),((gzzoj(j,l),j=0,mj),l=1,leqmax), &
				 ((grrup(j,l),j=0,mj),l=1,leqmax),((grtup(j,l),j=0,mj),l=1,leqmax),((grzup(j,l),j=0,mj),l=1,leqmax), &
				 ((gttup(j,l),j=0,mj),l=1,leqmax),((gtzup(j,l),j=0,mj),l=1,leqmax),((gzzup(j,l),j=0,mj),l=1,leqmax), &				 
				 ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
				 ((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &	
				 ((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &				 
				 ((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &	
				 ((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax),((dbsjtbj(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttr(j,l),j=0,mj),l=1,leqmax),((dgrrt(j,l),j=0,mj),l=1,leqmax),((dgrtt(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttt(j,l),j=0,mj),l=1,leqmax),((dgrrz(j,l),j=0,mj),l=1,leqmax),((dgrtz(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttz(j,l),j=0,mj),l=1,leqmax),((dgrtp(j,l),j=0,mj),l=1,leqmax),((dgttp(j,l),j=0,mj),l=1,leqmax), &
				 ((jsq(j,l),j=0,mj),l=1,leqmax),((bsgrt(j,l),j=0,mj),l=1,leqmax),((bsgtt(j,l),j=0,mj),l=1,leqmax), &
				 ((bsq(j,l),j=0,mj),l=1,leqmax),((bsqgtt(j,l),j=0,mj),l=1,leqmax),((lplrr(j,l),j=0,mj),l=1,leqmax), &
				 ((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
				 ((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax), &
				 ((lplzz(j,l),j=0,mj),l=1,leqmax),((bmod(j,l),j=0,mj),l=1,leqmax), &				 				 
				 ((sqgdroj(j,l),j=0,mj),l=1,leqmax),((sqgdthoj(j,l),j=0,mj),l=1,leqmax),((sqgdztoj(j,l),j=0,mj),l=1,leqmax), &
				 ((sqgdthojbst(j,l),j=0,mj),l=1,leqmax),((sqgdztojbst(j,l),j=0,mj),l=1,leqmax), &
				 ((sqgibmodith(j,l),j=0,mj),l=1,leqmax),((sqgibmodizt(j,l),j=0,mj),l=1,leqmax) 		
				 
        if(ieldamp_on .eq. 1) then
		  read(8)  ((eildrr(j,l),j=0,mj),l=1,leqmax),((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &					 
				   ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax),((eildzz(j,l),j=0,mj),l=1,leqmax), &		
				   ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax),((eildz(j,l),j=0,mj),l=1,leqmax)		
		end if	

        if(Trapped_on .eq. 1) then
		  read(8)  ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax),((omdzprp(j,l),j=0,mj),l=1,leqmax)
		end if	
				 
		read(8) (eta(j),j=0,mj),dt,time

		read(8) ((psi(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((phi(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((pr(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((nf(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((vprlf(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((vthprlf(j,lmap(l)),j=0,mj),l=1,lmaxo)		
        if(alpha_on .eq. 1) then 
		  read(8) ((nalp(j,lmap(l)),j=0,mj),l=1,lmaxo)	
		  read(8) ((vprlalp(j,lmap(l)),j=0,mj),l=1,lmaxo)
        end if		 
 
		read(8) etascl,reta,eta0,etalmb,ietaeq,stdifp,stdifu,stdifv,stdifnf,stdifvf,stdifnalp,stdifvalp, &
		        s,gamma,ipert,(widthio(l),l=1,lmaxo),(gammaio(l),l=1,lmaxo),pertsclo,norm_eildump

		rewind(8)
		close(unit=8)

	end subroutine rddump

	subroutine vmec

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none
 
		integer :: j,k,l,l1,lq,mjeqp,mjeqm1,mpol,ntor,lb0,m,mb,m2,mm1,mm2,mm12,mm22,mm2p1,mm2p12,mp2,mp22,isg,idummy,mband,n,nnn,mmm,lbm2, &
			   ir,jp,lp,nbst,lth,lzt,ithzt
		real(IDP) :: bigrn,pror,dummy,bet0in,twopi,twophip,one,gc,iota,iotap,p,pprimel,pprime,qmin,qmax,xl
		integer, dimension(:), allocatable :: llc,lls,lbst	
		real(IDP), dimension(:), allocatable :: rbinv,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,rsb,Abst,Wbst,phip_c
		real(IDP), dimension(:,:), allocatable :: rmnb,zmnb,pmnb,sqgib,sqgb,grrb,grtb,gttb,bmodb,grrojb,grtojb,gttojb,jbgrrb,jbgrtb,jbgttb
		real(IDP), dimension(:,:), allocatable :: rmneq,zmneq,pmneq

		character(len=32) :: formatt='("r",127(a1,i4,"/",i4))'
		character(len=32) :: formatv='(1pe13.6,127(a1,1pe15.8))'
		character(len=1) :: t
		character(len=16) :: confil
		integer :: m_sel,ltop

		interface
			subroutine mmblims
			end subroutine mmblims
			subroutine eqsplns(yeq,yeqr,yfar,m,m2,choice,smoo)
				use param
				implicit none
				integer :: m,m2
				character(len=6) :: choice
				real(IDP) :: smoo
				real(IDP), dimension(0:) :: yeq,yeqr
				real(IDP), dimension(0:) :: yfar
			end subroutine eqsplns
			subroutine fitter(yeq,yeqr,yfar,ir,m,mb)
				use param
				implicit none
				integer :: ir,m,mb
				real(IDP), dimension(0:) :: yeq,yeqr,yfar
			end subroutine fitter	
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine root(t,ft,b,c,relerr,abserr,iflag)
				use param
				implicit none
				integer :: iflag
				real(IDP) :: t,ft,b,c,relerr,abserr
			end subroutine root
		end interface

!  		reads in high beta equilibrium data from file eq_name on unit 25

		bet0in=bet0

		allocate (rmn(0:mj,0:leqmax),zmn(0:mj,0:leqmax),pmn(0:mj,0:leqmax))	  

!		Here we open the equilibria file and read the equilibria parameters		
 		open(unit=25,file=eq_name,status='old',convert='big_endian',form='unformatted')	

		read(25) nfp,lbmax,mjeqp,lasym

		mjeq=mjeqp-1
		mjeqm1=mjeq-1
		lbm2=2*lbmax-1

		if (lasym) then
			allocate (mmb(lbm2),nnb(lbm2),llc(lbmax),lls(lbmax))
		else
			allocate (mmb(lbmax),nnb(lbmax))
		end if	

		allocate (rfar(0:mjeq),rbinv(0:mjeq),qfar(0:mjeq),pfar(0:mjeq),phip(0:mjeq),curfar(0:mjeq),ffar(0:mjeq), &
			  sfar1(0:mjeq),sfar2(0:mjeq),sfar3(0:mjeq),sfar4(0:mjeq))
		if (lasym) then
			allocate (rmnb(0:mjeq,0:lbm2),zmnb(0:mjeq,0:lbm2),pmnb(0:mjeq,0:lbm2), &
                                  sqgb(0:mjeq,0:lbm2),sqgib(0:mjeq,0:lbm2),bmodb(0:mjeq,0:lbm2), &
				  grrb(0:mjeq,0:lbm2),grtb(0:mjeq,0:lbm2),gttb(0:mjeq,0:lbm2), &
				  grrojb(0:mjeq,0:lbm2),grtojb(0:mjeq,0:lbm2),gttojb(0:mjeq,0:lbm2), &
				  jbgrrb(0:mjeq,0:lbm2),jbgrtb(0:mjeq,0:lbm2),jbgttb(0:mjeq,0:lbm2))
		else
			allocate (rmnb(0:mjeq,0:lbmax),zmnb(0:mjeq,0:lbmax),pmnb(0:mjeq,0:lbmax), &
                                  sqgb(0:mjeq,0:lbmax),sqgib(0:mjeq,0:lbmax),bmodb(0:mjeq,0:lbmax), &
				  grrb(0:mjeq,0:lbmax),grtb(0:mjeq,0:lbmax),gttb(0:mjeq,0:lbmax), &
				  grrojb(0:mjeq,0:lbmax),grtojb(0:mjeq,0:lbmax),gttojb(0:mjeq,0:lbmax), &
				  jbgrrb(0:mjeq,0:lbmax),jbgrtb(0:mjeq,0:lbmax),jbgttb(0:mjeq,0:lbmax))
		end if

		rmnb(:,0)=0.0_IDP
		zmnb(:,0)=0.0_IDP
		pmnb(:,0)=0.0_IDP
		sqgib(:,0)=0.0_IDP
		sqgb(:,0)=0.0_IDP
		grrb(:,0)=0.0_IDP
		grtb(:,0)=0.0_IDP
		gttb(:,0)=0.0_IDP
		grrojb(:,0)=0.0_IDP
		grtojb(:,0)=0.0_IDP
		gttojb(:,0)=0.0_IDP

		read(25) (mmb(l),nnb(l),l=1,lbmax)
		read(25) (phip(j),qfar(j),curfar(j),ffar(j),pfar(j),j=1,mjeq)
		do j=1,mjeq
			read(25) (rmnb(j,l),zmnb(j,l),pmnb(j,l),bmodb(j,l),sqgb(j,l),sqgib(j,l),l=1,lbmax)
		end do
		do j=1,mjeq
			read(25) (grrb(j,l),grtb(j,l),gttb(j,l),grrojb(j,l),grtojb(j,l),gttojb(j,l),jbgrrb(j,l),jbgrtb(j,l),jbgttb(j,l),l=1,lbmax)
		end do

		lb0=0
		do l=1,lbmax
			if (mmb(l) == 0 .and. nnb(l) == 0) exit
		end do
		lb0=l
		if (lb0 == 0 .or. lb0 > lbmax) then
			write (6,'("  vmec: lb0=0")')
			stop
		end if

		if (lasym) then
			l1=lbmax
			do l=1,lbmax
				if (l == lb0) cycle
				l1=l1+1
				llc(l)=l1
				lls(l)=l1
			end do
			llc(lb0)=lb0
			lls(lb0)=0

			do j=1,mjeq
				read(25) (rmnb(j,lls(l)),zmnb(j,lls(l)),pmnb(j,lls(l)),bmodb(j,lls(l)),sqgb(j,lls(l)),sqgib(j,lls(l)),l=1,lbmax)
			end do
			do j=1,mjeq
				read(25) (grrb(j,lls(l)),grtb(j,llc(l)),gttb(j,lls(l)),grrojb(j,lls(l)),grtojb(j,llc(l)),gttojb(j,lls(l)), &
					  jbgrrb(j,lls(l)),jbgrtb(j,llc(l)),jbgttb(j,lls(l)),l=1,lbmax)
			end do
		end if

		do l=1,lbmax
			if (nnb(l) == 0 .and. mmb(l) < 0) then
				mmb(l)=-mmb(l)
				zmnb(:,l)=-zmnb(:,l)
				pmnb(:,l)=-pmnb(:,l)
				if (lasym) rmnb(:,lls(l))=-rmnb(:,lls(l))
			end if
			if (nnb(l) < 0) then
				mmb(l)=-mmb(l)
				nnb(l)=-nnb(l)
				zmnb(:,l)=-zmnb(:,l)
				pmnb(:,l)=-pmnb(:,l)
				if (lasym) rmnb(:,lls(l))=-rmnb(:,lls(l))
			end if
		end do

		if (lasym) then
			l1=lbmax
			do l=1,lbmax
				if (l == lb0) cycle
				l1=l1+1
				mmb(l1)=-mmb(l)
				nnb(l1)=-nnb(l)
			end do
		else
			lbm2=lbmax
		end if

		do j=0,mjeq
			rfar(j)=sqrt(1.0_IDP*j/mjeq)
		end do

		rbinv(0)=0.0_IDP
		do j=1,mjeq
			rbinv(j)=1.0_IDP/rfar(j)
		end do

		rmnb(0,:)=0.0_IDP
		zmnb(0,:)=0.0_IDP
		pmnb(0,:)=0.0_IDP

		do l=1,lbm2
			if (mmb(l) == 0) then
				rmnb(0,l)=1.5*rmnb(1,l)-0.5*rmnb(2,l)
				zmnb(0,l)=1.5*zmnb(1,l)-0.5*zmnb(2,l)
				pmnb(0,l)=1.5*pmnb(1,l)-0.5*pmnb(2,l)
			endif
		end do

!		Here we use splines to transfer the equilibria parameters to code parameters
                allocate (rmneq(0:mj_r,0:lbm2),zmneq(0:mj_r,0:lbm2),pmneq(0:mj_r,0:lbm2))	

		rmneq=0.0_IDP
		zmneq=0.0_IDP
		pmneq=0.0_IDP

		call grid

!		t=char(9)
!		open(unit=92,file="rmnb.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=1,mjeq
!			write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-0.5_IDP)/mjeq),(t,rmnb(j,l),l=1,lbm2)
!		end do
!		close(unit=92)
!		open(unit=92,file="zmnb.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=1,mjeq
!			write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-0.5_IDP)/mjeq),(t,zmnb(j,l),l=1,lbm2)
!		end do
!		close(unit=92)
!		open(unit=92,file="pmnb.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=1,mjeq
!			write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-0.5_IDP)/mjeq),(t,pmnb(j,l),l=1,lbm2)
!		end do
!		close(unit=92)

		ir=max(15,mjeq/10)
		do l=1,lbm2
			mb=abs(mmb(l))
			m=min(mb,30)
			dummy=1.5*rmnb(mjeq,l)-0.5*rmnb(mjeqm1,l)
			do j=1,mjeqm1
				rmnb(j,l)=0.5*(rmnb(j+1,l)+rmnb(j,l))
			end do
			rmnb(mjeq,l)=dummy
			call fitter(rmneq(:,l),sd1_r,rmnb(:,l),ir,m,mb)
			dummy=1.5*zmnb(mjeq,l)-0.5*zmnb(mjeqm1,l)
			do j=1,mjeqm1
				zmnb(j,l)=0.5*(zmnb(j+1,l)+zmnb(j,l))
			end do
			zmnb(mjeq,l)=dummy
			call fitter(zmneq(:,l),sd1_r,zmnb(:,l),ir,m,mb)
			dummy=1.5*pmnb(mjeq,l)-0.5*pmnb(mjeqm1,l)
			do j=1,mjeqm1
				pmnb(j,l)=0.5*(pmnb(j+1,l)+pmnb(j,l))
			end do
			pmnb(mjeq,l)=dummy
			call fitter(pmneq(:,l),sd1_r,pmnb(:,l),ir,m,mb)
		end do

!		open(unit=92,file="rmnfit.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=0,mj_r
!			write(92,'(f7.5,500(a1,1pe13.6))') r_red(j),(t,rmneq(j,l),l=1,lbm2)
!		end do
!		close(unit=92)
!		open(unit=92,file="zmnfit.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=0,mj_r
!			write(92,'(f7.5,500(a1,1pe13.6))') r_red(j),(t,zmneq(j,l),l=1,lbm2)
!		end do
!		close(unit=92)
!		open(unit=92,file="pmnfit.txt",recl=8192)
!		write(92,'("r",500(a1,i4,"/",i4))') (t,mmb(l),nnb(l),l=1,lbm2)
!		do j=0,mj_r
!			write(92,'(f7.5,500(a1,1pe13.6))') r_red(j),(t,pmneq(j,l),l=1,lbm2)
!		end do
!		close(unit=92)

!	neta is the least common multiple of lth_r and lzt_r
		neta=max(lth_r,lzt_r)
		do
			if (mod(neta,lth_r) == 0 .and. mod(neta,lzt_r) == 0) exit
			neta=neta+1
		end do
		ith=neta/lth_r
		izt=neta/lzt_r
		twopi=8.0_IDP*atan(1.0_IDP)
		allocate (cs(0:neta),sn(0:neta))
		do k=0,neta-1
			cs(k)=cos(k*twopi/neta)
			sn(k)=sin(k*twopi/neta)
		end do
		cs(neta)=1.0
		sn(neta)=0.0

		do lzt=0,lzt_r
			do lth=0,lth_r
				rmn_c(:,lth,lzt)=0.0
				zmn_c(:,lth,lzt)=0.0
				pmn_c(:,lth,lzt)=lzt*twopi/lzt_r
				do l=1,lbmax
					ithzt=modulo(mmb(l)*ith*lth+nnb(l)*izt*lzt,neta)
					rmn_c(:,lth,lzt)=rmn_c(:,lth,lzt)+rmneq(:,l)*cs(ithzt)
					zmn_c(:,lth,lzt)=zmn_c(:,lth,lzt)+zmneq(:,l)*sn(ithzt)
					pmn_c(:,lth,lzt)=pmn_c(:,lth,lzt)+pmneq(:,l)*sn(ithzt)
					if (lasym) then
						rmn_c(:,lth,lzt)=rmn_c(:,lth,lzt)+rmneq(:,lls(l))*sn(ithzt)
						zmn_c(:,lth,lzt)=zmn_c(:,lth,lzt)+zmneq(:,llc(l))*cs(ithzt)
						pmn_c(:,lth,lzt)=pmn_c(:,lth,lzt)+pmneq(:,llc(l))*cs(ithzt)
					end if
				end do
			end do
		end do
		lbmax=lbm2

		deallocate (rmneq)
		deallocate (zmneq)
		deallocate (pmneq)

		deallocate (mmb)
		deallocate (nnb)
		deallocate (rfar)
		deallocate (rbinv)
		deallocate (qfar)
		deallocate (pfar)
		deallocate (phip)
		deallocate (curfar)
		deallocate (ffar)
		deallocate (sqgib)
		deallocate (sqgb)
		deallocate (grrb)
		deallocate (grtb)
		deallocate (gttb)
		deallocate (bmodb)
		deallocate (grrojb)
		deallocate (grtojb)
		deallocate (gttojb)
		deallocate (jbgrrb)
		deallocate (jbgrtb)
		deallocate (jbgttb)
		deallocate (rmnb)
		deallocate (zmnb)
		deallocate (pmnb)

        end subroutine vmec

	subroutine sum_var

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: lth,lzt,ithzt
		real(IDP) :: scmn,scmx,sclmn,sclmx,scnorm
		character(len=1) :: t
		real(IDP) :: pi = 3.1415926
		integer :: i,l,j,iend,n,m,m2

		character(len=32) :: formatv='(1pe15.8,3(a1,1pe15.8))'
		character(len=32) :: formatt='("r",127(a1,i4,"/",i4))'
		character(len=32) :: formatx='(1pe13.6,127(a1,1pe15.8))'
		character(len=12) :: confil

		real(IDP), dimension(:,:), allocatable :: psi_r,phi_r,pr_r,uzt_r,nf_r,vprlf_r,vthprlf_r, &
                                                          Vr_r,Vth_r,Br_r,Bth_r,sqgi_r,nalp_r,vprlalp_r
		real(IDP), dimension(:), allocatable :: preqred,nfeqred,qqinvred,rred,nalpeqred

		interface
			subroutine dbydth(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydth
			subroutine dbydr(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydr
			subroutine dbydzt(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzt
			subroutine eqtodyn(adyn,aeq,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: adyn,aeq
			end subroutine eqtodyn
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine mult
			subroutine splns(yeq,yeqr,yfar,m,m2,choice,smoo)
				use param
				implicit none
				integer :: m,m2
				character(len=6) :: choice
				real(IDP) :: smoo
				real(IDP), dimension(0:) :: yeq,yeqr
				real(IDP), dimension(0:) :: yfar
			end subroutine splns
		end interface

                psi_b(:,:,:) = 0.0
                phi_b(:,:,:) = 0.0
                pr_b(:,:,:) = 0.0
                uzt_b(:,:,:) = 0.0
                nf_b(:,:,:) = 0.0
                vprlf_b(:,:,:) = 0.0
                vthprlf_b(:,:,:) = 0.0
                Vr_b(:,:,:) = 0.0
                Vth_b(:,:,:) = 0.0
                Br_b(:,:,:) = 0.0
                Bth_b(:,:,:) = 0.0
                Bzt_b(:,:,:) = 0.0
                sqgi_b(:,:,:) = 0.0
                if(alpha_on .eq. 1) then
                    nalp_b(:,:,:) = 0.0
                    vprlalp_b(:,:,:) = 0.0
                end if
                xx1(:,:,:) = 0.0
                xx2(:,:,:) = 0.0
                xx3(:,:,:) = 0.0

                sc1(:,:) = 0.0
                sc2(:,:) = 0.0
                sc3(:,:) = 0.0
                sc4(:,:) = 0.0
                sc5(:,:) = 0.0
                sc6(:,:) = 0.0
                sc7(:,:) = 0.0

!  vr up
		call dbydth(sc1,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr(sc2,phi,0.0_IDP,1.0_IDP,0)

		scmn=0.0_IDP
		scmx=0.0_IDP
		do l=1,lmaxn
			sclmn=minval(sc1(:,l))
			sclmx=maxval(sc1(:,l))
			scmn=min(scmn,sclmn)
			scmx=max(scmx,sclmx)
		end do
		if(scmx > abs(scmn)) then
			scnorm=scmx
		else
			scnorm=scmn
		end if

		if (scnorm /= 0.0_IDP) then
			do l=1,lmaxn
				psi(:,l)=psi(:,l)/scnorm
!				where (abs(psi(:,l)) < 1.e-50_IDP) psi(:,l)=0
				phi(:,l)=phi(:,l)/scnorm
!				where (abs(phi(:,l)) < 1.e-50_IDP) phi(:,l)=0
				pr(:,l)=pr(:,l)/scnorm
!				where (abs(pr(:,l)) < 1.e-50_IDP) pr(:,l)=0
				uzt(:,l)=uzt(:,l)/scnorm
!				where (abs(uzt(:,l)) < 1.e-50_IDP) uzt(:,l)=0
				sc1(:,l)=sc1(:,l)/scnorm
!				where (abs(sc1(:,l)) < 1.e-50_IDP) sc1(:,l)=0
				sc2(:,l)=sc2(:,l)/scnorm
!				where (abs(sc2(:,l)) < 1.e-50_IDP) sc2(:,l)=0
				nf(:,l)=nf(:,l)/scnorm
!				where (abs(nf(:,l)) < 1.e-50_IDP) nf(:,l)=0
				vprlf(:,l)=vprlf(:,l)/scnorm
!				where (abs(vprlf(:,l)) < 1.e-50_IDP) vprlf(:,l)=0
				vthprlf(:,l)=vthprlf(:,l)/scnorm
!				where (abs(vthprlf(:,l)) < 1.e-50_IDP) vthprlf(:,l)=0
				if(alpha_on .eq. 1) then
				  nalp(:,l)=nalp(:,l)/scnorm
!				  where (abs(nalp(:,l)) < 1.e-50_IDP) nalp(:,l)=0
				  vprlalp(:,l)=vprlalp(:,l)/scnorm
!				  where (abs(vprlalp(:,l)) < 1.e-50_IDP) vprlalp(:,l)=0
				end if
			end do
		end if

!  br up
		call dbydth(sc3,psi,1,0.0_IDP,-s,0)
!  bth up
		call dbydr(sc4,psi,0.0_IDP,s,0)

		call eqtodyn(sc7,sqgi,0.0_IDP,1.0_IDP)

		call mult(sc5,sc7,1,sc3,-1,0.0_IDP,1.0_IDP)
		call mult(sc6,sc7,1,sc4,1,0.0_IDP,1.0_IDP)

		allocate (psi_r(0:mj_r,0:ldim),phi_r(0:mj_r,0:ldim),pr_r(0:mj_r,0:ldim),uzt_r(0:mj_r,0:ldim), &
                          nf_r(0:mj_r,0:ldim),vprlf_r(0:mj_r,0:ldim),vthprlf_r(0:mj_r,0:ldim),Vr_r(0:mj_r,0:ldim), &
                          Vth_r(0:mj_r,0:ldim),Br_r(0:mj_r,0:ldim),Bth_r(0:mj_r,0:ldim),sqgi_r(0:mj_r,0:ldim), &
                          nalp_r(0:mj_r,0:ldim),vprlalp_r(0:mj_r,0:ldim))

		do l=1,lmaxn
			m=abs(mm(l))
			m2=m
			call splns(psi_r(:,l),sd1_r,psi(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(phi_r(:,l),sd1,phi(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(pr_r(:,l),sd1,pr(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(uzt_r(:,l),sd1,uzt(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(nf_r(:,l),sd1,nf(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(vprlf_r(:,l),sd1,vprlf(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(vthprlf_r(:,l),sd1,vthprlf(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(Vr_r(:,l),sd1,sc1(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(Vth_r(:,l),sd1,sc2(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(Br_r(:,l),sd1,sc5(:,l),m,m2,"icsscu",1.e-07_IDP)
			call splns(Bth_r(:,l),sd1,sc6(:,l),m,m2,"icsscu",1.e-07_IDP)
!			call splns(sqgi_r(:,l),sd1,sqgi(:,l),m,m2,"icsscu",1.e-07_IDP)
                       if(alpha_on .eq. 1) then
			    call splns(nalp_r(:,l),sd1,nalp(:,l),m,m2,"icsscu",1.e-07_IDP)
			    call splns(vprlalp_r(:,l),sd1,vprlalp(:,l),m,m2,"icsscu",1.e-07_IDP)
                        end if
		end do

!               Check interpolation from rrdump variables to vtk file dimensions

!		write(confil,'("test_pr")')
!		open(unit=200,file=confil,recl=2048)
!		write(200,formatt) (t,mm(l),nn(l),l=1,lmaxn)
!		do j=1,mj
!			write(200,formatx) r(j),(t,pr(j,l),l=1,lmaxn)
!		end do
!		close(200)
!		write(confil,'("test_pr_r")')
!		open(unit=201,file=confil,recl=2048)
!		write(201,formatt) (t,mm(l),nn(l),l=1,lmaxn)
!		do j=1,mj_r
!			write(201,formatx) r_red(j),(t,pr_r(j,l),l=1,lmaxn)
!		end do
!		close(201)

		do lth = 0,lth_r
			do lzt = 0,lzt_r
				do l=1,lmaxn
					if (signl(lln(l)) == 0) then  !(0,0) part
						psi_b(:,lth,lzt)=psi_b(:,lth,lzt)+psi_r(:,lln(l))
						phi_b(:,lth,lzt)=phi_b(:,lth,lzt)+phi_r(:,lln(l))
						pr_b(:,lth,lzt)=pr_b(:,lth,lzt)+pr_r(:,lln(l))
						uzt_b(:,lth,lzt)=uzt_b(:,lth,lzt)+uzt_r(:,lln(l))
						nf_b(:,lth,lzt)=nf_b(:,lth,lzt)+nf_r(:,lln(l))
						vprlf_b(:,lth,lzt)=vprlf_b(:,lth,lzt)+vprlf_r(:,lln(l))
						vthprlf_b(:,lth,lzt)=vthprlf_b(:,lth,lzt)+vthprlf_r(:,lln(l))
						Vr_b(:,lth,lzt)=Vr_b(:,lth,lzt)+Vr_r(:,lln(l))
						Vth_b(:,lth,lzt)=Vth_b(:,lth,lzt)+Vth_r(:,lln(l))
						Br_b(:,lth,lzt)=Br_b(:,lth,lzt)+Br_r(:,lln(l))
						Bth_b(:,lth,lzt)=Bth_b(:,lth,lzt)+Bth_r(:,lln(l))
!						sqgi_b(:,lth,lzt)=sqgi_b(:,lth,lzt)+sqgi_r(:,lln(l))
						if(alpha_on .eq. 1) then
							nalp_b(:,lth,lzt)=nalp_b(:,lth,lzt)+nalp_r(:,lln(l))
							vprlalp_b(:,lth,lzt)=vprlalp_b(:,lth,lzt)+vprlalp_r(:,lln(l))
						end if
					else if (signl(lln(l)) > 0) then  !sin part
						ithzt=modulo(mm(lln(l))*ith*lth+nn(lln(l))*izt*lzt,neta)
						psi_b(:,lth,lzt)=psi_b(:,lth,lzt)+psi_r(:,lln(l))*cs(ithzt)
						phi_b(:,lth,lzt)=phi_b(:,lth,lzt)+phi_r(:,lln(l))*sn(ithzt)
						pr_b(:,lth,lzt)=pr_b(:,lth,lzt)+pr_r(:,lln(l))*cs(ithzt)
						uzt_b(:,lth,lzt)=uzt_b(:,lth,lzt)+uzt_r(:,lln(l))*sn(ithzt)
						nf_b(:,lth,lzt)=nf_b(:,lth,lzt)+nf_r(:,lln(l))*cs(ithzt)
						vprlf_b(:,lth,lzt)=vprlf_b(:,lth,lzt)+vprlf_r(:,lln(l))*sn(ithzt)
						vthprlf_b(:,lth,lzt)=vthprlf_b(:,lth,lzt)+vthprlf_r(:,lln(l))*sn(ithzt)
						Vr_b(:,lth,lzt)=Vr_b(:,lth,lzt)+Vr_r(:,lln(l))*cs(ithzt)
						Vth_b(:,lth,lzt)=Vth_b(:,lth,lzt)+Vth_r(:,lln(l))*sn(ithzt)
						Br_b(:,lth,lzt)=Br_b(:,lth,lzt)+Br_r(:,lln(l))*sn(ithzt)
						Bth_b(:,lth,lzt)=Bth_b(:,lth,lzt)+Bth_r(:,lln(l))*cs(ithzt)
!						sqgi_b(:,lth,lzt)=sqgi_b(:,lth,lzt)+sqgi_r(:,lln(l))*cs(ithzt)
						if(alpha_on .eq. 1) then
							nalp_b(:,lth,lzt)=nalp_b(:,lth,lzt)+nalp_r(:,lln(l))*cs(ithzt)
							vprlalp_b(:,lth,lzt)=vprlalp_b(:,lth,lzt)+vprlalp_r(:,lln(l))*sn(ithzt)
						end if
					else  !cos part
						ithzt=modulo(-mm(lln(l))*ith*lth-nn(lln(l))*izt*lzt,neta)
						psi_b(:,lth,lzt)=psi_b(:,lth,lzt)+psi_r(:,lln(l))*sn(ithzt)
						phi_b(:,lth,lzt)=phi_b(:,lth,lzt)+phi_r(:,lln(l))*cs(ithzt)
						pr_b(:,lth,lzt)=pr_b(:,lth,lzt)+pr_r(:,lln(l))*sn(ithzt)
						uzt_b(:,lth,lzt)=uzt_b(:,lth,lzt)+uzt_r(:,lln(l))*cs(ithzt)
						nf_b(:,lth,lzt)=nf_b(:,lth,lzt)+nf_r(:,lln(l))*sn(ithzt)
						vprlf_b(:,lth,lzt)=vprlf_b(:,lth,lzt)+vprlf_r(:,lln(l))*cs(ithzt)
						vthprlf_b(:,lth,lzt)=vthprlf_b(:,lth,lzt)+vthprlf_r(:,lln(l))*cs(ithzt)
						Vr_b(:,lth,lzt)=Vr_b(:,lth,lzt)+Vr_r(:,lln(l))*sn(ithzt)
						Vth_b(:,lth,lzt)=Vth_b(:,lth,lzt)+Vth_r(:,lln(l))*cs(ithzt)
						Br_b(:,lth,lzt)=Br_b(:,lth,lzt)+Br_r(:,lln(l))*cs(ithzt)
						Bth_b(:,lth,lzt)=Bth_b(:,lth,lzt)+Bth_r(:,lln(l))*sn(ithzt)
!						sqgi_b(:,lth,lzt)=sqgi_b(:,lth,lzt)+sqgi_r(:,lln(l))*sn(ithzt)
						if(alpha_on .eq. 1) then
							nalp_b(:,lth,lzt)=nalp_b(:,lth,lzt)+nalp_r(:,lln(l))*sn(ithzt)
							vprlalp_b(:,lth,lzt)=vprlalp_b(:,lth,lzt)+vprlalp_r(:,lln(l))*cs(ithzt)
						end if
					end if
				end do
			end do
		end do

!       The variables mapped to the cylindrical coordinates are:
!                        ( A_b(j,lth,lzt),rmn_c(j,lth,lzt),pmn_c(j,lth,lzt),zmn_c(j,lth,lzt) )

!       Variables transformed to cartesian coordinates:

                do j = 0,mj_r
			 do lth = 0,lth_r
			    do lzt = 0,lzt_r
                              xx1(j,lth,lzt) = rmn_c(j,lth,lzt)*cos(pmn_c(j,lth,lzt))
                              xx2(j,lth,lzt) = -rmn_c(j,lth,lzt)*sin(pmn_c(j,lth,lzt))
                              xx3(j,lth,lzt) = zmn_c(j,lth,lzt)
		            end do
		        end do
                end do  
  
!               Check transformation to cartesian

!		write(confil,'("test_pr_c")')
!		open(unit=200,file=confil,recl=2048)
!		write(200,'("R",a1,"Theta",a1,"Z",a1,"Pr")') t,t,t
!		do j=1,mj_r
!			write(200,formatv) xx1(j,1,1),t,xx2(j,1,1),t,xx3(j,1,1),t,pr_b(j,1,1)                   
!		end do
!		close(200) 

!       The variables mapped to the cartesian coordinates are:
!                        ( A_b(j,lth,lzt),xx1(j,lth,lzt),xx2(j,lth,lzt),xx3(j,lth,lzt) )

		deallocate (rmn)
		deallocate (zmn)
		deallocate (pmn)
		deallocate (rmn_c)
		deallocate (zmn_c)
		deallocate (pmn_c)

		deallocate (psi_r)
		deallocate (phi_r)
		deallocate (pr_r)
		deallocate (uzt_r)
		deallocate (nf_r)
		deallocate (vprlf_r)
		deallocate (vthprlf_r)
		deallocate (Vr_r)
		deallocate (Vth_r)
		deallocate (Br_r)
		deallocate (Bth_r)
		deallocate (sqgi_r)
                if(alpha_on .eq. 1) then
		    deallocate (nalp_r)
		    deallocate (vprlalp_r)
                end if
  
        end subroutine sum_var

	subroutine vtk_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Write data in VTK format
!
!  Collection of basic functions to write legacy VTK binary files.
!  Files can be written in serial or parallel mode and consists of the
!  basic five parts :

!  1. File version and identifier
!  2. Header, consisting of a string
!  3. File format
!  4. Dataset structure: describes the geometry and topology of the dataset     
!  5. Dataset attributes. This section is used to write the actual binary data as a vector or scalar data   

!  The mesh topology and the variable datasets are written using single 
!  precision (4 bytes) binary format.
!  VTK file are written following big endian order. 

!  The WriteVTK_Header() function provides the basic functionality for 
!  steps 1, 2, 3 and 4. The default grid topology is a structured grid.
   
!  The data is written as scalars: 

!  Reference:
!  http://www.vtk.org/VTK/img/file-formats.pdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: j,lth,lzt,nx1,nx2,nx3,points
		character(len=12) :: confil
		character(len=24) :: text
		character(len=1) :: t
		real(IDP), dimension(0:jdim,0:2) :: node_coord
		real(IDP), dimension(0:jdim) :: data_in

!!!!!!!!!!!!!!  Get global domain sizes

                nx1 = mj_r+1;
                nx2 = lth_r+1;
                nx3 = lzt_r+1;

!!!!!!!!!!!!!!  Open the file

		write(confil,'("pr_",2a2,".vtk")') numrun(1),numrun(2)
		open(unit=100,file=confil,status='new',recl=2048)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!   1. File version and identifier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

		write(100,'("# vtk DataFile Version 2.0")')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   2. Header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(100,'("FAR3d VTK Data")')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   3. File format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(100,'("ASCII")')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   4. Dataset structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(100,'("DATASET STRUCTURED_GRID")')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   5. Mesh dimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(100,'("DIMENSIONS ",i4," ",i2," ",i2))') nx1,nx2,nx3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   6. Allocate memory and define node coordinates 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                points = nx1*nx2*nx3
		write(100,'("POINTS ",i10," float")') points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   7. Write structured grid information 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		do j=0,mj_r
                    do lth = 0,lth_r
                        do lzt = 0,lzt_r
			    write(100,'(1pe15.8," ",1pe15.8," ",1pe15.8)') xx1(j,lth,lzt),xx2(j,lth,lzt),xx3(j,lth,lzt)
		        end do
		    end do
		end do

		write(100,'("POINT_DATA ",i10)') points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   8. Write scalar data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		write(100,'("SCALARS Pr float")')
		write(100,'("LOOKUP_TABLE default")')

		do j=0,mj_r
                    do lth = 0,lth_r
                        do lzt = 0,lzt_r
			write(100,'(1pe15.8)') pr_b(j,lth,lzt)
		        end do
		    end do
		end do

                close(100)

!                 psi_b(j,lth,lzt)
!                 phi_b(j,lth,lzt)
!                 pr_b(j,lth,lzt)
!                 uzt_b(j,lth,lzt)
!                 nf_b(j,lth,lzt)
!                 vprlf_b(j,lth,lzt)
!                 vthprlf_b(j,lth,lzt)
!                 Vr_b(j,lth,lzt)
!                 Vth_b(j,lth,lzt)
!                 Br_b(j,lth,lzt)
!                 Bth_b(j,lth,lzt)
!                 Bzt_b(j,lth,lzt)

	end subroutine vtk_file

	subroutine dbydth(d,a,ltype,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,ltype,l,m
		real(IDP) :: c1,c2,temp
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in theta 	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*rinv*d(a)/dth		
		
		a(:,0)=0.0_IDP
		do l=1,lmax
			dold=d(:,l) 
			temp=-mm(l)*ltype
			d(:,l)=temp*rinv*a(:,l) 
			d(0,l)=0.0_IDP
			m=abs(mm(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
			d(:,l)=c1*dold+c2*d(:,l) 
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydth

	subroutine dbydr(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in rho 	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr		
		
		a(:,0)=0.0_IDP
		do l=1,lmax
			dold=d(:,l) 
			do j=1,mjm1
				d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
			end do
			d(0,l)=0.0_IDP
			m=abs(mm(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=rinv(1)*a(1,l)
			d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
			d(:,l)=c1*dold+c2*d(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydr

	subroutine dbydr0(d,a,c1,c2,k)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: k,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in rho (variables with only radial dependency)	
!		d = term name
!		a = element to be derivated
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr				
		
		dold=d
		do j=1,mjm1
			d(j)=dc1m(j)*(a(j-1)-a(j))+dc1p(j)*(a(j+1)-a(j))
		end do
		d(0)=0.
		m=0
		if (k > 1) m=abs(abs(m-1)-1)
		if (k == 1 .or. k == 3) m=m+1
		if (m == 1) d(0)=rinv(1)*a(1)
		d(mj)=(a(mj)-a(mjm1))/(r(mj)-r(mjm1))
		d=c1*dold+c2*d

	end subroutine dbydr0

	subroutine eqtodyn(adyn,aeq,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: adyn,aeq
		integer :: le,l
		
!		Equilibrium variables are transformed to dynamic variables		
		
		adyn=c1*adyn
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			adyn(:,l)=adyn(:,l)+c2*aeq(:,le) 
		end do
		adyn(:,0)=0. 
		aeq(:,0)=0. 

	end subroutine eqtodyn

	subroutine mult (f, g, gtype, h, htype, c1, c2)

		use param
		use domain
		implicit none

		integer :: gtype, htype,index

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f

		interface
			subroutine multc11 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc11
			subroutine multc12 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc12
			subroutine multc22 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc22
			subroutine multr11 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr11
			subroutine multr12 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr12
			subroutine multr22 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr22
		end interface

		index=abs(gtype*htype)
		if (nmin < 0) then
			select case (index)
				case (1)
					call multc11 (f, g, h, gtype, htype, c1, c2)
				case (2)
					if (abs(gtype) == 2) then
						call multc12 (f, g, h, gtype, htype, c1, c2)
					else
						call multc12 (f, h, g, htype, gtype, c1, c2)
					end if
				case (4)
					call multc22 (f, g, h, gtype, htype, c1, c2)
			end select
		else
			select case (index)
				case (1)
					call multr11 (f, g, h, gtype, htype, c1, c2)
				case (2)
					if (abs(gtype) == 2) then
						call multr12 (f, g, h, gtype, htype, c1, c2)
					else
						call multr12 (f, h, g, htype, gtype, c1, c2)
					end if
				case (4)
					call multr22 (f, g, h, gtype, htype, c1, c2)
			end select
		end if

	end subroutine mult

	subroutine multc11 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.


		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmax
			fstart = mmstart (frow)
			fend = mmend (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = frow-nmax, nmax

				!   Loop through all row pairs whose sum is frow

				hrow = frow - grow

				!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
				!   are [fstart-gend, fend-gstart] and [hstart,hend]

				gstart = mmstart (grow)
				gend = mmend (grow)
				hstart = mmstart (hrow)
				hend = mmend (hrow)

				if (gend < gstart .or. hend < hstart) cycle
				if (fend-gstart < hstart .or. fstart-gend > hend) cycle

				!   Compute all contributions to the current frow from the current grow-hrow
				!   pair

				do midx = fstart,fend

					!   calculate the loop limits for this particular midx. The ranges are
					!   [midx-gend,midx-gstart] and [hstart,hend]

					mstart = hstart
					if (midx-gend > hstart) mstart = midx-gend
					mend = hend
					if (midx-gstart < hend) mend = midx-gstart

					if (mend < mstart) cycle

					!   Do some convolving! Since only the positive bands of the complex
					!   arrays are stored, negative band values are extracted from corresponding
					!   positive bands. Also, the bands are not stored directly at their indices,
					!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

					!   It is possible that the following loops may be speeded up somewhat by
					!   moving the address calculations into the loop limits: i.e. let hidx
					!   be the do loop index directly.

					if (grow >= 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
						end do

					else if (grow < 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = -g1idx - mmstart(-grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
						end do

					else
!  						(grow >= 0 .and. hrow < 0)

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = -h1idx - mmstart(-hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
						end do

					end if

				end do

			end do

			!   The current row of f has now been completed. Store it back in its
			!   packed form and loop to the next result row

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc11

	subroutine multc12 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.


		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmax
			fstart = mmstart (frow)
			fend = mmend (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

				!   Loop through all row pairs whose sum is frow

				hrow = frow - grow

				!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
				!   are [fstart-gend, fend-gstart] and [hstart,hend]

				gstart = mmstart (grow)
				gend = mmend (grow)
				hstart = mmstart (hrow)
				hend = mmend (hrow)

				if (gend < gstart .or. hend < hstart) cycle
				if (fend-gstart < hstart .or. fstart-gend > hend) cycle

				!   Compute all contributions to the current frow from the current grow-hrow
				!   pair

				do midx = fstart,fend

					!   calculate the loop limits for this particular midx. The ranges are
					!   [midx-gend,midx-gstart] and [hstart,hend]

					mstart = hstart
					if (midx-gend > hstart) mstart = midx-gend
					mend = hend
					if (midx-gstart < hend) mend = midx-gstart

					if (mend < mstart) cycle

					!   Do some convolving! Since only the positive bands of the complex
					!   arrays are stored, negative band values are extracted from corresponding
					!   positive bands. Also, the bands are not stored directly at their indices,
					!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

					!   It is possible that the following loops may be speeded up somewhat by
					!   moving the address calculations into the loop limits: i.e. let hidx
					!   be the do loop index directly.

					if (grow >= 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
						end do

					else if (grow < 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = -g1idx - mmstart(-grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
						end do

					else
!  						(grow >= 0 .and. hrow < 0)

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = -h1idx - mmstart(-hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
						end do

					end if

				end do

			end do

			!   The current row of f has now been completed. Store it back in its
			!   packed form and loop to the next result row

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc12

	subroutine multc22 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.


		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmaxeq, nfp
			fstart = mmstart (frow)
			fend = mmend (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = frow-nmaxeq, nmaxeq, nfp

				!   Loop through all row pairs whose sum is frow

				hrow = frow - grow

				!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
				!   are [fstart-gend, fend-gstart] and [hstart,hend]

				gstart = mmstart (grow)
				gend = mmend (grow)
				hstart = mmstart (hrow)
				hend = mmend (hrow)

				if (gend < gstart .or. hend < hstart) cycle
				if (fend-gstart < hstart .or. fstart-gend > hend) cycle

				!   Compute all contributions to the current frow from the current grow-hrow
				!   pair

				do midx = fstart,fend

					!   calculate the loop limits for this particular midx. The ranges are
					!   [midx-gend,midx-gstart] and [hstart,hend]

					mstart = hstart
					if (midx-gend > hstart) mstart = midx-gend
					mend = hend
					if (midx-gstart < hend) mend = midx-gstart

					if (mend < mstart) cycle

					!   Do some convolving! Since only the positive bands of the complex
					!   arrays are stored, negative band values are extracted from corresponding
					!   positive bands. Also, the bands are not stored directly at their indices,
					!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

					!   It is possible that the following loops may be speeded up somewhat by
					!   moving the address calculations into the loop limits: i.e. let hidx
					!   be the do loop index directly.

					if (grow >= 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
						end do

					else if (grow < 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstart(hrow)
							gidx = -g1idx - mmstart(-grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
						end do

					else
!  						(grow >= 0 .and. hrow < 0)

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = -h1idx - mmstart(-hrow)
							gidx = g1idx - mmstart(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
						end do

					end if

				end do

			end do

			!   The current row of f has now been completed. Store it back in its
			!   packed form and loop to the next result row

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc22

	subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmax:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstart(frow)
			fend = mmend(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,ll(midx,frow)) =   c1 * f(:,ll(midx,frow))   + c2 * 2.0 * fr(:,midx)
					f(:,ll(-midx,-frow)) = c1 * f(:,ll(-midx,-frow)) - c2 * 2.0 * fi(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,ll(midx,frow))   =  c1 * f(:,ll(midx,frow))   - c2 * 2.0 * fi(:,midx)
					f(:,ll(-midx,-frow)) =  c1 * f(:,ll(-midx,-frow)) + c2 * 2.0 * fr(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,ll(0,0)) = c1 * f(:,ll(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmend(0)
				if (ftype > 0) then
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0))  + c2 * 2.0 * fr(:,midx)
					f(:,ll(-midx,0)) = c1 * f(:,ll(-midx,0)) - c2 * 2.0 * fi(:,midx)
				else
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0))  - c2 * 2.0 * fi(:,midx)
					f(:,ll(-midx,0)) = c1 * f(:,ll(-midx,0)) + c2 * 2.0 * fr(:,midx)
				end if
			end do

		end if

	end subroutine cpx2cs

	subroutine cs2cpx (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!      	[fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!      	[fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!        [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!        [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,lr,li

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

		!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		if (abs(ftype) == 1) then
		
			do frow = 1,nmax

				mptr = 0
				do midx = mmstart(frow),mmend(frow)
					if (ftype == 1) then
						lr = ll(midx,frow)
						li = ll(-midx,-frow)
					else
						lr = ll(-midx,-frow)
						li = ll(midx,frow)
					end if

					fr(:,mptr,frow) =  0.5*f(:,lr)
					fi(:,mptr,frow) = -0.5*f(:,li)

					mptr = mptr + 1
				end do

			end do
			
		end if

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstart(0), mmend(0)
			if ((ftype > 0 .and. midx >= 0) .or. (ftype < 0 .and. midx < 0)) then
				lr = ll(midx,0)
				li = ll(-midx,0)
			else
				lr = ll(-midx,0)
				li = ll(midx,0)
			end if

			if (midx > 0) then
				fr(:,mptr,0) =  0.5*f(:,lr)
				fi(:,mptr,0) = -0.5*f(:,li)
			else if (midx < 0) then
				fr(:,mptr,0) = 0.5*f(:,lr)
				fi(:,mptr,0) = 0.5*f(:,li)
			else
				fr(:,mptr,0) = facfl0 * f(:,lr)
				fi(:,mptr,0) = 0.0
		  	end if
			mptr = mptr + 1
		end do

	end subroutine cs2cpx

	subroutine multr11 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.

		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmax, nmax

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmax, nmax

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmax, nmax

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmax, nmax

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr11

	subroutine multr12 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.

		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart,fend
		integer :: h1idx,g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr12

	subroutine multr22 (f, g, h, gtype, htype, c1, c2)

	!   form the convolution of the two functions g and h, and store it in f.
	!   The logic is based on the complex exponential form of the functions,
	!   for which the convolution is

	!        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

	!   where m', n' range over all values for which both (m-m',n-n') and (m',n')
	!   are in the range of the data.

	!   The input data, however is in a cos-sin rather than complex exponential
	!   format, so it is converted to the other format as needed. The basic
	!   algorithm is: for each row (constant n) of exponential form output,
	!   find each row pair (n-n' and n') that contributes to that output row,
	!   then build that row-pair from the input functions and perform all
	!   calculations in the complex domain. When all row pairs have been
	!   considered the output row is done, and the results are mapped back into
	!   the cos-sin form in the result.

	!   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
	!   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
	!   For type -1 functions this is reversed.

	!   Since the functions represented are real, the complex exponential
	!   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
	!   So only positive values of n need be computed. (negative n rows of
	!   the input functions still contribute to the convolution, however.)

	!   This is a banded convolution, since the input functions are banded.
	!   The arrays mmstart and mmend give, for each row n in the complex
	!   representation, the lower and upper bounds of m values for which
	!   data is defined. Outside these bounds the data is assumed to be zero.

		use param
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart,fend
		integer :: h1idx,g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstart (grow)
					gend = mmend (grow)
					hstart = mmstart (hrow)
					hend = mmend (hrow)

					if (gend < gstart .or. hend < hstart) cycle
					if (fend-gstart < hstart .or. fstart-gend > hend) cycle

	!   Compute all contributions to the current frow from the current grow-hrow
	!   pair

					do midx = fstart,fend

	!   calculate the loop limits for this particular midx. The ranges are
	!   [midx-gend,midx-gstart] and [hstart,hend]

						mstart = hstart
						if (midx-gend > hstart) mstart = midx-gend
						mend = hend
						if (midx-gstart < hend) mend = midx-gstart

						if (mend < mstart) cycle

	!   Do some convolving! Since only the positive bands of the complex
	!   arrays are stored, negative band values are extracted from corresponding
	!   positive bands. Also, the bands are not stored directly at their indices,
	!   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

	!   It is possible that the following loops may be speeded up somewhat by
	!   moving the address calculations into the loop limits: i.e. let hidx
	!   be the do loop index directly.

						if (grow >= 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr22

	subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmaxx:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstart(frow)
			fend = mmend(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,ll(midx,frow)) =   c1 * f(:,ll(midx,frow)) + c2 * 2.0 * fr(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,ll(midx,frow))  =  c1 * f(:,ll(midx,frow)) - c2 * 2.0 * fi(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,ll(0,0)) = c1 * f(:,ll(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmend(0)
				if (ftype > 0) then
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0)) + c2 * 2.0 * fr(:,midx)
				else
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0)) - c2 * 2.0 * fi(:,midx)
				end if
			end do

		end if

	end subroutine rli2cs

	subroutine cs2rli (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!         [fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!         [fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!         [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!         [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,l

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

	!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		if (abs(ftype) == 1) then

			do frow = 1,nmax

				mptr = 0
				do midx = mmstart(frow),mmend(frow)
					l = ll(midx,frow)

					if (ftype == 1) then
						fr(:,mptr,frow) =  0.5*f(:,l)
						fi(:,mptr,frow) =  0.0
					else
						fr(:,mptr,frow) =  0.0
						fi(:,mptr,frow) = -0.5*f(:,l)
					end if

					mptr = mptr + 1
				end do

			end do
			
		end if

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstart(0), mmend(0)
			if (midx >= 0) then
				l = ll(midx,0)
			else
				l = ll(-midx,0)
			endif

			if (ftype > 0) then
				if (midx > 0) then
					fr(:,mptr,0) =  0.5*f(:,l)
					fi(:,mptr,0) =  0.0
				else if (midx < 0) then
					fr(:,mptr,0) = 0.5*f(:,l)
					fi(:,mptr,0) = 0.0
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			else
				if (midx > 0) then
					fr(:,mptr,0) =  0.0
					fi(:,mptr,0) = -0.5*f(:,l)
				else if (midx < 0) then
					fr(:,mptr,0) = 0.0
					fi(:,mptr,0) = 0.5*f(:,l)
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			end if
			mptr = mptr + 1
		end do

	end subroutine cs2rli

	subroutine eqsplns(yeq,yeqr,yfar,m,m2,choice,smoo)

		use param
		use domain
		implicit none

		integer :: m,m2,mjeqp,mjeqm,ier,j,jj
		character(len=6) :: choice
		real(IDP) :: smoo,fac1,fac2,dspl,rjm,rjmp,fac,facp
		real(IDP), dimension(0:) :: yeq,yeqr
		real(IDP), dimension(0:) :: yfar
		real(IDP), dimension(:), allocatable :: xspl,dfspl,fnspl,yspl
		real(IDP), dimension(:,:), allocatable :: cspl
		real(IDP), dimension(:,:), allocatable :: wkspl

		interface
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
			subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
				use param
				implicit none
				integer :: nx,ic,ier
				real(IDP) :: sm
				real(IDP), dimension(:) :: x,f,df,y
				real(IDP), dimension(:,:) :: c
				real(IDP), dimension(:,:) :: wk
			end subroutine icsscu
		end interface

!		Calculates splines for the equilibrium variables in vmec		
	
		allocate (xspl(0:mjeq))
		allocate (dfspl(0:mjeq))
		allocate (fnspl(0:mjeq))
		allocate (yspl(0:mjeq))
		allocate (cspl(0:mjeq,3))
		allocate (wkspl(7,mjeq+3))

		mjeqp=mjeq+1
		mjeqm=mjeq-1
		xspl=rfar**2
		fac1=xspl(2)/(xspl(2)-xspl(1))
		fac2=-xspl(1)/(xspl(2)-xspl(1))
		dfspl(1:mjeq)=1./rfar(1:mjeq)**m2
		fnspl(1:mjeq)=yfar(1:mjeq)/rfar(1:mjeq)**m
		dfspl(0)=dfspl(1)
		fnspl(0)=fac1*fnspl(1)+fac2*fnspl(2)
		if (m == 0) then
			dfspl(0)=1.0_IDP
			fnspl(0)=yfar(0)
		end if
		yspl=fnspl
		ier=0

		if (choice == "spline") call spline(mjeqp,xspl,yspl,cspl(:,1),cspl(:,2),cspl(:,3))
		if (choice == "icsscu") call icsscu(xspl,fnspl,dfspl,mjeqp,smoo,yspl,cspl,mjeqp,wkspl,ier)

		do j=0,mj_r
			do jj=1,mjeqm
				if (rfar(jj) > r_red(j)) exit
			end do
			jj=jj-1
			dspl=r_red(j)**2-xspl(jj)
			rjm=1.0_IDP
			if (m /= 0) rjm=0.
			if (r_red(j) /= 0.) rjm=r_red(j)**m
			rjmp=0.
			if (m == 1) rjmp=1.0_IDP
			if (m /= 0 .and. r_red(j) /= 0.) rjmp=m*r_red(j)**(m-1)
			fac=((cspl(jj,3)*dspl+cspl(jj,2))*dspl+cspl(jj,1))*dspl+yspl(jj)
			facp=r_red(j)*((6.*cspl(jj,3)*dspl+4.*cspl(jj,2))*dspl+2.*cspl(jj,1))
			yeq(j)=rjm*fac
			yeqr(j)=rjmp*fac+rjm*facp
		end do

		deallocate (wkspl)
		deallocate (cspl)
		deallocate (yspl)
		deallocate (fnspl)
		deallocate (dfspl)
		deallocate (xspl)

	end subroutine eqsplns

	subroutine spline (n, x, y, b, c, d)

		use param
		implicit none

		integer :: n
		real(IDP), dimension(:) :: x, y, b, c, d

!  	the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  	for a cubic interpolating spline

!  	s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

!  	for  x(i) .le. x .le. x(i+1)

!  	input..

!  	n = the number of data points or knots (n.ge.2)
!  	x = the abscissas of the knots in strictly increasing order
!  	y = the ordinates of the knots

!  	output..

!  	b, c, d  = arrays of spline coefficients as defined above.

!  	using  p  to denote differentiation,

!  	y(i) = s(x(i))
!  	b(i) = sp(x(i))
!  	c(i) = spp(x(i))/2
!  	d(i) = sppp(x(i))/6  (derivative from the right)

!  	the accompanying function subprogram  seval  can be used
!  	to evaluate the spline.


		integer :: nm1, ib, i
		real(IDP) :: t

		nm1 = n-1

		if ( n > 2 ) then

!  	set up tridiagonal system

!  	b = diagonal, d = offdiagonal, c = right hand side.

			d(1) = x(2) - x(1)
			c(2) = (y(2) - y(1))/d(1)
			do i = 2, nm1
				d(i) = x(i+1) - x(i)
				b(i) = 2.*(d(i-1) + d(i))
				c(i+1) = (y(i+1) - y(i))/d(i)
				c(i) = c(i+1) - c(i)
			end do

!  	end conditions.  third derivatives at  x(1)  and  x(n)
!  	obtained from divided differences

			b(1) = -d(1)
			b(n) = -d(n-1)
			c(1) = 0.
			c(n) = 0.
			if ( n > 3 ) then
				c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
				c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
				c(1) = c(1)*d(1)**2/(x(4)-x(1))
				c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
			end if

!  	forward elimination

			do i = 2, n
				t = d(i-1)/b(i-1)
				b(i) = b(i) - t*d(i-1)
				c(i) = c(i) - t*c(i-1)
			end do

!  	back substitution

			c(n) = c(n)/b(n)
			do ib = 1, nm1
				i = n-ib
				c(i) = (c(i) - d(i)*c(i+1))/b(i)
			end do

!  	c(i) is now the sigma(i) of the text

!  	compute polynomial coefficients

			b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
			do i = 1, nm1
				b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
				d(i) = (c(i+1) - c(i))/d(i)
				c(i) = 3.*c(i)
			end do
			c(n) = 3.*c(n)
			d(n) = d(n-1)

		else if (n == 2) then

			b(1) = (y(2)-y(1))/(x(2)-x(1))
			c(1) = 0.
			d(1) = 0.
			b(2) = b(1)
			c(2) = 0.
			d(2) = 0.

		end if

	end subroutine spline

	subroutine root(t,ft,b,c,relerr,abserr,iflag)

		use param
		implicit none

!		root computes a root of the nonlinear equation f(x)=0
!		where f(x) is a continuous real function of a single real
!		variable x.  the method used is a combination of bisection
!		and the secant rule.

!		normal input consists of a continuous function f and an
!		interval (b,c) such that f(b)*f(c).le.0.0.  each iteration
!		finds new values of b and c such that the interval (b,c) is
!		shrunk and f(b)*f(c).le.0.0.  the stopping criterion is

!				abs(b-c).le.2.0*(relerr*abs(b)+abserr)

!		where relerr=relative error and abserr=absolute error are
!		input quantities.  set the flag, iflag, positive to initialize
!		the computation.  as b,c and iflag are used for both input and
!		output, they must be variables in the calling program.

!		if 0 is a possible root, one should not choose abserr=0.0.

!		the output value of b is the better approximation to a root
!		as b and c are always redefined so that abs(f(b)).le.abs(f(c)).

!		to solve the equation, root must evaluate f(x) repeatedly. this
!		is done in the calling program.  when an evaluation of f is
!		needed at t, root returns with iflag negative.  evaluate ft=f(t)
!		and call root again.  do not alter iflag.

!		when the computation is complete, root returns to the calling
!		program with iflag positive:

!		 iflag=1  if f(b)*f(c).lt.0 and the stopping criterion is met.

!				=2  if a value b is found such that the computed value
!				    f(b) is exactly zero.  the interval (b,c) may not
!				    satisfy the stopping criterion.

!				=3  if abs(f(b)) exceeds the input values abs(f(b)),
!				    abs(f(c)).   in this case it is likely that b is close
!				    to a pole of f.

!				=4  if no odd order root was found in the interval.  a
!				    local minimum may have been obtained.

!				=5  if too many function evaluations were made.
!				    (as programmed, 500 are allowed.)

!		this code is a modification of the code  zeroin  which is completely
!		explained and documented in the text,  numerical computing:  an
!		introduction  by l. f. shampine and r. c. allen.

		integer :: iflag
		integer, save :: ic,kount
		real(IDP) :: t,ft,b,c,relerr,abserr,fb,cmb,acmb,tol,p,q
		real(IDP), save :: re,ae,a,fa,fc,fx,u,acbs
		integer, save :: mentry=1


		if (mentry == 1) then

!		if first entry then compute u

			u=epsilon(t)
			mentry = 2
		end if

		if (iflag < -1) then
			fb=ft
			if (iflag == -2) then
				fc=fa
				kount=2
				fx=max(abs(fb),abs(fc))
			else if (iflag == -3 .and. fb /= 0.0_IDP) then
				kount=kount+1
				if (sign(1.0_IDP,fb) == sign(1.0_IDP,fc)) then
					c=a
					fc=fa
				end if
			else
				iflag=2
				return
			end if
			if (abs(fc) < abs(fb)) then

!		interchange b and c so that abs(f(b)).le.abs(f(c)).

				a=b
				fa=fb
				b=c
				fb=fc
				c=a
				fc=fa
			end if
			cmb=0.5*(c-b)
			acmb=abs(cmb)
			tol=re*abs(b)+ae

!		test stopping criterion and function count.

			if (acmb <= tol) then
				if (sign(1.0_IDP,fb) == sign(1.0_IDP,fc)) then
					iflag=4
				else if (abs(fb) > fx) then
					iflag=3
				else
					iflag=1
				end if
				return
			end if
			if (kount >= 500) then
				iflag=5
				return
			end if

!		calculate new iterate implicitly as b+p/q
!		where we arrange p.ge.0.  the implicit
!		form is used to prevent overflow.

			p=(b-a)*fb
			q=fa-fb
			if (p < 0.0_IDP) then
				p=-p
				q=-q
			end if

!		update a, check if reduction in the size of bracketing
!		interval is satisfactory.  if not, bisect until it is.

			a=b
			fa=fb
			ic=ic+1
			if (ic > 3 .and. 8.0*acmb < acbs) then
				ic=0
				acbs=acmb
			end if
			if (ic < 4) then

!		test for too small a change.

				if (p > abs(q)*tol) then

!		root ought to be between b and (c+b)/2.

					if (p < cmb*q) then

!		use secant rule.

						b=b+p/q
					else

!		use bisection.

						b=0.5*(c+b)
					end if
					
				else

!		increment by tolerance.

					b=b+sign(tol,cmb)
				end if
			else

!		use bisection.

				b=0.5*(c+b)
			end if

!		have completed computation for new iterate b.

			t=b
			iflag=-3
		else if (iflag == -1) then
			fa=ft
			t=b
			iflag=-2
		else
			re=max(relerr,u)
			ae=max(abserr,0.0_IDP)
			ic=0
			acbs=abs(b-c)
			a=c
			t=a
			iflag=-1
		end if

	end subroutine root

	subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)

!  icsscu-------s-------library 2---------------------------------------

!  function            - cubic spline data smoothing
!  usage               - call icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
!  parameters   x      - vector of length nx containing the abscissae
!  							 of the nx data points (x(i),f(i)) i=1,...,
!  							 nx (input). x must be ordered so that
!  							 x(i) .lt. x(i+1).
!  				 f      - vector of length nx containing the ordinates
!  							 (or function values) of the nx data points
!  							 (input).
!  				 df     - vector of length nx (input).
!  							 df(i) is the relative weight of data
!  							 point i (see parameter sm below).
!  				 nx     - number of elements in x, f, df, and y (input).
!  							 nx must be .ge. 2.
!  				 sm     - a non-negative number which controls the
!  							 extent of smoothing (input). the spline
!  							 function s is determined such that the
!  							 sum from 1 to nx of
!  							 ((s(x(i))-f(i))/df(i))**2 .le. sm,
!  							 where equality holds unless s describes
!  							 a straight line.
!  				 y,c    - spline coefficients (output). y is a vector
!  							 of length nx. c is an nx-1 by 3 matrix.
!  							 the value of the spline approximation
!  							 at t is
!  							 s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!  							 where x(i) .le. t .lt. x(i+1) and
!  							 d = t-x(i).
!  				 ic     - row dimension of matrix c in the calling
!  							 program (input). ic must be .ge. nx-1.
!  				 wk     - work area vector of length greater than or
!  							 equal to 7*nx+14.
!  				 ier    - error parameter
!  							 terminal error
!  							 	ier = 129. ic is less than nx-1.
!  							 	ier = 130, nx is less than 2.
!  							 	ier = 131, input abscissae are not ordered
!  							 				  so that x(1) .lt. x(2) ... .lt. x(nx).
!  precision           - single
!  req'd imsl routines - uertst
!  language            - fortran

!  latest revision     - february 16, 1976
!  							 dec

		use param
		implicit none

		integer :: nx,ic,ier,m2,np1,i,np3,j
		real(IDP) :: sm,p,h,f2,ff,g,onedh,e,hmg
		real(IDP), dimension(:) :: x,f,df,y
		real(IDP), dimension(:,:) :: c
		real(IDP), dimension(:,:) :: wk

		interface
			subroutine uertst(ier,name)
				implicit none
				integer :: ier
				character(len=6) :: name
			end subroutine uertst
		end interface

!  	check error conditions
		ier = 0
		if (ic < nx-1) then
			ier = 129
			call uertst(ier,"icsscu")
      else if (nx < 2) then
			ier = 130
			call uertst(ier,"icsscu")
		else
!  	set up working areas
      	m2 = nx+2
      	np1 = nx+1
      	wk(1,1) = 0.0
      	wk(1,2) = 0.0
      	wk(2,np1) = 0.0
      	wk(3,m2) = 0.0
      	wk(3,np1) = 0.0
      	wk(6,1) = 0.0
      	wk(6,2) = 0.0
      	wk(6,m2) = 0.0
      	wk(6,np1) = 0.0
      	p = 0.0
      	h = x(2)-x(1)
      	if (h <= 0.0) then
				ier = 131
				call uertst(ier,"icsscu")
			else
				f2 = -sm
				ff = (f(2)-f(1))/h
				if (nx > 2) then
					do i=3,nx
						g = h
						h = x(i)-x(i-1)
      				if (h <= 0.0) then
							ier = 131
							call uertst(ier,"icsscu")
							return
						end if
						onedh = 1.0/h
						e = ff
						ff = (f(i)-f(i-1))*onedh
						y(i) = ff-e
						wk(4,i) = 2.*(g+h)/3.
						wk(5,i) = h/3.0
						wk(3,i) = df(i-2)/g
						wk(1,i) = df(i)*onedh
						wk(2,i) = -df(i-1)/g-df(i-1)*onedh
					end do
					do i=3,nx
						c(i-1,1) = wk(1,i)*wk(1,i)+wk(2,i)*wk(2,i)+wk(3,i)*wk(3,i)
						c(i-1,2) = wk(1,i)*wk(2,i+1)+wk(2,i)*wk(3,i+1)
						c(i-1,3) = wk(1,i)*wk(3,i+2)
					end do
				end if
!  	next iteration
				do while (.true.)
					if (nx > 2) then
						do i=3,nx
							wk(2,i-1) = ff*wk(1,i-1)
							wk(3,i-2) = g*wk(1,i-2)
							wk(1,i) = 1.0/(p*c(i-1,1)+wk(4,i)-ff*wk(2,i-1)-g*wk(3,i-2))
							wk(6,i) = y(i)-wk(2,i-1)*wk(6,i-1)-wk(3,i-2)*wk(6,i-2)
							ff = p*c(i-1,2)+wk(5,i)-h*wk(2,i-1)
							g = h
							h = c(i-1,3)*p
						end do
						np3 = nx+3
						do i=3,nx
							j = np3-i
							wk(6,j) = wk(1,j)*wk(6,j)-wk(2,j)*wk(6,j+1)-wk(3,j)*wk(6,j+2)
						end do
					end if
					e = 0.0
					h = 0.0
!  	compute u and accumulate e
					do i=2,nx
						g = h
						h = (wk(6,i+1)-wk(6,i))/(x(i)-x(i-1))
						hmg = h-g
						wk(7,i) = hmg*df(i-1)*df(i-1)
						e = e+wk(7,i)*hmg
					end do
					g = -h*df(nx)*df(nx)
					wk(7,np1) = g
					e = e-g*h
					g = f2
					f2 = e*p*p
					if (f2 >= sm .or. f2 <= g) exit
					ff = 0.0
					h = (wk(7,3)-wk(7,2))/(x(2)-x(1))
					if (nx > 2) then
						do i=3,nx
							g = h
							h = (wk(7,i+1)-wk(7,i))/(x(i)-x(i-1))
							g = h-g-wk(2,i-1)*wk(1,i-1)-wk(3,i-2)*wk(1,i-2)
							ff = ff+g*wk(1,i)*g
							wk(1,i) = g
						end do
					end if
					h = e-p*ff
					if (h <= 0.0) exit
!  	update the lagrange multiplier p for the next iteration
					p = p+(sm-f2)/((sqrt(sm/e)+p)*h)
				end do
!  	if e less than or equal to s, compute the coefficients and return.
				np1 = nx-1
				do i=1,np1
					y(i) = f(i)-p*wk(7,i+1)
					c(i,2) = wk(6,i+1)
					wk(1,i) = y(i)
				end do
				wk(1,nx) = f(nx)-p*wk(7,nx+1)
				y(nx) = wk(1,nx)
				do i=2,nx
					h = x(i)-x(i-1)
					c(i-1,3) = (wk(6,i+1)-c(i-1,2))/(h+h+h)
					c(i-1,1) = (wk(1,i)-y(i-1))/h-(h*c(i-1,3)+c(i-1,2))*h
				end do
			end if
		end if

	end subroutine icsscu

	subroutine uertst(ier,name)

!  uertst---------------library 2---------------------------------------

!  function            - error message generation
!  usage               - call uertst(ier,name)
!  parameters   ier    - error parameter. type + n  where
!  								 type= 128 implies terminal error
!  										  64 implies warning with fix
!  										  32 implies warning
!  								 n   = error code relevant to calling routine
!  				 name   - input scalar (double precision on dec)
!  							 containing the name of the calling routine
!  							 as a 6-character literal string.
!  language            - fortran

!  latest revision     - october 1,1975
!  							 dec

		implicit none
		
		integer :: ier,ier1,ier2,i
		character(len=6) :: name
		integer, dimension(4) :: ibit=(/32,64,128,0/)
		character(len=5), dimension(4,4) :: ityp
!  	character(len=5), dimension(4,4) :: ityp= (/"warni","ng   ","     ","     ","warni","ng(wi","th fi","x)   ", &
!  															  "termi","nal  ","     ","     ","non-d","efine","d    ","     "/)

		ityp(:,1)=(/"warni","ng   ","     ","     "/)
		ityp(:,2)=(/"warni","ng(wi","th fi","x)   "/)
		ityp(:,3)=(/"termi","nal  ","     ","     "/)
		ityp(:,4)=(/"non-d","efine","d    ","     "/)

 		ier2=ier
		if (ier2 < ibit(1)) then
!  	non-defined
			ier1=4
		else if (ier2 < ibit(2)) then
!  	warning
			ier1=1
		else if (ier2 < ibit(3)) then
!  	warning(with fix)
			ier1=2
		else
!  	terminal
			ier1=3
		end if
!  	extract 'n'
		ier2=ier2-ibit(ier1)
!  	print error message
      write (6,'(" *** i m s l(uertst) ***  ",4a5,2x,a6,2x,i2," (ier = ",i3,")")') (ityp(i,ier1),i=1,4),name,ier2,ier

	end subroutine uertst

	subroutine mmblims

	!   Calculate the minimum and maximum m values for each n value
	!   in the complex exponential representation.

		use param
		use domain
		implicit none

		integer :: n,l

		allocate (mmstartb(-nmaxb:nmaxb))
		allocate (mmendb(-nmaxb:nmaxb))

		do n = 0,nmaxb
			mmstartb(n) = 1000
			mmendb(n) = -1000
		end do

		do l = 1,lbmax

			n = nnb(l)
			if (n >= 0) then
				if (mmb(l) < mmstartb(n)) mmstartb(n) = mmb(l)
				if (mmb(l) > mmendb(n)) mmendb(n) = mmb(l)
			else
				if (-mmb(l) < mmstartb(-n)) mmstartb(-n) = -mmb(l)
				if (-mmb(l) > mmendb(-n)) mmendb(-n) = -mmb(l)
			end if

		end do

		!   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
		!   so limits for negative n rows are just negatives of corresponding
		!   positive n limits.

		do n = 0,nmaxb

			mmstartb(-n) = -mmendb(n)
			mmendb(-n) = -mmstartb(n)

		end do

	end subroutine mmblims

	subroutine splns(yeq,yeqr,yfar,m,m2,choice,smoo)

		use param
		use domain
		implicit none

		integer :: m,m2,mjp,mjm,ier,j,jj
		character(len=6) :: choice
		real(IDP) :: smoo,fac1,fac2,dspl,rjm,rjmp,fac,facp
		real(IDP), dimension(0:) :: yeq,yeqr
		real(IDP), dimension(0:) :: yfar
		real(IDP), dimension(:), allocatable :: xspl,dfspl,fnspl,yspl
		real(IDP), dimension(:,:), allocatable :: cspl
		real(IDP), dimension(:,:), allocatable :: wkspl

		interface
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
			subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
				use param
				implicit none
				integer :: nx,ic,ier
				real(IDP) :: sm
				real(IDP), dimension(:) :: x,f,df,y
				real(IDP), dimension(:,:) :: c
				real(IDP), dimension(:,:) :: wk
			end subroutine icsscu
		end interface

!		Calculates splines for the dynamic variables in var_sum		
		
		allocate (xspl(0:mj))
		allocate (dfspl(0:mj))
		allocate (fnspl(0:mj))
		allocate (yspl(0:mj))
		allocate (cspl(0:mj,3))
		allocate (wkspl(7,mj+3))

		mjp=mj+1
		mjm=mj-1
		xspl=r**2
		fac1=xspl(2)/(xspl(2)-xspl(1))
		fac2=-xspl(1)/(xspl(2)-xspl(1))
		dfspl(1:mj)=1./r(1:mj)**m2
		fnspl(1:mj)=yfar(1:mj)/r(1:mj)**m
		dfspl(0)=dfspl(1)
		fnspl(0)=fac1*fnspl(1)+fac2*fnspl(2)

		if (m == 0) then
			dfspl(0)=1.0_IDP
			fnspl(0)=yfar(0)
		end if
		yspl=fnspl
		ier=0

		if (choice == "spline") call spline(mjp,xspl,yspl,cspl(:,1),cspl(:,2),cspl(:,3))
		if (choice == "icsscu") call icsscu(xspl,fnspl,dfspl,mjp,smoo,yspl,cspl,mjp,wkspl,ier)

		do j=0,mj_r
			do jj=1,mjm
				if (r(jj) > r_red(j)) exit
			end do
			jj=jj-1
			dspl=r_red(j)**2-xspl(jj)
			rjm=1.0_IDP
			if (m /= 0) rjm=0.
			if (r_red(j) /= 0.) rjm=r_red(j)**m
			rjmp=0.
			if (m == 1) rjmp=1.0_IDP
			if (m /= 0 .and. r_red(j) /= 0.) rjmp=m*r_red(j)**(m-1)
			fac=((cspl(jj,3)*dspl+cspl(jj,2))*dspl+cspl(jj,1))*dspl+yspl(jj)
			facp=r(j)*((6.*cspl(jj,3)*dspl+4.*cspl(jj,2))*dspl+2.*cspl(jj,1))
			yeq(j)=rjm*fac
			yeqr(j)=rjmp*fac+rjm*facp
		end do

		deallocate (wkspl)
		deallocate (cspl)
		deallocate (yspl)
		deallocate (fnspl)
		deallocate (dfspl)
		deallocate (xspl)

	end subroutine splns

	subroutine grid

!  	this sub sets up the r grid
!  	and calculates the time independent geometric scale factors

		use param
		use domain
		implicit none

		integer :: mjp, j, jm, jp

		interface
			subroutine findr
			end subroutine findr
		end interface

		mjp=mj_r
		mjm1=mj_r-1
		mjm2=mj_r-2
		
!  	compute the r values

		call findr
		
	end subroutine grid

	subroutine findr

!  	this sub sets up the r grid

!  	requires
!  	nis==number of points in the island
!  	ni==number of points interior to the island
!  	ne==number of points exterior to the island
!  	mj==total number of points=ni+nis+ne
!  	delta==width of the uniform fine grid (island)
!  	rc==center of the fine grid (island)
!  	fti==fraction of interior point in transition region
!  	fte==fraction of exterior points in transition region

		use param
		use domain
		use findrf
		implicit none

		integer :: ni0,ne0,nis0,nset,nte,nti,maxj,j,jt,jj,minj
		real(IDP) :: fe,fi,de,di,dr,ss,b

		interface
			function zeroin(ax,bx,f,tol)
				use param
				implicit none
				real(IDP) :: ax,bx,tol,zeroin
				interface
					function f(x)
						use param
						implicit none
						real(IDP) :: f,x
					end function f
				end interface
			end function zeroin
			function findf(f)
				use param
				implicit none
				real(IDP) :: f,findf
			end function findf
		end interface

		ni =  (mj_r/2) - 1
		nis = (mj_r/4) + 1
		ne =  mj_r/4
		fe=0
		fi=0
		ni0=ni
		ne0=ne
		nis0=nis
		nset=0
		if (ni < 2) then
			nset=2
			ni=2
		end if
		if (ne < 3) then
			nset=2
			ne=3
		end if
		if(nset > 1) nis=mj_r-ne-ni
		if(nis < 2) then
			write(6,'(" grid: nis .lt. 2")')
			stop
		end if
		ne=ne-1
		nte=fte*ne
		nti=fti*ni
		if(nte < 1) nte=1
		if(nti < 1) nti=1
		dx=delta/(nis-1.)
		de=1.-(rc+0.5*delta)
		di=1.-delta-de
		if (de <= 0.) then
			write(6,'(" grid: de .le. 0")')
			stop
		end if
		if (di <= 0.) then
			write(6,'(" grid: di .le. 0")')
			stop
		end if
		n0=ne
		nt=nte
		d=de
		b=10.0_IDP
		if (nt > range(d)) b=10.0_IDP**(1.0_IDP*range(d)/nt)  
		fe=zeroin(0.0_IDP,b,findf,1.e-10_IDP)
		n0=ni
		nt=nti
		d=di
		b=10.0_IDP
		if (nt > range(d)) b=10.0_IDP**(1.0_IDP*range(d)/nt)
		fi=zeroin(0.0_IDP,b,findf,1.e-10_IDP)
		dr=dx*fi**nti
		ss=0.
		maxj=ni-nti+1  
		do j=1,maxj
			ss=ss+dr
			r_red(j)=ss
		end do
		do j=1,nti
			jt=nti-j+1
			dr=dx*fi**jt
			ss=ss+dr
			jj=maxj+j
			r_red(jj)=ss
		end do
		minj=ni+2
		maxj=ni+nis
		do j=minj,maxj
			ss=ss+dx
			r_red(j)=ss
		end do
		do j=1,nte
			dr=dx*fe**j
			ss=ss+dr
			jj=maxj+j
			r_red(jj)=ss
		end do
		minj=ni+nis+nte+1
		maxj=ni+nis+ne
		do j=minj,maxj
			ss=ss+dr
			r_red(j)=ss
		end do
		maxj=maxj+1
		if (mj_r /= maxj) then
			write(6,'(" grid: mj_r .ne. ni+ne+nis")')
			stop
		end if
		r_red(mj_r)=1.0_IDP
		ne=ne0
		ni=ni0
		nis=nis0
		r_red(0)=0.

	end subroutine findr

	function findf(f)

		use param
		use findrf
		implicit none
		real(IDP) :: f,findf,xn,xndx,fnt,sum

		xn=n0-nt+1
		xndx=xn*dx
		fnt=f**nt
		if (f /= 1.) sum=f*(1.-fnt)/(1.-f)
		if (f == 1.) sum=nt
		findf=dx*sum+xndx*fnt-d

	end function findf

	function zeroin(ax,bx,f,tol)

		use param
		implicit none
		
		interface
			function f(x)
				use param
				implicit none
				real(IDP) :: f,x
			end function f
		end interface
		
		real(IDP) :: ax,bx,tol,zeroin

!  	a zero of the function  f(x)  is computed in the interval ax,bx .

!  input..

!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!  		 the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!  		 final result ( .ge. 0.0)


!  output..

!  zeroin abcissa approximating a zero of  f  in the interval ax,bx


!  	it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!  	this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).


		real(IDP) :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s

!  compute eps, the relative machine precision

		eps = 1.0
		tol1 = 2.0
		do while (tol1 > 1.0)
			eps = eps/2.0
			tol1 = 1.0 + eps
		end do

!  initialization

		a = ax
		b = bx
		fa = f(a)
		fb = f(b)

!  begin step

		c = a
		fc = fa
		d = b - a
		e = d
		do while (.true.)
			if (abs(fc) < abs(fb)) then
      		a = b
      		b = c
      		c = a
      		fa = fb
      		fb = fc
      		fc = fa
			end if

!  convergence test

			tol1 = 2.0*eps*abs(b) + 0.5*tol
			xm = .5*(c - b)
			if (abs(xm) <= tol1) exit
			if (fb == 0.0) exit

!  is bisection necessary

      	if (abs(e) < tol1 .or. abs(fa) <= abs(fb)) then

!  bisection

				d = xm
				e = d
				
			else

!  is quadratic interpolation possible

				if (a == c) then

!  linear interpolation

					s = fb/fa
					p = 2.0*xm*s
					q = 1.0 - s

				else

!  inverse quadratic interpolation

					q = fa/fc
					r = fb/fc
					s = fb/fa
					p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
					q = (q - 1.0)*(r - 1.0)*(s - 1.0)

				end if

!  adjust signs

				if (p > 0.0) q = -q
				p = abs(p)

!  is interpolation acceptable

      		if ((2.0*p) >= (3.0*xm*q - abs(tol1*q)) .or. p >= abs(0.5*e*q)) then
					d = xm
					e = d
				else
					e = d
      			d = p/q
				end if
				
			end if

!  complete step

			a = b
			fa = fb
			if (abs(d) > tol1) then
				b = b + d
			else
				b = b + sign(tol1, xm)
			end if
			fb = f(b)
			if ((fb*(fc/abs(fc))) > 0.0) then
				c = a
				fc = fa
				d = b - a
				e = d
			end if
		end do

!  done

		zeroin = b

	end function zeroin

	subroutine fitter(yeq,yeqr,yfar,ir,m,mb)

		use param
		use domain
		implicit none

		integer :: ir,m,mb,ior
		real(IDP), dimension(0:) :: yeq,yeqr,yfar
		integer :: i,i1,j,jj,js,je,l,nr,irp1,kplus1,nrows,mf,lyf,np1,lwrk,liwrk,ifail
		real(IDP) :: pmin,pmax,ave,x,aux,dspl
		integer, dimension(1) :: ip
		integer, dimension(30) :: iwrk
		real(IDP), dimension(ir+1) :: flux2,w,temp
		real(IDP), dimension(mjeq-ir+3) :: xspl,yspl
		real(IDP), dimension(mjeq-ir+3,3) :: cspl
		real(IDP), dimension(1) :: xf
		real(IDP), dimension(3) :: yf
		real(IDP), dimension(5,5) :: a
		real(IDP), dimension(5) :: s
		real(IDP), dimension(0:4) :: xhfit
		real(IDP), dimension(600) :: wrk

		interface
			subroutine e02age(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, kplus1, liwrk, lwrk, lyf, m, mf, np1, nrows
				real(IDP), dimension(nrows,kplus1) :: a
				real(IDP), dimension(kplus1) :: s
				real(IDP), dimension(m) :: w, x, y
				real(IDP), dimension(lwrk) :: wrk
				real(IDP), dimension(mf) :: xf
				real(IDP), dimension(lyf) :: yf
				integer, dimension(mf) :: ip
				integer, dimension(liwrk) :: iwrk
			end subroutine e02age
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
		end interface

		ior=(3*ir)/10
		if (mb > m) ior=mb*ior/m
		irp1=ir+1
		flux2(1)=0.0_IDP
		do i=1,ir
			flux2(i+1)=rfar(i)*rfar(i)
		end do
!		write(*,'("irp1=",i3)') irp1
!		write(*,'("flux2=",1p10e13.6)') (flux2(i),i=1,irp1)

		do i=1,ir
			temp(i+1)=yfar(i)/rfar(i)**m
		end do
		temp(1)=5.*temp(5)-4.*temp(6)
!		write(*,'("temp=",1p10e13.6)') (temp(i),i=1,irp1)
		ave = 0.0_IDP
	 	do i=1,irp1
			ave = ave + abs(temp(i))
		end do
		ave = ave/irp1
		do i=1,irp1
			w(i) = 10.0_IDP
			if (abs(temp(i)) > 10.*ave) w(i) = ave/abs(temp(i))
		end do
		do i=1,ior
			w(i)=1.e-10_IDP
		end do
!		write(*,'("w=",1p10e13.6)') (w(i),i=1,irp1)
		pmin   = flux2(1)
		pmax   = flux2(irp1)
		xf(1)  = flux2(irp1)
		yf(1)  = temp(irp1)
		aux    = yfar(irp1)/rfar(irp1)**m
		yf(2)  = (aux-temp(ir))/(2.*flux2(2))
		yf(3)  = (aux-2.*temp(irp1)+temp(ir))/(flux2(2)*flux2(2))
!		write(*,'("yf=",1p10e13.6)') (yf(i),i=1,3)
		w(irp1)= 0.0_IDP
		kplus1 = 5
		nrows  = 5
		lyf    = 3
		mf     = 1
		ip(1)  = 2
		np1    = 4
		lwrk   = 600
		liwrk  = 30
		ifail  = 0

		call e02age(ir,kplus1,nrows,pmin,pmax,flux2,temp,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

		xhfit(0) = 0.5*a(5,1) - a(5,3) + a(5,5)
		xhfit(1) = a(5,2) - 3.*a(5,4)
		xhfit(2) = 2.*a(5,3) - 8.*a(5,5)
		xhfit(3) = 4.*a(5,4)
		xhfit(4) = 8.*a(5,5)

		do j=0,mj_r
			if (r_red(j) > rfar(ir)) exit
			x = ((r_red(j)*r_red(j)-pmin)-(pmax-r_red(j)*r_red(j)))/(pmax-pmin)
			yeq(j) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
			yeqr(j) = ((4.*xhfit(4)*x+3.*xhfit(3))*x+2.*xhfit(2))*x+xhfit(1)
			yeqr(j) = 2.0*yeqr(j)/(pmax-pmin)
		end do

		je=j-1
!		write(*,'("r_red=",1p10e13.6)') (r_red(j),j=0,je)
!		write(*,'("yeq=",1p10e13.6)') (yeq(j),j=0,je)
		xspl(1)=rfar(ir-2)
		x = ((xspl(1)*xspl(1)-pmin)-(pmax-xspl(1)*xspl(1)))/(pmax-pmin)
		yspl(1) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
		xspl(2)=rfar(ir-1)
		x = ((xspl(2)*xspl(2)-pmin)-(pmax-xspl(2)*xspl(2)))/(pmax-pmin)
		yspl(2) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
		if (m == 0) then
			do j=0,je
				yeqr(j)=2.*r_red(j)*yeqr(j)
			end do
		else if (m == 1) then
			do j=0,je
				yeqr(j)=yeq(j)+2.*r_red(j)*r_red(j)*yeqr(j)
				yeq(j)=r_red(j)*yeq(j)
			end do
			yspl(1)=xspl(1)*yspl(1)
			yspl(2)=xspl(2)*yspl(2)
		else
			do j=0,je
				yeqr(j)=m*yeq(j)*r_red(j)**(m-1)+2.*yeqr(j)*r_red(j)**(m+1)
				yeq(j)=yeq(j)*r_red(j)**m
			end do
			yspl(1)=yspl(1)*xspl(1)**m
			yspl(2)=yspl(2)*xspl(2)**m
		end if

	 	do i=ir,mjeq
			xspl(i-ir+3)=rfar(i)
			yspl(i-ir+3)=yfar(i)
		end do
		nr=mjeq-ir+3
!		write(*,'("xspl=",1p10e13.6)') (xspl(i),i=1,nr)
!		write(*,'("yspl=",1p10e13.6)') (yspl(i),i=1,nr)
		call spline(nr,xspl,yspl,cspl(:,1),cspl(:,2),cspl(:,3))

		js=3
		do j=je+1,mj_r
			do jj=js,nr
				if (xspl(jj) > r_red(j)) exit
			end do
			js=jj
			jj=jj-1
			dspl=r_red(j)-xspl(jj)
			yeq(j)=((cspl(jj,3)*dspl+cspl(jj,2))*dspl+cspl(jj,1))*dspl+yspl(jj)
			yeqr(j)=(3.*cspl(jj,3)*dspl+2.*cspl(jj,2))*dspl+cspl(jj,1)
		end do
!		write(*,'("r_red=",1p10e13.6)') (r_red(j),j=je+1,mj)
!		write(*,'("yeq=",1p10e13.6)') (yeq(j),j=je+1,mj)

	end subroutine fitter	

	subroutine e02ade(m,kplus1,nrows,x,y,w,work1,work2,a,s,ifail)

		use param
		implicit none

!  	nag library subroutine  e02ade

!  	e02ade  computes weighted least-squares polynomial
!  	approximations to an arbitrary set of data points.

!  	forsythe-clenshaw method with modifications due to
!  	reinsch and gentleman.

!  	uses nag library routine  p01aae.
!  	uses basic external function  sqrt.

!  	started - 1973.
!  	completed - 1976.
!  	author - mgc and jgh.

!  	work1  and  work2  are workspace areas.
!  	work1(1, r)  contains the value of the  r th  weighted
!  	residual for the current degree  i.
!  	work1(2, r)  contains the value of  x(r)  transformed
!  	to the range  -1  to  +1.
!  	work1(3, r)  contains the weighted value of the current
!  	orthogonal polynomial (of degree  i)  at the  r th
!  	data point.
!  	work2(1, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the current orthogonal polynomial
!  	(of degree  i).
!  	work2(2, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the previous orthogonal polynomial
!  	(of degree  i - 1).

!  	nag copyright 1975
!  	mark 5 release
!  	mark 6 revised  ier-84
!  	mark 11.5(f77) revised. (sept 1985.)

!  	check that the values of  m  and  kplus1  are reasonable


!  	.. parameters ..
		character(len=6), parameter :: srname="e02ade"
!  	.. scalar arguments ..
		integer :: ifail, kplus1, m, nrows
!  	.. array arguments ..
		real(IDP), dimension(nrows,kplus1) :: a
		real(IDP), dimension(kplus1) :: s
		real(IDP), dimension(m) :: w, x, y
		real(IDP), dimension(3,m) :: work1
		real(IDP), dimension(2,kplus1) :: work2
!  	.. local scalars ..
		real(IDP) :: alpip1, betai, bj, bjp1, bjp2, ci, d, df, di, dim1, dj, epsr, factor, pij, sigmai, wrpr, wrprsq, x1, xcapr, xm
		integer :: i, ierror, iplus1, iplus2, j, jplus1, jplus2, jrev, k, mdist, r
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic sqrt
!  	.. executable statements ..
		if (kplus1 < 1 .or. m < kplus1) then
			ierror = 4
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		k = kplus1 - 1

!  	test the validity of the data.

!  	check that the weights are strictly positive.

		do r = 1, m
			if (w(r) <= 0.0_IDP) then
				ierror = 1
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
		end do

!  	check that the values of  x(r)  are non-decreasing and
!  	determine
!  	the number  (mdist)  of distinct values of  x(r).

		mdist = 1
		do r = 2, m
			if (x(r) < x(r-1)) then
				ierror = 2
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			if (x(r) == x(r-1)) cycle
			mdist = mdist + 1
		end do

!  	if the  x(r)  all have the same value, i.e.  mdist = 1,
!  	the normalization of the independent variable is not
!  	possible.

		if (mdist == 1) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

!  	if the number of distinct values of  x(r)  fails to exceed
!  	the maximum degree  k  there is no unique polynomial
!  	approximation of that degree.

		if (mdist <= k) then
			ierror = 4
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

!  	check that  nrows  has been set sufficiently large.

		if (nrows < kplus1) then
			ierror = 5
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

		x1 = x(1)
		xm = x(m)
		d = xm - x1

!  	the initial values  work1(1, r)  (r = 1, 2, ..., m)  of the
!  	weighted residuals and the values  work1(2, r)  (r = 1, 2,
!  	..., m)
!  	of the normalized independent variable are computed.  note
!  	that
!  	work1(2, r)  is computed from the expression below rather
!  	than the
!  	more natural form

!  	(2.0*x(r) - x1 - xm)/d,

!  	since the former guarantees the computed value to differ from
!  	the true value by at most  4.0*machine accuracy,  whereas the
!  	latter has no such guarantee.

		do r = 1, m
			work1(1,r) = w(r)*y(r)
			work1(2,r) = ((x(r)-x1)-(xm-x(r)))/d
		end do
		i = 1
		betai = 0.0_IDP
		do iplus1 = 1, kplus1

!  		set starting values for degree  i.

			iplus2 = iplus1 + 1
			if (iplus1 < kplus1) then
				do jplus1 = iplus2, kplus1
					a(iplus1,jplus1) = 0.0_IDP
				end do
				work2(1,iplus2) = 0.0_IDP
				work2(2,iplus2) = 0.0_IDP
			end if
			alpip1 = 0.0_IDP
			ci = 0.0_IDP
			di = 0.0_IDP
			a(i,iplus1) = 0.0_IDP
			work2(1,iplus1) = 1.0_IDP
			if (kplus1 > 1) work2(2,1) = work2(1,2)
			do r = 1, m
				xcapr = work1(2,r)

!  			the weighted value  work1(3, r)  of the orthogonal
!  			polynomial of
!  			degree  i  at  x = x(r)  is computed by recurrence from its
!  			chebyshev-series representation.

				if (iplus1 == 1) then
					wrpr = w(r)*0.5*work2(1,1)
					work1(3,r) = wrpr
				else
					j = iplus2
					if (xcapr > 0.5) then

!  			reinsch*s modified recurrence.

						factor = 2.0*(1.0-xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							dj = work2(1,j) + dj - factor*bj
							bj = bj + dj
						end do
						wrpr = w(r)*(0.5*work2(1,1)+dj-0.5*factor*bj)
						work1(3,r) = wrpr
					else if (xcapr >= -0.5) then

!  			clenshaw*s original recurrence.

						factor = 2.0*xcapr
						bjp1 = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							bjp2 = bjp1
							bjp1 = bj
							bj = work2(1,j) - bjp2 + factor*bjp1
						end do
						wrpr = w(r)*(0.5*work2(1,1)-bjp1+0.5*factor*bj)
						work1(3,r) = wrpr
					else

!  			gentleman*s modified recurrence.

						factor = 2.0*(1.0+xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							dj = work2(1,j) - dj + factor*bj
							bj = dj - bj
						end do
						wrpr = w(r)*(0.5*work2(1,1)-dj+0.5*factor*bj)
						work1(3,r) = wrpr
					end if
				end if

!  			the coefficient  ci  of the  i th  orthogonal polynomial and
!  			the
!  			coefficients  alpip1  and  beta i  in the
!  			three-term recurrence relation for the orthogonal
!  			polynomials are computed.

				wrprsq = wrpr**2
				di = di + wrprsq
				ci = ci + wrpr*work1(1,r)
				alpip1 = alpip1 + wrprsq*xcapr
			end do
			ci = ci/di
			if (iplus1 /= 1) betai = di/dim1
			alpip1 = 2.0*alpip1/di

!  		the weighted residuals  work1(1, r)  (r = 1, 2, ..., m)  for
!  		degree  i  are computed, together with their sum of squares,
!  		sigmai.

			sigmai = 0.0
			do r = 1, m
				epsr = work1(1,r) - ci*work1(3,r)
				work1(1,r) = epsr
				sigmai = sigmai + epsr**2
			end do

!  		the root mean square residual  s(i + 1)  for degree  i  is
!  		theoretically undefined if  m = i + 1  (the condition for the
!  		polynomial to pass exactly through the data points).  should
!  		this
!  		case arise the r.m.s. residual is set to zero.

			if (iplus1 < m) then
				df = m - iplus1
				s(iplus1) = sqrt(sigmai/df)
			else
				s(iplus1) = 0.0_IDP
			end if

!  		the chebyshev coefficients  a(i + 1, 1), a(i + 1, 2), ...,
!  		a(i + 1, i + 1)  together with the coefficients
!  		work2(1, 1), work2(1, 2), ..., work2(1, i + 1),   in the
!  		chebyshev-series representation of the  i th  orthogonal
!  		polynomial are computed.

			do jplus1 = 1, iplus1
				jplus2 = jplus1 + 1
				pij = work2(1,jplus1)
				a(iplus1,jplus1) = a(i,jplus1) + ci*pij
				if (jplus2 > kplus1) exit
				work2(1,jplus1) = work2(1,jplus2) + work2(2,jplus1) - alpip1*pij - betai*work2(2,jplus2)
				work2(2,jplus2) = pij
			end do
			if (iplus1 < kplus1) then
				dim1 = di
				i = iplus1
			end if
		end do
		ifail = 0

	end subroutine e02ade

	subroutine e02afe(nplus1,f,a,ifail)

		use param
		implicit none

!  	nag library subroutine  e02afe

!  	e02afe  computes the coefficients of a polynomial,
!  	in its chebyshev-series form, which interpolates
!  	(passes exactly through) data at a special set of
!  	points.  least-squares polynomial approximations
!  	can also be obtained.

!  	clenshaw method with modifications due to reinsch
!  	and gentleman.

!  	uses nag library routines  p01aae  and  x01aae.
!  	uses basic external function  sin.

!  	started - 1973.
!  	completed - 1976.
!  	author - mgc and jgh.

!  	nag copyright 1975
!  	mark 5 release
!  	mark 5b revised  ier-73
!  	mark 11.5(f77) revised. (sept 1985.)


!  	.. parameters ..
		character(len=6), parameter :: srname="e02afe"
!  	.. scalar arguments ..
		integer :: ifail, nplus1
!  	.. array arguments ..
		real(IDP), dimension(nplus1) :: a, f
!  	.. local scalars ..
		real(IDP) :: bk, bkp1, bkp2, dk, f0, factor, fli, fln, halffn, piby2n
		integer :: i, ierror, iplus1, j, k, krev, n, n2, nless1
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec
!  	.. data statements ..
		real(IDP) :: pi=3.14159265358979323846_IDP

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic sin
!  	.. executable statements ..
		ierror = 0
		if (nplus1 < 2) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
		else if (nplus1 == 2) then
			a(1) = f(1) + f(2)
			a(2) = 0.5*(f(1)-f(2))
			ifail = 0
		else
			n = nplus1 - 1
			fln = n
			n2 = 2*n
			nless1 = n - 1
			piby2n = 0.5*pi/fln
			f0 = f(1)
			halffn = 0.5*f(nplus1)
			do iplus1 = 1, nplus1
				i = iplus1 - 1
				k = nplus1
				j = 3*i
				if (j > n2) then

!  			gentleman*s modified recurrence.

					fli = n - i
					factor = 4.0*(sin(piby2n*fli))**2
					dk = halffn
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						dk = f(k) - dk + factor*bk
						bk = dk - bk
					end do
					a(iplus1) = (f0-2.0*dk+factor*bk)/fln
				else if (j >= n) then

!  			clenshaw*s original recurrence.

					fli = n - 2*i
					factor = 2.0*sin(piby2n*fli)
					bkp1 = 0.0_IDP
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						bkp2 = bkp1
						bkp1 = bk
						bk = f(k) - bkp2 + factor*bkp1
					end do
					a(iplus1) = (f0-2.0*bkp1+factor*bk)/fln
				else

!  			reinsch*s modified recurrence.

					fli = i
					factor = 4.0*(sin(piby2n*fli))**2
					dk = halffn
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						dk = f(k) + dk - factor*bk
						bk = bk + dk
					end do
 					a(iplus1) = (f0+2.0*dk-factor*bk)/fln
				end if
			end do
			a(nplus1) = 0.5*a(nplus1)
			ifail = 0
		end if

	end subroutine e02afe

	subroutine e02age(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 9 revised. ier-314 (sep 1981).
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 14c revised. ier-878 (nov 1990).

!  	**************************************************************

!  	npl algorithms library routine confit

!  	created  20/8/1979      updated  16/5/80     release no. 00/05

!  	author... gerald t anthony.
!  	national physical laboratory,
!  	teddington, middlesex tw11 0lw,
!  	england.

!  	**************************************************************

!  	e02age calls aewe01 to determine a polynomial mu(x) which
!  	interpolates the given constraints and a polynomial nu(x)
!  	which has value zero where a constrained value is specified.
!  	it then calls adze02 to fit y-mu(x) as a polynomial in x
!  	with factor nu(x). finally the coefficients of mu are added
!  	to these fits to give the coefficients of the constrained
!  	fits to y. all polynomials are expressed in chebyshev
!  	series form

!  	input parameters
!  		m        the number of data points to be fitted
!  		kplus1   fits with up to kplus1 coefficients are required
!  		nrows    first dimension of array a where coefficients are
!  					to be stored
!  		xmin,    end points of the range of the
!  		xmax     independent variable
!  		x, y, w  arrays of data values of the independent variable,
!  					dependent variable and weight, respectively
!  		mf       number of x values at which a constraint
!  					is specified
!  		xf       array of values of the independent
!  					variable at which constraints are
!  					specified
!  		yf       array of specified values and derivatives of the
!  					dependent variable in the order
!  						y1, y1 deriv, y1 2nd deriv,...., y2,....
!  		lyf      dimension of array yf
!  		ip       integer array of degrees of derivatives
!  					specified at each point xf

!  	output parameters
!  		a        on exit, 2 parameter array containing the
!  					coefficients of the chebyshev series
!  					representation of the fits, a(i+1, j+1)
!  					contains the coefficient of tj in the fit
!  					of degree i, i = n,n+1,...,k,  j =
!  					0,1,...,i  where n = np1 - 1
!  		s        on exit, array containing the r.m.s. residual for
!  					each degree of fit from n to k
!  		np1      on exit, contains n + 1, where n is the
!  					total number of interpolation conditions

!  		ifail    failure indicator
!  					 0 - successful termination
!  					 1 - at least one of the following conditions
!  						  has been violated
!  						  lyf    at least n
!  						  lwrk   at least 2*n + 2 + the larger of
!  									4*m + 3*kplus1 and 8*np1 +
!  									5*imax + mf - 3 where imax =
!  									1 + max(ip(i))
!  						  liwrk  at least 2*mf + 2
!  						  kplus1 at least np1
!  						  m      at least 1
!  						  nrows  at least kplus1
!  						  mf     at least 1
!  					 2 - for some i, ip(i) is less than 0
!  					 3 - xmin is not strictly less than xmax
!  						  or for some i, xf(i) is not in range
!  						  xmin to xmax or the xf(i) are not
!  						  distinct
!  					 4 - for some i, x(i) is not in range xmin to xmax
!  					 5 - the x(i) are not non-decreasing
!  					 6 - the number of distinct values of x(i) with
!  						  non-zero weight is less than kplus1 - np1
!  					 7 - aewe01 has failed to converge, ie
!  						  the constraint cannot be satisfied
!  						  with sufficient accuracy

!  	workspace parameters
!  		wrk      real workspace array
!  		lwrk     dimension of wrk.   lwrk must be at least
!  					2*n + 2 + the larger of
!  					4*m + 3*kplus1 and 8*np1 + 5*imax + mf - 3
!  					where imax = 1 + max(ip(i))
!  		iwrk     integer workspace array
!  		liwrk    dimension of iwrk.   liwrk must be at least
!  					2*mf + 2

!  	.. parameters ..
		character(len=6), parameter :: srname="e02age"
!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, kplus1, liwrk, lwrk, lyf, m, mf, np1, nrows
!  	.. array arguments ..
		real(IDP), dimension(nrows,kplus1) :: a
		real(IDP), dimension(kplus1) :: s
		real(IDP), dimension(m) :: w, x, y
		real(IDP), dimension(lwrk) :: wrk
		real(IDP), dimension(mf) :: xf
		real(IDP), dimension(lyf) :: yf
		integer, dimension(mf) :: ip
		integer, dimension(liwrk) :: iwrk
!  	.. local scalars ..
		real(IDP) :: amuj, xi, xmu
		integer :: i, ierror, im1, imax, ipi, iymux, j, lw, mdist, n, nanu, neps, nser, nwrk, nwrk1, nwrk2
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
!  	.. external subroutines ..
			subroutine aewe01(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, itmax, itmin, liwrk, lwrk, m, n, np1
				real(IDP), dimension(n) :: a, y
				real(IDP), dimension(np1) :: b
				real(IDP), dimension(lwrk) :: wrk
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
				integer, dimension(liwrk) :: iwrk
			end subroutine aewe01
			subroutine adze02(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, inup1, kall, kplus1, mfirst, mlast, mtot, ndv, nrows
				real(IDP), dimension(ndv,nrows,kplus1) :: a
				real(IDP), dimension(ndv,mlast) :: eps, y
				real(IDP), dimension(inup1) :: nu
				real(IDP), dimension(ndv,kplus1) :: s
				real(IDP), dimension(kplus1) :: serr
				real(IDP), dimension(mlast) :: w, x
				real(IDP), dimension(2,mtot) :: work1
				real(IDP), dimension(2,kplus1) :: work2
			end subroutine adze02
			subroutine e02ake(np1,xmin,xmax,a,ia1,la,x,result,ifail)
				use param
				implicit none
				real(IDP) :: result, x, xmax, xmin
				integer :: ia1, ifail, la, np1
				real(IDP), dimension(la) :: a
			end subroutine e02ake
		end interface

!  	.. data statements ..
		real(IDP) :: one=1.0_IDP, zero=0.0_IDP
!  	.. executable statements ..
		if (mf < 1) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		imax = 0
		np1 = 1
		do i = 1, mf
			ipi = ip(i) + 1
			if (ipi < 1) then
				ierror = 2
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			if (ipi > imax) imax = ipi
			np1 = np1 + ipi
		end do
		n = np1 - 1
		if (lyf < n .or. liwrk < 2*mf+2) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		i = 4*m + 3*kplus1
		lw = 8*np1 + 5*imax + mf - 3
		if (lw < i) lw = i
		lw = lw + 2*np1
		if (lw > lwrk .or. np1 > kplus1 .or. m < 1) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		nanu = np1 + 1
		nwrk = nanu + np1
		if (xmax <= xmin) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (xf(1) > xmax .or. xf(1) < xmin) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (mf > 1) then
			do i = 2, mf
				xi = xf(i)
				if (xi > xmax .or. xi < xmin) then
					ierror = 3
					ifail = p01abe(ifail,ierror,srname,0,p01rec)
					return
				end if
				im1 = i - 1
				do j = 1, im1
					if (xi == xf(j)) then
						ierror = 3
						ifail = p01abe(ifail,ierror,srname,0,p01rec)
						return
					end if
				end do
			end do
		end if
		xmu = xmin
		if (x(1) == xmin) xmu = xmax
		mdist = 0
		do i = 1, m
			xi = x(i)
			if (xi == xmu .or. w(i) == zero) cycle
			do j = 1, mf
				if (xi == xf(j)) exit
			end do
			if (j > mf) mdist = mdist + 1
			xmu = xi
		end do
		if (mdist < kplus1-n) then
			ierror = 6
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		iwrk(1) = 1
		ierror = 1
		call aewe01(mf,xmin,xmax,xf,yf,ip,n,np1,5,20,wrk,wrk(nanu),wrk(nwrk),lwrk-nwrk+1,iwrk,liwrk,ierror)
		if (ierror /= 0) then
			ierror = 7
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		do i = 1, m
			ierror = 1
			call e02ake(n,xmin,xmax,wrk,1,np1,x(i),xmu,ierror)
			if (ierror /= 0) then
				ierror = 4
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			iymux = nwrk + i - 1

!  		store y - mu(x) at ith data point

			wrk(iymux) = y(i) - xmu
		end do
		nwrk1 = nwrk + m
		nwrk2 = nwrk1 + 2*m
		nser = nwrk2 + 2*kplus1
		neps = nser + kplus1
		ierror = 1
		call adze02(1,m,m,kplus1,nrows,1,1,x,wrk(nwrk),w,xmin,xmax,np1,wrk(nanu),wrk(nwrk1),wrk(nwrk2),a,s,wrk(nser),wrk(neps),ierror)
		if (ierror == 0) then
			do j = 1, n
				amuj = wrk(j)
				do  i = np1, kplus1
					a(i,j) = a(i,j) + amuj
				end do
			end do
		else
			ierror = ierror + 3
			if (ierror == 8) ierror = 1
		end if
		ifail = p01abe(ifail,ierror,srname,0,p01rec)

	end subroutine e02age

	subroutine e02ake(np1,xmin,xmax,a,ia1,la,x,result,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine tval1c

!  	created 19/3/79    updated 6/7/79    release no. 00/04.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1      np1 = n + 1. n is the degree of the
!  					chebyshev series
!  		xmin     minimum value of x
!  		xmax     maximum value of x
!  		a        the array where the coefficients are stored
!  		ia1      the address increment of a
!  		la       dimension of a
!  		x        unnormalized argument in the range (xmin, xmax)

!  	output parameters
!  		result   value of the summation
!  		ifail    error indicator

!  	np1 chebyshev coefficients a0, a1, ..., an, are
!  	stored in the array a in positions 1, 1+ia1, 1+2*ia1, ...,
!  	1+n*ia1, where n = np1 - 1.
!  	ia1 must not be negative.
!  	la must be at least equal to 1 + n*ia1.
!  	the value of the polynomial of degree n
!  	a0t0(xcap)/2 + a1t1(xcap) + a2t2(xcap) + + ... + antn(xcap),
!  	is calculated for the argument xcap, where xcap is
!  	the normalized value of x in the range (xmin, xmax),
!  	storing it in result.
!  	unless the routine detects an error, ifail contains
!  	zero on exit.
!  	ifail = 1 indicates at least one of the restrictions on
!  		input parameters is violated - ie
!  	np1  >  0
!  	ia1  >=  0
!  	la  >=  1 + n * ia1
!  	xmin  <  xmax
!  	ifail = 2 indicates that
!  	x does not satisfy the restriction xmin  <=  x  <=  xmax.
!  	the recurrence relation by clenshaw, modified by reinsch
!  	and gentleman, is used.

!  	.. parameters ..
		character(len=6), parameter :: srname="e02ake"
!  	.. scalar arguments ..
		real(IDP) :: result, x, xmax, xmin
		integer :: ia1, ifail, la, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
!  	.. local scalars ..
		real(IDP) :: xcap
		integer :: ierror
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
			subroutine akye02(xmin,xmax,x,xcap)
				use param
				implicit none
				real(IDP) :: x, xcap, xmax, xmin
			end subroutine akye02
			subroutine akze02(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine akze02
		end interface

!  	.. executable statements ..
		if (np1 < 1 .or. ia1 < 1 .or. la < 1+(np1-1)*ia1 .or. xmax <= xmin) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (x > xmax .or. x < xmin) then
			ierror = 2
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		ierror = 0
		call akye02(xmin,xmax,x,xcap)
		call akze02(np1,a,ia1,la,xcap,result)
		ifail = p01abe(ifail,ierror,srname,0,p01rec)

	end subroutine e02ake

	subroutine adze02(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)

		use param
		implicit none

!  	mark 7 release. nag copyright 1978.
!  	mark 8 revised. ier-228 (apr 1980).
!  	mark 11.5(f77) revised. (sept 1985.)

!  	adze02  computes weighted least-squares polynomial
!  	approximations to an arbitrary set of data points,
!  	with, if required, several sets of values of the
!  	dependent variable.

!  	forsythe-clenshaw method with modifications due to
!  	reinsch and gentleman.

!  	started - 1973.
!  	completed - 1978.
!  	authors - mgc and gta.

!  	work1  and  work2  are workspace areas.
!  	work1(1, r)  contains the value of  x(r)  transformed
!  	to the range  -1  to  +1.
!  	work1(2, r)  contains the weighted value of the current
!  	orthogonal polynomial (of degree  i)  at the  r th
!  	data point.
!  	work2(1, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the current orthogonal polynomial
!  	(of degree  i).
!  	work2(2, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the previous orthogonal polynomial
!  	(of degree  i - 1).

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, inup1, kall, kplus1, mfirst, mlast, mtot, ndv, nrows
!  	.. array arguments ..
		real(IDP), dimension(ndv,nrows,kplus1) :: a
		real(IDP), dimension(ndv,mlast) :: eps, y
		real(IDP), dimension(inup1) :: nu
		real(IDP), dimension(ndv,kplus1) :: s
		real(IDP), dimension(kplus1) :: serr
		real(IDP), dimension(mlast) :: w, x
		real(IDP), dimension(2,mtot) :: work1
		real(IDP), dimension(2,kplus1) :: work2
!  	.. local scalars ..
		real(IDP) :: alpip1, betai, bj, bjp1, bjp2, cil, d, df, di, dim1, dj, epslr, factor, pij, sigmai, wr, wrpr, wrprsq, x1, &
						 xcapr, xm
		integer :: i, ii, im1, inu, iplus1, iplus2, j, jplus1, jplus2, jrev, k, l, m, mdist, mr, r
		logical :: wnz
!  	.. local arrays ..
		real(IDP), dimension(10) :: ci
!  	.. intrinsic functions ..
!  	intrinsic sqrt
!  	.. executable statements ..
		k = kplus1 - 1
		inu = inup1 - 1

!  	test the validity of the data.

!  	check input parameters.

		m = mlast - mfirst + 1
		i = kplus1 - inu
		if (mfirst < 1 .or. inup1 < 1 .or. kplus1 < inup1 .or. m < i .or. ndv < 1 .or. (kall /= 1 .and. kall /= 0)) then
			ifail = 5
			return
		end if

!  	check that the values of x(r) are non-decreasing and
!  	determine the number (mdist) of distinct values of x(r)
!  	with non-zero weight

		mdist = 1
		if (w(mfirst) == 0.0_IDP) mdist = 0
		l = mfirst + 1
		if (l <= mlast) then
			wnz = w(mfirst)  /=  0.0_IDP
			do r = l, mlast
				if (x(r) < x(r-1)) then
					ifail = 2
					return
				end if
				if (x(r) > x(r-1)) wnz = .false.
				if (w(r) == 0.0_IDP .or. wnz) cycle
				mdist = mdist + 1
				wnz = .true.
			end do
		end if

!  	check that xmin < xmax and that xmin and xmax span the data
!  	x values.

		if (xmin > x(mfirst) .or. xmax < x(mlast) .or. xmin >= xmax) then
			ifail = 1
			return
		end if

!  	if the number of distinct values of  x(r)  with non-zero
!  	weight is less than the number of independent coefficients
!  	in the fit of maximum degree  k  there is no unique
!  	polynomial
!  	approximation of that degree.

		l = k - inu
		if (mdist <= l) then
			ifail = 3
			return
		end if

!  	check that  nrows  has been set sufficiently large.

		if (kall == 1 .and. nrows < kplus1) then
			ifail = 5
			return
		end if
		if (inup1 > 1) then

!  	normalize the forcing factor so that its leading coefficient
!  	is unity, checking that this coefficient was not zero.

			di = nu(inup1)
			if (di == 0.0_IDP) then
				ifail = 4
				return
			end if
			do i = 1, inup1
				work2(1,i) = nu(i)/di
				work2(2,i) = 0.0_IDP
			end do
		end if

		x1 = xmin
		xm = xmax
		d = xm - x1

!  	the initial values of eps(l,r) (l = 1,2,....ndv and r =
!  	mfirst, mfirst+1,....mlast) of the weighted residuals and
!  	the values work1(1,r)(r=1,2...m) of the normalized
!  	independent variable are computed. n.b. work1(1,r) is
!  	computed from the expression below rather than the more
!  	natural form   (2.0*x(r) - x1 - xm)/d
!  	since the former guarantees the computed value to differ from
!  	the true value by at most  4.0*machine accuracy,  whereas the
!  	latter has no such guarantee.

!  	mdist is now used to record the total number of data points
!  	with non-zero weight.

		mdist = 0
		do r = mfirst, mlast
			wr = w(r)
			if (wr /= 0.0_IDP) mdist = mdist + 1
			mr = r - mfirst + 1
			do l = 1, ndv
				eps(l,r) = wr*y(l,r)
			end do
			work1(1,mr) = ((x(r)-x1)-(xm-x(r)))/d
		end do
		im1 = inu*kall + 1
		betai = 0.0_IDP
		do jplus1 = 1, kplus1
			serr(jplus1) = 0.0_IDP
			do l = 1, ndv
				a(l,im1,jplus1) = 0.0_IDP
			end do
		end do
		do iplus1 = inup1, kplus1

!  		set starting values for degree  i.

			ii = (iplus1-1)*kall + 1
			iplus2 = iplus1 + 1
			if (iplus1 < kplus1) then
				if (kall == 1) then
					do jplus1 = iplus2, kplus1
						do l = 1, ndv
							a(l,ii,jplus1) = 0.0_IDP
						end do
					end do
				end if
				work2(1,iplus2) = 0.0_IDP
				work2(2,iplus2) = 0.0_IDP
			end if
			alpip1 = 0.0_IDP
			di = 0.0_IDP
			do l = 1, ndv
				ci(l) = 0.0_IDP
			end do
			work2(1,iplus1) = 1.0_IDP
			if (kplus1 > 1) work2(2,1) = work2(1,2)
			do r = mfirst, mlast
				if (w(r) == 0.0_IDP) cycle
				mr = r - mfirst + 1
				xcapr = work1(1,mr)

!  			the weighted value work1(2, r)  of the orthogonal polynomial
!  			of degree i at x = x(r) is computed by recurrence from its
!  			chebyshev-series representation.

				if (iplus1 == 1) then
					wrpr = w(r)*0.5*work2(1,1)
					work1(2,mr) = wrpr
				else
					j = iplus2
					if (xcapr > 0.5_IDP) then

!  					reinsch*s modified recurrence.

						factor = 2.0*(1.0-xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							dj = work2(1,j) + dj - factor*bj
							bj = bj + dj
						end do
						wrpr = w(r)*(0.5*work2(1,1)+dj-0.5*factor*bj)
						work1(2,mr) = wrpr
					else if (xcapr >= -0.5_IDP) then

!  					clenshaw*s original recurrence.

						factor = 2.0*xcapr
						bjp1 = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							bjp2 = bjp1
							bjp1 = bj
							bj = work2(1,j) - bjp2 + factor*bjp1
						end do
						wrpr = w(r)*(0.5*work2(1,1)-bjp1+0.5*factor*bj)
						work1(2,mr) = wrpr
					else

!  					gentleman*s modified recurrence.

						factor = 2.0*(1.0+xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							dj = work2(1,j) - dj + factor*bj
							bj = dj - bj
						end do
						wrpr = w(r)*(0.5*work2(1,1)-dj+0.5*factor*bj)
						work1(2,mr) = wrpr
					end if
				end if

!  			the coefficients ci(l) of the i th orthogonal polynomial
!  			l=1,2....ndv and the coefficients alpip1 and betai in the
!  			three-term recurrence relation for the orthogonal
!  			polynomials are computed.

				wrprsq = wrpr**2
				di = di + wrprsq
				do l = 1, ndv
					ci(l) = ci(l) + wrpr*eps(l,r)
				end do
				alpip1 = alpip1 + wrprsq*xcapr
			end do
			do l = 1, ndv
				ci(l) = ci(l)/di
			end do
			if (iplus1 /= inup1) betai = di/dim1
			alpip1 = 2.0*alpip1/di

!  		the weighted residuals eps(l,r)(l=1,2....ndv and r=mfirst,
!  		mfirst+1....mlast) for degree i are computed, together
!  		with their sum of squares, sigmai

			df = mdist - (iplus1-inu)
			do l = 1, ndv
				cil = ci(l)
				sigmai = 0.0_IDP
				do r = mfirst, mlast
					if (w(r) == 0.0_IDP) cycle
					mr = r - mfirst + 1
					epslr = eps(l,r) - cil*work1(2,mr)
					eps(l,r) = epslr
					sigmai = sigmai + epslr**2
				end do

!  			the root mean square residual  s(l, i + 1)  for degree  i
!  			is theoretically undefined if  m = i + 1 - inu  (the
!  			condition for the polynomial to pass exactly through the
!  			data points). should this case arise the r.m.s. residual
!  			is set to zero.

				if (df <= 0.0_IDP) s(l,iplus1) = 0.0_IDP
				if (df > 0.0_IDP) s(l,iplus1) = sqrt(sigmai/df)
			end do

!  		the chebyshev coefficients a(l, i+1, 1), a(l, i+1, 2)....
!  		a(l, i+1, i+1) in the polynomial approximation of degree i
!  		to each set of values of the independent variable
!  		(l=1,2,...,ndv) together with the coefficients
!  		work2(1, 1), work2(1, 2), ..., work2(1, i + 1),   in the
!  		chebyshev-series representation of the  (i + 1) th
!  		orthogonal polynomial are computed.

			do jplus1 = 1, iplus1
				jplus2 = jplus1 + 1
				pij = work2(1,jplus1)
				serr(jplus1) = serr(jplus1) + pij**2/di
				do l = 1, ndv
					a(l,ii,jplus1) = a(l,im1,jplus1) + ci(l)*pij
				end do
				if (jplus1 == kplus1) exit
				work2(1,jplus1) = work2(1,jplus2) + work2(2,jplus1) - alpip1*pij - betai*work2(2,jplus2)
				work2(2,jplus2) = pij
			end do
  			if (iplus1 < kplus1) then
				dim1 = di
				im1 = ii
			end if
		end do
		do iplus1 = 1, kplus1
			serr(iplus1) = 1.0/sqrt(serr(iplus1))
		end do
		ifail = 0

	end subroutine adze02

	subroutine aeue01(m,x,ip,np1,b,w)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 13 revised. use of mark 12 x02 functions (apr 1988).

!  	*******************************************************

!  	npl algorithms library routine q0poly

!  	created 02 05 80.  updated 13 05 80.  release 00/08

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england.

!  	*******************************************************

!  	aeue01.  an algorithm to determine the chebyshev series
!  	representation of a zeroizing polynomial  q0(x),
!  	i.e. a polynomial which takes on zero values (and
!  	possibly zero derivative values) at specified points

!  	input parameters
!  		m        number of distinct x-values.
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		ip       highest order of derivative at each x-value
!  		np1      n + 1,  where  n = number of zeros (including
!  						those of derivatives) to be taken on by
!  						q0(x).  n = m + ip(1) + ip(2) + ... + ip(m).

!  	output parameters
!  		b        chebyshev coefficients of  q0(x)

!  	workspace parameters
!  		w        workspace

!  	.. scalar arguments ..
		integer :: m, np1
!  	.. array arguments ..
		real(IDP), dimension(np1) :: b, w
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
!  	.. local scalars ..
		real(IDP) :: ai, anbig, ct, eps, factor, ovfl, ri, sfac, test, unfl, xi, xtrema
		integer :: i, i2, ifail, ip1, k, l, n, nu

		interface
!  	.. external functions ..
			function x02ame(x)
				use param
				implicit none
				real(IDP) :: x, x02ame
			end function x02ame
!  	.. external subroutines ..
			subroutine e02afe(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine e02afe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, log, sin
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, sxtnth=0.0625_IDP, one=1.0_IDP, two=2.0_IDP, sxteen=16.0_IDP, pi=3.14159265358979323846_IDP
!  	.. executable statements ..

!  	evaluate  q0(x)  at the extrema of the chebyshev polynomial
!  	of degree  n,  employing dynamic scaling to avoid the
!  	possibility of overflow or underflow

		ovfl = sxtnth/x02ame(sxtnth)
		unfl = sxtnth*ovfl
		eps = epsilon(sxtnth)
		n = np1 - 1
		factor = 2*n
		factor = pi/factor
		do i = 1, np1
			w(i) = one
		end do
		do k = 1, m
			ip1 = ip(k) + 1
			xi = x(k)
			do l = 1, ip1
				anbig = zero
				i2 = n + 2
				do i = 1, np1
					i2 = i2 - 2
					ri = i2
					xtrema = sin(factor*ri)
					ai = w(i)*(xtrema-xi)
					w(i) = ai
					if (abs(ai) > anbig) anbig = abs(ai)
				end do
				do
					sfac = one
					if (anbig > ovfl) sfac = sxtnth
					if (anbig < unfl) sfac = sxteen
					if (sfac == one) exit
					anbig = anbig*sfac
					do i = 1, np1
						w(i) = w(i)*sfac
					end do
				end do
			end do
		end do
		ct = ovfl
		do i = 1, np1
			test = abs(w(i))/anbig
			if (test <= eps) w(i) = zero
			if (test > eps .and. test < ct) ct = test
		end do
		ct = ct*anbig
		sfac = one
		do
		 	if (ct < one) exit
			ct = ct*sxtnth
			sfac = sfac*sxtnth
		end do
		do i = 1, np1
			w(i) = w(i)*sfac
		end do

!  	determine the chebyshev representation of  q0(x)

		call e02afe(np1,w,b,ifail)

!  	set the leading coefficient of  q0(x)  to an
!  	exact power of  2

		ai = b(np1)
		ct = log(abs(ai))/log(two)
		nu = ct
		sfac = two**nu
		b(np1) = sfac
		sfac = sfac/ai
		do i = 1, n
			b(i) = b(i)*sfac
		end do

	end subroutine aeue01

	subroutine aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 13 revised. use of mark 12 x02 functions (apr 1988).

!  	*******************************************************

!  	npl algorithms library routine presid

!  	created 18 02 80.  updated 14 05 80.  release 00/28

!  	author ... maurice g. cox.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aeve01.  forms performance indices and residuals
!  	for a polynomial approximation  p(x)  to a set of data
!  	which may contain derivative values.  also indicates
!  	whether the polynomial defined by the specified
!  	coefficients is a better approximation than the best
!  	so far.

!  	input parameters
!  		withq0   true if zeroizing polynomial, else false
!  		withpi   if true, performance indices produced under
!  						all circumstances.  otherwise, produced
!  						only if  p(x)  is an improvement
!  		m        number of x-values.  all distinct
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        x-values.  normalized to  (-1, 1)
!  		n        number of y-values
!  		y        values and derivatives of dependent variable
!  		ip       highest order of derivative at each x-value
!  		imax     one greater than largest value of  ip
!  		a        chebyshev coefficients of  p(x)
!  		la       dimension of  a  and  adif.
!  						 >=  n + 1 if withq0 is true,
!  						 >=  n     otherwise
!  		it       iteration number

!  	input/output parameters
!  		rnmbst   2-norms of residuals corresponding to the
!  						best polynomial so far and its derivatives

!  	output parameters
!  		rnm      2-norms of residuals corresponding to
!  						p(x)  and its derivatives
!  		improv   true if  p(x)  is an improvement, else false
!  		adif     chebyshev coefficients of  (imax - 1)-st
!  						derivative of  p(x)
!  		res      residuals corresponding to y-values
!  	 * pmax     largest performance index
!  	 * pindex   performance indices

!  	note.  the parameters marked  *  are provided only
!  			 if  improv  is true

!  	key local variables
!  		asumax   for current value of  l,  maximum value
!  						over  j = 0, 1, ..., l  of sum of moduli
!  						of chebyshev coefficients of derivative
!  						of order  j  of  p(x)
!  		eps      relative machine precision
!  		iy       location in  y  of current specified
!  						derivative value of order  l
!  		nl       number of specified derivative values
!  						of order  l
!  		nterms   number of terms in chebyshev series
!  						representation of current derivative
!  						of  p(x)
!  		res2nm   2-norm of the residuals corresponding to the
!  						specified derivative values of order  l - 1

!  	.. scalar arguments ..
		real(IDP) :: pmax, xmax, xmin
		integer :: imax, it, la, m, n
		logical :: improv, withpi, withq0
!  	.. array arguments ..
		real(IDP), dimension(la) :: a, adif
		real(IDP), dimension(imax) :: pindex, rnm, rnmbst
		real(IDP), dimension(n) :: res, y
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
!  	.. local scalars ..
		real(IDP) :: absres, asum, asumax, eps, p, plarge, resmax, rscale, t
		integer :: i, ia, iy, l, lm1, nl, nterms
		logical :: imp

		interface
!  	.. external subroutines ..
			subroutine ahze02(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)
				use param
				implicit none
				real(IDP) :: patm1, xmax, xmin
				integer :: ia1, iadif1, la, ladif, np1
				real(IDP), dimension(la) :: a
				real(IDP), dimension(ladif) :: adif
			end subroutine ahze02
			subroutine akze02(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine akze02
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, sqrt
!  	.. data statements ..
		real(IDP) :: half=0.5_IDP, zero=0.0_IDP, one=1.0_IDP, mltplr=8.0_IDP
!  	.. executable statements ..
		pmax = zero
		eps = epsilon(pmax)
		nterms = n
		if (withq0) nterms = nterms + 1
		asumax = zero
		do i = 1, nterms
			adif(i) = a(i)
		end do
		nterms = nterms + 1
		do l = 1, imax
			nterms = nterms - 1

!  		determine sum of moduli  asum  of chebyshev coefficients
!  		of derivative of order  l - 1  of  p(x),  and update
!  		asumax

			asum = half*abs(adif(1))
			if (nterms > 1) then
				do ia = 2, nterms
					asum = asum + abs(adif(ia))
				end do
			end if
			if (asum > asumax) asumax = asum

!  		pindex(l)  is used temporarily to hold the
!  		value of  asumax

			pindex(l) = asumax
			iy = l
			nl = 0

!  		to reduce the possibility of underflow and overflow,
!  		compute  res2nm  as  resmax*sqrt(rscale),  where
!  		resmax  and  rscale  are updated for each residual
!  		corresponding to a specified derivative value of
!  		order  l - 1.  at any stage,  resmax  holds the
!  		modulus of the residual of maximum magnitude so far,
!  		and  rscale  the sum so far of the squares of the
!  		residuals, each scaled by  resmax.

			resmax = zero
			rscale = one
			do i = 1, m

!  			skip if no derivative value of order  l - 1  is
!  			specified at  i-th  x-value

				if (ip(i)+1 >= l) then
					nl = nl + 1

!  				evaluate  p,  the  (l - 1)-st  derivative of  p(x)
!  				at  x = x(i)

					call akze02(nterms,adif,1,la,x(i),p)

!  				save residual corresponding to this value

					if (withq0) res(iy) = -p
					if ( .not. withq0) res(iy) = y(iy) - p

!  				update  resmax  and  rscale

					absres = abs(res(iy))
					if (absres /= zero .and. absres <= resmax) rscale = rscale + (absres/resmax)**2
					if (absres /= zero .and. absres > resmax) rscale = rscale*(resmax/absres)**2 + one
					if (absres /= zero .and. absres > resmax) resmax = absres
				end if
				iy = iy + ip(i) + 1
			end do
			rnm(l) = resmax*sqrt(rscale)

!  		form coefficients in chebyshev series representation
!  		of  l-th  derivative of  p(x)  from those of
!  		(l - 1)-st derivative

			if (l < imax) call ahze02(nterms,xmin,xmax,adif,1,la,t,adif,1,la)
		end do

!  	if not on zero-th iteration,
!  	detect whether there has been an improvement ...

		imp = (it == 0)
		if (it /= 0) then
			do l = 1, imax
				if (rnm(l) < rnmbst(l)) imp = .true.
			end do
		end if
		if (imp .or. withpi) then

!  	... and, if so, or if zero-th iteration, or if
!  	they are required anyway, form the performance
!  	indices corresponding to the improved
!  	approximation

			plarge = zero
			do l = 1, imax

!  			nl  is the number of derivative values of order  l - 1

				nl = 0
				lm1 = l - 1
				do i = 1, m
					if (ip(i) >= lm1) nl = nl + 1
				end do
				t = nl

!  			retrieve the value of asumax

				asumax = pindex(l)
				if (asumax /= zero) pindex(l) = rnm(l)/(mltplr*eps*asumax*sqrt(t))
				if (asumax == zero) pindex(l) = zero
				if (pindex(l) > plarge) plarge = pindex(l)
			end do
			pmax = plarge
		end if
		improv = imp

	end subroutine aeve01

	subroutine aewe01(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine pntrpa

!  	created 20 12 79.  updated 13 05 80.  release 00/47

!  	authors ... gerald t. anthony, maurice g. cox,
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aewe01. a routine, without checks, which determines and
!  	refines a polynomial interpolant  q(x)  to data which
!  	may contain derivatives, and an associated zeroizing
!  	polynomial  q0(x).

!  		m        number of distinct x-values
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        independent variable values (distinct)
!  		y        values and derivatives of
!  						dependent variable.
!  		ip       highest order of derivative at each x-value.
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		np1      value of  n + 1
!  		itmin,
!  		itmax    minimum and maximum number of iterations to be
!  						performed.
!  		iwrk(1)  see workspace (and associated
!  						dimension) parameters

!  	output parameters
!  		a        chebyshev coefficients of  q(x)
!  		b        chebyshev coefficients of  q0(x)

!  	workspace (and associated dimension) parameters
!  		wrk      real workspace array.  the first imax elements
!  						contain, on exit, performance indices for
!  						the interpolating polynomial, and the next
!  						n  elements the computed residuals
!  		lwrk     dimension of wrk. lwrk must be at least
!  						8*n + 5*imax + m + 5, where
!  						imax is one more than the largest element
!  						of the array ip.
!  		iwrk     integer workspace array.
!  						the first element of this array is used
!  						as an input parameter (which is destroyed
!  						on exit).  the zeroizing polynomial  q0(x)
!  						is constructed or not according to whether
!  						iwrk(1)  is non-zero or zero.
!  		liwrk    dimension of iwrk.  at least 2*m + 2.

!  	failure indicator parameter
!  		ifail    failure indicator
!  						0 - successful termination
!  						1 - iteration limit in deriving  q(x)
!  						2 - divergent iteration in deriving  q(x)
!  						3 - iteration limit in deriving  q0(x)
!  						4 - divergent iteration in deriving  q0(x)

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, itmax, itmin, liwrk, lwrk, m, n, np1
!  	.. array arguments ..
		real(IDP), dimension(n) :: a, y
		real(IDP), dimension(np1) :: b
		real(IDP), dimension(lwrk) :: wrk
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
		integer, dimension(liwrk) :: iwrk
!  	.. local scalars ..
		real(IDP) :: pmax
		integer :: i, iadif, iatrl, ibdif, ibtrl, ic, id, ida, idb, ierror, iftau, ilocx, ilocy, imax, initq, initq0, ipiq, ipiq0, &
					  iptrl, ires, irnm, irtrnm, iw, ix
		logical :: withq0

		interface
!  	.. external subroutines ..
			subroutine aeye01(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w, &
									adif,da,rnm,rtrlnm,locx,locy,ifail)
				use param
				implicit none
				real(IDP) :: pmax, xmax, xmin
				integer :: ifail, imax, itmax, itmin, la, m, n, nit, np1
				logical :: withq0
				real(IDP), dimension(la) :: a, adif, atrial, da
				real(IDP), dimension(n) :: c, d, ftau, res, y
				real(IDP), dimension(imax) :: pindex, ptrial, rnm, rtrlnm
				real(IDP), dimension(np1) :: w
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx, locy
			end subroutine aeye01
			subroutine akye02(xmin,xmax,x,xcap)
				use param
				implicit none
				real(IDP) :: x, xcap, xmax, xmin
			end subroutine akye02
		end interface

!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP
!  	.. executable statements ..
		ierror = 0
		imax = 0
		do i = 1, m
			if (ip(i) > imax) imax = ip(i)
		end do
		imax = imax + 1

!  	indicate whether  q0(x)  is required

		withq0 = (iwrk(1) /= 0)

!  	treat the case n = 1 separately.

		if (n == 1) then
			a(1) = two*y(1)
			if (withq0) then
				call akye02(xmin,xmax,x(1),b(1))
				b(1) = -two*b(1)
				b(2) = one
			end if

!  	set to zero the numbers of iterations taken, the (sole)
!  	residual and the values of the performance indices
!  	in this special case, and finish

			iwrk(1) = 0
			if (withq0) iwrk(2) = 0
			wrk(1) = zero
			wrk(2) = zero
			if (withq0) wrk(3) = zero
		end if

!  	transform the  x*s  to  (-1, 1)

		ix = 2*imax + n
		do i = 1, m
			ix = ix + 1
			call akye02(xmin,xmax,x(i),wrk(ix))
		end do
		if (n == 1) then
			ifail = ierror
			return
		end if
		if (withq0) then

!  		workspace allocation for call to  aeye01  ...

			ires = imax + 1
			ipiq0 = ires + n
			ix = ipiq0 + imax
			ibtrl = ix + m
			iptrl = ibtrl + np1
			iftau = iptrl + imax
			ic = iftau + n
			id = ic + n
			iw = id + np1
			ibdif = iw + np1
			idb = ibdif + np1
			irnm = idb + np1
			irtrnm = irnm + imax
			initq0 = 2
			ilocx = initq0 + 1
			ilocy = ilocx + m

!  		... to determine  q0(x)

			call aeye01(.true.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,b,np1,wrk(ires),pmax,wrk(ipiq0),iwrk(initq0), &
							wrk(ibtrl),wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(ibdif),wrk(idb),wrk(irnm),wrk(irtrnm), &
							iwrk(ilocx),iwrk(ilocy),ierror)

!  		re-code failure indicator

			if (ierror /= 0) ierror = ierror + 2
			if (ierror /= 0) then
				ifail = ierror
				return
			end if
		end if

!  	workspace allocation for call to  aeye01  ...

		ipiq = 1
		ires = ipiq + imax
		ix = ires + n + imax
		iatrl = ix + m
		iptrl = iatrl + n
		iftau = iptrl + imax
		ic = iftau + n
		id = ic + n
		iw = id + n
		iadif = iw
		ida = iadif + np1
		irnm = ida + np1
		irtrnm = irnm + imax
		initq = 1
		ilocx = initq + 2
		ilocy = ilocx + m

!  	... to determine  q(x)

		call aeye01(.false.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,n,wrk(ires),pmax,wrk(ipiq),iwrk(initq),wrk(iatrl), &
						wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(iadif),wrk(ida),wrk(irnm),wrk(irtrnm),iwrk(ilocx),iwrk(ilocy), &
						ierror)

		ifail = ierror

	end subroutine aewe01

	subroutine aexe01(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine divdif

!  	created 17 07 79.  updated 14 05 80.  release 00/08

!  	authors ... maurice g. cox and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aexe01.  an algorithm to determine the next coefficient
!  	in the newton form of an interpolating polynomial

!  	input parameters
!  		m        number of distinct x-values.
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		ip       highest order of derivative at each x-value
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		locx     pointers to x-values in constructing
!  						newton form of polynomial
!  		!  		newton coefficients determined so far
!  		n!  	  number of newton coefficients determined so far
!  		xnew     element of  x  associated with new
!  						newton coefficient
!  		ixnext   number of x-values so far incorporated
!  						(including  xnew)
!  		ynew     scaled derivative value corresponding to
!  						xnew  and  nordp1
!  		nordp1   one plus order of derivative
!  						associated with  ynew

!  	input/output parameters
!  		d        elements in previous, and then new, upward
!  						sloping diagonal of divided difference table

!  	output parameters
!  		cnew     new newton coefficient generated

!  	.. scalar arguments ..
		real(IDP) :: cnew, xnew, ynew
		integer :: ixnext, m, n, nc, nordp1
!  	.. array arguments ..
		real(IDP), dimension(n) :: c, d
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx
!  	.. local scalars ..
		real(IDP) :: dif
		integer :: ic, is, ix, k, locxi
!  	.. executable statements ..
		ic = nc - nordp1 + 1
		d(1) = ynew
		if (ixnext /= 1) then
			is = 0
			ix = 0
			do k = 1, ic
				if (k > is) then
					ix = ix + 1
					locxi = locx(ix)
					is = is + ip(locxi) + 1
					dif = x(locxi) - xnew
				end if
				if (nordp1 == 1) d(k+1) = (c(k)-d(k))/dif
				if (nordp1 > 1) d(k+1) = (d(k+1)-d(k))/dif
			end do
		end if
		cnew = d(ic+1)

	end subroutine aexe01

	subroutine aeye01(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w,adif,da, &
							rnm,rtrlnm,locx,locy,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine refh

!  	created 17 07 79.  updated 14 05 80.  release 00/47

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aeye01.  a routine to approximate and then refine an
!  	interpolating polynomial  q(x)  or a zeroizing
!  	polynomial  q0(x)  in its chebyshev representation

!  	input parameters
!  		withq0   true if zeroizing polynomial, else false
!  		m        the number of distinct data points.
!  		x        array containing the distinct x-values,
!  						normalized if necessary to (-1, 1).
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  	 * y        array containing values and derivatives of
!  						the dependent variable.
!  		ip       array specifying the highest order of
!  						derivative at each x-value.
!  		imax     one more than the largest element of the
!  						array ip.
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		np1      value of  n + 1
!  		itmin,
!  		itmax    the lower and upper limits on the iterative
!  						process.

!  	output (and associated dimension) parameters
!  		a        chebyshev coefficients of polynomial
!  		la       dimension of  a.
!  						 >=  n  if interpolating polynomial
!  						 >=  n + 1  if zeroizing polynomial
!  		res      residuals of polynomial
!  		pmax     largest performance index
!  		pindex   performance indices
!  		nit      number of iterations taken

!  	workspace parameters
!  		atrial   trial values of the chebyshev coefficients
!  		ptrial   performance indices corresponding to  atrial
!  	 * ftau     scaled values of  y.  if  y(i)  is the
!  						value of an  r-th  derivative, then
!  						((xmax - xmin)/2)**r/(factorial r)
!  						times  y(i)  is the value of  ftau(i)
!  	 * !  		coefficients in newton form of polynomial
!  	 * d        intermediate divided difference values
!  	 **w        values of correction polynomial at
!  						chebyshev extrema
!  		adif     chebyshev coefficients of a derivative of
!  						an approximation to the polynomial
!  		da       chebyshev coefficients of a
!  						correction polynomial
!  		rnm      residual norms corresponding to  a
!  		rtrlnm   residual norms corresponding to  atrial
!  	 * locx     pointers to x-values in constructing
!  						newton form of polynomial
!  	 * locy     pointers to y-values corresponding to x-values

!  	failure indicator parameter
!  		ifail    failure indicator
!  						0 - successful termination
!  						1 - iteration limit exceeded
!  						2 - iteration divergent

!  			notes.  (1) the elements of the arrays marked  *  are
!  							not accessed if  withq0  is  true.
!  					  (2) the elements of the array marked  **  is
!  							not accessed if  withq0  is  false.

!  	.. scalar arguments ..
		real(IDP) :: pmax, xmax, xmin
		integer :: ifail, imax, itmax, itmin, la, m, n, nit, np1
		logical :: withq0
!  	.. array arguments ..
		real(IDP), dimension(la) :: a, adif, atrial, da
		real(IDP), dimension(n) :: c, d, ftau, res, y
		real(IDP), dimension(imax) :: pindex, ptrial, rnm, rtrlnm
		real(IDP), dimension(np1) :: w
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx, locy
!  	.. local scalars ..
		real(IDP) :: amax, atrnrm, danrm, pmxtrl, scale
		integer :: i, ierror, it, itemp, itmxp1, itp1, l, nfref, npilt1, nterms
		logical :: improv, withpi, zeroda

		interface
!  	.. external subroutines ..
			subroutine aeue01(m,x,ip,np1,b,w)
				use param
				implicit none
				integer :: m, np1
				real(IDP), dimension(np1) :: b, w
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
			end subroutine aeue01
			subroutine aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)
				use param
				implicit none
				real(IDP) :: pmax, xmax, xmin
				integer :: imax, it, la, m, n
				logical :: improv, withpi, withq0
				real(IDP), dimension(la) :: a, adif
				real(IDP), dimension(imax) :: pindex, rnm, rnmbst
				real(IDP), dimension(n) :: res, y
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
			end subroutine aeve01
			subroutine aeze01(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: m, n
				real(IDP), dimension(n) :: a, c, d, ftau, y
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx, locy
			end subroutine aeze01
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, log
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, one=1.0_IDP, sxteen=16.0_IDP
!  	.. executable statements ..
		ierror = 0

!  	number of terms in polynomial

		nterms = n
		if (withq0) nterms = nterms + 1

!  	number of performance indices less than one

		npilt1 = 0

!  	indicate that performance indices are to be produced
!  	only if an improvement is obtained

		withpi = .false.

!  	indicate that the fine refinement stage has not yet started

		nfref = -2

!  	set residuals initially equal to specified y-values
!  	(if  q(x)  required) or zero (if  q0(x)  required)

		do i = 1, n
			if ( .not. withq0) res(i) = y(i)
			if (withq0) res(i) = zero
		end do

!  	initialize trial chebyshev coefficients

		do i = 1, nterms
			atrial(i) = zero
		end do

!  	commence iterative refinement

		itmxp1 = itmax + 1
		do itp1 = 1, itmxp1

!  		it  is the actual iteration number,  it = 0
!  		corresponding to the first estimate of the polynomial

			it = itp1 - 1

!  		determine chebyshev coefficients  da  of polynomial
!  		approximately satisfying the conditions in  res

			if (withq0 .and. it == 0) call aeue01(m,x,ip,np1,da,w)
			if ( .not. withq0 .or. (withq0 .and. it > 0)) call aeze01(m,xmin,xmax,x,res,ip,n,da,locx,locy,ftau,d,c)

!  		skip test for divergence if on zero-th iteration

			if (it > 0) then

!  		determine the norms of  da  and (the previous)  atrial

				danrm = half*abs(da(1))
				atrnrm = half*abs(atrial(1))
				if (n > 1) then
					do i = 2, n
						danrm = danrm + abs(da(i))
						atrnrm = atrnrm + abs(atrial(i))
					end do
				end if
				if (withq0) atrnrm = atrnrm + abs(atrial(np1))

!  		assume divergence if the norm of  da  is not
!  		less than that of  atrial  ...

				if (danrm >= atrnrm) ierror = 2
				if (danrm >= atrnrm) exit
			end if

!  		... otherwise determine new trial approximation

			zeroda = .true.
			do i = 1, n
				atrial(i) = atrial(i) + da(i)
				if (da(i) /= zero) zeroda = .false.
			end do
			if (withq0 .and. it == 0) then
				atrial(np1) = da(np1)
				if (da(np1) /= zero) zeroda = .false.
			end if

!  		determine residuals, performance indices and
!  		largest performance index corresponding to
!  		trial coefficients  atrial

			call aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,atrial,la,it,rnm,rtrlnm,improv,adif,res,pmxtrl,ptrial)

!  		set dummy, non-zero, value of  pmxtrl  if no
!  		improvement, otherwise it is undefined

			if ( .not. improv) pmxtrl = sxteen

!  		if on first iteration, or if the largest performance
!  		index is zero, or if all components of  da  are zero,
!  		take the trial set of coefficients and performance
!  		indices as the best (so far)

			if (it == 0 .or. pmxtrl == zero .or. zeroda) then
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
				pmax = pmxtrl
			end if

!  		finish if largest performance index is zero
!  		or if all components of  da  are zero
!  		(i.e. no further improvement is possible)

			if (pmxtrl == zero .or. zeroda) exit

!  		indicate whether the fine refinement stage has commenced
!  		(i.e. for the first time all performance indices are
!  		less than one)

			if (nfref == -2 .and. pmxtrl < one) nfref = -1

!  		branch according to whether the process is in the
!  		fine refinement stage  (nfref  >=  0)  or not
!  		(nfref  ==  -1)

			if (nfref < -1) then

!  		the process is in the course refinement phase.
!  		update the coefficients and the corresponding
!  		norms and performance indices if
!  			(i)  there has been an improvement in (at
!  				  least) one of the residual norms, and
!  			(ii) the number of performance indices
!  				  that are less than one has not
!  				  increased compared with those of
!  				  the best polynomial so far.

				if ( .not. improv) cycle
				itemp = 0
				do l = 1, imax
					if (ptrial(l) < one) itemp = itemp + 1
				end do
 				if (itemp < npilt1) cycle
				npilt1 = itemp
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
				pmax = pmxtrl
				cycle
			end if

!  		the process is in the fine refinement phase.
!  		update the coefficients and the corresponding
!  		norms and performance indices if
!  			(i)  there has been an improvement in (at
!  				  least) one of the residual norms, and
!  			(ii) the largest performance index is less
!  				  than the largest of that of the best
!  				  polynomial so far.
!  		increment the number of fine refinements (the number
!  		of refinements since the first occasion when all
!  		performance indices were less than unity), exiting
!  		if as many as  itmin  fine refinements have been
!  		performed

			if (improv .and. pmxtrl < pmax) then
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
 				pmax = pmxtrl
			end if
			nfref = nfref + 1
			if (nfref >= itmin) exit
		end do

!  	the process has not succeeded in reducing all the
!  	performance indices to less than unity

		if (itp1 > itmxp1) then
			ierror = 1
			it = itmax
		end if

!  	number of iterations actually taken

		nit = it
		if (withq0) then

!  	in the case of  q0(x),  scale its coefficients by
!  	an integral power of  16  such that the largest
!  	coefficient is of order unity

			amax = zero
			do i = 1, np1
				if (abs(a(i)) > amax) amax = abs(a(i))
			end do
			if (amax > zero) then
				i = log(amax)/log(sxteen)
				scale = sxteen**(-i)
				do i = 1, np1
					a(i) = scale*a(i)
				end do
			end if
		end if

!  	return residuals corresponding to selected coefficients

		withpi = .true.
		call aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnm,rtrlnm,improv,adif,res,pmax,pindex)
		ifail = ierror

	end subroutine aeye01

	subroutine aeze01(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine qpoly

!  	created 02 05 80.  updated 14 05 80.  release 00/09

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england.

!  	*******************************************************

!  	aeze01. an algorithm to determine the chebyshev series
!  	representation of a polynomial interpolant  q(x)  to
!  	arbitrary data points where derivative information may
!  	be given.

!  	input parameters
!  		m        number of distinct x-values.
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		y        values and derivatives of dependent variable
!  		ip       highest order of derivative at each x-value
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).

!  	output parameters
!  		a        chebyshev coefficients of  q(x)

!  	workspace parameters
!  		locx     pointers to x-values in constructing
!  						newton form of polynomial
!  		locy     pointers to y-values corresponding to x-values
!  		ftau     scaled values of  y.  if  y(i)  is the
!  						value of an  r-th  derivative, then
!  						((xmax - xmin)/2)**r/(factorial r)
!  						times  y(i)  is the value of  ftau(i)
!  		d        intermediate divided difference values
!  		!  		newton coefficients of  q(x)

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: m, n
!  	.. array arguments ..
		real(IDP), dimension(n) :: a, c, d, ftau, y
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx, locy
!  	.. local scalars ..
		real(IDP) :: cmin, cnew, factor, ri, rj, s, scale, sfac, v, xch
		integer :: i, i2, ic, icmin, ifail, iftau, ip1, isave, iy, j, jmax, k, krev, l, lmax, locxi, locxj, locxk, locyi, nc

		interface
!  	.. external subroutines ..
			subroutine aexe01(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)
				use param
				implicit none
				real(IDP) :: cnew, xnew, ynew
				integer :: ixnext, m, n, nc, nordp1
				real(IDP), dimension(n) :: c, d
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx
			end subroutine aexe01
			subroutine e02afe(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine e02afe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, sin
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP, pi=3.14159265358979323846_IDP
!  	.. executable statements ..
		scale = (xmax-xmin)/two

!  	initialize x- and y-pointers

		iy = 0
		do i = 1, m
			locx(i) = i
			iy = iy + 1
			locy(i) = iy
			ftau(iy) = y(iy)
			jmax = ip(i)
			if (jmax == 0) cycle

!  		form the scaled derivatives, i.e. an  r-th  derivative
!  		value is divided by  factorial r  and multiplied
!  		by the  r-th  power of  (xmax - xmin)/2

			sfac = one
			do j = 1, jmax
				iy = iy + 1
				rj = j
				sfac = sfac*scale/rj
				ftau(iy) = y(iy)*sfac
			end do
		end do

!  	form successive upward sloping diagonals of
!  	the divided difference table

		nc = 0
	outer:do j = 1, m

!  			choose each x-value in turn to make the corresponding
!  			newton coefficient as small in magnitude as possible

				do i = j, m
					locxi = locx(i)
					locyi = locy(locxi)
					call aexe01(m,x,ip,n,locx,c,nc,x(locxi),j,ftau(locyi),1,cnew,d)
					if (i > j .and. abs(cnew) >= abs(cmin)) cycle
					cmin = cnew
					icmin = i
				end do
				c(nc+1) = cmin
				isave = locx(j)
				locxj = locx(icmin)
				locx(icmin) = isave
				locx(j) = locxj

!  			calculate the resulting newton coefficient (i.e.
!  			repeat the above computation, but only in the case
!  			leading to the smallest new coefficient)

				iftau = locy(locxj) - 1
				ip1 = ip(locxj) + 1
	inner:	do i = 1, ip1
					iftau = iftau + 1
					call aexe01(m,x,ip,n,locx,c,nc,x(locxj),j,ftau(iftau),i,c(nc+1),d)
					nc = nc + 1
					if (nc == n) exit outer
				end do inner
			end do outer

!  	evaluate  q(x)  (from its newton form) at the extrema
!  	of the chebyshev polynomial of degree  n - 1  ...

		factor = 2*n - 2
		factor = pi/factor
		i2 = n + 1
		do i = 1, n
			i2 = i2 - 2
			ri = i2
			xch = sin(factor*ri)
			s = c(n)
			ic = n
			k = m + 1
			do krev = 1, m
				k = k - 1
				locxk = locx(k)
				lmax = ip(locxk) + 1
				if (k == m) lmax = lmax - 1
				if (lmax <= 0) cycle
				v = xch - x(locxk)
				do l = 1, lmax
					ic = ic - 1
					s = s*v + c(ic)
				end do
			end do
			d(i) = s
		end do

!  	... in order to determine the coefficients in its
!  	chebyshev representation

		call e02afe(n,d,a,ifail)

	end subroutine aeze01

	subroutine ahze02(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine auxcdf

!  	created 1/5/79    updated 23/1/80     release no. 00/03

!  	authors.. gerald t anthony, maurice g cox, j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1    = n+1 where n is degree of given polynomial
!  		xmin   lower limit of range of x
!  		xmax   upper limit of range of x
!  		a      coefficients a0, a1,...an of the given polynomial
!  		ia1       are stored in array a in positions 1, 1+ia1,...
!  					 1+n*ia1, respectively
!  		la     the declared dimension of array a

!  	output parameters
!  		patm1  the value of the given polynomial at xmin
!  		adif   the coefficients of the derivative polynomial
!  		iadif1    are returned in array adif in positions
!  					 1, 1+iadif1,...1+(n-1)*iadif1
!  		ladif  the declared dimension of array adif

!  	differentiate the series with coefficients a of degree n
!  	(i.e. np1 coefficients) to obtain the series with coefficients
!  	adif of degree n-1. also set next higher coefficient to zero.

!  	.. scalar arguments ..
		real(IDP) :: patm1, xmax, xmin
		integer :: ia1, iadif1, la, ladif, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
		real(IDP), dimension(ladif) :: adif
!  	.. local scalars ..
		real(IDP) :: ptemp, r, sclftr, u, v, w
		integer :: i, n, na, nadif
!  	.. data statements ..
		real(IDP) :: two=2.0_IDP
!  	.. executable statements ..
		u = 0.0_IDP
		v = u
		sclftr = two/(xmax-xmin)
		n = np1 - 1
		nadif = n*iadif1 + 1
		ptemp = u
		if (n > 0) then
			na = n*ia1 + 1
			do i = 1, n
				r = np1 - i
				w = u + two*r*a(na)
				ptemp = a(na) - ptemp

!  			store coeff formed previous time round. first time round
!  			store zero as coeff of degree n.

				adif(nadif) = sclftr*v
				u = v
				v = w
				na = na - ia1
				nadif = nadif - iadif1
			end do
		end if
		adif(nadif) = sclftr*v
		patm1 = a(1)/two - ptemp

	end subroutine ahze02

	subroutine akye02(xmin,xmax,x,xcap)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine nrmlze

!  	created 5/5/78    updated 11/12/78    release no. 00/02.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		xmin     minimum value of x
!  		xmax     maximum value of x
!  		x        unnormalized argument in range (xmin, xmax)

!  	output parameter
!  		xcap     normalized value of x

!  	a value of x is given, such that
!  	xmin  <=  x  <=  xmax.
!  	xcap is calculated so that -1  <=  x  <=  +1.

!  	this form for xcap ensures that the computed value has a
!  	very small absolute error.

!  	.. scalar arguments ..
		real(IDP) :: x, xcap, xmax, xmin
!  	.. executable statements ..
		xcap = ((x-xmin)-(xmax-x))/(xmax-xmin)

	end subroutine akye02

	subroutine akze02(np1,a,ia1,la,xcap,result)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine tval1

!  	created 9/5/78    updated 6/4/79    release no. 00/07.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1      np1 = n + 1. n is the degree of the
!  					chebyshev series
!  		a        the array where the coefficients are stored
!  		ia1      the address increment of a
!  		la       dimension of a
!  		xcap     normalized argument of the polynomial

!  	output parameter
!  		result   value of the summation

!  	np1 chebyshev coefficients a0, a1, ..., an, are
!  	stored in the array a in positions 1, 1+ia1, 1+2*ia1, ...,
!  	1+n*ia1, where n = np1 - 1.
!  	ia1 must not be negative.
!  	la must be at least equal to 1 + n*ia1.
!  	the argument xcap is assumed to lie in the range
!  	-1  <=  xcap  <=  +1.
!  	the value of the polynomial of degree n
!  	a0t0(xcap)/2 + a1t1(xcap) + a2t2(xcap) + + ... + antn(xcap),
!  	is calculated for the argument xcap storing it in result.
!  	the recurrence relation by clenshaw, modified by reinsch
!  	and gentleman, is used.

!  	.. scalar arguments ..
		real(IDP) :: result, xcap
		integer :: ia1, la, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
!  	.. local scalars ..
		real(IDP) :: aj, bj, cj, factor, sum
		integer :: j, jrev, n
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, two=2.0_IDP
!  	.. executable statements ..
		if (np1 == 1) then
			sum = half*a(1)
		else
			n = np1 - 1
			aj = zero
			bj = zero
			j = 1 + np1*ia1
			if (xcap > half) then

!  		reinschs modified recurrence.

				factor = two - (xcap+xcap)

!  		bracketing necessary in order to avoid errors

				do jrev = 1, n
					j = j - ia1
					aj = a(j) + aj - bj*factor
					bj = aj + bj
				end do
				sum = half*a(1) + aj - half*factor*bj
			else if (xcap >= -half) then

!  		clenshaws original recurrence.

				factor = xcap + xcap
				do jrev = 1, n
					j = j - ia1
					cj = bj
					bj = aj
					aj = a(j) - cj + bj*factor
				end do
				sum = half*a(1) - bj + half*factor*aj
			else

!  		gentlemans modified recurrence.

				factor = two + (xcap+xcap)

!  		bracketing necessary so as to avoid errors

				do jrev = 1, n
					j = j - ia1
					aj = a(j) - aj + bj*factor
					bj = aj - bj
				end do
				sum = half*a(1) - aj + half*factor*bj
			end if
		end if
		result = sum

	end subroutine akze02

	integer function p01abe(ifail,ierror,srname,nrec,rec)

		use param
		implicit none

!  	mark 11.5(f77) release. nag copyright 1986.
!  	mark 13 revised. ier-621 (apr 1988).
!  	mark 13b revised. ier-668 (aug 1988).

!  	p01abe is the error-handling routine for the nag library.

!  	p01abe either returns the value of ierror through the routine
!  	name (soft failure), or terminates execution of the program
!  	(hard failure). diagnostic messages may be output.

!  	if ierror = 0 (successful exit from the calling routine),
!  	the value 0 is returned through the routine name, and no
!  	message is output

!  	if ierror is non-zero (abnormal exit from the calling routine),
!  	the action taken depends on the value of ifail.

!  	ifail =  1: soft failure, silent exit (i.e. no messages are
!  					output)
!  	ifail = -1: soft failure, noisy exit (i.e. messages are output)
!  	ifail =-13: soft failure, noisy exit but standard messages from
!  					p01abe are suppressed
!  	ifail =  0: hard failure, noisy exit

!  	for compatibility with certain routines included before mark 12
!  	p01abe also allows an alternative specification of ifail in which
!  	it is regarded as a decimal integer with least significant digits
!  	cba. then

!  	a = 0: hard failure  a = 1: soft failure
!  	b = 0: silent exit   b = 1: noisy exit

!  	except that hard failure now always implies a noisy exit.

!  	s.hammarling, m.p.hooper and j.j.du croz, nag central office.

!  	.. scalar arguments ..
		integer :: ierror, ifail, nrec
		character(*) :: srname
!  	.. array arguments ..
		character(*), dimension(:) :: rec
!  	.. local scalars ..
		integer :: i, nerr
		character(len=72) :: mess
!  	.. intrinsic functions ..
!  	intrinsic  				 abs, mod

		interface
!  	.. external subroutines ..
			subroutine abzp01
			end subroutine abzp01
			subroutine x04aae(i,nerr)
				implicit none
				integer :: i, nerr
			end subroutine x04aae
			subroutine x04bae(nout,rec)
				implicit none
				integer :: nout
				character(*) :: rec
			end subroutine x04bae
		end interface

!  	.. executable statements ..
		if (ierror /= 0) then
!  	   abnormal exit from calling routine
			if (ifail == -1 .or. ifail == 0 .or. ifail == -13 .or. (ifail > 0 .and. mod(ifail/10,10) /= 0)) then
!  			noisy exit
				call x04aae(0,nerr)
				do i = 1, nrec
					call x04bae(nerr,rec(i))
				end do
				if (ifail /= -13) then
					write (mess,'(" ** abnormal exit from nag library routine ",a,": ifail"," =",i6)') srname, ierror
					call x04bae(nerr,mess)
					if (abs(mod(ifail,10)) /= 1) then
!  					hard failure
						call x04bae(nerr," ** nag hard failure - execution terminated")
						call abzp01
					else
!  					soft failure
						call x04bae(nerr," ** nag soft failure - control returned")
					end if
				end if
			end if
		end if
		p01abe = ierror

	end function p01abe

	subroutine abzp01

!  	mark 11.5(f77) release. nag copyright 1986.

!  	terminates execution when a hard failure occurs.

!  	******************** implementation note ********************
!  	the following stop statement may be replaced by a call to an
!  	implementation-dependent routine to display a message and/or
!  	to abort the program.
!  	*************************************************************
!  	.. executable statements ..
		stop

	end subroutine abzp01

	subroutine x04aae(i,nerr)

		implicit none

!  	mark 7 release. nag copyright 1978
!  	mark 7c revised ier-190 (may 1979)
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 14 revised. ier-829 (dec 1989).
!  	if i = 0, sets nerr to current error message unit number
!  	(stored in nerr1).
!  	if i = 1, changes current error message unit number to
!  	value specified by nerr.

!  	.. scalar arguments ..
		integer :: i, nerr
!  	.. local scalars ..
		integer, save :: nerr1=0
!  	.. executable statements ..
		if (i == 0) nerr = nerr1
		if (i == 1) nerr1 = nerr

	end subroutine x04aae

	subroutine x04bae(nout,rec)

		implicit none

!  	mark 11.5(f77) release. nag copyright 1986.

!  	x04bae writes the contents of rec to the unit defined by nout.

!  	trailing blanks are not output, except that if rec is entirely
!  	blank, a single blank character is output.
!  	if nout < 0, i.e. if nout is not a valid fortran unit identifier,
!  	then no output occurs.

!  	.. scalar arguments ..
		integer :: nout
		character(*) :: rec
!  	.. local scalars ..
		integer :: i
!  	.. intrinsic functions ..
!  	intrinsic  	    len
!  	.. executable statements ..
		if (nout >= 0) then
!  	   remove trailing blanks
			do i = len(rec), 2, -1
				if (rec(i:i) /= " ") exit
			end do
!  	   write record to external file
			write (nout,'(a)') rec(1:i)
		end if

	end subroutine x04bae

	function x02ame(x)

		use param
		implicit none

!  	mark 12 release. nag copyright 1986.

!  	returns the 'safe range' parameter
!  	i.e. the smallest positive model number z such that
!  	for any x which satisfies x >= z and x <= 1/z
!  	the following can be computed without overflow, underflow or other
!  	error

!  	   -x
!  	   1.0/x
!  	   sqrt(x)
!  	   log(x)
!  	   exp(log(x))
!  	   y**(log(x)/log(y)) for any y

!  	x is a dummy argument

		real(IDP) :: x, x02ame

!  	.. intrinsic functions ..
!  	intrinsic  			 sqrt
!  	.. executable statements ..
		x02ame = sqrt(2.0_IDP)*tiny(x)

	end function x02ame
