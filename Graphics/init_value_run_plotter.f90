!  This code is formed by taking pieces (modules, rddump, part of main code from
!   stability code, so that dump file ca be read and plots made (initial value version).
!   May need modification for new versions of FAR3d if new variables are added.
!
	module param
		implicit none
		integer, parameter :: IDP = kind(1.d0)
	end module param												
													
	module cotrol
		use param
		implicit none
		save
		character(len=8), dimension(:), allocatable :: numhist
		character(len=2), dimension(3) :: numrun,numruno,numruns
		character(len=5) :: numvac
		
		real(IDP) :: stdifp,stdifu,stdifnf,stdifvf,stdifv,stdifnalp,stdifvalp,dt0,Adens,Bdens,Adensalp,Bdensalp, &
		             LcA0,LcA1,LcA2,LcA3,LcA0alp,LcA1alp,LcA2alp,LcA3alp, &
                     ext_prof,omegar,iflr,r_epflr,r_epflralp,dpres,DIIID_u,betath_factor,etascl,reta,eta0,etalmb,deltaq,deltaiota		
		integer :: ihist,nocpl,maxstp,nstep,ndump,nprint,lplots,itime,nstep1,nonlin,noeqn,edge_p, &
                   iflr_on,epflr_on,ieldamp_on,twofl_on,alpha_on,inalp,ivalp,iq,iw,ix1,ix2,iwa,ix1a,ix2a, &
				   EP_dens_on,Alpha_dens_on,EP_vel_on,Alpha_vel_on,Trapped_on,ext_prof_len,q_prof_on,Eq_vel_on, & 
                                                           Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,Auto_grid_on,Edge_on
                                                                         
		character(len=40) :: eq_name,ext_prof_name
		character(len=8) :: rdump_file_name		
		logical :: matrix_out
		integer :: leqdim,ldim,jdim,nstres,variable
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				   
			!  numhist is the simulation name you want to continue
			!  numrun adn numruns are the run number													
			!  numruno is the previous run number
			!  numvac is the equilibrium number
			!  stdifp is the pressure diffusivity
			!  stdifu is the vorticity diffusivity
		                !  stdifv is the parallel velocity diffusivity		
			!  stdifnf is the fast particle density diffusivity												
			!  stdifvf is the fast particle parallel velocity diffusivity		
			!  stdifnfalp is the fast particle second species density diffusivity									
			!  stdifvfalp is the fast particle second apscies parallel velocity diffusivity					
			!  dt0 is the time step		
			!  Adens is the fast particle density flatness (if no external profiles)					
			!  Bdens is the location of the fast particle density gradient (if no external profiles)														
			!  Adensalp is the fast particle 2nd species density flatness (if no external profiles)					
			!  Bdensalp is the location of the fast particle 2nd species density gradient (if no external profiles)		
			!  LcA0 Landau closure 1
			!  LcA1 Landau closure 2
			!  LcA2 correction to the fast particle beta
			!  LcA3 correction to the ratio between fast particle thermal velocity and Alfven velocity
			!  LcA0alp Landau closure 1 2nd fast particle species
			!  LcA1alp Landau closure 2 2nd fast particle species
			!  LcA2alp correction to the fast particle beta 2nd fast particle species
            !  LcA3alp correction to the ratio between fast particle thermal velocity and Alfven velocity 2nd fast particle species
			!  ext_prof introduce the experimental profiles	using an external text file		
			!  omegar eigenmode frequency (required to introduce the damping effects)		
			!  iflr thermal ions Larmour radius normalized to the device minor radius
			!  r_epflr fast particle Larmour radius normalized to the device minor radius	
			!  r_epflr 2nd species fast particle Larmour radius normalized to the device minor radius		
			!  dpres electron pressure normalized to the total pressure	(required in the two fluid approximation)	
			!  DIIID_u activates external profiles using TRANSP code units (cm)	
			!  betath_factor thermal beta factor (if 1 the same value than the equilibria)
			!  etascl user defined constant value for the resistivity (if ietaeq=2)
			!  deltaq user defined safety factor profile displacement
			!  deltaiota user defined iota profile displacement
			!  reta,eta0,etalmb user defined resistivity profile (if ietaeq=3)
			!  ihist previous run number
			!  nocpl if zero the run doesn't include couplings
			!  maxstp number of time steps in the run
			!  nstep and nstep1 internal loop time step number
			!  ndump indicates the number of step per eigenfunctions output 
			!  nprint indicates the number of step per output in farprt output
			!  lplots number of modes in the eigenfunctions output 
			!  itime time step normalization option 
			!  noeqn number of equations of the model
			!  edge_p grid point from where the VMEC data is extrapolated
			!  inalp,ivalp,iq,iw,ix1,ix2,iwa,ix1a,ix2a internal variables
			!  iflr_on activates thermal ion FLR effects if 1
			!  epflr_on activates energetic particle FLR effects if 1													
			!  ieldamp_on activates electron-ion Landau damping if 1
			!  twofl_on activates two fluids effects if 1
			!  alpha_on activates a 2nd energetic particle species if 1		
			!  EP_dens_on activates user defined density profiles of energetic particle species if 1
			!  Alpha_dens_on activates user defined density profiles of 2nd species energetic particle species if 1
			!  EP_vel_on user defined fast particle vth/vA0 profile
			!  Alpha_vel_on user defined 2nd species fast particle vth/vA0 profile
			!  eq_name name of the equilibria
			!  ext_prof_name external profile file name
			!  ext_prof_len number of lines in the external profile file
			!  deltaq safety factor displace
			!  deltaiota iota displace
			!  matrix_out create the input of the eigensolver (not available yet)		
			!  leqdim is the number of equilibrium modes
			!  ldim is the total mumber of modes (equilibrium + dynamic)
			!  jdim is the number of radial point
			!  nstres indicates if the run is new (if 0) or a continuation (if 1)
                        !  Eq_vel_on activates external profiles for the equilibrium toroidal velocity if 1
                        !  Eq_velp_on activates external profiles for the equilibrium poloidal velocity if 1
                        !  Eq_Presseq_on activates external profiles for the thermal pressure if 1
                        !  Eq_Presstot_on activates external profiles for the thermal and EP pressures if 1
                        !  Auto_grid_on activates a default grid spacing
                        !  Edge_on activates the VMEC data extrapolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	
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
		real(IDP), dimension(:), allocatable :: r,rinv
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
                                                        pthermal,pep,ptot,vthermalep,vAlfven,vtherm_ionP, &
							dnalpha,dnalphann,talpha,talphann,vzt_eqp,vth_eqp
												                                            												
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			!  mjeq number of radial grids in the equilibrium file		
			!  nfp indicates the number of field periods of the equilibrium
			!  rfar normalized radius in the equilibrium file
			!  qfar safety factor in the equilibrium file
			!  mj last point of the radial grid
			!  lmax total number of modes (input)
			!  leqmax total number of modes (input)
			!  lbmax total number of modes in the equilibrium file
			!  mmin smallest poloidal mode for each toroidal mode family
			!  mmax largest poloidal mode for each toroidal mode family
			!  nmin smallest toroidal mode
			!  nmax largest toroidal mode
			!  mmineq smallest equilibrium poloidal mode for each equilibrium toroidal mode family
			!  mmaxeq largest equilibrium poloidal mode for each equilibrium toroidal mode family
			!  nmineq smallest equilibrium toroidal mode
			!  nmaxeq largest equilibrium toroidal mode
			!  mjm1 is mj-1
			!  mjm2 is mj-2
			!  lmax0 total number of equilibrium poloidal modes
			!  lmaxn total number of dynamic poloidal modes
			!  l0 indicates the (0,0) mode
			!  leq0 indicates the (0,0) equilibrium mode
			!  lhmax indicate the helicity number of dynamic modes
			!  lheqmx indicate the helicity number of equilibrium modes
			!  nnum indicates the toroidal mode order number 
			!  m0dy indicates the number of evolving equilibrium modes (not available, only non linear)
			!  mxmband number of poloidal modes for each toroidal family including both parities
			!  nnd number of n-values with one parity
			!  nst number of n-values with both parity		
			!  mminb smallest poloidal mode in the equilibrium file
			!  mmaxb largest poloidal mode in the equilibrium file
			!  nminb smallest toroidal mode in the equilibrium file
			!  nmaxb smallest poloidal mode in the equilibrium file
			!  mxmbandb number of poloidal modes for each toroidal family including both parities in the equilibrium file
			!  mmaxx absolute value of the largest poloidal mode
			!  mmaxxb absolute value of the largest poloidal mode in the equilibrium file
			!  mbandeq number of equilibrium poloidal modes for each toroidal family including both parities
			!  mm Vector with the dynamic polidal mode numbers
			!  nn Vector with the dynamic toroidal mode numbers
			!  mh Vector with the helicities dynamic polidal mode numbers
			!  nh Vector with the helicities dynamic toroidal mode numbers															
			!  mmeq Vector with the equilibrium polidal mode numbers
			!  nneq Vector with the equilibrium toroidal mode numbers
			!  mheq Vector with the helicities equilibrium polidal mode numbers
			!  nheq Vector with the helicities equilibrium toroidal mode numbers	
			!  mmb Vector with the dynamic polidal mode numbers in the equilibrium file
			!  nnb Vector with the dynamic toroidal mode numbers in the equilibrium file													
			!  mmstart lower bound of m values
			!  mmend upper bounds of m values
			!  mmstartb lower bound of m values in the equilibrium file
			!  mmendb upper bounds of m values in the equilibrium file	
			!  r is the normalized minor radius
			!  rinv is the inverse of the normalized minor radius
			!  dc1m,dc1p,dc2m,dc2p,del2cm,del2cp are rerivate weigths to be used in grid subroutine
			!  wt1m,wt10,wt1p,wt2m,wt20,wt2p are derivative weights to be used in blockj, block0, b2lx and b2lx0 subroutines
			!  ll matrix with the dynamic totoidal and poloidal mode index
			!  llb matrix with the dynamic totoidal and poloidal mode index	in the equilibrium file
			!  lleq matrix with the equilibrium totoidal and poloidal mode index		
			!  m1n lower poloidal mode in the mode vector
			!  mrang rank of the mode vector
			!  ll0 index of the equilibrium poloidal mode 
			!  lln index of the equilibrium and dynamic poloidal mode 
			!  lo index of the equilibrium poloidal mode  if no couplings
			!  llno index of the dynamic poloidal mode in linstart subroutine
			!  signl vector with the sign of the dynamic modes
			!  signleq vector with the sign of the equilibrium modes	
			!  lnumn indicates the total number of equilibrium and dynamic poloidal modes  	
			!  mnumn indicates the total number of equilibrium and dynamic poloidal modes 
			!  ni number of points interior to the island
			!  nis number of points in the island
			!  ne number of points exterior to the island
			!  delta width of the uniform fine grid (island)
			!  rc center of the fine grid (island)
			!  fti fraction of interior point in transition region
			!  fte fraction of exterior points in transition region
			!  bmodn normalized module of the magnetic field
			!  mu0 vacuum magnetic permeability (MKS)
			!  uion fast particle Z number
			!  vthi normalized ion thermal velocity at the magnetic axis
			!  vthe normalized electron thermal velocity at the magnetic axis													
			!  xnuelc0 electron-ion collision FR axis (MKS)
			!  coul_log Coulomb logarithm (Te > 10 eV)								
			!!!!!!!!!!!!!!!!!!!!!! External profiles variable: !!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!  dnnbi normalized fast particle density
			!  dne normalized thermal electron density
			!  dni normalized thermal ion density
			!  pthermal normalized fast particle temperature
			!  ti  normalized thermal ion temperature
			!  te  normalized thermal electron temperature
			!  vzt_eqp normalized equilibrium toroidal rotation (to the Alfven velocity)
			!  vth_eqp normalized equilibrium poloidal rotation (to the Alfven velocity)
			!  qprofile safety factor
			!  dnnbinn fast particle density (MKS)
			!  dnenn thermal electron density (MKS)
			!  dninn thermal ion density (MKS)
			!  tinn thermal ion temperature (eV)
			!  tenn thermal electron temperature (eV)
			!  tbn normalized fast particle temperature
			!  tbnnn fast particle temperature (eV)
			!  pthermalnn normalized thermal pressure
			!  vthermalep thermal velocity of the fast particles (MKS)
			!  vAlfven normalized Alfven velocity
			!  vtherm_ionP normalized thermal ions velocity
			!  dnalpha normalized density of the 2nd fast particle species (MKS)
			!  dnalphann density of the 2nd fast particle species  (MKS)
			!  talpha normalized temperature of the 2nd fast particle species
			!  talphann temperature of the 2nd fast particle species (eV)
			!  pep fast particle pressure
			!  pepnn normalized fast particle pressure
			!  ptot equilibrium + fast particle pressure
			!  ptot normalized equilibrium + fast particle pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
 	end module domain

	module equil
		use param
		implicit none
		save
		real(IDP), dimension(:), allocatable :: qq,qqinv,qqinvp,denseq,denseqr,psieq,chieq,preq,feq,cureq,teeq,nfeq,dnfeqdr,vfova,vfova2, &
		                                vtherm_ion,vzt_eq,vth_eq,vth_eq1,nalpeq,valphaova,valphaova2,dnalpeqdr
		character(len=8), dimension(2) :: ndevice
		integer :: ngeneq,leq
		real(IDP) :: eps,bet0
		real(IDP) :: omcy,bet0_f,bet0_alp,omcyalp,omcyb,rbound,norm_eildump
		real(IDP), dimension(:), allocatable :: rs
		integer, dimension(1:10) :: nstep_count		
		real(IDP), dimension(0:10) :: cnep,ctep,cnfp,cvfp,cvep,cnfpalp,cvfpalp,eqvt,eqvp
		real(IDP), dimension(:,:), allocatable :: grr,grt,gtt,grz,gtz,gzz,sqgi,sqg, &
							  grroj,grtoj,gttoj,grzoj,gtzoj,gzzoj,grrup,grtup,grzup,gttup,gtzup,gzzup, &
							  bmod,bst,jbgrr,jbgrt,jbgtt,jbgrz,jbgtz,omdrprp,omdtprp,omdzprp, &
							  omdr,omdt,omdz,djroj,djtoj,djzoj,dbsjtoj,dbsjzoj,dbsjtbj,dgttr,dgrrt,dgrtt,dgttt,dgrrz, &
							  dgrtz,dgttz,dgrtp,dgttp,jsq,bsgrt,bsgtt,bsq,bsqgtt,lplrr,lplrt,lplrz,lplr,lpltt, &
							  lpltz,lplt,lplz,lplzz,eildr,eildt,eildz,eildrr,eildrt,eildrz,eildtt,eildtz,eildzz, &
							  sqgdroj,sqgdthoj,sqgdztoj,sqgdthojbst,sqgdztojbst,sqgibmodith,sqgibmodizt, &
                                                          test,testr,testt,testrr,testtt,testrt
							   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
			!  qq safety factor
			!  qqinv inve3rse of the safety factor (iota)
			!  qqinvp radial derivate of the safety factor inverse
			!  denseq normalized thermal plasma density
			!  denseqr radial derivate of the normalized thermal plasma density
			!  psieq normalized equilibrium poloidal flux
			!  chieq normalized equilibrium toroidal flux
			!  preq normalized equilibrium thermal plasma pressure
			!  feq normalized equilibrium current density
			!  cureq normalized equilibrium electric currect
			!  teeq normalized equilibrium plasma temperature
			!  nfeq normalized equilibrium fast particle density
			!  dnfeqdr radial derivate of the normalized fast particle density
			!  vfova fast particle thermal velocity normalized to the Alfven velocity
			!  vfova2 is vfova^2
			!  vtherm_ion user defined normalized thermal ion velocity
			!  vzt_eq user defined normalized equilibrium toroidal rotation (to the Alfven velocity)
			!  vth_eq user defined normalized equilibrium poloidal rotation (to the Alfven velocity)
			!  norm_eildump factor in the terms that include the electron-ion Landau damping effects
			!  nalpeq normalized equilibrium density of the second fast particle species	
			!  valphaova Ffst particle thermal velocity normalized to the Alfven velocity of the second fast particle species
			!  valphaova2 is vfalphaova^2
			!  dnalpeqdr radial derivate of the normalized density of the second fast particle species
			!  ndevice type of device
			!  ngeneq type of equilibria input (if 1 VMEC)
			!  leq dummy variable 
			!  ndat indicates the number of radial points in several subroutines
			!  eps ratio of the minor and major radius
			!  bet0 thermal beta
			!  omcy cyclotron frequency of the fast particles (normalized to the Alfven time) 
			!  omcyalp cyclotron frequency of the 2nd species of fast particles (normalized to the Alfven time) 
			!  omcyb bounce frequency of the helically trapped energetic particles (normalized Alfven time)
			!  rbound bounce distance of the helically trapped energetic particles (normalized to the major radius)
			!  bet0_f fast particle beta 
			!  bet0_alp 2nd species fast particle beta 
			!  rs eigenfunction radial width
			!  xr dummy value for the radial point in several subroutines
			!  xdat dummy value
			!  nstep_count time step screen out
			!!!!!!!!!!!!!!!!!!!!!!! Expansion elements of the user defined profiles  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			!  cnep thermal plasma density
			!  ctep thermal plasma temperature
			!  cnfp fast particle density
			!  cvfp fast particle parallel velocity
			!  cvep thermal ions parallel velocity
			!  cnfpalp 2nd species fast particle density
			!  cvfpalp 2nd species fast particle parallel velocity
			!  eqvt normalized equilibrium toroidal rotation (to the Alfven velocity)
			!  eqvp normalized equilibrium poloidal rotation (to the Alfven velocity)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! metric elements  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!  grr,grt,gtt,grz,gtz,gzz covariant metric elements with r=radial, t=poloidal, z=toroidal
			!  sqgi inverse of the metric Jacobian
			!  sqg  metric Jacobian 
			!  grroj,grtoj,gttoj,grzoj,gtzoj,gzzoj  covariant metric elements divided by the Jacobian
			!  grrup,grtup,grzup,gttup,gtzup,gzzup contravariant metric elements
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! variable defined in VMEC subroutine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			!  bmod normalized magnetic field module (to the value in the mag. axis)
			!  bst = -rinv*(d(preq)/dr)*bet0*sd2*sqg(:,l)/(2*(mmb(l)*qqinv-nnb(l)))
			!  jbgrr = sqg*grr
			!  jbgrt = sqg*grt
			!  jbgtt = sqg*gtt
			!  jbgrr = sqg*grr
			!  jbgrz = sqg*grz
			!  jbgtz = sqg*gtz
			!  omdr,omdt,omdz Omega d operator (r, th and zt components)
			!  omdrprp,omdt[prp,omdzprp Omega d operator helically trapped particles (r, th and zt components)
			!  djroj = sqgi*[ d/dr (sqg) ]
			!  djtoj= (sqgi/rho)*[ d/dth (sqg) ]
			!  djzoj= sqgi*[ d/dzt (sqg) ]
			!  dbsjtoj= sqgi*[ d/dth (sqg*bst ) ]
			!  dbsjzoj= sqgi*[ d/dzt (sqg*bst ) ]
			!  dbsjtbj= sqg*[ d/dzt (sqg*bst ) ]
			!  dgttr= sqg*[ d/dr (gtt) ]
			!  dgrrt= (sqg/rho)*[ d/dth (grr) ]
			!  dgrtt= (sqg/rho)*[ d/dth (grt) ]
			!  dgttt= (sqg/rho)*[ d/dth (gtt) ]
			!  dgrrz= sqg*[ d/dzt (grr) ]
			!  dgrtz= sqg*[ d/dzt (grt) ]
			!  dgttz= sqg*[ d/dzt (gtt) ]
			!  dgrtp= sqgi*[ (d/dzt - qqinv*d/dth)(grt) ]
			!  dgttp= sqgi*[ (d/dzt - qqinv*d/dth)(gtt) ]
			!  jsq= sqg*sqg
			!  bsgrt= bst*grt*sqgi
			!  bsgtt= bst*gtt*sqgi
			!  bsq= bst*bst
			!  bsqgtt= bst* bst*gtt*sqgi
			!  sqgdroj = sqgi*d(sqg)/dr
			!  sqgdthoj = (sqgi/rho)*d(sqg)/dth
			!  sqgdztoj = sqgi*d(sqg)/dzt
			!  sqgdztojbst = bst*sqgdztoj
			!  sqgdthojbst = bst*sqgdthoj
			!  sqgibmodith = bmod*sqg*[ d/dth (sqgi) ] / rho 
			!  sqgibmodizt = bmod*sqg*[ d/dtz (sqgi) ] 
			!  eildrr,eildrt,eildrz,eildtt,eildtz,eildzz,eildr,eildt,eildz electron-ion Landau damping terms
			!  lplrr,lplrt,lplrz,lplr,lpltt,lpltz,lplt,lplz,lplzz perpendicular gradient operator terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
	end module equil	

	module dynamo
		use param
		implicit none
		save
		integer :: ietaeq,ipert
		real(IDP) :: dt,dtd2,time,s,gamma,pertscl		
		real(IDP), dimension(:), allocatable :: eta,etann
		real(IDP), dimension(:), allocatable :: widthi,gammai		
		real(IDP), dimension(:,:), allocatable :: psi,phi,pr,uzt,nf,vprlf,vthprlf,nalp,vprlalp,plot_var
		real(IDP), dimension(:,:,:), allocatable :: cmamm,cmamp,cmapm,cmapp
		real(IDP), dimension(:,:,:), allocatable :: amat,bmat,cmat,amatw,bmatw,cmatw,amatwalp,bmatwalp,cmatwalp
		real(IDP), dimension(:,:), allocatable :: xt,yt,xw
		integer, dimension(:,:), allocatable :: ipc,ipcw,ipcwalp
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
			!  ietaeq index to select the type of eta (if 1 the electron temperature is used)
			!  ipert index of the perturbation type
			!  dt time step
			!  dtd2 = dt/2
			!  time indicates the simulation time
			!  s magnetic Lundquist number
			!  gamma polytropic index
			!  pertscl perturbation scale factor
			!  eta magnetic diffusivity
			!  widthi size of the perturbation
			!  gammai perturbation		
			!  cmamm,cmamp,cmapm,cmapp sign matrix 
			!  amat,bmat,cmat tridiagonal matrix where the model terms are stored using block subroutine
                        !  amatw,bmatw,cmatw,amatwalp,bmatwalp,cmatwalp tridiagonal matrix where the model terms are stored using block subroutine
			!  xt,yt,xw matrix where the model terms are stored using b2lx subroutine
			!  ipc,ipcw,ipcwalp dummy matrix used in solbt subroutine
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Normalized fluctuating variable  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
			!  psi poloidal flux (norm minor radius^2 * mag. axis magnetic field)
			!  phi electrostatic potential (norm minor radius^2 * mag. axis magnetic field / resistive time)
			!  pr pressure (norm pressure in the mag. axis)
			!  uzt vorticity (toroidal component, norm Alfven time)
			!  nf fast particle density (norm fast particle density in the mag. axis)
			!  vprlf fast particle parallel velocity (norm to the Alfven time / major radius)
			!  vthprlf thermal parallel velocity (norm to the Alfven time / minor radius)
			!  nfalp density of the second fast particle species (norm 2ndfast particle density in the mag. axis)
			!  vprlfalp parallel velocity of the second fast particle species (norm to the Alfven time / major radius)		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
		
	end module dynamo

	module scratch
		use param
		implicit none
		save
		real(IDP), dimension(:,:), allocatable :: sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10,sc11
 		real(IDP), dimension(:,:), allocatable :: sceq1,sceq2,sceq3,sceq4,sceq5,sceq6,sceq7
		real(IDP), dimension(:), allocatable :: sd1,sd2,sd3,sd4,sd5,sd6
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			!  sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8 dummy dynamic variables with radial and angular dependency
			!  sceq1,sceq2,sceq3,sceq4,sceq5,sceq6,sceq7 dummy equilibrium variables with radial and angular dependency
			!  sd1,sd2,sd3,sd4,sd5,sd6 dummy variables with only radial dependency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
	end module scratch

	module findrf
		use param
		implicit none
		save
		integer :: n0,nt
		real(IDP) :: d,dx
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			!  n0,nt,d,dx parameters of the subroutine findr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
		
	end module findrf
	program far

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none
		character(len=10) :: timew,datew,zone_h
		integer, dimension(8) :: values_s,values_e
		character(len=8) :: confil
		integer :: i,l,j,iend,count
		real(IDP) :: d,dk
		character*1 :: cdum0
		character(len=132) :: char5			
                character*8 :: arg1
                character*2 :: arg2

		interface
			subroutine inputlist
			end subroutine inputlist		
			subroutine rddump
			end subroutine rddump
		end interface
		rdump_file_name="fs0000z"
		variable = 2
		count = command_argument_count()
		if(count .eq. 1) then
		 call getarg(1, arg1)
                 rdump_file_name=trim(adjustl(arg1))
		else if(count .eq. 2) then
		 call getarg(1, arg1)
                 rdump_file_name=trim(adjustl(arg1))
                 call getarg(2, arg2)
                 read(arg2,'(i2)') variable
		end if
!		write(0,'(/)')
		write(0,'(" READING FROM ",a8," dump-file")')rdump_file_name
		if(variable .eq. 1) write(0,'("  PLOTTING 2D psi contours")')
		if(variable .eq. 2) write(0,'("  PLOTTING 2D phi contours")')
		if(variable .eq. 3) write(0,'("  PLOTTING 2D pressure contours")')
		if(variable .eq. 4) write(0,'("  PLOTTING 2D n_fast contours")')
		if(variable .eq. 5) write(0,'("  PLOTTING 2D v_fast_parallel contours")')
		if(variable .eq. 6) write(0,'("  PLOTTING 2D vth_parallel contours")')
!		Input read	 
!		write(0,'(" ====> Checking input list ... ")')		
	
		call inputlist	

!		write(0,'(" ====> Input list check DONE !! ")')		
		
		lmax=ldim
		leqmax=leqdim	
		mj=jdim
		write(0,*) mj
		call rddump	
        end program far


	subroutine inputlist
	
		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

!		namelist/nam/maxstp,ndump,nprint,lplots,itime,dt0,nonlin,mj,lmax,leqmax,mm,nn,mmeq,nneq,Adens,Bdens,ni,nis,ne, &
!			     delta,rc,fti,fte,stdifp,stdifu,ngeneq,eps,bet0,rs,leq,tfl_Un,tfl_Psin,tfl_pn, &
!			     ndat,xr,xdat,etascl,reta,eta0,etalmb,ietaeq,LcA0,LcA1,LcA2,LcA3,ext_prof,omegar, &
!				 iflr,r_epflr,iflr_on,s,gamma,ipert,widthi,gammai,pertscl,m0dy,nocpl,cnep,ctep,cvep,omcy,bet0_f,cnfp,cvfp,stdifnf,stdifvf,dpres, &
!			     ieldamp_on,epflr_on,twofl_on,LcA0alp,LcA1alp,LcA2alp,LcA3alp,alpha_on,stdifnfalp,stdifvfalp, &
!				 Adensalp,Bdensalp,bet0_falp,cnfpalp,cvfpalp,omcyalp,DIIID_u,r_epflralp,EP_dens_on,Alpha_dens_on,ext_prof_len,betath_factor, &
!                EP_vel_on,Alpha_vel_on,deltaq,deltaiota,matrix_out,leqdim,ldim,jdim,nstres,etascl,reta,eta0,etalmb, &	
!                                            q_prof_on, Trapped_on,omcyb,rbound			 
		character*1 :: cdum0
		character(len=132) :: char5		
		integer :: i
		character(len=8) :: confil	
		real(IDP) :: widthix	

!		Read the input files
		
		open (unit=5,file="Input_Model",status="old",form="formatted")
!		open (unit=6,file="farprt",status="unknown")
		
!		write(6,'(" copy of Input_Model:")')
		
		do
			read(5,'(a)',end=20) char5
!			write(6,*) char5
		end do
20    rewind(5)		
		
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
		
		if (nstres /= 0) then

			confil="fs"//numruno(1)//numruno(2)//numruno(3)
			open(unit=8,file=confil,status='unknown',form='unformatted')

			read(8) ihist
			rewind(8)

		end if	
		
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
		allocate (plot_var(0:jdim,0:ldim))		!	Normalized fluctuating parallel velocity of the second fast particle species
		
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
		read(5,'(a1)') cdum0		
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
		
!		close(6)
		
	end subroutine inputlist

!*************************************
	subroutine rddump

		use param
		use cotrol
		use domain
		use equil
		use dynamo
!DIR$ DEFINE PLPLOT_WORKING=1
!DIR$ IF(PLPLOT_WORKING .eq. 1)
                USE plplot
!DIR$ ENDIF
		use bspline
		implicit none
                include "/Users/dsp/fortran_code_development/silo-4.10.2/include/silo.inc"

		integer :: i,l,j,lmaxo,m,n,ihisto,l1,icycl,icol,ic,mjbz
		integer, dimension(ldim) :: mmo,nno
		integer, dimension(ldim) :: lmap
		integer, allocatable, dimension(:) :: color_set
		real*8, allocatable, dimension(:) :: yy,phi_all,rr,temp_all
		real(IDP) :: pertsclo
		real(IDP), dimension(ldim) :: widthio,gammaio
		real(IDP) :: ymin, ymax, scale
		character(len=2), dimension(3) :: numrunp
		character(len=8) :: confil
				
            real*8, dimension(:,:), allocatable :: rmnc, zmns, phimns
            real*8, dimension(:,:), allocatable :: rmns, zmnc, phimnc
            integer, dimension(:), allocatable :: mmbz,nnbz
            integer :: icc, lp, itht, ierr, ndims, npts, dbfile, err
            real, dimension(:), allocatable :: x, z, tht_p, phi_plot, rad, xouter, zouter
            real, dimension(:), allocatable :: flux_surfs,yp,temp,qqz,qcontour
	    logical :: i_exist, lasym
	    real*8 :: twopi, arg, zeta
            real*8, dimension(:), allocatable :: bcoef,xfit,yfit,xknot
            real*8, dimension(:,:), allocatable :: plot_var_intrp
            real*8 :: nw_yfit, xval
            integer :: kxord
	    twopi = 8.*atan(1.0)
      
!

!		The data of the previous run is read to continue the simulation		

		confil="fs0000z"		
		open(unit=8,file=rdump_file_name,status='unknown',form='unformatted',POSITION="REWIND")

		read(8) ihist
		rewind(8)
		ihist=ihist+1
		
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
		allocate (jbgrz(0:mj,0:leqmax))
		allocate (jbgtz(0:mj,0:leqmax))
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
		jbgrz=0.0_IDP
		jbgtz=0.0_IDP
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
                                 ((jbgrz(j,l),j=0,mj),l=1,leqmax),((jbgtz(j,l),j=0,mj),l=1,leqmax), &
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
		
!		write(6,'("  rddump: have read fs",2a2,a1)') (numruno(i),i=1,3)
		
		do l=1,leqmax
			m=mmeq(l)
			n=nneq(l)
			sgnleq(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) sgnleq(l)=-1
			if (n == 0 .and. m == 0) sgnleq(l)=0
		end do
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			signl(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) signl(l)=-1
			if (n == 0 .and. m == 0) signl(l)=0
		end do
!
!    Select field for subsequent plots:
!		
		 do j = 0,mj
		 do l = 1,lmaxo
		 if(variable .eq. 1) then
		   plot_var(j,l) = psi(j,l)
		 else if(variable .eq. 2) then
		   plot_var(j,l) = phi(j,l)
		 else if(variable .eq. 3) then
		   plot_var(j,l) = pr(j,l)
		 else if(variable .eq. 4) then
		   plot_var(j,l) = nf(j,l)
		 else if(variable .eq. 5) then
		   plot_var(j,l) = vprlf(j,l)
		 else if(variable .eq. 6) then
		   plot_var(j,l) = vthprlf(j,l)
		 end if
		 end do
		 end do
		 
!DIR$ IF(PLPLOT_WORKING .eq. 1)		
      allocate(color_set(9))
      allocate(yy(mj))
      allocate(rr(mj))
      allocate(phi_all(mj*lmaxo))
      ic = 0
      do l = 1,lmaxo
       do j=1,mj
        ic = ic + 1
        phi_all(ic) = plot_var(j,lmap(l))
       end do
      end do 
      ymax = maxval(phi_all); ymin = minval(phi_all)
      scale = max(abs(ymax),abs(ymin))
      ymax = 1.1*ymax/scale; ymin = 1.1*ymin/scale
      do j=1,mj
       rr(j) = r(j)
      end do     
      color_set(1) = 15;color_set(2) = 1;color_set(3) = 3
      color_set(4) = 4;color_set(5) = 8;color_set(6) = 9
      color_set(7) = 10;color_set(8) = 12;color_set(9) = 13
      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plsfam(1,1,200000)
      call plssub(1,1)
      call plinit
      call plscolor(1)
      call plcol0(15)
      call plfont(2)
      call plenv(0.d0, 1.d0, ymin, ymax, 0, 0)
      call plschr(0.d0,0.7d0)
      call pllab('<r>/<a>','Phi_mn','FAR3D Phi eigenmode')
      do l = 1,lmaxo
       if(nno(l) .ne. 0) then
!	if(l .eq.lmaxo/2) write(0,*) mm(l),nn(l)
!	if(l .eq.lmaxo/2) write(0,*) mm(lmap(l)),nn(lmap(l))
       do j=1,mj
        yy(j) = plot_var(j,lmap(l))/scale
!	if(l .eq.lmaxo/2) write(0,*) rr(j),yy(j)
       end do
       icycl = mod(l,8) + 1
       icol = color_set(icycl)
       call plcol0(icol)
       call plwidth(3.d0)
       call plline(rr,yy)
       call plwidth(1.d0)
       call plcol0(15)
      end if
      end do
      call plend
!DIR$ ENDIF

!      Generate Visit 2D plot
!
!      Check to see if geom_bzr_to_cyl file is available, and if so read it in.
!      If not present, then stop.
!	
        i_exist = .false.
	inquire(file='geom_bzr_to_cyl.dat',exist=i_exist)
	if(.not.i_exist) then
	 write(0,'("Need geom_bzr_to_cyl.dat file to be present")')
	 stop 88
	end if
	open(unit=30,file='geom_bzr_to_cyl.dat',status='old')
	read(30,*) mjbz, lbmax, lasym
	 allocate (rmnc(lbmax,mjbz))
	 allocate (zmns(lbmax,mjbz))
	 allocate (phimns(lbmax,mjbz))
	 allocate (qqz(mjbz))
	if(lasym) then
	 allocate (rmns(lbmax,mjbz))
	 allocate (zmnc(lbmax,mjbz))
	 allocate (phimnc(lbmax,mjbz))
	endif
	allocate (rad(mjbz))
	allocate (mmbz(lbmax))
	allocate (nnbz(lbmax))
	do l=1,lbmax
          read(30,*) mmbz(l),nnbz(l)
	end do
	do j=1,mjbz
	  read(30,*) rad(j), qqz(j)
	  write(0,*) rad(j), qqz(j)
	end do
	do l=1,lbmax
	 do j=1,mjbz
	if(lasym) then
          read(30,'(e15.8,5(2x,e15.8))') rmnc(l,j),zmns(l,j),phimns(l,j),  &
     	 rmns(l,j), zmnc(l,j),phimnc(l,j)
	else
          read(30,'(e15.8,2(2x,e15.8))') rmnc(l,j),  &
                      zmns(l,j),phimns(l,j)
	endif
	 end do		 
	end do
	close(unit=30)	
	rmnc(1,1) = rmnc(1,2)
		write(0,*) mj,mjbz
!
!       Generate flux surfaces in 2D cylindrical geometry and calculate field values
!
!       Normalize
		 ic = 0
		 allocate(temp_all(mj*lmaxo))
		 do j = 1,mj
		 do l = 1,lmaxo
		  ic = ic + 1
		  temp_all(ic) = plot_var(j,l)
		 end do
		 end do
		 
                 ymax = maxval(temp_all); ymin = minval(temp_all)
                 scale = max(abs(ymax),abs(ymin))
		 
		 do j = 1,mj
		 do l = 1,lmaxo
		  plot_var(j,l) = plot_var(j,l)/scale
		 end do
		 end do

!      Interpolate FAR3d data in plot_var onto the same radial grid that geom_bzr_to_cyl.dat
!       was calculated on. i.e., this allows
!	
       kxord = 3
       allocate(plot_var_intrp(mjbz,lmaxo))
       allocate(xfit(mj),yfit(mj))
       allocate(xknot(mj+kxord),bcoef(mj))
       do l=1,lmaxo
	do j=1,mj
	 yfit(j) = plot_var(j,l)
	 xfit(j) = r(j)**2
	end do
	call dbsnak(mj,xfit,kxord,xknot)
	call dbsint(mj,xfit,yfit,kxord,xknot,bcoef)
	do j=1,mjbz
	 xval = rad(j)
	 if(xval .gt. xfit(1) .and. xval .lt. xfit(mj)) then
	  nw_yfit = dbsval(xval,kxord,xknot,mj,bcoef)
	 else
	  nw_yfit = 0.
	 end if
	 plot_var_intrp(j,l) = nw_yfit
	end do
       end do

		write(0,*) mj,mjbz
		itht = 600; nfp = 1
		allocate(xouter(itht), zouter(itht))
		zeta = 0.
		allocate(x(itht*mjbz),z(itht*mjbz),phi_plot(itht*mjbz),flux_surfs(itht*mjbz))
		allocate(qcontour(itht*mjbz))
		allocate(tht_p(itht))
		x(:) = 0.
		z(:) = 0.
		phi_plot(:) = 0.
		icc = 0
		do j = 1,mjbz
		 do i = 1,itht
		  icc = icc + 1
		  tht_p(i) = twopi*real(i-1)/real(itht-1)
                  x(icc) = 0.; z(icc) = 0.
		  
		  do l = 1,lbmax
		  
		   if(lasym) then
		    arg = tht_p(i)*real(mmbz(l)) - zeta*real(nnbz(l))
		    x(icc) = x(icc) + rmnc(l,j)*cos(arg) + rmns(l,j)*sin(arg)
		    z(icc) = z(icc) + zmns(l,j)*sin(arg) + zmnc(l,j)*cos(arg)

     		   else
		    arg = tht_p(i)*real(mmbz(l)) - zeta*real(nnbz(l))
		    x(icc) = x(icc) + rmnc(l,j)*cos(arg)
		    z(icc) = z(icc) + zmns(l,j)*sin(arg)
     
                   endif
		   
		  end do
		  
		  if(j .eq. mjbz) then
		    xouter(i) = x(icc)
		    zouter(i) = z(icc)
		  endif
		  
		 phi_plot(icc) = 0.
		 		   		   
		 do l = 1,lmaxo
		  arg = real(mmo(l))*tht_p(i)-zeta*real(nno(l))
		   if(signl(l) .eq. 1) then
		     phi_plot(icc) = phi_plot(icc) + plot_var_intrp(j,lmap(l))*sin(arg)
		    else if(signl(l) .eq. -1) then
		     phi_plot(icc) = phi_plot(icc) + plot_var_intrp(j,lmap(l))*cos(arg)
		    endif
		   end do
		   flux_surfs(icc) = real(rad(j))/real(rad(mjbz))
		   qcontour(icc) = qqz(j)
		   
		 end do
		end do

!    Write data out to Visit Silo data file
!
      if(variable .eq. 1) then
       err = dbcreate("2D_psi.silo", 11, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
      else if(variable .eq. 2) then
       err = dbcreate("2D_phi.silo", 11, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
      else if(variable .eq. 3) then
       err = dbcreate("2D_pressure.silo", 16, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
      else if(variable .eq. 4) then
       err = dbcreate("2D_nfast.silo", 13, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
     else if(variable .eq. 5) then
       err = dbcreate("2D_vprlf.silo", 13, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
      else if(variable .eq. 6) then
       err = dbcreate("2D_vthprl.silo", 14, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
      end if
	if(err .ne. 0) stop 27
        ndims = 2; npts = itht*mjbz
       err = dbputpm(dbfile, "pointmesh",9, ndims, x, z, DB_F77NULL,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 28
	
       if(variable .eq. 1) then
        err=dbputpv1(dbfile,"psi",3,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       else if(variable .eq. 2) then
        err=dbputpv1(dbfile,"phi",3,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       else if(variable .eq. 3) then
        err=dbputpv1(dbfile,"pres",4,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       else if(variable .eq. 4) then
        err=dbputpv1(dbfile,"nfast",5,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       else if(variable .eq. 5) then
        err=dbputpv1(dbfile,"vprl",4,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       else if(variable .eq. 6) then
        err=dbputpv1(dbfile,"vthprl",6,"pointmesh",9,phi_plot,npts,DB_FLOAT,DB_F77NULL,ierr)
       end if
       
	if(ierr .ne. 0) stop 29
       err = dbputpv1(dbfile, "flux", 4, "pointmesh", 9, flux_surfs,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 30
       err = dbputpv1(dbfile, "q", 1, "pointmesh", 9, qcontour,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 30
      err = dbclose(dbfile)
      
      open(unit=8,file="bound.curve",status="unknown")
      do i=1,itht
       write(8,*) xouter(i), zouter(i)
      end do
      close(8)

	end subroutine rddump
