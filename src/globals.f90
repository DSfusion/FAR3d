module param
  implicit none
  integer, parameter :: IDP = kind(1.d0)
end module param

module processor
  implicit none
  save
  integer :: myPE,numPEs,numPEsm1,numPElm1,numdev
end module processor

module var_para
  use param
  implicit none
  save
  integer :: n_start,n_end
  integer :: mj_start,mj_end
  integer, dimension(:), allocatable :: n_st,n_nd,lmaxPE
  integer, dimension(:), allocatable :: mj_br,mj_inc,mj_st,mj_dl
  real(IDP), dimension(:,:), allocatable :: scp,scu
end module var_para

module cotrol
  use param
  implicit none
  save
  character(len=8), dimension(:), allocatable :: numhist
  character(len=2), dimension(3) :: numrun,numruns
  character(len=2), dimension(3):: numruno
  real(IDP) :: stdifp,stdifu,stdifnf,stdifvf,stdifv,stdifnalp,stdifvalp,dt0,LcA0,LcA1,LcA2,LcA3,dpres,omegar, &
       iflr,r_epflr,LcA0alp,LcA1alp,LcA2alp,LcA3alp,r_epflralp
  real(IDP), dimension(:), allocatable :: stdifpn,stdifun,stdifnfn,stdifvfn,stdifvn,stdifnalpn,stdifvalpn, &
       fctr_dif,dfctr_difdr
  integer :: noeqn,nvar,ihist,nocpl,maxstp,nstep,ndump,nprint,ndiag,lplots,itime,ext_prof, &
       nstep1,nonlin,twofl_on,epflr_on,alpha_on,iflr_on,ieldamp_on,ivalp, &
       iq,iw,ix1,ix2,iwa,ix1a,ix2a,nstres,trapped_on,DIIID_u,difnr_on, &
       nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on,src_sink_th_on,src_sink_EP1_on,src_sink_EP2_on, &
       src_sink_DIIID_on,src_sink_ITER_on, &
       edge_p,inalp,EP_dens_on,Alpha_dens_on,EP_vel_on,Alpha_vel_on,q_prof_on,Eq_vel_on, &                     
       Eq_velp_on,Eq_Presseq_on,Eq_Presstot_on,Auto_grid_on,Edge_on,B_par_on
  real(IDP) :: rsrc,wsrc,asrc,rsrc_EP1,wsrc_EP1,asrc_EP1,rsrc_EP2,wsrc_EP2,asrc_EP2,AWfctr,Nfctr, &
       AWfctr_dif,Rfctr,Wfctr
  character(len=40) :: eq_name,ext_prof_name
  logical :: matrix_out, old_rd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  numhist is the number of runs in the series prior to this run
  !  numruns is the first run number in the series
  !  numrun is the present run number
  !  numruno is the previous run number
  !  old_rd = .true. when the fs##### dump file for numruno is older than this version of the code
  !  stdifp is the pressure diffusivity
  !  stdifu is the vorticity diffusivity
  !  stdifv is the parallel velocity diffusivity  
  !  stdifnf is the fast particle density diffusivity
  !  stdifvf is the fast particle parallel velocity diffusivity  
  !  stdifnalp is the fast particle second species density diffusivity
  !  stdifvalp is the fast particle second apscies parallel velocity diffusivity
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
  !  ext_prof introduce the experimental profiles using an external text file  
  !  omegar eigenmode frequency (required to introduce the damping effects)  
  !  iflr thermal ions Larmor radius normalized to the device minor radius
  !  r_epflr fast particle Larmour radius normalized to the device minor radius
  !  r_epflr 2nd species fast particle Larmour radius normalized to the device minor radius  
  !  dpres electron pressure normalized to the total pressure (required in the two fluid approximation)
  !  DIIID_u activates external profiles using TRANSP code units (cm)
  !  betath_factor thermal beta factor (if 1 the same value than the equilibria)
  !  etascl user defined constant value for the resistivity (if ietaeq=2)
  !  deltaq user defined safety factor profile displacement
  !  deltaiota user defined iota profile displacement
  !  reta,eta0,etalmb user defined resistivity profile (if ietaeq=3)
  !  ihist previous run number
  !  nocpl if one cylindrical geometry (no toroidal or helical couplings)
  !  maxstp number of time steps in the run
  !  nstep and nstep1 internal loop time step number
  !  ndump indicates the number of step per eigenfunctions output 
  !  nprint indicates the number of step per output in farprt output
  !  lplots number of modes in the eigenfunctions output 
  !  itime time step calculation option 
  !  noeqn number of equations of the model
  !  edge_p grid point from where the VMEC data is extrapolated
  !  inalp,ivalp,iq,iw,ix1,ix2,iwa,ix1a,ix2a internal variables
  !  difnr_on enables diffusion dependency on n and/or rho if 1
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
  !  Auto_grid_on activates grid with evenly spaced points
  !  Edge_on activates the VMEC data extrapolation
  !  vtk_on activates the vtk format in the output
  !  spe1 is the ion mass to proton mass ratio for the first EP population
  !  spe2 is the ion mass to proton mass ratio for the second EP population
  !  nopsievol_on neglects psi nonlinear terms if 1
  !  noprevol_on neglects thermal pressure nonlinear terms if 1
  !  nonfevol_on neglects EP density nonlinear terms if 1
  !  nonalpevol_on neglects second population EP density nonlinear terms if 1
!!!!!!!!!!!!!!!!!!!!!! Source and sinks: !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  src_sink_th_on activate source and sinks in the thermal plasma if 1
  !  src_sink_EP1_on activate source and sinks in the first EP population if 1 
  !  src_sink_EP2_on activate source and sinks in the second EP population if 1 
  !  src_sink_DIIID_on activate source and sinks DIII-D default case if 1
  !  src_sink_ITER_on activate source and sinks ITER default case if 1
  !  rsrc radial location of the source in thermal plasma
  !  wsrc width of the source in thermal plasma
  !  asrc intensity of the source in thermal plasma
  !  rsrc_EP1 radial location of the source in the first EP population
  !  wsrc_EP1 width of the source in the first EP population
  !  asrc_EP1 intensity of the source in first EP population
  !  rsrc_EP2 radial location of the source in the sceond EP population
  !  wsrc_EP2 width of the source in the second EP population
  !  asrc_EP2 intensity of the source in second EP population

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module cotrol

module domain
  use param
  implicit none
  save
  integer :: mjeq,nfp
  real(IDP), dimension(:), allocatable :: rfar,qfar
  integer :: mj,lmax,lmaxo,leqmax,lbmax,mmin,mmax,nmin,nmax,mmineq,mmaxeq,nmineq,nmaxeq,mjm1,mjm2,mjp1,lmax0,lmaxn,lmx,l0, &
       leq0,lhmax,lheqmx,nnum,m0dy,mxmband,mbandeq,mnxxx,mminb,mmaxb,nminb,nmaxb,mxmbandb,mmaxx,mmaxxb,mmaxxeq
  integer, dimension(:), allocatable :: mm,nn,mh,nh
  integer, dimension(:), allocatable :: mmeq,nneq,mheq,nheq
  integer, dimension(:), allocatable :: mmb,nnb
  integer, dimension(:), allocatable :: mmstart,mmend,mmstreq,mmendeq
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
  integer, dimension(:), allocatable :: signl,jsl
  integer, dimension(:), allocatable :: sgnleq
  integer, dimension(:), allocatable :: lnumn
  integer, dimension(:), allocatable :: mnumn,nskpn,nskpxn,nskpin,nskp2n
  integer :: ni,nis,ne
  real(IDP) :: pi,delta,rc,fti,fte
  LOGICAL :: lasym

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
  !  llb matrix with the dynamic totoidal and poloidal mode index in the equilibrium file
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module domain

module equil
  use param
  implicit none
  save
  real(IDP), dimension(:), allocatable :: qq,qqinv,qqinvp,denseq,denseqr,preq,feq,cureq,teeq,tieq,nfeq,dnfeqdr, &
       dpreqdr,vfova,vfova2,vzt_eq,vth_eq,vtherm_elc,nalpeq,valphaova,valphaova2,dnalpeqdr,vtherm_ion
  character(len=8), dimension(2) :: ndevice
  integer :: spe1,spe2
  real(IDP) :: bigrn,bmodn,eps,bet0,xle,omcy,omcyalp,bet0_f,bet0_alp,Adens,Bdens,Adensalp,Bdensalp, &
       uion,vthi,vthe,xnuelc0,coul_log,omcyb,rbound,norm_eildump,betath_factor,deltaiota,deltaq  
  real(IDP), dimension(:), allocatable :: dnnbi,dne,dni,ti,te,temp_ep,qprofile, &
       dnnbinn,dnenn,dninn,pthermalnn,tinn,tenn,tbn,tbnnn,pepnn,ptotnn, &
       pthermal,pep,ptot,vthermalep,vAlfven,vtherm_ionP,vtherm_elecP, &
       dnalpha,dnalphann,talpha,talphann,vzt_eqp,vth_eqp,temp_epnn
  real(IDP), dimension(:), allocatable :: rs
  real(IDP), dimension(0:10) :: cnep,ctep,cvep,cnfp,cvfp,cnfpalp,cvfpalp,eqvt,eqvp,srcsinkth,srcsinkEP1,srcsinkEP2
  real(IDP), dimension(:,:), allocatable :: R_eq,Z_eq,Phimns_eq,grr,grt,gtt,sqgi,sqg,bmod,grroj,grtoj,gttoj,bst, &
       jbgrr,jbgrt,jbgtt,lplr,lplt,lplz,djroj,djtoj,djzoj,omdr,omdt,omdz,dbsjtoj,dbsjzoj,dbsjtbj, &
       dgttr,dgrrt,dgrtt,dgttt,dgrrz,dgrtz,dgttz,dgrtp,dgttp,jsq,bsgrt,bsgtt,bsq,bsqgtt,lplrr,lplrt,lplrz, &
       lpltr,lpltt,lpltz,lplzr,lplzt,lplzz,eildr,eildt,eildz,eildrr,eildrt,eildrz,eildtt,eildtz,eildzz, &
       omdrprp,omdtprp,omdzprp

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
  !  eps ratio of the minor and major radius
  !  bet0 thermal beta
  !  bmodn normalized module of the magnetic field
  !  uion fast particle Z number
  !  vthi normalized ion thermal velocity at the magnetic axis
  !  vthe normalized electron thermal velocity at the magnetic axis
  !  xnuelc0 electron-ion collision FR axis (MKS)
  !  coul_log Coulomb logarithm (Te > 10 eV)
  !  omcy cyclotron frequency of the fast particles (normalized to the Alfven time) 
  !  omcyalp cyclotron frequency of the 2nd species of fast particles (normalized to the Alfven time) 
  !  omcyb bounce frequency of the helically trapped energetic particles (normalized Alfven time)
  !  rbound bounce distance of the helically trapped energetic particles (normalized to the major radius)
  !  bet0_f fast particle beta 
  !  bet0_alp 2nd species fast particle beta 
  !  rs eigenfunction radial width
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
  !  srcsinkth sink and source thermal plasma
  !  srcsinkEP1 sink and source first EP population
  !  srcsinkEP2 sink and source second EP population
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
  !  bsqgtt= bst*bst*gtt*sqgi
  !  bsjgrt= bst*sqg*grt
  !  bsjgtt= bst*sqg*gtt
  !  bstg= bst*sqg*sqg
  !  bsqg= bst*bst*sqg*sqg
  !  eildrr,eildrt,eildrz,eildtt,eildtz,eildzz,eildr,eildt,eildz electron-ion Landau damping terms
  !  lplrr,lplrt,lplrz,lplr,lpltt,lpltz,lplt,lplz,lplzz,lplr_r,lplt_r perpendicular Laplacian terms

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

end module equil

module dynamo
  use param
  implicit none
  save
  integer :: ietaeq,ipert,idump,idiag
  real(IDP) :: dt,dtd2,time,etascl,reta,eta0,etalmb,s,gamma,pertscl
  real(IDP), dimension(:), allocatable :: eta,etann
  real(IDP), dimension(:,:), allocatable :: psi,phi,pr,uzt,nf,vprlf,vthprlf,nalp,vprlalp,psip,prp,uztp,nfpp,vprlfp,vthprlfp, &
       nalpp,vprlalpp
  real(IDP), dimension(:,:,:), allocatable :: psi_nl,pr_nl,uzt_nl,nf_nl,vprlf_nl,vthprlf_nl,nalp_nl,vprlalp_nl
  real(IDP), dimension(:), allocatable :: widthi,gammai
  real(IDP), dimension(:,:), allocatable :: cmamm,cmamp,cmapm,cmapp
  real(IDP), dimension(:), allocatable :: amat,bmat,cmat
  real(IDP), dimension(:,:,:), allocatable :: amatw,bmatw,cmatw,amatwalp,bmatwalp,cmatwalp
  real(IDP), dimension(:), allocatable :: xt,yt
  real(IDP), dimension(:,:), allocatable :: xw,eilnd
  integer, dimension(:), allocatable :: ipc
  integer, dimension(:,:), allocatable :: ipcw,ipcwalp
  real(IDP), dimension(:,:), allocatable :: epsi,ephi,epr,eprnc,ekenc,eke,emenc,eme,ealp,ealpnc
  real(IDP), dimension(:), allocatable :: src,src_nf,src_nalpha

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
  !  eilnd electron-ion Landau damping terms
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
  !  src source and sinks thermal plasma
  !  src_nf source and sinks first EP population
  !  src_nalpha source and sinks second EP population     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dynamo
