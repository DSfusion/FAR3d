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
		
		open (unit=7,file="Input_Model",status="old",form="formatted")
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

		write(6,'(" copy of Input_Model:")')
		
		do
			read(7,'(a)',end=20) char5
			write(6,*) char5
		end do
20    rewind(7)		
		
		read(7,'(a1)') cdum0										 
		read(7,'(a1)') cdum0
		read(7,*) nstres
		read(7,'(a1)') cdum0
		read(7,'(2a2)') (numrun(i),i=1,2)
		read(7,'(a1)') cdum0
		read(7,'(2a2,a1)') (numruno(i),i=1,3)
		read(7,'(a1)') cdum0
		read(7,'(a5)') numvac	
        read(7,'(a1)') cdum0
        read(7,*) nonlin	
        read(7,'(a1)') cdum0
        read(7,*) ngeneq
        read(7,'(a1)') cdum0
        read(7,'(a40)') eq_name			
        read(7,'(a1)') cdum0					
        read(7,*) maxstp	
        read(7,'(a1)') cdum0
        read(7,*) dt0	
        read(7,'(a1)') cdum0		
        read(7,*) ldim		
        read(7,'(a1)') cdum0
        read(7,*) leqdim	
        read(7,'(a1)') cdum0	
        read(7,*) jdim		
		
		if (nstres /= 0) then

			confil="fs"//numruno(1)//numruno(2)//numruno(3)
			open(unit=8,file=confil,status='unknown',form='unformatted')

			read(8) ihist
			rewind(8)

		end if	
		
        close (7)		
		
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
		allocate (sceq1(0:jdim,0:leqdim))		! dummy variable for equilibrium modes
		allocate (sceq2(0:jdim,0:leqdim))	
		allocate (sceq3(0:jdim,0:leqdim))	
		allocate (sceq4(0:jdim,0:leqdim))	
		allocate (sceq5(0:jdim,0:leqdim))	
		allocate (sceq6(0:jdim,0:leqdim))	
		allocate (sceq7(0:jdim,0:leqdim))		

		call dfault	

		open (unit=7,file="Input_Model",status="old",form="formatted")		

		read(7,'(a1)') cdum0										 
		read(7,'(a1)') cdum0
		read(7,*) nstres
		read(7,'(a1)') cdum0
		read(7,'(2a2)') (numrun(i),i=1,2)
		read(7,'(a1)') cdum0
		read(7,'(2a2,a1)') (numruno(i),i=1,3)
		read(7,'(a1)') cdum0
		read(7,'(a5)') numvac	
		read(7,'(a1)') cdum0
		read(7,*) nonlin	
		read(7,'(a1)') cdum0
		read(7,*) ngeneq	
		read(7,'(a1)') cdum0
		read(7,'(a40)') eq_name			
		read(7,'(a1)') cdum0	
		read(7,*) maxstp	
		read(7,'(a1)') cdum0
		read(7,*) dt0	
		read(7,'(a1)') cdum0		
		read(7,*) ldim		
		read(7,'(a1)') cdum0
		read(7,*) leqdim	
		read(7,'(a1)') cdum0	
		read(7,*) jdim	
		read(7,'(a1)') cdum0	
		read(7,*) ext_prof	
		read(7,'(a1)') cdum0	
		read(7,'(a40)') ext_prof_name	
		read(7,'(a1)') cdum0	
		read(7,*) ext_prof_len		
		read(7,'(a1)') cdum0	
		read(7,*) iflr_on	
		read(7,'(a1)') cdum0		
		read(7,*) epflr_on	
		read(7,'(a1)') cdum0	
		read(7,*) ieldamp_on	
		read(7,'(a1)') cdum0	
		read(7,*) twofl_on	
		read(7,'(a1)') cdum0
		read(7,*) alpha_on	
		read(7,'(a1)') cdum0	
		read(7,*) Trapped_on	
		read(7,'(a1)') cdum0	
		read(7,*) matrix_out
		read(7,'(a1)') cdum0
		read(7,*) m0dy
		read(7,'(a1)') cdum0	
					
		read(7,'(a1)') cdum0										 
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0		
		read(7,*) (mm(i),i=1,ldim)
		read(7,'(a1)') cdum0		
		read(7,*) (nn(i),i=1,ldim)				
		read(7,'(a1)') cdum0		
		read(7,*) (mmeq(i),i=1,leqdim)				
		read(7,'(a1)') cdum0		
		read(7,*) (nneq(i),i=1,leqdim)	  
		read(7,'(a1)') cdum0	
		read(7,'(a1)') cdum0		
		read(7,*) ipert	
		read(7,'(a1)') cdum0		
		read(7,*) widthix
		read(7,'(a1)') cdum0		
		read(7,*) Auto_grid_on		
		read(7,'(a1)') cdum0		
		read(7,*) ni
		read(7,'(a1)') cdum0		
		read(7,*) nis
		read(7,'(a1)') cdum0		
		read(7,*) ne
		read(7,'(a1)') cdum0		
		read(7,*) delta
		read(7,'(a1)') cdum0		
		read(7,*) rc	
		read(7,'(a1)') cdum0			
		read(7,'(a1)') cdum0		
		read(7,*) gamma			
		read(7,'(a1)') cdum0		
		read(7,*) s	
		read(7,'(a1)') cdum0		
		read(7,*) betath_factor
		read(7,'(a1)') cdum0		
		read(7,*) ietaeq		
		read(7,'(a1)') cdum0		
		read(7,*) bet0_f	
		read(7,'(a1)') cdum0		
		read(7,*) bet0_alp		
		read(7,'(a1)') cdum0		
		read(7,*) omcy			
		read(7,'(a1)') cdum0		
		read(7,*) omcyb			
		read(7,'(a1)') cdum0
		read(7,*) rbound			
		read(7,'(a1)') cdum0		
		read(7,*) omcyalp		
		read(7,'(a1)') cdum0		
		read(7,*) itime	
		read(7,'(a1)') cdum0		
		read(7,*) dpres		
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0		
		read(7,*) stdifp			
		read(7,'(a1)') cdum0		
		read(7,*) stdifu	
		read(7,'(a1)') cdum0		
		read(7,*) stdifv	
		read(7,'(a1)') cdum0	
		read(7,*) stdifnf		
		read(7,'(a1)') cdum0		
		read(7,*) stdifvf	
		read(7,'(a1)') cdum0		
		read(7,*) stdifnalp	
		read(7,'(a1)') cdum0		
		read(7,*) stdifvalp
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0		
		read(7,*) LcA0			
		read(7,'(a1)') cdum0		
		read(7,*) LcA1	
		read(7,'(a1)') cdum0		
		read(7,*) LcA2		
		read(7,'(a1)') cdum0		
		read(7,*) LcA3	
		read(7,'(a1)') cdum0		
		read(7,*) LcA0alp	
		read(7,'(a1)') cdum0		
		read(7,*) LcA1alp
		read(7,'(a1)') cdum0		
		read(7,*) LcA2alp	
		read(7,'(a1)') cdum0		
		read(7,*) LcA3alp		
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0		
		read(7,*) omegar			
		read(7,'(a1)') cdum0		
		read(7,*) iflr	
		read(7,'(a1)') cdum0		
		read(7,*) r_epflr		
		read(7,'(a1)') cdum0		
		read(7,*) r_epflralp	
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0
		read(7,*) lplots			
		read(7,'(a1)') cdum0		
		read(7,*) nprint	
		read(7,'(a1)') cdum0		
		read(7,*) ndump	
		read(7,'(a1)') cdum0
		read(5,*) vtk_on
		read(5,'(a1)') cdum0
		read(7,'(a1)') cdum0
		read(7,*) DIIID_u	
		
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0
		read(7,'(a1)') cdum0		
		read(7,*) EP_dens_on			
		read(7,'(a1)') cdum0		
		read(7,*) Adens			
		read(7,'(a1)') cdum0		
		read(7,*) Bdens	
		read(7,'(a1)') cdum0		
		read(7,*) Alpha_dens_on		
		read(7,'(a1)') cdum0		
		read(7,*) Adensalp		
		read(7,'(a1)') cdum0		
		read(7,*) Bdensalp		
		read(7,'(a1)') cdum0		
		read(7,*) EP_vel_on	
		read(7,'(a1)') cdum0	
		read(7,*) Alpha_vel_on	
		read(7,'(a1)') cdum0		
		read(7,*) q_prof_on
		read(7,'(a1)') cdum0		
		read(7,*) Eq_vel_on
		read(7,'(a1)') cdum0	
		read(7,*) Eq_velp_on				
		read(7,'(a1)') cdum0	
		read(7,*) Eq_Presseq_on			
		read(7,'(a1)') cdum0
		read(7,*) Eq_Presstot_on			
		read(7,'(a1)') cdum0
		read(7,*) deltaq			
		read(7,'(a1)') cdum0	
		read(7,*) deltaiota	
		read(7,'(a1)') cdum0	
		read(7,*) etascl
		read(7,'(a1)') cdum0	
		read(7,*) eta0
		read(7,'(a1)') cdum0	
		read(7,*) reta
		read(7,'(a1)') cdum0	
		read(7,*) etalmb		
		read(7,'(a1)') cdum0			
		read(7,*) (cnep(i),i=0,10)		
		read(7,'(a1)') cdum0		
		read(7,*) (ctep(i),i=0,10)	
		read(7,'(a1)') cdum0	
		read(7,*) (cnfp(i),i=0,10)			
		read(7,'(a1)') cdum0		
		read(7,*) (cvep(i),i=0,10)
		read(7,'(a1)') cdum0		
		read(7,*) (cvfp(i),i=0,10)	
		read(7,'(a1)') cdum0		
		read(7,*) (cnfpalp(i),i=0,10)	
		read(7,'(a1)') cdum0		
		read(7,*) (cvfpalp(i),i=0,10)	
		read(7,'(a1)') cdum0		
		read(7,*) (eqvt(i),i=0,10)
		read(7,'(a1)') cdum0		
		read(7,*) (eqvp(i),i=0,10)	

		widthi(:)=widthix

		if (Auto_grid_on == 1) then

		  ni =  (jdim/2) - 1
		  nis = (jdim/4) + 1
		  ne =  jdim/4
                  
		end if	
		
		close (7)
		
		close(6)
		
	end subroutine inputlist
