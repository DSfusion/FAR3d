	subroutine ae_profiles

!   	Set up profiles using an external source

#ifdef IMAS
   use transp_eq
   use far3d_wout
   use dcon_ids_mod
#endif

        use param
        use cotrol
        use domain
        use equil
        use dynamo
        use scratch
        use ffunctions
        implicit none
		
		interface		
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
		end interface		
	
		integer :: 	j,ic,i,l
		real(IDP) ::  xpi,B0_e,R0_e,a_e,uion_e,kappa_e,delt_e,beta0_e,rmax_e,xx
		real(IDP) ::  etactte,qe,me,epsil,kblotz
        integer :: ns0, nunit, je	
		real(IDP), dimension(:), allocatable ::  rho_e,qprof,den_beam_e,den_ion_e,den_elec_e,den_imp_e,temp_beam_e,temp_ion_e, &
                                      temp_elec_e,pres_beam_e,pres_thermal_e,pres_equil_e,zeff_e,tor_rot_freq_e,tor_rot_vel_e, &
									  pol_rot_vel_e,den_alpha_e,temp_alpha_e
		real(IDP) ::  B_inboard,B_outboard,Rin,Rout,btor,bpol,bnrm,tecnt,ticnt,dnecnt,dnicnt
		real(IDP) ::  bt0,rmajr,rminr,mu0,mu0prof,va0,omgcya,betalf,bet00,vf0,vtor0		
		real(IDP) :: dnnbinn_max,dnenn_max,dninn_max,tinn_max,tenn_max,tbnnn_max, &
		             pthermalnn_max,ptotnn_max,dnalphann_max,talphann_max, &
			     ddm1, ddm2		
		real(IDP), dimension(:), allocatable :: bspl,cspl,dspl
		character*1 :: tb,cdum,cdum2
#ifdef IMAS
  TYPE(bdstruc) :: bd
  integer        :: ierr, rhovar
  write(*,*) 'TRANSP2FAR3D Profiles'
  call ids2dcon(6, ext_prof_name, 'd3d', ierr)
  if (ierr /= 0) then
     write(*,*) 'ids2dcon returned ierr = ',ierr
     stop
  endif
  call ids_beamGet(bd, 6, ierr) ! set equilibrium data
  if (ierr /= 0) then
     write(*,*) 'ids_beamGet returned ierr = ',ierr
     stop
  endif

  alpha_on = bd%lalpha

  ext_prof_len = bd%nr
  B0_e = bd%bt
  R0_e = bd%rgcen
  a_e = bd%minrad
  kappa_e = bd%elong
  delt_e = bd%triavg
  uion_e = bd%ionmass
  beta0_e = bd%beta0
  rmax_e = bd%rm
#endif
		tb = char(9)

		allocate(rho_e(ext_prof_len))		
		allocate(qprof(ext_prof_len))	
		allocate(den_beam_e(ext_prof_len))	
		allocate(den_ion_e(ext_prof_len))	
		allocate(den_elec_e(ext_prof_len))		
		allocate(den_imp_e(ext_prof_len))	
		allocate(temp_beam_e(ext_prof_len))	
		allocate(temp_ion_e(ext_prof_len))			
		allocate(temp_elec_e(ext_prof_len))	
		allocate(pres_beam_e(ext_prof_len))	
		allocate(pres_thermal_e(ext_prof_len))	
		allocate(pres_equil_e(ext_prof_len))	
		allocate(zeff_e(ext_prof_len))	
		allocate(pol_rot_vel_e(ext_prof_len))	
		allocate(tor_rot_vel_e(ext_prof_len))	
		allocate(tor_rot_freq_e(ext_prof_len))	
	           if(alpha_on .eq. 1) then
		    allocate(den_alpha_e(ext_prof_len))	
		    allocate(temp_alpha_e(ext_prof_len))	
                                end if
					
	     if(alpha_on .eq. 0) then
		 
!  		Option to read parameters and profiles directly from a data file
!  		assumed to be ext_prof_name. The profiles are interpolated. We use IS units.
!		Only one energetic particle species	

	      xpi = 4.*atan(1.)
	
#ifdef IMAS
  rho_e = bd%rho
  qprof = bd%q
  den_beam_e = bd%bdens
  den_ion_e = bd%idens
  den_elec_e = bd%edens
  den_imp_e = bd%zdens
  temp_beam_e = bd%tbeam
  temp_ion_e = bd%ti
  temp_elec_e = bd%te
  pres_beam_e = bd%pbeam
  pres_thermal_e = bd%pt
  pres_equil_e = bd%pe
  tor_rot_vel_e = bd%trot
  pol_rot_vel_e = bd%prot
  ns0 = ext_prof_len  ! number of rows
  write(*,*) 'END TRANSP2FAR3D'
#else
          nunit = 17		    
          open(unit=nunit,file=ext_prof_name,status="old",form="formatted")

          ns0 = ext_prof_len  ! number of rows 
  	        
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,*) B0_e
          read(nunit,'(a1)') cdum
          read(nunit,*) R0_e
          read(nunit,'(a1)') cdum
          read(nunit,*) a_e
          read(nunit,'(a1)') cdum
          read(nunit,*) kappa_e
          read(nunit,'(a1)') cdum
          read(nunit,*) delt_e
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,*) uion_e
          read(nunit,'(a25,f10.7,a7,f7.5)') cdum,beta0_e,cdum2,rmax_e
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum

		  if(DIIID_u .eq. 0) then		  
            do i=1,ns0
	          read(nunit,*) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),pres_beam_e(i), &
                            pres_thermal_e(i),pres_equil_e(i),tor_rot_vel_e(i),pol_rot_vel_e(i)       
  	        end do    
	        close(unit=nunit)
		  end if
		  
		  if(DIIID_u .eq. 1) then
            do i=1,ns0
	          read(nunit,*) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),pres_beam_e(i), pres_thermal_e(i), &
                            pres_equil_e(i),zeff_e(i),tor_rot_freq_e(i),tor_rot_vel_e(i)
			    pol_rot_vel_e(i) = 0.
  	        end do    
	        close(unit=nunit)
		  end if
		  	  
		  if(DIIID_u .eq. 2) then		  
            do i=1,ns0
	          read(nunit,*) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),pres_beam_e(i), pres_thermal_e(i), &
                            pres_equil_e(i),zeff_e(i),tor_rot_freq_e(i),tor_rot_vel_e(i), &
			    ddm1,ddm2
			    pol_rot_vel_e(i) = 0.
  	        end do    
	        close(unit=nunit)
		  end if	  
#endif
          uion = uion_e                                                                   !! ion species
          bt0 = B0_e                                                                      !! magnetic field axis
          rmajr = R0_e                                                                    !! major radius
          rminr = a_e                                                                     !! minor radius
          mu0prof = 4.*3.1415926e-7                                                       !! magnetic permeability vacuum
		  mu0=mu0prof
          va0 = B0_e/sqrt(uion_e*mu0prof*den_ion_e(1)*1.e+20*1.672e-27)                   !! Alfven velocity axis  
	  if(DIIID_u .eq. 1 .or. DIIID_u .eq. 2) &
	    va0 = B0_e/sqrt(uion_e*mu0prof*den_ion_e(1)*1.e+19*1.672e-27)		  
          omgcya = 1.602e-19*B0_e*R0_e/(1.672e-27*uion_e*va0)                             !! cyclotron FR
          betalf = 2.*mu0prof*(pres_beam_e(1)*1.e+3)/(B0_e**2)                            !! beta EP axis
          bet00 = betath_factor*2.*mu0prof*(pres_thermal_e(1)*1.e+3)/(B0_e**2)            !! beta thermal plasma axis
          vf0 = sqrt(2.)*(sqrt(1000.*temp_beam_e(1)*1.6e-19/(uion_e*1.672e-27)))/va0      !! normalized EP thermal velocity axis
          vthi = sqrt(1000.*temp_ion_e(1)*1.6e-19/(uion_e*1.672e-27))/va0                 !! normalized ion thermal velocity axis
          vthe = sqrt(1000.*temp_elec_e(1)*1.6e-19/(9.109e-31))/va0                              !! normalized electron thermal velocity axis
          vtor0 = tor_rot_vel_e(1)*1.e+5/va0                                              !! normalized toroidal velocity axis
	  qe = 1.602e-19								  !! electron charge (C)
	  me = 9.109e-31								  !! electron mass (Kg)		
	  epsil = 8.85e-12								  !! vacuum permittivity (F/m)
	  kblotz = 1.602e-19								  !! Conversion factor eV to J
  
          write(0,'("External Profiles: omgcya = ",f8.4," betalf = ",f8.4, &
                    " bet00 = ",f8.4,/,"vthi = ",f8.4," vthe = ",f8.4, &
                    " vtor0 = ",f8.4)') omgcya, betalf, bet00, vthi, vthe, vtor0

          write(0,'("bt0 = ",e12.5," rmajr = ",f8.4, &
                    " rminr = ",f8.4,/,"va0 = ",e15.7, &
                    " pres_beam0 = ",e15.7)') B0_e, rmajr,rminr,va0, pres_beam_e(1)
          write(0,'("Conversion factor from code freqeuncy to kHz: ",e15.7)') va0/(2000.*xpi*R0_e)		  

          allocate(bspl(ns0),cspl(ns0),dspl(ns0))
		  allocate(qprofile(0:mj))
		  allocate(dnnbi(0:mj))
		  allocate(dnnbinn(0:mj))		  
		  allocate(dne(0:mj))
		  allocate(dnenn(0:mj))		  
		  allocate(dni(0:mj))
		  allocate(dninn(0:mj))			  		  
		  allocate(ti(0:mj))
		  allocate(tinn(0:mj))	
	      allocate(te(0:mj))
	      allocate(tenn(0:mj))	
	      allocate(tbn(0:mj))	
	      allocate(tbnnn(0:mj))		  
		  allocate(vzt_eqp(0:mj))
		  allocate(vth_eqp(0:mj))
          allocate(pthermal(0:mj))
          allocate(pthermalnn(0:mj))
          allocate(pep(0:mj))
          allocate(pepnn(0:mj))	
          allocate(ptot(0:mj))
          allocate(ptotnn(0:mj))	
          allocate(vthermalep(0:mj))	
          allocate(vAlfven(0:mj))
          allocate(vtherm_ionP(0:mj))
          allocate(vtherm_elecP(0:mj))		  

		  call spline(ns0,rho_e,qprof,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            qprofile(je) = seval(ns0,xx,rho_e,qprof,bspl,cspl,dspl)
          end do  
		  
		  call spline(ns0,rho_e,den_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            dnnbinn(je) = 1.e+20*seval(ns0,xx,rho_e,den_beam_e,bspl,cspl,dspl)
			dnnbinn(je) = max(dnnbinn(je),1e15)
          end do  
                                    dnnbinn_max=maxval(dnnbinn(0:mj))
		  dnnbi(:)=dnnbinn(:)/dnnbinn_max
		  
          call spline(ns0,rho_e,den_elec_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
			dnenn(je) = 1.e+20*seval(ns0,xx,rho_e,den_elec_e,bspl,cspl,dspl)
          end do   	
		  dnenn_max=maxval(dnenn(0:mj))	
		  dne(:)=dnenn(:)/dnenn_max		  
		  
          call spline(ns0,rho_e,den_ion_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            dninn(je) = 1.e+20*seval(ns0,xx,rho_e,den_ion_e,bspl,cspl,dspl)
			vAlfven(je) = (B0_e/sqrt(uion_e*mu0prof*dninn(je)*1.672e-27_IDP))/va0
          end do  	
		  dninn_max=maxval(dninn(0:mj))	
		  dni(:)=dninn(:)/dninn_max		  
		  
          call spline(ns0,rho_e,temp_ion_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tinn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_ion_e,bspl,cspl,dspl)	
		    vtherm_ionP(je) = sqrt(tinn(je)*1.6e-19_IDP/(uion_e*1.672e-27_IDP))/va0			
          end do	
		  tinn_max=maxval(tinn(0:mj))	
		  ti(:)=tinn(:)/tinn_max			  
		  
          call spline(ns0,rho_e,temp_elec_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tenn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_elec_e,bspl,cspl,dspl)
            vtherm_elecP(je) = sqrt(tenn(je)*1.6e-19/(uion_e*9.1094e-31))/va0	
          end do	
		  tenn_max=maxval(tenn(0:mj))	
		  te(:)=tenn(:)/tenn_max		  
		  
          call spline(ns0,rho_e,temp_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tbnnn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_beam_e,bspl,cspl,dspl)
          end do	
		  tbnnn_max=maxval(tbnnn(0:mj))	
		  tbn(:)=tbnnn(:)/tbnnn_max			  

          call spline(ns0,rho_e,tor_rot_vel_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            vzt_eqp(je) = seval(ns0,xx,rho_e,tor_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
          end do 				  

          call spline(ns0,rho_e,pol_rot_vel_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            vth_eqp(je) = seval(ns0,xx,rho_e,pol_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
          end do 

          call spline(ns0,rho_e,pres_thermal_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            pthermalnn(je) = seval(ns0,xx,rho_e,pres_thermal_e,bspl,cspl,dspl)
          end do 
		  pthermalnn_max=maxval(pthermalnn(0:mj))	
		  pthermal(:)=pthermalnn(:)/pthermalnn_max	

          call spline(ns0,rho_e,pres_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            pepnn(je) = seval(ns0,xx,rho_e,pres_beam_e,bspl,cspl,dspl)
                        ptotnn(je) = pepnn(je) + pthermalnn(je)
          end do 

		  ptotnn_max=maxval(ptotnn(0:mj))	
		  ptot(:)=ptotnn(:)/ptotnn_max

		  if(DIIID_u .eq. 1 .or. DIIID_u .eq. 2) then
		 
            do je=0,mj
			  if(EP_dens_on .eq. 0) then
			    dnnbinn(je) = dnnbinn(je)/10._IDP
			    dnnbi(je) = dnnbinn(je)/dnnbinn(1)
			  end if
			  dnenn(je) = dnenn(je)/10._IDP
              dne(je) = dnenn(je)/dnenn(1)	
              dninn(je) = dninn(je)/10._IDP
			  dni(je) = dninn(je)/dninn(1)			
		  	  vAlfven(je) = (B0_e/sqrt(uion_e*mu0prof*dninn(je)*1.672e-27_IDP))/va0		
              vzt_eqp(je) = 100._IDP*vzt_eqp(je)			  
              vth_eqp(je) = 100._IDP*vth_eqp(je)			  
            end do 			  
 
		  end if		  
 
	      dnnbi(0) = dnnbi(1)		  
          dne(0) = dne(1)
          dni(0) = dni(1)	  
          te(0) = te(1)		  
          ti(0) = ti(1)
          tbnnn(0) = tbnnn(1)	  
          vzt_eqp(0) = vzt_eqp (1)	
          pthermal(0) = pthermal(1)  

          coul_log = 24. - log(1.e-3*sqrt(dnenn_max)/tenn_max)                                      !! Coulomb logarithm (Te > 10 eV)
	  etactte = 0.02116*uion_e*uion_e*coul_log*qe*qe*sqrt(me)/(epsil*epsil*sqrt(kblotz*kblotz*kblotz))  !! ctte resistivity

          do je=0,mj		
		    denseq(je)=dne(je)
		    teeq(je)=te(je)	
		    if(EP_dens_on .eq. 0) then
              nfeq(je)=dnnbi(je) 
			end if  
			if(EP_vel_on .eq. 0) then
			  vfova(je)=sqrt(tbnnn(je)*1.602e-19/(uion_e*1.672e-27))/(va0*LcA3)
			end if  
            vthermalep(je)=sqrt(tbnnn(je)*1.602e-19/(uion_e*1.672e-27))	

		    if(q_prof_on .eq. 1) then	
			qq(je)=qprofile(je)
		    end if  
		        qq(je)=qq(je)+deltaq		
		    qqinv(je)=(1.0_IDP/qq(je))+deltaiota
		    etann(je)=etactte*dninn(je)/(dnenn(je)*sqrt(tenn(je)*tenn(je)*tenn(je)))

                    vzt_eq(je)=0._IDP; vth_eq(je)=0._IDP   !defaults if Eq_vel_on = 0 and Eq_velp_on
		    if(Eq_vel_on .eq. 1) then
                                        vzt_eq(je)=vzt_eqp(je)
		    end if  
		    if(Eq_velp_on .eq. 1) then
                                        vth_eq(je)=vth_eqp(je)
		    end if 	
		    if(Eq_Presseq_on .eq. 1) then
                                        preq(je)=pthermal(je)
                                                 		    if(Eq_Presstot_on .eq. 1) then
                                                                     preq(je)=ptot(je)
		                  end if 
		    end if  
		
          end do 	
 
          if(Edge_on .eq. 1) then      
            do je=edge_p,mj
	      etann(je)=etann(je-2) + (r(je)-r(je-2))*(etann(je-1)-etann(je-2))/(r(je-1)-r(je-2))
            end do
          end if

	   eta(:)=etann(:)/etann(1)	
  
          if(ieldamp_on .eq. 1) then			
             xnuelc0 = 2.89e-12*dninn_max*coul_log/(tenn_max*sqrt(tenn_max))         !! electron-ion collision FR axis (MKS) 			 
             xnuelc0 = rmajr*xnuelc0/va0                                       			        !! normalized electron-ion collision FR axis	 	
          end if		  
	  
		 end if	

		 close(17)

		 if(alpha_on .eq. 1) then		
	
!  		Option to read parameters and profiles directly from a data file
!  		assumed to be "ext_prof_name.txt". The profiles are interpolated. We use IS units.
!		Two energetic particle species		

	      xpi = 4.*atan(1.)
	
#ifdef IMAS
  rho_e = bd%rho
  qprof = bd%q
  den_beam_e = bd%bdens
  den_ion_e = bd%idens
  den_elec_e = bd%edens
  den_alpha_e = bd%adens
  den_imp_e = bd%zdens
  temp_beam_e = bd%tbeam
  temp_ion_e = bd%ti
  temp_elec_e = bd%te
  temp_alpha_e = bd%ta
  pres_beam_e = bd%pbeam
  pres_thermal_e = bd%pt
  pres_equil_e = bd%pe
  tor_rot_vel_e = bd%trot
  pol_rot_vel_e = bd%prot
  write(*,*) 'END TRANSP2FAR3D'
#else
          nunit = 17		  
          open(unit=nunit,file=ext_prof_name,status="old",form="formatted")

          ns0 = ext_prof_len ! number of rows 
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,*) B0_e
          read(nunit,'(a1)') cdum
          read(nunit,*) R0_e
          read(nunit,'(a1)') cdum
          read(nunit,*) a_e
          read(nunit,'(a1)') cdum
          read(nunit,*) kappa_e
          read(nunit,'(a1)') cdum
          read(nunit,*) delt_e
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
          read(nunit,*) uion_e
          read(nunit,'(a25,f10.7,a7,f7.5)') cdum,beta0_e,cdum2,rmax_e
          read(nunit,'(a1)') cdum
          read(nunit,'(a1)') cdum
	  
           do i=1,ns0
	         read(nunit,*) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_alpha_e(i),den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),temp_alpha_e(i),pres_beam_e(i), &
                            pres_thermal_e(i),pres_equil_e(i),tor_rot_vel_e(i),pol_rot_vel_e(i)         
  	       end do    
	       close(unit=nunit)	  	  
#endif

          uion = uion_e                                                                   !! ion species
          bt0 = B0_e                                                                      !! magnetic field axis
          rmajr = R0_e                                                                    !! major radius
          rminr = a_e                                                                     !! minor radius
          mu0prof = 4.*3.1415926e-7                                                       !! magnetic permeability vacuum
		  mu0=mu0prof
          va0 = B0_e/sqrt(uion_e*mu0prof*den_ion_e(1)*1.e+20*1.672e-27)                   !! Alfven velocity axis
          omgcya = 1.602e-19*B0_e*R0_e/(1.672e-27*uion_e*va0)                             !! cyclotron FR
          betalf = 2.*mu0prof*(pres_beam_e(1)*1.e+3)/(B0_e**2)                            !! beta EP axis
          bet00 = betath_factor*2.*mu0prof*(pres_thermal_e(1)*1.e+3)/(B0_e**2)            !! beta thermal plasma axis
          vf0 = sqrt(2.)*(sqrt(1000.*temp_beam_e(1)*1.6e-19/(uion_e*1.672e-27)))/va0      !! normalized EP thermal velocity axis
          vthi = sqrt(1000.*temp_ion_e(1)*1.6e-19/(uion_e*1.672e-27))/va0                 !! normalized ion thermal velocity axis
          vthe = sqrt(1000.*temp_elec_e(1)*1.6e-19/(9.109e-31))/va0                              !! normalized electron thermal velocity axis
          vtor0 = tor_rot_vel_e(1)*1.e+5_IDP/va0                                              !! normalized toroidal velocity axis
	  qe = 1.602e-19_IDP								  !! electron charge (C)
	  me = 9.109e-31_IDP								  !! electron mass (Kg)		
	  epsil = 8.85e-12_IDP								  !! vacuum permittivity (F/m)
	  kblotz = 1.602e-19_IDP								  !! Conversion factor eV to J
		  
          write(0,'("External Profiles: omgcya = ",f8.4," betalf = ",f8.4, &
                    " bet00 = ",f8.4,/,"vthi = ",f8.4," vthe = ",f8.4, &
                    " vtor0 = ",f8.4)') omgcya, betalf, bet00, vthi, vthe, vtor0

          write(0,'("bt0 = ",e12.5," rmajr = ",f8.4, &
                    " rminr = ",f8.4,/,"va0 = ",e15.7, &
                    " pres_beam0 = ",e15.7)') B0_e, rmajr,rminr,va0, pres_beam_e(1)
          write(0,'("Conversion factor from code freqeuncy to kHz: ",e15.7)') va0/(2000.*xpi*R0_e)		  
  
          allocate(bspl(ns0),cspl(ns0),dspl(ns0))
		  allocate(qprofile(0:mj))
		  allocate(dnnbi(0:mj))
		  allocate(dnnbinn(0:mj))		  
		  allocate(dne(0:mj))
		  allocate(dnenn(0:mj))		  
		  allocate(dni(0:mj))
		  allocate(dninn(0:mj))		
		  allocate(dnalpha(0:mj))
		  allocate(dnalphann(0:mj))		  		  
		  allocate(ti(0:mj))
		  allocate(tinn(0:mj))	
	      allocate(te(0:mj))
	      allocate(tenn(0:mj))	
	      allocate(tbn(0:mj))	
	      allocate(tbnnn(0:mj))
		  allocate(talpha(0:mj))
		  allocate(talphann(0:mj))			  
		  allocate(vzt_eqp(0:mj))
		  allocate(vth_eqp(0:mj))
          allocate(pthermal(0:mj))
          allocate(pthermalnn(0:mj))
          allocate(pep(0:mj))
          allocate(pepnn(0:mj))	
          allocate(ptot(0:mj))
          allocate(ptotnn(0:mj))		
          allocate(vthermalep(0:mj))	
          allocate(vAlfven(0:mj))
          allocate(vtherm_ionP(0:mj))	
          allocate(vtherm_elecP(0:mj))	  

		  call spline(ns0,rho_e,qprof,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            qprofile(je) = seval(ns0,xx,rho_e,qprof,bspl,cspl,dspl)
          end do  
		  
		  call spline(ns0,rho_e,den_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            dnnbinn(je) = 1.e+20*seval(ns0,xx,rho_e,den_beam_e,bspl,cspl,dspl)
			dnnbinn(je) = max(dnnbinn(je),1e15)			
          end do  
		  dnnbinn_max=maxval(dnnbinn(0:mj))
		  dnnbi(:)=dnnbinn(:)/dnnbinn_max
		  
          call spline(ns0,rho_e,den_elec_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
			dnenn(je) = 1.e+20*seval(ns0,xx,rho_e,den_elec_e,bspl,cspl,dspl)
          end do   	
		  dnenn_max=maxval(dnenn(0:mj))	
		  dne(:)=dnenn(:)/dnenn_max				  
		  
          call spline(ns0,rho_e,den_ion_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            dninn(je) = 1.e+20*seval(ns0,xx,rho_e,den_ion_e,bspl,cspl,dspl)
			vAlfven(je) = (B0_e/sqrt(uion_e*mu0prof*dninn(je)*1.672e-27_IDP))/va0
          end do  
		  dninn_max=maxval(dninn(0:mj))	
		  dni(:)=dninn(:)/dninn_max			  

		  call spline(ns0,rho_e,den_alpha_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            dnalphann(je) = 1.e+20*seval(ns0,xx,rho_e,den_alpha_e,bspl,cspl,dspl)
			dnalphann(je) = max(dnalphann(je),1e15)			
          end do  	
		  dnalphann_max=maxval(dnalphann(0:mj))	
		  dnalpha(:)=dnalphann(:)/dnalphann_max			  
		  
          call spline(ns0,rho_e,temp_ion_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tinn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_ion_e,bspl,cspl,dspl)	
		    vtherm_ionP(je) = sqrt(tinn(je)*1.6e-19_IDP/(uion_e*1.672e-27_IDP))/va0			
          end do	
		  tinn_max=maxval(tinn(0:mj))	
		  ti(:)=tinn(:)/tinn_max				  
		  
          call spline(ns0,rho_e,temp_elec_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tenn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_elec_e,bspl,cspl,dspl)
            vtherm_elecP(je) = sqrt(tenn(je)*1.6e-19/(uion_e*9.1094e-31))/va0
          end do	
		  tenn_max=maxval(tenn(0:mj))	
		  te(:)=tenn(:)/tenn_max				  
		  
          call spline(ns0,rho_e,temp_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            tbnnn(je) = 1.e+3*seval(ns0,xx,rho_e,temp_beam_e,bspl,cspl,dspl)
          end do	
		  tbnnn_max=maxval(tbnnn(0:mj))	
		  tbn(:)=tbnnn(:)/tbnnn_max			  

          call spline(ns0,rho_e,temp_alpha_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            talphann(je) = 1.e+3*seval(ns0,xx,rho_e,temp_alpha_e,bspl,cspl,dspl)
          end do	
		  talphann_max=maxval(talphann(0:mj))	
		  talpha(:)=talphann(:)/talphann_max			  

          call spline(ns0,rho_e,tor_rot_vel_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            vzt_eqp(je) = seval(ns0,xx,rho_e,tor_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
          end do 	

          call spline(ns0,rho_e,pol_rot_vel_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            vth_eqp(je) = seval(ns0,xx,rho_e,pol_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
          end do 			  

          call spline(ns0,rho_e,pres_thermal_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            pthermalnn(je) = seval(ns0,xx,rho_e,pres_thermal_e,bspl,cspl,dspl)
          end do 	 
		  pthermalnn_max=maxval(pthermalnn(0:mj))	
		  pthermal(:)=pthermalnn(:)/pthermalnn_max	

          call spline(ns0,rho_e,pres_beam_e,bspl,cspl,dspl)
          do je=0,mj
            xx = r(je)
            pepnn(je) = seval(ns0,xx,rho_e,pres_beam_e,bspl,cspl,dspl)
                        ptotnn(je) = pepnn(je) + pthermalnn(je)
          end do 

		  ptotnn_max=maxval(ptotnn(0:mj))	
		  ptot(:)=ptotnn(:)/ptotnn_max	
		  
	      dnnbi(0) = dnnbi(1)		  
          dne(0) = dne(1)
          dni(0) = dni(1)
	      dnalpha(0) = dnalpha(1)		  
          te(0) = te(1)		  
          ti(0) = ti(1)
          tbnnn(0) = tbnnn(1)
          talphann(0) = talphann(1)		  
          vzt_eqp(0) = vzt_eqp(1)	
          pthermal(0) = pthermal(1)  

          coul_log = 24. - log(1.e-3*sqrt(dnenn_max)/tenn_max)                                          !! Coulomb logarithm (Te > 10 eV)
	  etactte = 0.02116*uion_e*uion_e*coul_log*qe*qe*sqrt(me)/(epsil*epsil*sqrt(kblotz*kblotz*kblotz))  !! ctte resistivity
	
          do je=0,mj		
		    denseq(je)=dne(je)
		    teeq(je)=te(je)	
		    if(EP_dens_on .eq. 0) then			
              nfeq(je)=dnnbi(je) 
			end if  
			if(EP_vel_on .eq. 0) then
			  vfova(je)=sqrt(tbnnn(je)*1.602e-19/(uion_e*1.672e-27))/(va0*LcA3)
			end if
            vthermalep(je)=sqrt(tbnnn(je)*1.602e-19/(uion_e*1.672e-27))		
                          if(q_prof_on .eq. 1) then	
			qq(je)=qprofile(je)
		    end if  
		        qq(je)=qq(je)+deltaq		
		    qqinv(je)=(1.0_IDP/qq(je))+deltaiota
		    if(Alpha_dens_on .eq. 0) then			
              nalpeq(je)=dnalpha(je) 	
			end if  
			if(Alpha_vel_on .eq. 0) then
			  valphaova(je)=sqrt(talphann(je)*1.602e-19/(2*1.672e-27))/(va0*LcA3alp)
                                               end if
		    etann(je)=etactte*dninn(je)/(dnenn(je)*sqrt(tenn(je)*tenn(je)*tenn(je)))		  

                    vzt_eq(je)=0._IDP; vth_eq(je)=0._IDP   !defaults if Eq_vel_on = 0 and Eq_velp_on
		    if(Eq_vel_on .eq. 1) then
                                    vzt_eq(je)=vzt_eqp(je)
			end if  
		    if(Eq_velp_on .eq. 1) then
                                        vth_eq(je)=vth_eqp(je)
		    end if 	
		    if(Eq_Presseq_on .eq. 1) then
                                        preq(je)=pthermal(je)
                                                 		    if(Eq_Presstot_on .eq. 1) then
                                                                     preq(je)=ptot(je)
		                  end if 
		    end if  

          end do 	

	   eta(:)=etann(:)/etann(1)	  		  
	  
          if(ieldamp_on .eq. 1) then			
             xnuelc0 = 2.89e-12*dninn_max*coul_log/(tenn_max*sqrt(tenn_max))        !! electron-ion collision FR axis (MKS) 			 
             xnuelc0 = rmajr*xnuelc0/va0                                                   	   !! normalized electron-ion collision FR axis
          end if	
		  
		  close(17)
		  
	     end if
		 
	end subroutine ae_profiles
