!    This utility code reads the FAR3d profile data stored in Data.txt and then writes it back out to a
!    file called Data_mod.txt. This is useful in case one wants to make a modification to the initial
!    profiles or parameters that are in Data.txt.
!
	module param
		implicit none
		integer, parameter :: IDP = kind(1.d0)
	end module param												
        program read_modify_write
	use param
	  real(IDP) ::  xpi,B0_e,R0_e,a_e,uion_e,kappa_e,delt_e,beta0_e,rmax_e
	  real(IDP), dimension(:), allocatable :: rho_e,qprof,den_beam_e,den_ion_e,den_elec_e, &
                            den_imp_e,temp_beam_e,temp_ion_e, &
                            temp_elec_e,pres_beam_e, pres_thermal_e, &
                            pres_equil_e,zeff_e,tor_rot_freq_e,tor_rot_vel_e, &
			    ddm1,ddm2, pol_rot_vel_e,den_alpha_e,temp_alpha_e
	  character*1 :: tb,cdum,cdum2
	  integer :: 	i,nunit,ext_prof_len,ns0,DIIID_u
          nunit = 17
	  ext_prof_len = 101
	  DIIID_u = 2		  
          open(unit=nunit,file="Data.txt",status="old",form="formatted")

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
	  allocate(rho_e(ns0),qprof(ns0),den_beam_e(ns0),den_ion_e(ns0),den_elec_e(ns0), &
                            den_imp_e(ns0),temp_beam_e(ns0),temp_ion_e(ns0), &
                            temp_elec_e(ns0),pres_beam_e(ns0), pres_thermal_e(ns0), &
                            pres_equil_e(ns0),zeff_e(ns0),tor_rot_freq_e(ns0),tor_rot_vel_e(ns0), &
			    ddm1(ns0),ddm2(ns0),pol_rot_vel_e(ns0),den_alpha_e(ns0),temp_alpha_e(ns0))

	   if(alpha_on .eq. 0) then
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
	         read(nunit,*)  rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),pres_beam_e(i), pres_thermal_e(i), &
                            pres_equil_e(i),zeff_e(i),tor_rot_freq_e(i),tor_rot_vel_e(i), &
			    ddm1(i),ddm2(i)
		 
  	       end do    
	       close(unit=nunit)
             end if
            end if
	   
	   if(alpha_on .eq. 1) then
           do i=1,ns0
	         read(nunit,*) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_alpha_e(i),den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),temp_alpha_e(i),pres_beam_e(i), &
                            pres_thermal_e(i),pres_equil_e(i),tor_rot_vel_e(i),pol_rot_vel_e(i)         
  	       end do    
	       close(unit=nunit)	  	  
	   end if

           open(unit=nunit,file="Data_mod.txt",status="unknown",form="formatted")
           write(nunit,'("PLASMA GEOMETRY")')
           write(nunit,'("Vacuum Toroidal magnetic field at R=1.69550m [Tesla]")')
           write(nunit,'(f8.5)') B0_e
           write(nunit,'("Geometric Center Major radius [m]")')
           write(nunit,'(f8.5)') R0_e
           write(nunit,'("Minor radius [m]")')
           write(nunit,'(f8.5)') a_e
           write(nunit,'("Avg. Elongation")')
           write(nunit,'(f8.5)') kappa_e
           write(nunit,'("Avg. Top/Bottom Triangularity")')
           write(nunit,'(f8.5)') delt_e
           write(nunit,'("Main Contaminant Species")')
           write(nunit,'("	12C")')
           write(nunit,'("Main Ion Species mass/proton mass")')
           write(nunit,'(f3.1)') uion_e
           write(nunit,'("TRYING TO GET TO BETA(0)=",f9.7," , Rmax=",f7.5,/)') beta0_e,rmax_e
           write(nunit,'("  Rho(norml_sqrt_toroid_flux), q, BeamIonDensity(10^13cm^-3), IonDensity(10^13cm^-3),&
	    ElecDensity(10^13cm^-3), ImpurityDensity(10^13cm^-3), BeamIonEffectiveTemp(keV), IonTemp(keV),&
	     ElectronTemp(keV), BeamPressure(kPa), ThermalPressure(kPa), EquilPressure(kPa), Zeff, TorRot(kHz),&
	      TorRot(10^5m/s), ClassicalBeamPressure(kPa), ClassicalBeamDensity(10^13cm^-3)")')
           do i=1,ns0
	         write(nunit,'(f9.5,16(1x,f9.5))') rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                            den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                            temp_elec_e(i),pres_beam_e(i), pres_thermal_e(i), &
                            pres_equil_e(i),zeff_e(i),tor_rot_freq_e(i),tor_rot_vel_e(i), &
			    ddm1(i),ddm2(i)         
  	   end do    
	   close(unit=nunit)


	   stop
	   end program read_modify_write	  	  
