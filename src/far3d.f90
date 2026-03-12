!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Open-source MIT LICENSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    FAR3d gyrofluid code ver 2.0 (2 EP real nonlinear version)                !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Copyright (C) 2024  D. Spong, L. Garcia, J. Varela and Y. Ghai            !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  						  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is free software: you can redistribute it and/or modify      !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    it under the terms of the open-source MIT license.                        !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!										  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    This program is distributed in the hope that it will be useful,           !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of            !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   		          !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program far3d

  use mpi
  use param
  use processor
  use cotrol
  use domain
  use equil
  use dynamo
  use initialize
  use diagnostics
  use output_mod
  use transfer
  use tools
  use scratch
  use openacc
  
  implicit none

  namelist/nam_par/maxstp,ndump,nprint,ndiag,lplots,itime,dt0,nonlin,mj,lmax,leqmax,ni,nis,ne,delta,rc,fti,fte, &
       difnr_on,stdifp,stdifu,stdifv,eps,bet0,etascl,reta,eta0,etalmb,ietaeq, &
       s,gamma,ipert,pertscl,m0dy,nocpl,xle,omcy,bet0_f,stdifnf,stdifvf,stdifnalp,stdifvalp, &
       ext_prof,ext_prof_name,epflr_on,r_epflr,alpha_on,iflr_on,iflr,LcA0,LcA1,LcA2,LcA3, &
       Adens,Bdens,twofl_on,dpres,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp, &
       omegar,ieldamp_on,omcyb,rbound,trapped_on,betath_factor,spe1,spe2,EP_dens_on, &
       EP_vel_on,Alpha_dens_on,Alpha_vel_on,DIIID_u,Eq_vel_on,Eq_velp_on,q_prof_on,deltaq,deltaiota,Eq_Presseq_on, &
       Eq_Presstot_on,Edge_on,edge_p,Auto_grid_on,nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on, &
       src_sink_th_on,src_sink_EP1_on,src_sink_EP2_on,src_sink_DIIID_on,src_sink_ITER_on,rsrc,wsrc,asrc, &
       rsrc_EP1,wsrc_EP1,asrc_EP1,rsrc_EP2,wsrc_EP2,asrc_EP2,AWfctr,Nfctr,AWfctr_dif,Rfctr,Wfctr,B_par_on,old_rd
  namelist/nam_arr/mm,nn,mmeq,nneq,widthi,gammai,cnep,ctep,cvep,cnfp,cvfp,cnfpalp,cvfpalp,eqvt,eqvp, &
       srcsinkth,srcsinkEP1,srcsinkEP2

  integer :: i,l,j,iend,l2t,mnum,idum,mjo,leqmaxo,ierror
  integer, dimension(8) :: values_s,values_e
  real(IDP) :: d,dk,scmn,scmx,sclmn,sclmx,scnorm,dum
  character(len=132) :: char5
  character(len=10) :: timew,datew,zone_h
  character(len=16) :: confil
  character(len=1) :: t
  character(len=10) :: hifrq_out
  character(len=14) :: phi_hifrq_out
  character(len=2), dimension(3) :: numrunp
  character(len=5) :: numvac

!  Input read 
                                                   
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

! Initialize MPI
  call MPI_INIT(ierror)
! Find number of PEs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numPEs, ierror)
! Find the ID of this PE
  call MPI_COMM_RANK(MPI_COMM_WORLD, myPE  , ierror)
  numPEsm1 = numPEs-1
  numdev = acc_get_num_devices(ACC_DEVICE_NVIDIA)
  call acc_set_device_num(myPE,ACC_DEVICE_NVIDIA)
  call acc_init(ACC_DEVICE_NVIDIA)
  
  ! if (myPE == 0) write(0,'(" ====> Checking input list ... ")')

  open (unit=5,file="farin",status="old")
  
  read(5,'(i1,2x,2a2,2x,2a2,a1,2x,a40)') nstres,(numrun(i),i=1,2),(numruno(i),i=1,3),eq_name
  rewind(5)

  t=char(9)

  if (myPE == 0) then

     confil="farprt"//numrun(1)//numrun(2)
     open (unit=6,file=confil)

     write(6,'(" copy of in:")')

     do while (.true.)
        read(5,'(a)',end=20) char5
        write(6,*) char5
     end do
20   rewind(5)

! Date and time

     call date_and_time(datew,timew,zone_h,values_s)
     write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
                                                                                  datew(7:8),datew(5:6),datew(1:4)
  end if

  call dfault

  read(5,'(i1,2x,2a2,2x,2a2,a1,2x,a40)') nstres,(numrun(i),i=1,2),(numruno(i),i=1,3),eq_name

  if (nstres /= 0) then

     confil="fs"//numruno(1)//numruno(2)//numruno(3)
     open(unit=8,file=confil,status='old',convert='big_endian',form='unformatted')

     read(8) ihist
     rewind(8)

  end if

  ihist=ihist+1

  allocate (numhist(ihist))

  if (nstres /= 0) then

     read(8) idum,numrunp,numrunp,numrunp,idum,numvac,(numhist(i),i=1,ihist-1),(idum,i=1,10),dum,(idum,i=1,6), &
             (dum,i=1,8),ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,idum,ieldamp_on,(dum,i=1,9),r_epflralp
     read(8) mjo,lmaxo,leqmaxo
     rewind(8)

     mj=mjo
     lmax=lmaxo
     leqmax=leqmaxo

  end if

  read(5,nam_par)

  if (nstres /= 0) then

     if (mj /= mjo) then
        mj=mjo
        if (myPE == 0) write(6,'(/"Grid points cannot be changed"/)')
     end if
     if (leqmax /= leqmaxo) then
        leqmax=leqmaxo
        if (myPE == 0) write(6,'(/"Equilibrium modes cannot be changed"/)')
     end if

  end if

! Identify the number of equations in the model

  noeqn=7
  ivalp=0
  iq=0
  iw=0
  ix1=0
  ix2=0
  iwa=0
  ix1a=0
  ix2a=0
  if (alpha_on == 1) then
     noeqn=noeqn+2
     ivalp=noeqn
  end if
  nvar=noeqn
  if (ext_prof == 0 .and. iflr_on == 2) iflr_on = 1
  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) noeqn=noeqn+1
  iq=noeqn
  if (epflr_on > 0) then
     if (r_epflr > 0.0 .or. epflr_on == 2) then
        iw=noeqn+1
        ix1=noeqn+2
        noeqn=noeqn+3
        ix2=noeqn
     end if
     if (alpha_on == 1 .and. (r_epflralp > 0.0 .or. epflr_on == 2)) then
        iwa=noeqn+1
        ix1a=noeqn+2
        noeqn=noeqn+3
        ix2a=noeqn
     end if
  end if

  if (ext_prof == 0) ieldamp_on = 0

  allocate (mm(lmax),nn(lmax),mh(lmax),nh(lmax),mmeq(leqmax),nneq(leqmax),mheq(leqmax),nheq(leqmax))
  allocate (r(0:mj),rinv(0:mj),dc1m(mj),dc1p(mj),dc2m(mj),dc2p(mj),del2cm(mj),del2cp(mj),rs(lmax))
  allocate (wt1m(mj,2),wt10(mj,2),wt1p(mj,2),wt2m(mj,2),wt20(mj,2),wt2p(mj,2))
  allocate (signl(lmax),jsl(lmax),sgnleq(leqmax))
  allocate (qq(0:mj),qqinv(0:mj),qqinvp(0:mj),denseq(0:mj),denseqr(0:mj),preq(0:mj),feq(0:mj), &
            cureq(0:mj),teeq(0:mj),tieq(0:mj),nfeq(0:mj),dnfeqdr(0:mj),dpreqdr(0:mj),vfova(0:mj),vzt_eq(0:mj), &
            vth_eq(0:mj),vfova2(0:mj),vtherm_elc(0:mj),nalpeq(0:mj),valphaova(0:mj),valphaova2(0:mj),dnalpeqdr(0:mj))
  allocate (eta(0:mj),widthi(lmax),gammai(lmax))
  allocate (sceq1(0:mj,0:leqmax),sceq2(0:mj,0:leqmax),sceq3(0:mj,0:leqmax),sceq4(0:mj,0:leqmax),sceq5(0:mj,0:leqmax), &
            sceq6(0:mj,0:leqmax))
  allocate (sd1(0:mj),sd2(0:mj),sd3(0:mj),sd4(0:mj),sd5(0:mj),sd6(0:mj),sd7(0:mj))
  allocate (epsi(lmax,2),ephi(lmax,2),epr(lmax,2),eprnc(lmax,2),ekenc(lmax,2),eke(lmax,2),emenc(lmax,2),eme(lmax,2), &
            ealp(lmax,2),ealpnc(lmax,2))

  sd1=0.0_IDP;sd2=0.0_IDP;sd3=0.0_IDP;sd4=0.0_IDP;sd5=0.0_IDP;sd6=0.0_IDP;sd7=0.0_IDP
  sceq1=0.0_IDP;sceq2=0.0_IDP;sceq3=0.0_IDP;sceq4=0.0_IDP;sceq5=0.0_IDP;sceq6=0.0_IDP
  r=0.0_IDP;rinv=0.0_IDP;dc1m=0.0_IDP;dc1p=0.0_IDP;dc2m=0.0_IDP
  dc2p=0.0_IDP;del2cm=0.0_IDP;del2cp=0.0_IDP;rs=0.0_IDP
  wt1m=0.0_IDP;wt10=0.0_IDP;wt1p=0.0_IDP;wt2m=0.0_IDP;wt20=0.0_IDP;wt2p=0.0_IDP
  signl=0.0_IDP;jsl=0.0_IDP;sgnleq=0.0_IDP
  widthi=0.0
  gammai=0.0
  read(5,nam_arr)

!  nstres indicates if the run is a new run or a continuation
!  If it is a new run the inital subroutine is called
!  If it is a continuation the resume subroutine is called

  if (nstres == 0) call inital

  if (nstres == 0) write(0,'(" ====> Preparing new run ...")')

  if (nstres /= 0) call resume

  if (nstres /= 0) write(0,'(" ====> Preparing run continuation ...")')

!  Modification of the equilibrium thermal beta
  if (betath_factor /= 1.0_IDP) write(0,'(" WARNING THERMAL BETA FACTOR ACTIVE = ",1pe13.6)') betath_factor 

!  Modification of the safety factor / iota profiles
  if (deltaq /= 0) write(0,'(" WARNING SAFETY FACTOR PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaq
  if (deltaiota /= 0) write(0,'(" WARNING IOTA PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaiota

  if (itime == 1) then
     dk=0.
     do l=1,lmax
        do j=0,mj
           d=nn(l)-mm(l)*qqinv(j)
           dk=max(dk,abs(d))
        end do
     end do
     dt=dt0*2./dk
     dtd2=dt*.5
  else
     dt=dt0
     dtd2=dt*.5
  end if
  if (myPE == 0) then

     confil="fs"//numrun(1)//numrun(2)
     numhist(ihist)=confil

!  The model parameter are included in farprt output file

     call output

     if(nonlin .ne. 0 .and. myPE == 0) then
        hifrq_out="hifrq_"//numrun(1)//numrun(2)
        phi_hifrq_out="phi_hifrq_"//numrun(1)//numrun(2)
        open(unit=82,file=hifrq_out,status='unknown',form='formatted')
        open(unit=77,file=phi_hifrq_out,status='unknown',form='formatted')
        if(nstres .eq. 0) write(82, &
             '("time",a1,"bth_r.2",a1,"bth_r.4",a1,"bth_r.6",a1,"bth_r.8", &
             a1,"bth_r1",a1,"avg(phi)",a1,"avg(psi)",a1,"avg(denf)",a1,"avg(vprlf)")') t,t,t,t,t,t,t,t,t 
        if(nstres .eq. 0) write(77, &
             '("time",a1,"phi_r.1",a1,"phi_r.2",a1,"phi_r.3",a1,"phi_r.4",a1,"phi_r.5",  &
             a1,"phi00_r.1",a1,"phi00_r.2",a1,"phi00_r.3",a1,"phi00_r.4",   &
             a1,"phi00_r.5")') t,t,t,t,t,t,t,t,t,t
        close(unit=82)
        close(unit=77)
     endif

     !  The model parameter are saved in fs####z file

     if (nstres == 0) call wrdump(0,.false.)
     ! call numinc
  end if

!  From here the model is advanced in time

  ! write(0,'(" ====> Time stepping begins ...")')  

  ! do i = 1,10
  !   nstep_count(i)=maxstp*i/10
  ! end do
  ! if (nstres /= 0) nstep_count(:)=nstep_count(:)+nstep
  ! i=1

  if (maxstp > 0) then

     nstep1=nstep

!  Subroutine linstart creates the tridiagonal matrix where the right and left side of the
!  model equations are added  

     call linstart

     iend = 0
     do while (nstep < nstep1+maxstp)

        nstep=nstep+1

        ! if (nstep_count(i) == nstep) then
        !    write(0,'(" ====> Time step = ",i8)') nstep
        !    i=i+1
        ! end if

        time=time+dt
        if (nstep == nstep1+maxstp .and. mod(nstep,nprint) /= 0) iend=1

!  Subroutine listep calculates perturvation variable time advance
!  Subroutine energy calculates the radial and poloidal component of the magentic and velocity fields
!  as well as kinetic and magnetic energy of the system.

        if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(1)
        call solve
        if (mod(nstep,nprint) == 0 .or. iend /= 0) call energy(2)

        call trnsfr0(psi,1)
        call trnsfr0(phi,-1)
        call trnsfr0(nf,1)
        call trnsfr0(vprlf,-1)

        if (mod(nstep,ndump) == 0 .and. nstep /= nstep1+maxstp) then
           call trnsfr0(pr,1)
           call trnsfr0(vthprlf,-1)
           if (alpha_on == 1) then
              call trnsfr0(nalp,1)
              call trnsfr0(vprlalp,-1)
           end if
           !  Subroutine wrdump writes an output file to continue the run
           if (myPE == 0) then
              call wrdump(idump,.false.)
              idump = idump + 1
              ! call numinc
           end if
        end if

        !   Write high-frequency (every time step) data:
        if(nonlin .ne. 0) then
           if(nstep .ge. nstep1+1 .and. myPE .eq. 0) then
              open(unit=82,file=hifrq_out,position='append',status='unknown')
              open(unit=77,file=phi_hifrq_out,position='append',status='unknown')
           endif
           if(nstep .lt. nstep1+maxstp) call hifreq(.false.)
           if(nstep .eq. nstep1+maxstp) call hifreq(.true.)
           close(unit=82)
           close(unit=77)
        endif

     end do

!  stepping done.
!  Subroutine lincheck calculates the growth rate and the frequency of the instability

     if (nonlin == 0) call lincheck

     call trnsfr0(psi,1)
     call trnsfr0(phi,-1)
     call trnsfr0(pr,1)
     call trnsfr0(nf,1)
     call trnsfr0(vprlf,-1)
     call trnsfr0(vthprlf,-1)
     if (alpha_on == 1) then
        call trnsfr0(nalp,1)
        call trnsfr0(vprlalp,-1)
     end if

     numrun(3)="z"
     if (myPE == 0) call wrdump(idump,.true.)
     call endrun

  end if

  write(0,'(" ====> Simulation DONE !! ")')

  if (myPE == 0) then
     call date_and_time(datew,timew,zone_h,values_e)
     write(6,'(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') timew(1:2),timew(3:4),timew(5:6), &
                   datew(7:8),datew(5:6),datew(1:4)

     call elapsed_time(values_s,values_e)
  end if

  CALL MPI_FINALIZE(ierror)

end program far3d
