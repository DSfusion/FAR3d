MODULE initialize

  use param
  use cotrol
  use domain
  use equil
  use dynamo

  PRIVATE :: etachi, mmlims, pert
  PUBLIC  :: inital, resume
  
CONTAINS

  subroutine inital

    use processor
    use var_para
    use equilibrium
    use scratch

    implicit none

    ! This subroutine is called if the simulation is new, not a continuation

    numruns=numrun

    rewind(5)
    if (rc == 0.) rc=0.5_IDP
    if (delta >= rc .or. delta >= (1.-rc)) delta=min(rc,1.-rc)

    ! Subroutine setmod set up the modes of the model
    call setmod

    ! n value in the complex exponential representation.
    call mmlims

    ! Subroutine grid set up the time independent geometric scale factors
    call grid

    if (myPE == 0) then
       allocate (psi(0:mj,0:lmax),phi(0:mj,0:lmax),uzt(0:mj,0:lmax),pr(0:mj,0:lmax),nf(0:mj,0:lmax), &
            vprlf(0:mj,0:lmax),vthprlf(0:mj,0:lmax))
    else
       allocate (psi(mj_start:mj_end,0:lmax),phi(mj_start:mj_end,0:lmax),uzt(mj_start:mj_end,0:lmax), &
            pr(mj_start:mj_end,0:lmax),nf(mj_start:mj_end,0:lmax),vprlf(mj_start:mj_end,0:lmax), &
            vthprlf(mj_start:mj_end,0:lmax))
    end if
    allocate (psip(mj_start:mj_end,0:lmax),uztp(mj_start:mj_end,0:lmax),prp(mj_start:mj_end,0:lmax), &
         nfpp(mj_start:mj_end,0:lmax),vprlfp(mj_start:mj_end,0:lmax),vthprlfp(mj_start:mj_end,0:lmax))
    allocate (psi_nl(mj_start:mj_end,0:lmax,3),uzt_nl(mj_start:mj_end,0:lmax,3),pr_nl(mj_start:mj_end,0:lmax,3), &
         nf_nl(mj_start:mj_end,0:lmax,3),vprlf_nl(mj_start:mj_end,0:lmax,3),vthprlf_nl(mj_start:mj_end,0:lmax,3))
    allocate (sc1(mj_start:mj_end,0:lmax),sc2(mj_start:mj_end,0:lmax),sc3(mj_start:mj_end,0:lmax),sc4(mj_start:mj_end,0:lmax), &
         sc5(mj_start:mj_end,0:lmax),sc6(mj_start:mj_end,0:lmax),sc7(mj_start:mj_end,0:lmax),sc8(mj_start:mj_end,0:lmax), &
         sc9(mj_start:mj_end,0:lmax),sc10(mj_start:mj_end,0:lmax))

    psi=0.0_IDP
    phi=0.0_IDP
    pr=0.0_IDP
    nf=0.0_IDP
    vprlf=0.0_IDP
    vthprlf=0.0_IDP
    psip=0.0_IDP
    uztp=0.0_IDP
    prp=0.0_IDP
    nfpp=0.0_IDP
    vprlfp=0.0_IDP
    vthprlfp=0.0_IDP
    psi_nl=0.0_IDP
    uzt_nl=0.0_IDP
    pr_nl=0.0_IDP
    nf_nl=0.0_IDP
    vprlf_nl=0.0_IDP
    vthprlf_nl=0.0_IDP

    if (alpha_on == 1) then
       if (myPE == 0) then
          allocate (nalp(0:mj,0:lmax),vprlalp(0:mj,0:lmax))
       else
          allocate (nalp(mj_start:mj_end,0:lmax),vprlalp(mj_start:mj_end,0:lmax))
       end if
       allocate (nalpp(mj_start:mj_end,0:lmax),vprlalpp(mj_start:mj_end,0:lmax))
       allocate (sc11(mj_start:mj_end,0:lmax),sc12(mj_start:mj_end,0:lmax))
       allocate (nalp_nl(mj_start:mj_end,0:lmax,3),vprlalp_nl(mj_start:mj_end,0:lmax,3))
       nalp=0.0_IDP
       vprlalp=0.0_IDP
       nalpp=0.0_IDP
       vprlalpp=0.0_IDP
       nalp_nl=0.0_IDP
       vprlalp_nl=0.0_IDP
    end if

    ! Subroutine seteq set up the equilibria
    call seteq

    ! Subroutine etachi set up the magnetic diffusivity of the model
    call etachi

    ! Subroutine pert set up the perturbation of the equilibria
    call pert

  end subroutine inital

  subroutine resume

    use mpi
    use processor
    use var_para
    use equilibrium
    use output_mod
    use scratch

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

    integer :: i,j,l,idum,ierr,iPE,tag
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(IDP) :: dum

    ! This subroutine is called if the simulation is a continuation run

    ! Subroutine rddump reads the data required to continue the run
    call rddump

    rewind(5)
    read(5,'(i1,2x,2a2,2x,2a2,a1)') nstres,(numrun(i),i=1,2),(numruno(i),i=1,3)
    read(5,nam_par)
    read(5,nam_arr)

    ! Subroutine setmod set up the modes of the model
    call setmod

    ! n value in the complex exponential representation.
    call mmlims

    ! Needed when iflr_on = 2 or ieldamp_on = 1
    if (iflr_on == 2 .or. ieldamp_on == 1) call ae_profiles

    ! Change the thermal beta
    bet0=bet0*betath_factor

    ! Displacement included to the safety factor / iota
    if (deltaq .ne. 0.0) then
       qq = 1.0_IDP/qqinv + deltaq
       qqinv = 1.0_IDP/qq
    else if (deltaiota .ne. 0.0) then
       qqinv = qqinv + deltaiota
       qq = 1.0_IDP/qqinv
    end if

    if (myPE == 0) then
       allocate (psi(0:mj,0:lmax),phi(0:mj,0:lmax),uzt(0:mj,0:lmax),pr(0:mj,0:lmax),nf(0:mj,0:lmax), &
            vprlf(0:mj,0:lmax),vthprlf(0:mj,0:lmax))
    else
       allocate (psi(mj_start:mj_end,0:lmax),phi(mj_start:mj_end,0:lmax),uzt(mj_start:mj_end,0:lmax), &
            pr(mj_start:mj_end,0:lmax),nf(mj_start:mj_end,0:lmax),vprlf(mj_start:mj_end,0:lmax), &
            vthprlf(mj_start:mj_end,0:lmax))
    end if
    allocate (psip(mj_start:mj_end,0:lmax),uztp(mj_start:mj_end,0:lmax),prp(mj_start:mj_end,0:lmax), &
         nfpp(mj_start:mj_end,0:lmax),vprlfp(mj_start:mj_end,0:lmax),vthprlfp(mj_start:mj_end,0:lmax))
    allocate (psi_nl(mj_start:mj_end,0:lmax,3),uzt_nl(mj_start:mj_end,0:lmax,3),pr_nl(mj_start:mj_end,0:lmax,3), &
         nf_nl(mj_start:mj_end,0:lmax,3),vprlf_nl(mj_start:mj_end,0:lmax,3),vthprlf_nl(mj_start:mj_end,0:lmax,3))
    allocate (sc1(mj_start:mj_end,0:lmax),sc2(mj_start:mj_end,0:lmax),sc3(mj_start:mj_end,0:lmax),sc4(mj_start:mj_end,0:lmax), &
         sc5(mj_start:mj_end,0:lmax),sc6(mj_start:mj_end,0:lmax),sc7(mj_start:mj_end,0:lmax),sc8(mj_start:mj_end,0:lmax), &
         sc9(mj_start:mj_end,0:lmax),sc10(mj_start:mj_end,0:lmax))

    psip=0.0_IDP
    uztp=0.0_IDP
    prp=0.0_IDP
    nfpp=0.0_IDP
    vprlfp=0.0_IDP
    vthprlfp=0.0_IDP
    psi_nl=0.0_IDP
    uzt_nl=0.0_IDP
    pr_nl=0.0_IDP
    nf_nl=0.0_IDP
    vprlf_nl=0.0_IDP
    vthprlf_nl=0.0_IDP

    if (alpha_on == 1) then
       if (myPE == 0) then
          allocate (nalp(0:mj,0:lmax),vprlalp(0:mj,0:lmax))
       else
          allocate (nalp(mj_start:mj_end,0:lmax),vprlalp(mj_start:mj_end,0:lmax))
       end if
       allocate (nalpp(mj_start:mj_end,0:lmax),vprlalpp(mj_start:mj_end,0:lmax))
       allocate (sc11(mj_start:mj_end,0:lmax),sc12(mj_start:mj_end,0:lmax))
       allocate (nalp_nl(mj_start:mj_end,0:lmax,3),vprlalp_nl(mj_start:mj_end,0:lmax,3))
       nalpp=0.0_IDP
       vprlalpp=0.0_IDP
       nalp_nl=0.0_IDP
       vprlalp_nl=0.0_IDP
    end if

    sc1=0.
    sc2=0.
    sc3=0.
    sc4=0.
    sc5=0.
    sc6=0.
    sc7=0.
    sc8=0.
    sc9=0.
    sc10=0.

    if (lmax < lmaxo) then
       if (myPE == 0) then
          read(8) ((psi(j,l),j=0,mj),l=1,lmax)
          read(8) ((phi(j,l),j=0,mj),l=1,lmax)
          read(8) ((pr(j,l),j=0,mj),l=1,lmax)
          read(8) ((nf(j,l),j=0,mj),l=1,lmax)
          read(8) ((vprlf(j,l),j=0,mj),l=1,lmax)
          read(8) ((vthprlf(j,l),j=0,mj),l=1,lmax)
          if (alpha_on == 1) then
             read(8) ((nalp(j,l),j=0,mj),l=1,lmax)
             read(8) ((vprlalp(j,l),j=0,mj),l=1,lmax)
          end if
       else
          read(8) ((dum,j=0,mj_start-1),(psi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          read(8) ((dum,j=0,mj_start-1),(phi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          read(8) ((dum,j=0,mj_start-1),(pr(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          read(8) ((dum,j=0,mj_start-1),(nf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          read(8) ((dum,j=0,mj_start-1),(vprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          read(8) ((dum,j=0,mj_start-1),(vthprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          if (alpha_on == 1) then
             read(8) ((dum,j=0,mj_start-1),(nalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
             read(8) ((dum,j=0,mj_start-1),(vprlalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmax)
          end if
       end if
    else
       if (myPE == 0) then
          read(8) ((psi(j,l),j=0,mj),l=1,lmaxo)
          read(8) ((phi(j,l),j=0,mj),l=1,lmaxo)
          read(8) ((pr(j,l),j=0,mj),l=1,lmaxo)
          read(8) ((nf(j,l),j=0,mj),l=1,lmaxo)
          read(8) ((vprlf(j,l),j=0,mj),l=1,lmaxo)
          read(8) ((vthprlf(j,l),j=0,mj),l=1,lmaxo)
          if (alpha_on == 1) then
             read(8) ((nalp(j,l),j=0,mj),l=1,lmaxo)
             read(8) ((vprlalp(j,l),j=0,mj),l=1,lmaxo)
          end if
       else
          read(8) ((dum,j=0,mj_start-1),(psi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          read(8) ((dum,j=0,mj_start-1),(phi(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          read(8) ((dum,j=0,mj_start-1),(pr(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          read(8) ((dum,j=0,mj_start-1),(nf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          read(8) ((dum,j=0,mj_start-1),(vprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          read(8) ((dum,j=0,mj_start-1),(vthprlf(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          if (alpha_on == 1) then
             read(8) ((dum,j=0,mj_start-1),(nalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
             read(8) ((dum,j=0,mj_start-1),(vprlalp(j,l),j=mj_start,mj_end),(dum,j=mj_end+1,mj),l=1,lmaxo)
          end if
       end if
    end if

    if (old_rd) then
       read(8) etascl,reta,eta0,etalmb,ietaeq,dum,idum,dum,stdifp,stdifu,stdifnf,stdifv,stdifvf, &
            stdifnalp,stdifvalp,s,gamma,xle,ipert,(dum,l=1,lmaxo),(dum,l=1,lmaxo),dum
    else
       read(8) etascl,reta,eta0,etalmb,ietaeq,stdifp,stdifu,stdifnf,stdifv,stdifvf, &
            stdifnalp,stdifvalp,s,gamma,xle,ipert,(dum,l=1,lmaxo),(dum,l=1,lmaxo),dum
       read(8) ext_prof_name,difnr_on,nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on, &
            src_sink_th_on,src_sink_EP1_on,src_sink_EP2_on,src_sink_DIIID_on,src_sink_ITER_on,rsrc,wsrc,asrc, &
            rsrc_EP1,wsrc_EP1,asrc_EP1,rsrc_EP2,wsrc_EP2,asrc_EP2,AWfctr,Nfctr,AWfctr_dif,Rfctr,Wfctr
       read(8) (srcsinkth(i),i=0,10),(srcsinkEP1(i),i=0,10),(srcsinkEP2(i),i=0,10)
    end if

    rewind(8)
    close(unit=8)

    ! Equilibrium profiles are modified to flux surface averages when nonlin=0

    if (nonlin == 0) then
       if (myPE == 0) then
          preq=preq+pr(:,l0)
          nfeq=nfeq+nf(:,l0)
          do iPE=1,numPEsm1
             tag=iPE
             call MPI_SEND(preq(0),mjp1,MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,ierr)
             tag=numPEs+iPE
             call MPI_SEND(nfeq(0),mjp1,MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,ierr)
          end do
       else
          tag=myPE
          call MPI_RECV(preq(0),mjp1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
          tag=numPEs+myPE
          call MPI_RECV(nfeq(0),mjp1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
       end if
    end if

    ! Subroutine pert set up the perturbation of the equilibria
    call pert

  end subroutine resume

  subroutine mmlims

    !   Calculate the minimum and maximum m values for each n value
    !   in the complex exponential representation.

    implicit none

    integer :: n,l

    allocate (mmstart(-nmax:nmax))
    allocate (mmend(-nmax:nmax))
    allocate (mmstreq(-nmaxeq:nmaxeq))
    allocate (mmendeq(-nmaxeq:nmaxeq))

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

    do n = 0,nmaxeq
       mmstreq(n) = 1000
       mmendeq(n) = -1000
    end do

    do l = 1,leqmax

       n = nneq(l)
       if (n >= 0) then
          if (mmeq(l) < mmstreq(n)) mmstreq(n) = mmeq(l)
          if (mmeq(l) > mmendeq(n)) mmendeq(n) = mmeq(l)
       else
          if (-mmeq(l) < mmstreq(-n)) mmstreq(-n) = -mmeq(l)
          if (-mmeq(l) > mmendeq(-n)) mmendeq(-n) = -mmeq(l)
       end if

    end do

    !   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
    !   so limits for negative n rows are just negatives of corresponding
    !   positive n limits.

    do n = 0,nmaxeq

       mmstreq(-n) = -mmendeq(n)
       mmendeq(-n) = -mmstreq(n)

    end do

  end subroutine mmlims

  subroutine etachi

    implicit none

    ! Set up the magnetic diffusivity radial profile	

    select case (ietaeq)
    case(1)
       if (ext_prof == 0) eta=1.0_IDP/teeq**1.5_IDP
    case(2)
       eta=etascl
    case(3)
       eta=eta0*(1.+(r/reta)**(2.*etalmb))**(1./etalmb) 
    end select

  end subroutine etachi

  subroutine pert

    use mpi
    use processor
    use var_para
    use tools

    implicit none

    integer :: l,l1,lp,iamm,iann,j,jj,jm,jp,ndivi,n,k,i,clock,ierr,iPE,tag,dl_max
    real(IDP) :: sq2,teps,xran,yran,p,sigwid,qp,psipp,beta,rbeta,awid,psitil,phitil,yprb,xx,psij,phij,r32,y,rr,dr,gam,rmrbar,rsl,dum
    integer, dimension(:), allocatable :: seed
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(IDP), dimension(:,:), allocatable :: sct

    !  ipert = 0  indicates perturbation performed the old way
    !  ipert = 1            perturbation performed the new way

    if (pertscl /= 1.0_IDP) then
       psi=psi*pertscl 
       phi=phi*pertscl 
       pr=pr*pertscl 
       nf=nf*pertscl 
       vprlf=vprlf*pertscl 
       vthprlf=vthprlf*pertscl 
       if (alpha_on == 1) then
          nalp=nalp*pertscl
          vprlalp=vprlalp*pertscl
       end if
    endif

    pi=4.0_IDP*atan(1.0_IDP)
    sq2=sqrt(2.0_IDP)
    teps = 0.025
    if (myPE == 0) then
       k=8
       call random_seed(size=k)
       allocate(seed(k))
       call system_clock(count=clock)
       seed=clock+37*(/(i-1,i=1,k)/)
       call random_seed(put=seed)
       deallocate(seed)
    end if

    do l1=1,lmaxn
       l=lln(l1)
       lp=l
       if (lo(l1) /= 0) lp=lln(lo(l1))
       if (myPE == 0 .and. (signl(l) > 0 .or. lo(l1) == 0)) then
          call random_number(xran)
          call random_number(yran)
          xran=1.0_IDP
          yran=0.125_IDP
          if (ipert == 0) then
             if (lo(l1) /= 0) gammai(lp)=gammai(lp)*xran*sin(2.*yran*pi)
             gammai(l)=gammai(l)*xran*cos(2.*yran*pi)
          else
             if (lo(l1) /= 0) widthi(lp)=widthi(l)*xran*sin(2.*yran*pi)
             widthi(l)=widthi(l)*xran*cos(2.*yran*pi)
          end if
       end if

       !    compute rs

       iamm=mm(l)*sign(1,nn(l))
       iann=abs(nn(l))
       jsl(l)=0
       p=0.
       if (iann /= 0) p=iamm*1.0_IDP/iann
       do j=1,mj
          if (((1./qqinv(j-1))-p)*((1./qqinv(j))-p) <= 0.) exit
       end do
       if (j > mj) then
          rs(l)=.5
          jsl(l)=mj/2
       else
          rs(l)=r(j-1)+(r(j)-r(j-1))*((1./qqinv(j-1))-p)/((1./qqinv(j-1))-(1./qqinv(j)))
          jsl(l)=j
       end if

       if (myPE > 0 .or. widthi(l) == 0) cycle
       sigwid=sign(1.0_IDP,widthi(l))


       select case (ipert)
       case (0)

          !       old perturbation

          do j=1,mj
             if (rs(l) < r(j)) exit
          end do
          if (j > mj) j = mj
          qp=((1./qqinv(j))-(1./qqinv(j-1)))/(r(j)-r(j-1))
          psipp=rs(l)*qp*qqinv(j)**2
          if(iamm == 1 .and. iann == 1) then
             !          m=n=1, so use exact linear perturbation.
             beta=-3.
             awid=abs(widthi(l))
             rbeta=beta*awid+rs(l)
             phitil=-abs(psipp)*s*sqrt(2.*pi)*widthi(l)**2/(32.*rs(l))
             psitil=-abs(psipp)*(widthi(l)**2/16.)*(exp(-beta*beta*.5)-.5*sqrt(2.*pi)*beta*(1.-erf(beta/sq2)))
             do j=1,mjm1
                if (rbeta < r(j)) exit
             end do
             yprb=rbeta*(-1.+1.*qqinv(j))
             do j=1,mjm1
                xx=(r(j)-rs(l))/awid
                phij=phitil*r(j)*(1.-erf(xx/sq2))
                if (r(j) < rbeta) then
                   psij=(r(j)*(-1.+qqinv(j))/yprb)*psitil
                else
                   psij=-abs(psipp)*(widthi(l)**2/16.)*(exp(-xx*xx*.5)-.5*sqrt(2.*pi)*xx*(1.-erf(xx/sq2)))
                end if
                phi(j,l)=phi(j,l)+phij*sigwid
                !                psi(j,l)=psi(j,l)+psij*sigwid
             end do
          else
             !          guess a perturbation.
             psitil=-(widthi(l)*.25)**2*psipp*sigwid/(rs(l)**iamm*(1.-rs(l)))
             do j=1,mjm1
                psij=psitil*r(j)**iamm*(1.-r(j))*2./(1.+exp(10.*(-1.+r(j)/rs(l))))
                phi(j,l)=phi(j,l)+psij*gammai(l)*(rs(l)-r(j))/(rs(l)+r(j))
                !                psi(j,l)=psi(j,l)+psij
             end do
          end if

       case (1)

          !          0/0 perturbation

          if (l == l0) then
             p=1.5
             do j=1,mj
                if(((1./qqinv(j-1))-p)*((1./qqinv(j))-p).le.0.) exit
             end do
             if (j > mj) then
                r32=.5
             else
                r32=r(j-1)+(r(j)-r(j-1))*((1./qqinv(j-1))-p)/((1./qqinv(j-1))-(1./qqinv(j)))
             end if
             phitil=0.
             phi(mj,l) = phi(mj,l) + phitil
             y=phitil
             ndivi = 1000
             do jj=1,mj
                jm=mj-jj
                jp=jm+1
                dr=(r(jp)-r(jm))/ndivi
                rr=r(jp)+.5*dr
                do n=1,ndivi
                   rr=rr-dr
                   y=y-dr*rr*tanh((rr-r32)/xle)
                end do
                phitil=widthi(l)*y/r32
                phi(jm,l) = phi(jm,l) + phitil
             end do
          else

             !          new perturbation

             rsl = rs(l)
             gam = teps
             if(iamm > 0) gam = teps*iamm**(-1./3.)
             do j = 1,mjm1
                rmrbar = ((r(j)-rsl)/gam)**2

                !                set perturbations

                phitil = widthi(l)*exp(-0.5*rmrbar)
                psitil = 0.0

                !                perturb phi and psi

                phi(j,l) = phi(j,l) + phitil
                !                psi(j,l) = psi(j,l) + psitil
             end do
          end if

       end select

    end do

    dl_max=maxval(mj_dl)
    allocate (sct(dl_max,lmax))
    if (myPE == 0) then
       do iPE=1,numPEsm1
          do l=1,lmax
             sct(1:mj_dl(iPE),l)=phi(mj_st(iPE):mj_st(iPE)+mj_dl(iPE)-1,l)
          end do
          tag=iPE
          call MPI_SEND(sct,dl_max*lmax,MPI_DOUBLE_PRECISION,iPE,tag,MPI_COMM_WORLD,ierr)
       end do
    else
       tag=myPE
       call MPI_RECV(sct,dl_max*lmax,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
       do l=1,lmax
          phi(mj_start:mj_end,l)=sct(1:mj_dl(myPE),l)
       end do
    end if
    deallocate (sct)

  end subroutine pert

END MODULE initialize
