subroutine linstart

  use mpi
  use param
  use processor
  use var_para
  use cotrol
  use domain
  use equil
  use dynamo
  use om_mod
  use block_mod
  use dbyd
  use gam_mod
  use matrix
  use transfer
  use tools
  use scratch

  implicit none

  integer :: i,j,l,m,n,ibnd,mcnt2,nbdim,nbxdim,nbidim,it1,it2,le,lp,ln,l1,l2,l3,l1t,l2t,l3t,mnum3,mnum6,loca,loci, &
             ier,it,imat,imt1,imt2,nend,lskp,it3,is1,is2,is3,m1,m2,m3,mm3,mp3,n1,n2,n3,np3,nm3,sgn12,l3p, &
             i1,ip,ierr,iPE,tag,tagst,ind,lsk,mnum2,lpn,lq,ann
  integer, dimension(MPI_STATUS_SIZE) :: status
  real(IDP) :: epsq,oneos,betfc,betfc_f,betfc_alp,coef,omcyd,omcydalp,beteom,betiom,rmrbar,acoef,bcoef,x,aux,xnuion0,abkprl, &
               xsii,xsie,xsi2,xsi3,xsi4,ztr,zti,zir,zii,zer,zei,vther,vthir,tauie,xnuelc,xnuion,fctr,rsq
  complex(IDP) :: zetai,zetae,zetai2,zetai3,zetai4,zi,y0i,y1i,y2i,ddii,zetae2,zetae3,zetae4,ze,y0e,y1e,y2e,ddee,rei,sei, &
                  stfe,stfi,reii,sei1,sei2,sei3,zetaiinv,zetaeinv,cmplx1
  real(IDP), dimension(:), allocatable :: xa
  real(IDP), dimension(:,:), allocatable :: scp2,scp3
  character*1 :: tb
  
  interface
     subroutine dlstar(ss,ff,itypf,wk1,wk2,wkeq1,wkeq2,c1,c2)
        use param
        use var_para
        implicit none
        integer :: itypf
        real(IDP) :: c1,c2
        real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2
        real(IDP), dimension(0:,0:) :: wkeq1,wkeq2
     end subroutine dlstar
     subroutine dlsq(ss,ff,itypf,wk1,wk2,c1,c2)
        use param
        use var_para
        implicit none
        integer :: itypf
        real(IDP) :: c1,c2
        real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2
     end subroutine dlsq
  end interface

  tb = char(9)
  epsq=eps*eps
  oneos=0.0_IDP
  if (s > 0.0_IDP) oneos=1.0_IDP/s
  betfc=bet0/(2.*epsq)
  betfc_f=LcA2*bet0_f/(2.*epsq)
  betfc_alp=LcA2alp*bet0_alp/(2.*epsq)
  omcyd = omcy
  omcydalp = omcyalp
  beteom=dpres*bet0/(2*epsq*omcyd)
  betiom=(1-dpres)*bet0/(2*epsq*omcyd)

! edge-factor for vfova, nfeq, denseq and preq
  ! if (nstep1 == 0) then
  !    vfova = vfova*(1.0-exp(-(1.0-r)**2/1.e-3))
  !    acoef=nfeq(925)-37.5*(nfeq(924)-nfeq(926))/6.0
  !    bcoef=500.0*(nfeq(924)-nfeq(926))/(6*0.075**5)
  !    do j=926,mj
  !       nfeq(j)=acoef+bcoef*(1.0-r(j))**6
  !    end do
  !    acoef=denseq(925)-37.5*(denseq(924)-denseq(926))/6.0
  !    bcoef=500.0*(denseq(924)-denseq(926))/(6*0.075**5)
  !    do j=926,mj
  !       denseq(j)=acoef+bcoef*(1.0-r(j))**6
  !    end do
  !    acoef=preq(925)-37.5*(preq(924)-preq(926))/6.0
  !    bcoef=500.0*(preq(924)-preq(926))/(6*0.075**5)
  !    do j=926,mj
  !       preq(j)=acoef+bcoef*(1.0-r(j))**6
  !    end do
  ! end if

  call dbydr0(qqinvp,qqinv,0.0_IDP,1.0_IDP,0)
  call dbydr0(denseqr,denseq,0.0_IDP,1.0_IDP,0)
  call dbydr0(dnfeqdr,nfeq,0.0_IDP,1.0_IDP,0)
  call dbydr0(dpreqdr,preq,0.0_IDP,1.0_IDP,0)

!  define viscosities

  if (difnr_on == 1) then

     allocate (stdifpn(nnum))
     allocate (stdifun(nnum))
     allocate (stdifnfn(nnum))
     allocate (stdifvfn(nnum))
     allocate (stdifvn(nnum))

     do i=1,nnum
        l=lnumn(i-1)+1
        n=nn(lln(l))
        ann=abs(n)
        fctr=1.0+0.5*(AWfctr-1.0)*(1.+tanh(ann-Nfctr))    ! 1 except at n > Nfctr where it jumps to AWfctr
        stdifpn(i)=stdifp*fctr
        stdifun(i)=stdifu*fctr
        stdifnfn(i)=stdifnf*fctr
        stdifvfn(i)=stdifvf*fctr
        stdifvn(i)=stdifv*fctr
        if (lln(l) == l0) then
            stdifpn(i)=stdifp
            stdifun(i)=stdifu
            stdifnfn(i)=stdifnf
            stdifvfn(i)=stdifvf
            stdifvn(i)=stdifv
        end if
     end do

     allocate (fctr_dif(0:mj),dfctr_difdr(0:mj))
     fctr_dif=1.0_IDP+0.5*(AWfctr_dif-1.0)*(1.0+tanh((r-Rfctr)/Wfctr))                          ! 1 except at the edge (Rfctr) where it jumps to AWfctr_dif
     dfctr_difdr=0.5*(AWfctr_dif-1.0)*(1.0-tanh((r-Rfctr)/Wfctr)*tanh((r-Rfctr)/Wfctr))/Wfctr   ! radial derivative of fctr_dif

  end if

!!!!!!!!!!!!!   Source and sinks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate (src(0:mj),src_nf(0:mj),src_nalpha(0:mj))

! thermal plasma
  if (src_sink_th_on == 1) then

     if (srcsinkth(0) > 0.0_IDP) then
        src=srcsinkth(0)
        do j=0,mj
           rsq=r(j)*r(j)
           do i=1,10
              src(j)=src(j)+srcsinkth(i)*rsq**i
              if (src(j) < 3.e-9_IDP) src(j)=3.e-9_IDP     !to keep fit .gt. 0 at edge
           end do
        end do
     else
        do j=0,mj
           rmrbar = (r(j)-rsrc)/wsrc
           src(j)=0.
           if (abs(rmrbar) < 15.0_IDP) src(j)=asrc*exp(-rmrbar**2)
        end do
     end if

  end if

! first EP population
  if (src_sink_EP1_on == 1) then

     if (srcsinkEP1(0) > 0.0_IDP) then
        src_nf=srcsinkEP1(0)
        do j=0,mj
           rsq=r(j)*r(j)
           do i=1,10
              src_nf(j)=src_nf(j)+srcsinkEP1(i)*rsq**i
              if (src_nf(j) < 3.e-9_IDP) src_nf(j)=3.e-9_IDP     !to keep fit .gt. 0 at edge
           end do
        end do
     else
        do j=0,mj
           rmrbar = (r(j)-rsrc_EP1)/wsrc_EP1
           src_nf(j)=0.
           if (abs(rmrbar) < 15.0_IDP) src_nf(j)=asrc_EP1*exp(-rmrbar**2)
        end do
     end if

  end if

! second EP population
  if (src_sink_EP2_on == 1) then

     if (srcsinkEP2(0) > 0.0_IDP) then
        src_nalpha=srcsinkEP2(0)
        do j=0,mj
           rsq=r(j)*r(j)
           do i=1,10
              src_nalpha(j)=src_nalpha(j)+srcsinkEP2(i)*rsq**i
              if (src_nalpha(j) < 3.e-9_IDP) src_nalpha(j)=3.e-9_IDP     !to keep fit .gt. 0 at edge
           end do
        end do
     else
        do j=0,mj
           rmrbar = (r(j)-rsrc_EP2)/wsrc_EP2
           src_nalpha(j)=0.
           if (abs(rmrbar) < 15.0_IDP) src_nalpha(j)=asrc_EP2*exp(-rmrbar**2)
        end do
     end if

  end if

!  Set default case DIIID (only one EP population):

  if (src_sink_DIIID_on == 1) then

     do j=0,mj
        src_nf(j) = 2.37e-7*(1. + 1.1073*r(j) - 16.78*r(j)**2 &
                  + 40.487*r(j)**3 - 39.507*r(j)**4 + 13.817*r(j)**5) &   !NBI source - DIII-D
                  - 0.9*2.37e-7*exp(-(1. - r(j))**2/0.05)                 !edge sink - DIII-D
        if (src_nf(j) < 3.e-9_IDP) src_nf(j)=3.e-9_IDP                    !to keep fit .gt. 0 at edge
     end do

  end if

!  Set default case ITER (NBI EP and alpha particles):

  if (src_sink_ITER_on == 1) then

     do j=0,mj
        src_nf(j)=0.042017995138940*(r(j)**12)-0.252754706166547*(r(j)**11)+0.663402385523094*(r(j)**10) &
                 -0.995693681922380*(r(j)**9)+0.940423102262354*(r(j)**8)-0.579060405577752*(r(j)**7) &
                 +0.232561832976523*(r(j)**6)-0.059074775228074*(r(j)**5)+0.008831371170105*(r(j)**4) &
                 -0.000671393343300*(r(j)**3)+ 0.000017935955667*(r(j)**2) &
                 +0.000000206593558*r(j)+0.000000134257161
        if (src_nf(j) < 3.e-9_IDP) src_nf(j)=3.e-9_IDP             !to keep fit .gt. 0 at edge
        src_nalpha(j)=1.0e-04*(0.056315398901465*(r(j)**6)-0.213402569730035*(r(j)**5) &
                     +0.281146465825509*(r(j)**4)-0.141647500548349*(r(j)**3) &
                     +0.014708082536526*(r(j)**2)-0.001987451900372*r(j)+0.004806440502696)
        if (src_nalpha(j) < 3.e-9_IDP) src_nalpha(j)=3.e-9_IDP     !to keep fit .gt. 0 at edge
     end do

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vfova2 = vfova*vfova
  if (alpha_on == 1) then
     valphaova2=valphaova*valphaova
     call dbydr0(dnalpeqdr,nalpeq,0.0_IDP,1.0_IDP,0)
  end if

  if (myPE == 0 .and. nstres == 0) then
     open(unit=44,file="profiles.dat",status="unknown")
     if (ext_prof == 0) then
        if (alpha_on == 0) then
           write(44,'("r",a1,"q",a1,"J",a1,"I",a1,"preq",a1,"denseq",a1,"teeq",a1,"tieq",a1,"nfeq",a1,"vfova",a1, &
                      "vzt_eq",a1,"vth_eq",a1,"vtherm_elc",a1,"eta",a1,"iota")') (tb,i=1,14)
           do j=0,mj
              write(44,'(0pf12.8,14(a1,1pe13.5))') r(j),tb,qq(j),tb,feq(j),tb,cureq(j),tb,preq(j),tb,denseq(j),tb, &
                   teeq(j),tb,tieq(j),tb,nfeq(j),tb,vfova(j),tb,vzt_eq(j),tb,vth_eq(j),tb,vtherm_elc(j),tb,eta(j), &
                   tb,qqinv(j)
           end do
        else
           write(44,'("r",a1,"q",a1,"J",a1,"I",a1,"preq",a1,"denseq",a1,"teeq",a1,"tieq",a1,"nfeq",a1,"vfova",a1, &
                      "nalpeq",a1,"valphaova",a1,"vzt_eq",a1,"vth_eq",a1,"vtherm_elc",a1,"eta",a1,"iota")') &
                (tb,i=1,16)
           do j=0,mj
              write(44,'(0pf12.8,16(a1,1pe13.5))') r(j),tb,qq(j),tb,feq(j),tb,cureq(j),tb,preq(j),tb,denseq(j),tb, &
                   teeq(j),tb,tieq(j),tb,nfeq(j),tb,vfova(j),tb,nalpeq(j),tb,valphaova(j),tb,vzt_eq(j),tb,vth_eq(j),tb, &
                   vtherm_elc(j),tb,eta(j),tb,qqinv(j)
           end do
        end if
     else
        if (alpha_on == 0) then
           write(44,'("r",a1,"q",a1,"J",a1,"I",a1,"preq",a1,"denseq",a1,"teeq",a1,"nfeq",a1,"vfova",a1, &
                      "vzt_eq",a1,"vth_eq",a1,"vtherm_elecP",a1,"vAlfven",a1,"eta",a1,"iota")') (tb,i=1,14)
           do j=0,mj
              write(44,'(0pf12.8,14(a1,1pe13.5))') r(j),tb,qq(j),tb,feq(j),tb,cureq(j),tb,preq(j),tb,denseq(j),tb, &
                   teeq(j),tb,nfeq(j),tb,vfova(j),tb,vzt_eq(j),tb,vth_eq(j),tb,vtherm_elecP(j),tb,vAlfven(j), &
                   tb,eta(j),tb,qqinv(j)
           end do
        else
           write(44,'("r",a1,"q",a1,"J",a1,"I",a1,"preq",a1,"denseq",a1,"teeq",a1,"nfeq",a1,"vfova",a1, &
                      "nalpeq",a1,"valphaova",a1,"vzt_eq",a1,"vth_eq",a1,"vtherm_elecP",a1,"vAlfven",a1, &
                      "eta",a1,"iota")') (tb,i=1,16)
           do j=0,mj
              write(44,'(0pf12.8,16(a1,1pe13.5))') r(j),tb,qq(j),tb,feq(j),tb,cureq(j),tb,preq(j),tb,denseq(j),tb, &
                   teeq(j),tb,nfeq(j),tb,vfova(j),tb,nalpeq(j),tb,valphaova(j),tb,vzt_eq(j),tb,vth_eq(j),tb, &
                   vtherm_elecP(j),tb,vAlfven(j),tb,eta(j),tb,qqinv(j)
           end do
        end if
     end if
     close(unit=44)

!     if(ext_prof .eq. 1) then
!        open(unit=45,file="profiles_ex.dat",status="unknown")
!        write(45,'("r",a1,"dnnbi",a1,"dne",a1,"dni",a1,"tbn",a1, &
!              "ti",a1,"te",a1,"vzt_eq",a1,"temp_ep",a1,"vfova")') &
!              tb,tb,tb,tb,tb,tb,tb,tb,tb
!        do j=0,mj
!           write(45,'(e15.7,9(a1,e15.7))') r(j),tb,dnnbi(j),tb,dne(j),tb,dni(j),tb,tbn(j), &
!              tb,ti(j),tb,te(j),tb,vzt_eq(j),tb,temp_ep(j),tb,vfova(j)
!        end do
!        close(unit=45)
!     end if
  end if
  
  allocate (nskpn(nnum))
  allocate (nskp2n(nnum))
  allocate (nskpxn(nnum))
  allocate (nskpin(nnum))

  nskp2n(1)=0
  nskpxn(1)=0
  do i=1,nnum-1
     nskpxn(i+1)=nskpxn(i)+mnumn(i)*noeqn*mj
     nskp2n(i+1)=nskp2n(i)+mnumn(i)*mnumn(i)
  end do
  nbxdim=lmaxn*noeqn*mj
  mcnt2=nskp2n(nnum)+mnumn(nnum)*mnumn(nnum)

  do ip=0,numPElm1
     nskpn(n_st(ip))=0
     nskpin(n_st(ip))=0
     do i=n_st(ip),n_nd(ip)-1
        nskpn(i+1)=nskpn(i)+mnumn(i)*mnumn(i)*noeqn*noeqn*mjm1
        nskpin(i+1)=nskpin(i)+mnumn(i)*noeqn*mjm1
     end do
  end do

  if (myPE <= numPElm1) then
     nbdim=nskpn(n_nd(myPE))+mnumn(n_nd(myPE))*mnumn(n_nd(myPE))*noeqn*noeqn*mjm1
     nbidim=nskpin(n_nd(myPE))+mnumn(n_nd(myPE))*noeqn*mjm1
  end if

  if (myPE == 0) then
     write (6,'("mcnt2=",i12)') mcnt2
     write (6,'("nbdim=",i12)') nbdim
     write (6,'("nbidim=",i12)') nbidim
     write (6,'("nbxdim=",i9)') nbxdim
  else if (myPE <= numPElm1) then
     write (0,'("iPE=",i3," nbdim=",i12)') myPE,nbdim
     write (0,'("iPE=",i3," nbidim=",i12)') myPE,nbidim
  end if
  
  n_start=n_st(myPE)
  n_end=n_nd(myPE)

  allocate (lmaxPE(0:numPEsm1))

  lmaxPE=0
  do ip=0,numPElm1
     do i=n_st(ip),n_nd(ip)
        lmaxPE(ip)=lmaxPE(ip)+mnumn(i)
     end do
  end do

  if (myPE <= numPElm1) then
     allocate (scp(lmaxPE(myPE),0:mj))
  else
     allocate (scp(1,0:mj))
  end if
  allocate (scu(maxval(lmaxPE)*maxval(mj_dl),0:numPEsm1))
  
  if (myPE == 0) then
     do ip=0,numPElm1
        write(6,'("iPE =",i3," n =",2i6," lmaxPE =",i4)') ip,n_st(ip),n_nd(ip),lmaxPE(ip)
     end do
     do ip=0,numPEsm1
        write(6,'("iPE =",i3," mj_br =",i4," mj_inc =",i4," mj_st =",i4," mj_dl =",i4)') &
             ip,mj_br(ip),mj_inc(ip),mj_st(ip),mj_dl(ip)
     end do
  end if

  if (myPE < numPElm1) then
     allocate (cmamm(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqmax))
     allocate (cmamp(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqmax))
     allocate (cmapm(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqmax))
     allocate (cmapp(nskp2n(n_st(myPE))+1:nskp2n(n_nd(myPE)+1),leqmax))
     allocate (amat(nbdim))
     allocate (bmat(nbdim))
     allocate (cmat(nbdim))
     allocate (ipc(nbidim))
     allocate (xt(nskpxn(n_st(myPE))+1:nskpxn(n_nd(myPE)+1)))
     allocate (yt(nskpxn(n_st(myPE))+1:nskpxn(n_nd(myPE)+1)))
  else if (myPE == numPElm1) then
     allocate (cmamm(nskp2n(n_st(myPE))+1:mcnt2,leqmax))
     allocate (cmamp(nskp2n(n_st(myPE))+1:mcnt2,leqmax))
     allocate (cmapm(nskp2n(n_st(myPE))+1:mcnt2,leqmax))
     allocate (cmapp(nskp2n(n_st(myPE))+1:mcnt2,leqmax))
     allocate (amat(nbdim))
     allocate (bmat(nbdim))
     allocate (cmat(nbdim))
     allocate (ipc(nbidim))
     allocate (xt(nskpxn(n_st(myPE))+1:nbxdim))
     allocate (yt(nskpxn(n_st(myPE))+1:nbxdim))
  else
     allocate (cmamm(1,1))
     allocate (cmamp(1,1))
     allocate (cmapm(1,1))
     allocate (cmapp(1,1))
     allocate (amat(1))
     allocate (bmat(1))
     allocate (cmat(1))
     allocate (ipc(1))
     allocate (xt(1))
     allocate (yt(1))
  end if

  cmamm=0.
  cmamp=0.
  cmapm=0.
  cmapp=0.

  amat=0.
  bmat=0.
  cmat=0.

  ipc=0
  xt=0.
  yt=0.

!  calculate derivative weights to be used in blockj, block0, b2lx and b2lx0

  wt1m(1,1)=0.
  wt2m(1,1)=0.
  wt1m(1,2)=0.
  wt2m(1,2)=0.
  wt10(1,1)=-dc1m(1)-dc1p(1)
  wt20(1,1)=-dc2m(1)-dc2p(1)
  wt10(1,2)=r(2)**2*dc1m(1)/(r(2)**2-r(1)**2)-dc1m(1)-dc1p(1)
  wt20(1,2)=r(2)**2*dc2m(1)/(r(2)**2-r(1)**2)-dc2m(1)-dc2p(1)
  wt1p(1,1)=dc1p(1)
  wt2p(1,1)=dc2p(1)
  wt1p(1,2)=-r(1)**2*dc1m(1)/(r(2)**2-r(1)**2)+dc1p(1)
  wt2p(1,2)=-r(1)**2*dc2m(1)/(r(2)**2-r(1)**2)+dc2p(1)
  wt1m(2:mjm1,1)=dc1m(2:mjm1)
  wt2m(2:mjm1,1)=dc2m(2:mjm1)
  wt10(2:mjm1,1)=-dc1m(2:mjm1)-dc1p(2:mjm1)
  wt20(2:mjm1,1)=-dc2m(2:mjm1)-dc2p(2:mjm1)
  wt1p(2:mjm2,1)=dc1p(2:mjm2)
  wt2p(2:mjm2,1)=dc2p(2:mjm2)
  wt1p(mjm1,1)=0.
  wt2p(mjm1,1)=0.
  wt1m(2:mjm2,2)=dc1m(2:mjm2)
  wt2m(2:mjm2,2)=dc2m(2:mjm2)
  wt1m(mjm1,2)=dc1m(mjm1)-dc1p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
  wt2m(mjm1,2)=dc2m(mjm1)-dc2p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
  wt10(2:mjm2,2)=-dc1m(2:mjm2)-dc1p(2:mjm2)
  wt20(2:mjm2,2)=-dc2m(2:mjm2)-dc2p(2:mjm2)
  wt10(mjm1,2)=-dc1m(mjm1)-dc1p(mjm1)+dc1p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
  wt20(mjm1,2)=-dc2m(mjm1)-dc2p(mjm1)+dc2p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
  wt1p(2:mjm2,2)=dc1p(2:mjm2)
  wt2p(2:mjm2,2)=dc2p(2:mjm2)
  wt1p(mjm1,2)=0.
  wt2p(mjm1,2)=0.

!  read start or continue values into proper arrays

  call dlstar(uzt,phi,-1,sc1,sc2,sceq1,sceq2,0.0_IDP,1.0_IDP)
  if (myPE == numPEsm1) uzt(mj,:)=0.0

  do ind=1,nvar
     select case (ind)
        case (1)
           call trnsfr(psi,1,2)
        case (2)
           call trnsfr(phi,-1,2)
        case (3)
           call trnsfr(pr,1,2)
        case (4)
           call trnsfr(uzt,-1,2)
        case (5)
           call trnsfr(nf,1,2)
        case (6)
           call trnsfr(vprlf,-1,2)
        case (7)
           call trnsfr(vthprlf,-1,2)
        case (8)
           call trnsfr(nalp,1,2)
        case (9)
           call trnsfr(vprlalp,-1,2)
     end select
     
     l1=0
     do i=n_start,n_end
        mnum3=noeqn*mnumn(i)
        lsk=(ind-1)*mnumn(i)+nskpxn(i)
        do l1t=1,mnumn(i)
           l1=l1+1
           lskp=l1t+lsk
           do j=1,mj
              imt1=lskp+mnum3*(j-1)
              xt(imt1)=scp(l1,j)
           end do
        end do
     end do
  end do

!  calculate "c" coefficients to be used in blockj, b2lx and b2lx0

  allocate (llno(lmax))
  do l=1,lmax
     do l1=1,lmaxn
        if(lln(l1) == l) exit
     end do
     if (l1 > lmaxn) then
        llno(l)=0
     else
        llno(l)=l1
     end if
  end do

  do it1=-1,1,2

     do it2=-1,1,2

        it3=it1*it2

        do le=1,leqmax
           m1=mmeq(le)*sgnleq(le)
           n1=nneq(le)*sgnleq(le)
           is1=it1*sgnleq(le)

           do i=n_start,n_end
              do l2t=1,mnumn(i)
                 l2=l2t+lnumn(i-1)
                 l=lln(l2)
                 m2=mm(l)*signl(l)
                 n2=nn(l)*signl(l)
                 is2=it2*signl(l)

                 if (is1 == 0) then
                    lp=nskp2n(i)+mnumn(i)*(l2t-1)+l2t
                    if (it1 == 1) then
                       if (it2 == 1) then
                          cmapp(lp,le)=1.0_IDP
                       else
                          cmapm(lp,le)=1.0_IDP
                       end if
                    else
                       l3=lo(l2)
                       if (l3 > 0) then
                          l3t=l3-lnumn(i-1)
                          if (l3t <= mnumn(i)) then
                             lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                             if (it2 == 1) then
                                cmamp(lp,le)=1.0_IDP
                             else
                                cmamm(lp,le)=1.0_IDP
                             end if
                          end if
                       end if
                    end if

                 else if (is2 == 0) then
                    m3=it2*mmeq(le)
                    n3=it2*nneq(le)
                    if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) then
                       l3p=ll(m3,n3)
                       if (l3p > 0) then
                          l3=llno(l3p)
                          if (l3 > 0) then
                             l3t=l3-lnumn(i-1)
                             if (l3t <= mnumn(i)) then
                                lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                                if (it2 == 1) then
                                   if (it1 == 1) then
                                      cmapp(lp,le)=1.0_IDP
                                   else
                                      cmamp(lp,le)=1.0_IDP
                                   end if
                                else
                                   if (it1 == 1) then
                                      cmapm(lp,le)=1.0_IDP
                                   else
                                      cmamm(lp,le)=1.0_IDP
                                   end if
                                end if
                             end if
                          end if
                       end if
                    end if

                 else
                    is3=is1*is2

                    if (n1 < n2) then
                       sgn12=-1
                       mm3=m2-m1
                       nm3=n2-n1
                    else if (n1 > n2) then
                       sgn12=1
                       mm3=m1-m2
                       nm3=n1-n2
                    else
                       sgn12=sign(1,m1-m2)
                       mm3=abs(m1-m2)
                       nm3=0
                    end if
                    m3=is3*it3*mm3
                    n3=is3*it3*nm3
                    l3p=0
                    if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) l3p=ll(m3,n3)
                    l3=0
                    if (l3p > 0) l3=llno(l3p)
                    if (l3 > 0) then
                       l3t=l3-lnumn(i-1)
                       lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                       coef=0.5_IDP
                       if (is3 == -1) then
                          if (is1*sgn12 == 1) coef=-0.5_IDP
                          if (mm3 == 0 .and. nm3 == 0) coef=0.0_IDP
                       end if
                       if (it2 == 1) then
                          if (it1 == 1) then
                             cmapp(lp,le)=coef
                          else
                             cmamp(lp,le)=coef
                          end if
                       else
                          if (it1 == 1) then
                             cmapm(lp,le)=coef
                          else
                             cmamm(lp,le)=coef
                          end if
                       end if
                    end if

                    mp3=m1+m2
                    np3=n1+n2
                    m3=is3*it3*mp3
                    n3=is3*it3*np3
                    l3p=0
                    if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) l3p=ll(m3,n3)
                    l3=0
                    if (l3p > 0) l3=llno(l3p)
                    if (l3 > 0) then
                       l3t=l3-lnumn(i-1)
                       lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                       coef=0.5_IDP
                       if (is1 == -1 .and. is2 == -1) coef=-0.5_IDP
                       if (it2 == 1) then
                          if (it1 == 1) then
                             cmapp(lp,le)=coef
                          else
                             cmamp(lp,le)=coef
                          end if
                       else
                          if (it1 == 1) then
                             cmapm(lp,le)=coef
                          else
                             cmamm(lp,le)=coef
                          end if
                       end if
                    end if

                 end if
              end do
           end do
        end do
     end do
  end do

!   Ion FLR effects auxiliary variable

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then
     call dlsq(sc1,psi,1,sc2,sc3,0.0_IDP,1.0_IDP)
     call trnsfr(sc1,1,2)
     l1=0
     do i=n_start,n_end
        mnum3=noeqn*mnumn(i)
        lsk=(iq-1)*mnumn(i)+nskpxn(i)
        do l1t=1,mnumn(i)
           l1=l1+1
           lskp=l1t+lsk
           do j=1,mj
              imt1=lskp+mnum3*(j-1)
              xt(imt1)=scp(l1,j)
           end do
        end do
     end do
  end if

!   EP FLR effects auxiliary variables

  if (epflr_on > 0) then

     call dlsq(sc1,phi,-1,sc2,sc3,0.0_IDP,1.0_IDP)
     call dbydzt_par(sc2,psi,1,0.0_IDP,1.0_IDP)
     call dbydth_par(sc4,psi,1,0.0_IDP,1.0_IDP,0)
     do l=1,lmax
        sc3(:,l)=r(mj_start:mj_end)*sc4(:,l)
     end do

     allocate (xa(0:mj))
     allocate (scp2(lmaxPE(myPE),0:mj),scp3(lmaxPE(myPE),0:mj))
     call trnsfr(sc2,-1,2)
     scp2=scp
     call trnsfr(sc3,-1,2)
     scp3=scp
     call trnsfr(sc1,-1,2)

     lsk=0
     do i=n_start,n_end
        allocate (amatw(mnumn(i),mnumn(i),mjm1),bmatw(mnumn(i),mnumn(i),mjm1),cmatw(mnumn(i),mnumn(i),mjm1))
        allocate (xw(mnumn(i),mjm1))
        allocate (ipcw(mnumn(i),mjm1))

        amatw=0.0_IDP
        bmatw=0.0_IDP
        cmatw=0.0_IDP
        mnum3=noeqn*mnumn(i)

        do l2t=1,mnumn(i)
           l2=l2t+lnumn(i-1)
           m=mm(lln(l2))
           n=nn(lln(l2))
           ibnd=1
           if (m == 0) ibnd=2
           do l3t=1,mnumn(i)
              do lq=1,leqmax
                 lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                 x=cmapm(lp,lq)
                 if (x == 0.0_IDP) cycle
                 xa=x
                 do j=1,mjm2
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt20(j,ibnd)*xa(j)*lplrr(j,lq)
                    cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt2m(j,ibnd)*xa(j)*lplrr(j,lq)
                    bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt2p(j,ibnd)*xa(j)*lplrr(j,lq)
                 end do
                 amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt20(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
                 cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt2m(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
                 bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt2p(mjm1,1)*xa(mjm1)*lplrr(mjm1,lq)
                 do j=1,mjm2
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplr(j,lq)
                    cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplr(j,lq)
                    bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplr(j,lq)
                 end do
                 amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
                 cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
                 bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplr(mjm1,lq)
                 xa=-x*m*m*rinv*rinv
                 do j=1,mjm1
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltt(j,lq)
                 end do
                 xa=-x*m*n*rinv
                 do j=1,mjm1
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lpltz(j,lq)
                 end do
                 xa=-x*n*n
                 do j=1,mjm1
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplzz(j,lq)
                 end do
              end do
              do lq=1,leqmax
                 lp=nskp2n(i)+mnumn(i)*(l2t-1)+l3t
                 x=cmamp(lp,lq)
                 if (x == 0.0_IDP) cycle
                 xa=x*m*rinv
                 do j=1,mjm2
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrt(j,lq)
                    cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrt(j,lq)
                    bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrt(j,lq)
                 end do
                 amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
                 cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
                 bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrt(mjm1,lq)
                 do j=1,mjm1
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplt(j,lq)
                 end do
                 xa=x*n
                 do j=1,mjm2
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+wt10(j,ibnd)*xa(j)*lplrz(j,lq)
                    cmatw(l3t,l2t,j)=cmatw(l3t,l2t,j)+wt1m(j,ibnd)*xa(j)*lplrz(j,lq)
                    bmatw(l3t,l2t,j)=bmatw(l3t,l2t,j)+wt1p(j,ibnd)*xa(j)*lplrz(j,lq)
                 end do
                 amatw(l3t,l2t,mjm1)=amatw(l3t,l2t,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
                 cmatw(l3t,l2t,mjm1)=cmatw(l3t,l2t,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
                 bmatw(l3t,l2t,mjm1)=bmatw(l3t,l2t,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrz(mjm1,lq)
                 do j=1,mjm1
                    amatw(l3t,l2t,j)=amatw(l3t,l2t,j)+xa(j)*lplz(j,lq)
                 end do
              end do
           end do
        end do

        if (alpha_on == 1 .and. (r_epflralp > 0.0 .or. epflr_on == 2)) then
           allocate (amatwalp(mnumn(i),mnumn(i),mjm1),bmatwalp(mnumn(i),mnumn(i),mjm1),cmatwalp(mnumn(i),mnumn(i),mjm1))
           allocate (ipcwalp(mnumn(i),mjm1))
           if (epflr_on == 1) then
              coef=-r_epflralp*r_epflralp
              amatwalp=coef*amatw
              bmatwalp=coef*bmatw
              cmatwalp=coef*cmatw
           else if (epflr_on == 2) then
              sd1=-valphaova2/(epsq*omcydalp*omcydalp)
              do j=1,mjm1
                 amatwalp(:,:,j)=sd1(j)*amatw(:,:,j)
                 bmatwalp(:,:,j)=sd1(j)*bmatw(:,:,j)
                 cmatwalp(:,:,j)=sd1(j)*cmatw(:,:,j)
              end do
           end if

           bmatwalp(:,:,mjm1)=0.0_IDP
           cmatwalp(:,:,1)=0.0_IDP
           call decbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,ipcwalp,ier)
           if (ier /= 0) stop 18

           if (epflr_on == 1) then
              do j=1,mjm1
                 xw(:,j)=coef*scp(lsk+1:lsk+mnumn(i),j)
              end do
           else if (epflr_on == 2) then
              do j=1,mjm1
                 xw(:,j)=sd1(j)*scp(lsk+1:lsk+mnumn(i),j)
              end do
           end if
           call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
           do l=1,mnumn(i)
              lskp=l+(iwa-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do

           do j=1,mjm1
              xw(:,j)=scp2(lsk+1:lsk+mnumn(i),j)
           end do
           call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
           do l=1,mnumn(i)
              lskp=l+(ix1a-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do

           do j=1,mjm1
              xw(:,j)=scp3(lsk+1:lsk+mnumn(i),j)
           end do
           call solbt(mnumn(i),mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
           do l=1,mnumn(i)
              lskp=l+(ix2a-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do
           deallocate (amatwalp,bmatwalp,cmatwalp,ipcwalp)
        end if

        if (r_epflr > 0.0 .or. epflr_on == 2) then
           if (epflr_on == 1) then
              coef=-r_epflr*r_epflr
              amatw=coef*amatw
              bmatw=coef*bmatw
              cmatw=coef*cmatw
           else if (epflr_on == 2) then
              sd1=-vfova2/(epsq*omcyd*omcyd)
              do j=1,mjm1
                 amatw(:,:,j)=sd1(j)*amatw(:,:,j)
                 bmatw(:,:,j)=sd1(j)*bmatw(:,:,j)
                 cmatw(:,:,j)=sd1(j)*cmatw(:,:,j)
              end do
           end if
           do l=1,mnumn(i)
              do j=1,mjm1
                 amatw(l,l,j)=1.0+amatw(l,l,j)
              end do
           end do

           bmatw(:,:,mjm1)=0.0_IDP
           cmatw(:,:,1)=0.0_IDP

           call decbt(mnumn(i),mjm1,amatw,bmatw,cmatw,ipcw,ier)
           if (ier /= 0) stop 17

           if (epflr_on == 1) then
              do j=1,mjm1
                 xw(:,j)=coef*scp(lsk+1:lsk+mnumn(i),j)
              end do
           else if (epflr_on == 2) then
              do j=1,mjm1
                 xw(:,j)=sd1(j)*scp(lsk+1:lsk+mnumn(i),j)
              end do
           end if
           call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
           do l=1,mnumn(i)
              lskp=l+(iw-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do

           do j=1,mjm1
              xw(:,j)=scp2(lsk+1:lsk+mnumn(i),j)
           end do
           call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
           do l=1,mnumn(i)
              lskp=l+(ix1-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do

           do j=1,mjm1
              xw(:,j)=scp3(lsk+1:lsk+mnumn(i),j)
           end do
           call solbt(mnumn(i),mjm1,amatw,bmatw,cmatw,xw,ipcw)
           do l=1,mnumn(i)
              lskp=l+(ix2-1)*mnumn(i)+nskpxn(i)
              do j=1,mjm1
                 imt1=lskp+mnum3*(j-1)
                 xt(imt1)=xw(l,j)
              end do
           end do
        end if
        lsk=lsk+mnumn(i)
        deallocate (amatw,bmatw,cmatw,xw,ipcw)
     end do
     deallocate (xa,scp2,scp3)

  end if

! put in r.h.s. of equations

! psi equation

  sd1=1.0_IDP
  call om0(sd1,1,2,0,0,0,1.0_IDP)
  call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*sceq2(:,l)
  end do
  call blockj(sceq1,-1,1,1,1,0,0,oneos)
  call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*sceq2(:,l)
  end do
  call blockj(sceq1,1,1,1,0,1,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*grroj(:,l)
  end do
  call blockj(sceq1,1,1,1,2,0,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=-2.*eta*feq*grtoj(:,l)
  end do
  call blockj(sceq1,-1,1,1,1,1,0,oneos)
  do l=1,leqmax
     sceq1(:,l)=eta*feq*gttoj(:,l)
  end do
  call blockj(sceq1,1,1,1,0,2,0,oneos)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call blockj(sceq2,1,1,1,0,1,0,oneos)

! Ion FLR effects

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then

     if (ext_prof == 1) then
        if (iflr_on == 1) then
           sd1=iflr*iflr
        else
           sd1=vthi*vthi*tieq/(epsq*omcyd*omcyd)
        end if
        do l=1,leqmax
           sceq5(:,l)=1.2533*sd1*vAlfven*vAlfven*bmod(:,l)/(vtherm_elecP*(feq-qqinv*cureq))
        end do
     else
        do l=1,leqmax
           sceq5(:,l)=1.2533*iflr*iflr*bmod(:,l)/(denseq*vtherm_elc*(feq-qqinv*cureq))
        end do
     end if
     call blockj_landau_grad_parallel(sceq5,1,1,iq,0,0,0,1.0_IDP)

!  Ion FLR effects auxiliary equation

     sd1=1.0_IDP
     call block_dlsq(iq,1,sd1)
     call block0(sd1,iq,iq,0,0,0,-1.0_IDP)

  end if

! Two fluid terms

  if (twofl_on == 1) then
     sd1=beteom/denseq
     call om0(sd1,1,3,0,0,0,-1.0_IDP)
  end if

! u-zeta equation

  call dbydreq(sceq1,sqg,0.0_IDP,-1.0_IDP,0)
  call blockj(sceq1,1,2,3,1,0,0,betfc)
  call dbydtheq(sceq2,sqg,1,0.0_IDP,1.0_IDP,0)
  call blockj(sceq2,-1,2,3,0,1,0,betfc)
  if (myPE == 0) then
     write(*,'(/"bet0 = ",1pe12.4,"   betf0_f = ",e12.4)') bet0, bet0_f
     write(*,'("eps = <a>/<R> = ",1pe12.4,"   S = tau_R/tau_A = ",e12.4)') eps, s
     write(*,'("betfc = ",1pe12.4,"   betfc_f = ",e12.4)') betfc, betfc_f
  end if

! fast ion coupling

  call blockj(sceq1,1,2,5,1,0,0,betfc_f)
  call blockj(sceq2,-1,2,5,0,1,0,betfc_f)

  if (alpha_on == 1) then
     call blockj(sceq1,1,2,8,1,0,0,betfc_alp)
     call blockj(sceq2,-1,2,8,0,1,0,betfc_alp)
  end if

! Shared equilibrium toroidal flow velocity for u-zeta equation

  call block0(vzt_eq,2,4,0,0,1,-1.0_IDP)

  call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
  call om(sceq2,-1,2,1,1,0,0,1.0_IDP)
  call grpareq(sceq1,sceq2,-1,0.0_IDP,1.0_IDP)
  call blockj(sceq1,1,2,1,1,0,0,1.0_IDP)
  call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
  call om(sceq2,1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,sceq2,1,0.0_IDP,1.0_IDP)
  call blockj(sceq1,-1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,grtoj,-1,0.0_IDP,-2.0_IDP)
  call blockj(sceq1,1,2,1,1,1,0,1.0_IDP)
  call grpareq(sceq1,gttoj,1,0.0_IDP,1.0_IDP)
  call blockj(sceq1,-1,2,1,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)
  call grpareq(sceq1,grroj,1,0.0_IDP,1.0_IDP)
  call blockj(sceq1,-1,2,1,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-grtoj(:,l)
  end do
  call om(sceq1,-1,2,1,1,1,0,2.0_IDP*1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=grroj(:,l)
  end do
  call om(sceq1,1,2,1,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=gttoj(:,l)
  end do
  call om(sceq1,1,2,1,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=rinv*sceq1(:,l)
  end do
  call om(sceq2,1,2,1,0,1,0,1.0_IDP)
  call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)
  sd2=rinv*sd1
  call d2bydr20(sd2,cureq,1.0_IDP,-1.0_IDP,0)
  sd1=rinv*sd2/epsq
  call block0(sd1,2,1,1,0,0,1.0_IDP)
  call dbydtheq(sceq2,bst,-1,0.0_IDP,-1.0_IDP,0)
  do l=1,leqmax
     sceq1(:,l)=r*sceq2(:,l)/epsq
  end do
  call dbydreq(sceq2,sceq1,0.0_IDP,-1.0_IDP,0)
  call blockj(sceq2,1,2,1,1,0,0,1.0_IDP)
  call dbydtheq(sceq2,sceq1,1,0.0_IDP,1.0_IDP,0)
  call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)

! Ion FLR effects

  ! if (iflr_on == 1) then
  !    sd1=omegar*iflr*iflr
  !    call block_dlsq(2,4,sd1)
  ! end if

! Electron-ion Landau damping

  if (ieldamp_on == 1) then   
     coef=-(1-dpres)*bet0/(2.*epsq*omegar)
     cmplx1 = cmplx(0.0_IDP,1.0_IDP)                                               !! complex number i 
     tauie = tinn(0)/tenn(0)
     xnuelc0 = xnuelc0/omegar                                                      !! collision FR normalized to omegar
     xnuion0 = 0.0165*xnuelc0/(tauie*sqrt(uion*tauie))                             !! ion-ion collision FR axis

     allocate (eilnd(mjm1,lmaxn))
     eilnd=0.0_IDP
     do l=1,lmaxn     
        if (lln(l) == l0) cycle
        do j=1,mjm1
           vther = vthe*sqrt(teeq(j))
           vthir = vthi*sqrt(tieq(j))
           if(vthir .le. 0.) vthir=.001
           if(vther .le. 0.) vther=.01
           if(teeq(j) .le. 0.) teeq(j)=.01
           tauie = tinn(j)/tenn(j)
           sd1(j) = tauie                                                !! sd1 is tau_i = T_i/T_e
           xnuelc = xnuelc0*dni(j)/(teeq(j)*sqrt(teeq(j)))               !! electron-ion collision FR
           xnuion = xnuion0*dni(j)/(tieq(j)*sqrt(tieq(j)))               !! ion-ion collision FR
           abkprl = abs(nn(lln(l)) - mm(lln(l))*qqinv(j))                !! parallel gradient
           if(abkprl .lt. 1.e-50_IDP) then
              zetaiinv = (sqrt(2.0_IDP)*vthir*abkprl)/(omegar*(1. + cmplx1*xnuion))     
              zetaeinv = (sqrt(2.0_IDP)*vther*abkprl)/(omegar*(1. + cmplx1*xnuelc))
              y0i = -(1./(1.+cmplx1*xnuion))*(1. + 0.5*(zetaiinv**2))
              y1i = -(1./(1.+cmplx1*xnuion))*(1. + (zetaiinv**2))
              y2i = -(1.75/(1.+cmplx1*xnuion))*(1. +(23./14.)*(zetaiinv**2))
              y0e = -(1./(1.+cmplx1*xnuelc))*(1. + 0.5*(zetaeinv**2))
              y1e = -(1./(1.+cmplx1*xnuelc))*(1. + (zetaeinv**2))
              y2e = -(1.75/(1.+cmplx1*xnuelc))*(1. +(23./14.)*(zetaeinv**2))
              ddii = 1. + cmplx1*xnuion*y0i
              ddee = 1. + cmplx1*xnuelc*y0e
           else
              xsii = omegar/(sqrt(2.)*vthir*abkprl)
              xsie = xsii*vthir/vther
              xsi2 = xsii*xsii
              xsi3 = xsii*xsi2
              xsi4 = xsii*xsi3
              zetai = (1.+cmplx1*xnuion)*xsii
              zetae = (1.+cmplx1*xnuelc)*xsie
              zetai2 = zetai*zetai
              zetai3 = zetai*zetai2
              zetai4 = zetai*zetai3
              if(abs(zetai) .le. 10.0_IDP) then
                 ztr = real(zetai)
                 zti = aimag(zetai)
                 call zzdisp(ztr,zti,zir,zii)
                 zi = cmplx(zir,zii)
                 y0i = xsii*zi
                 y1i = xsii*(zetai+(0.5+zetai2)*zi)
                 y2i = xsii*(1.5*zetai + zetai3 +(0.5 + zetai2 + zetai4)*zi)
              else
                 y0i = -(xsii/zetai)*(1. + (0.5/zetai2))
                 y1i = -(xsii/zetai)*(1. + (1./zetai2))
                 y2i = -1.75*(xsii/zetai)*(1. + (23./(14.*zetai2)))
              endif
              ddii = 1. + cmplx1*xnuion*y0i
              zetae2 = zetae*zetae
              zetae3 = zetae*zetae2
              zetae4 = zetae*zetae3
              if(abs(zetae) .le. 10.0_IDP) then
                 ztr = real(zetae)
                 zti = aimag(zetae)
                 call zzdisp(ztr,zti,zer,zei)
                 ze = cmplx(zer,zei)
                 y0e = xsie*ze
                 y1e = xsie*(zetae+(0.5+zetae2)*ze)
                 y2e = xsie*(1.5*zetae + zetae3 +(0.5 + zetae2 + zetae4)*ze)
              else
                 y0e = -(xsie/zetae)*(1. + (0.5/zetae2))
                 y1e = -(xsie/zetae)*(1. + (1./zetae2))
                 y2e = -1.75*(xsie/zetae)*(1. + (23./(14.*zetae2)))
              endif
              ddee = 1. + cmplx1*xnuelc*y0e
           endif
           stfe = y0e/ddee
           stfi = (1. + y0i/ddii)/tauie
           reii = 1. + stfe + stfi
           rei = 1./reii
           sei1 = y2e + y2i*tauie
           stfe = (y1e*y1e)/ddee
           stfi = (y1i*y1i)/ddii
           sei2 = -cmplx1*(xnuion*stfi + xnuelc*stfe)
           sei3 = rei*(((y1e/ddee) - (y1i/ddii))**2)
           sei = sei1 + sei2 + sei3 
           eilnd(j,l) = coef*aimag(sei)  
        end do
     end do
     call blockjl(eildrr,1,2,2,0,2,0,eilnd)
     call blockjl(eildrt,-1,2,2,1,1,0,eilnd)
     call blockjl(eildrz,-1,2,2,0,1,1,eilnd)
     call blockjl(eildtt,1,2,2,2,0,0,eilnd)
     call blockjl(eildtz,1,2,2,1,0,1,eilnd)
     call blockjl(eildzz,1,2,2,0,0,2,eilnd)
     call blockjl(eildr,1,2,2,0,1,0,eilnd)
     call blockjl(eildt,-1,2,2,1,0,0,eilnd)
     call blockjl(eildz,-1,2,2,0,0,1,eilnd)
  end if

! Two fluid terms

  if (twofl_on == 1) then
     sd2=feq-qqinv*cureq
     do l=1,leqmax
        sceq1(:,l)=feq*dpreqdr*jbgrr(:,l)/sd2
        sceq2(:,l)=feq*dpreqdr*jbgrt(:,l)/sd2
        sceq3(:,l)=feq*dpreqdr*jbgtt(:,l)/sd2
     end do
     call blockj(sceq1,1,2,2,3,0,0,-betiom)
     call blockj(sceq2,-1,2,2,2,1,0,2*betiom)
     call blockj(sceq3,1,2,2,1,2,0,-betiom)
     call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,3)
     call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,3)
     do l=1,leqmax
        sceq4(:,l)=sceq4(:,l)+rinv*sceq2(:,l)
     end do
     call blockj(sceq4,-1,2,2,2,0,0,-betiom)
     call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,3)
     call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,3)
     call blockj(sceq4,1,2,2,1,1,0,betiom)
     do l=1,leqmax
        sceq1(:,l)=rinv*cureq*dpreqdr*jbgrr(:,l)/sd2
        sceq2(:,l)=rinv*cureq*dpreqdr*jbgrt(:,l)/sd2
        sceq3(:,l)=rinv*cureq*dpreqdr*jbgtt(:,l)/sd2
     end do
     call blockj(sceq1,1,2,2,2,0,1,betiom)
     call blockj(sceq2,-1,2,2,1,1,1,-2*betiom)
     call blockj(sceq3,1,2,2,0,2,1,betiom)
     call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,4)
     call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,4)
     call blockj(sceq4,-1,2,2,1,0,1,betiom)
     call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,4)
     call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,4)
     do l=1,leqmax
        sceq4(:,l)=sceq4(:,l)-rinv*sceq3(:,l)
     end do
     call blockj(sceq4,1,2,2,0,1,1,-betiom)
     do l=1,leqmax
        sceq1(:,l)=rinv*cureq*dpreqdr*dgrrz(:,l)
        sceq2(:,l)=rinv*cureq*dpreqdr*dgrtz(:,l)
        sceq3(:,l)=rinv*cureq*dpreqdr*dgttz(:,l)
     end do
     call blockj(sceq1,-1,2,2,2,0,0,0.5*betiom)
     call blockj(sceq2,1,2,2,1,1,0,-betiom)
     call blockj(sceq3,-1,2,2,0,2,0,0.5*betiom)
     sd3=rinv*feq*dpreqdr
     sd3(0)=2*feq(0)*(preq(1)-preq(0))/(r(1)*r(1))
     do l=1,leqmax
        sceq1(:,l)=sd3*dgrrt(:,l)-rinv*cureq*dpreqdr*dgrrz(:,l)
        sceq2(:,l)=feq*dpreqdr*dgttr(:,l)+2*sd3*jbgtt(:,l)/sd2-rinv*cureq*dpreqdr*dgrtz(:,l)
     end do
     call dbydtheq(sceq3,sceq1,-1,0.0_IDP,1.0_IDP,2)
     call dbydreq(sceq3,sceq2,1.0_IDP,-1.0_IDP,2)
     call dbydr0(sd4,qqinv,0.0_IDP,1.0_IDP,0)
     call dbydr0(sd5,cureq,0.0_IDP,1.0_IDP,0)
     call dbydreq(sceq2,jsq,0.0_IDP,0.5_IDP,0)
     do l=1,leqmax
        sceq1(:,l)=cureq*dpreqdr*(qqinv*(rinv*dgrtt(:,l)-dgttr(:,l))-(sd4+2*rinv*qqinv)*jbgtt(:,l)/sd2+ &
                   rinv*(dbsjtbj(:,l)-rinv*(cureq*sceq2(:,l)+sd5*jsq(:,l))/sd2)/(eps*eps))
     end do
     call dbydreq(sceq3,sceq1,1.0_IDP,-1.0_IDP,4)
     call blockj(sceq3,1,2,2,1,0,0,-0.5*betiom)
     call dbydtheq(sceq3,sceq1,1,0.0_IDP,1.0_IDP,4)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*(feq*dgttt(:,l)-cureq*dgttz(:,l))
        sceq2(:,l)=2*sd3*(dgrtt(:,l)-jbgtt(:,l)/sd2)-dpreqdr*(feq*dgttr(:,l)+rinv*cureq*dgrtz(:,l))
     end do
     call dbydreq(sceq5,sceq1,0.0_IDP,1.0_IDP,3)
     do l=1,leqmax
        sceq3(:,l)=sceq3(:,l)+rinv*sceq5(:,l)
     end do
     call dbydtheq(sceq3,sceq2,1,1.0_IDP,-1.0_IDP,2)
     call blockj(sceq3,-1,2,2,0,1,0,-0.5*betiom)
  end if

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call block0_dlsq(sd1,2,4,stdifu)
  else
     call block0_dlsqnr(fctr_dif,2,4,stdifun,dfctr_difdr)
  end if

! p equation

  call block0(dpreqdr,3,2,1,0,0,1.0_IDP)

!  do l=1,leqmax
!     sceq1(:,l)=preq*djroj(:,l)
!  end do
!  call blockj(sceq1,1,3,2,1,0,0,gamma)
!  do l=1,leqmax
!     sceq1(:,l)=-preq*djtoj(:,l)
!  end do
!  call blockj(sceq1,-1,3,2,0,1,0,gamma)

  sd2=feq/(feq-qqinv*cureq)
  sd3=cureq/(feq-qqinv*cureq)
  do l=1,leqmax
     sceq1(:,l)=-preq*(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
  end do
  call blockj(sceq1,-1,3,2,0,1,0,gamma)
  do l=1,leqmax
     sceq1(:,l)=-preq*(r*dbsjzoj(:,l)-sd2*djroj(:,l))
  end do
  call blockj(sceq1,1,3,2,1,0,0,gamma)
  do l=1,leqmax
     sceq1(:,l)=-preq*(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
  end do
  call blockj(sceq1,1,3,2,0,0,1,gamma)
  call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
  sd2=preq*sd4
  call block0(sd2,3,2,1,0,0,gamma)
  call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
  sd2=-rinv*preq*sd4
  call block0(sd2,3,2,0,0,1,gamma)

! parallel thermal velocity term  

  do l=1,leqmax
     sceq1(:,l)=preq*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call om(sceq1,1,3,7,0,0,0,-gamma)
  call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=preq*sceq2(:,l)/(feq-qqinv*cureq)
  end do
  call blockj(sceq1,-1,3,7,0,0,0,gamma)

! Shared equilibrium toroidal flow velocity for pressure equation

  call block0(vzt_eq,3,3,0,0,1,-1.0_IDP)

! Two fluid terms

  if (twofl_on == 1) then
     sd2=gamma*betiom*preq/denseq
     sd3=sd2/(feq-qqinv*cureq)
     call dbydr0(sd4,feq,0.0_IDP,1.0_IDP,0)
     sd4=sd3*sd4
     call block0(sd4,3,3,1,0,0,1.0_IDP)
     call dbydr0(sd4,cureq,0.0_IDP,1.0_IDP,0)
     sd4=rinv*sd3*sd4
     call block0(sd4,3,3,0,0,1,-1.0_IDP)
     call dbydzteq(sceq1,bst,-1,0.0_IDP,1.0_IDP)
     call dbydtheq(sceq2,bst,-1,0.0_IDP,1.0_IDP,0)
     do l=1,leqmax
        sceq1(:,l)=r*sd3*sceq1(:,l)
        sceq2(:,l)=r*sd3*sceq2(:,l)
     end do
     call blockj(sceq1,1,3,3,1,0,0,-1.0_IDP)
     call blockj(sceq2,1,3,3,0,0,1,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=sd2*omdr(:,l)
        sceq2(:,l)=sd2*omdt(:,l)
        sceq3(:,l)=sd2*omdz(:,l)
     end do
     call blockj(sceq1,-1,3,3,0,1,0,-2.0_IDP)
     call blockj(sceq2,1,3,3,1,0,0,-2.0_IDP)
     call blockj(sceq3,1,3,3,0,0,1,-2.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*sd3*grtoj(:,l)
        sceq2(:,l)=dpreqdr*sd3*gttoj(:,l)
     end do
     call om(sceq1,-1,3,1,1,0,0,-epsq)
     call om(sceq2,1,3,1,0,1,0,epsq)
     sd4=rinv*dpreqdr*sd3*cureq
     call block0(sd4,3,1,1,1,0,-1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=r*dpreqdr*sd3*bst(:,l)
     end do
     call blockj(sceq1,-1,3,1,2,0,0,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=dpreqdr*sd2*dgrtp(:,l)
        sceq2(:,l)=dpreqdr*sd2*dgttp(:,l)
     end do
     call blockj(sceq1,1,3,1,1,0,0,-epsq)
     call blockj(sceq2,-1,3,1,0,1,0,epsq)
     do l=1,leqmax
        sceq1(:,l)=rinv*dpreqdr*sd3*cureq*djtoj(:,l)
        sceq2(:,l)=dpreqdr*sd2*dbsjtoj(:,l)
     end do
     call blockj(sceq1,-1,3,1,0,1,0,-1.0_IDP)
     call blockj(sceq2,1,3,1,1,0,0,1.0_IDP)
  end if

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call block0_dlsq(sd1,3,3,stdifp)
  else
     call block0_dlsqnr(fctr_dif,3,3,stdifpn,dfctr_difdr)
  end if

! u-zeta expression

  do l=1,leqmax
     sceq2(:,l)=denseq*jbgrt(:,l)
  end do
  call blockj(sceq2,-1,4,2,1,1,0,-2.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=denseq*jbgrr(:,l)
  end do
  call blockj(sceq2,1,4,2,2,0,0,1.0_IDP)
  do l=1,leqmax
     sceq2(:,l)=denseq*jbgtt(:,l)
  end do
  call blockj(sceq2,1,4,2,0,2,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=rinv*sceq2(:,l)
  end do
  call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
  call clgam(sceq2,1,1,2,0.0_IDP,-1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=denseq*sceq2(:,l)-denseqr*jbgrt(:,l)
  end do
  call blockj(sceq1,-1,4,2,1,0,0,1.0_IDP)
  call clgam(sceq2,1,2,2,0.0_IDP,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=denseq*sceq2(:,l)+denseqr*jbgtt(:,l)
  end do
  call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
  sd1=-1.0_IDP
  call block0(sd1,4,4,0,0,0,1.0_IDP)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   NBI particle effects   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! Load 1st Omega-d terms in fast ion density equation:
  if (myPE == 0) write(*,'(/"vfova2 = ",1pe12.4)') vfova2(1)

  if (trapped_on .eq. 0) then
     do l=1,leqmax
        sceq1(:,l)=vfova2*omdr(:,l)/(epsq*omcyd)
        sceq2(:,l)=vfova2*omdt(:,l)/(epsq*omcyd)
        sceq3(:,l)=vfova2*omdz(:,l)/(epsq*omcyd)
     end do
  else
     do l=1,leqmax
        sceq1(:,l)=vfova2*omdrprp(:,l)/(epsq*omcyd)
        sceq2(:,l)=vfova2*omdtprp(:,l)/(epsq*omcyd)
        sceq3(:,l)=vfova2*omdzprp(:,l)/(epsq*omcyd)
     end do
  end if
  
  call blockj(sceq1,-1,5,5,0,1,0,-1.0_IDP)
  call blockj(sceq2,1,5,5,1,0,0,-1.0_IDP)
  call blockj(sceq3,1,5,5,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion density equation

  call block0(vzt_eq,5,5,0,0,1,-1.0_IDP)

! Load 1st Omega-d terms in fast ion parallel velocity equation:

  call blockj(sceq1,-1,6,6,0,1,0,-1.0_IDP)
  call blockj(sceq2,1,6,6,1,0,0,-1.0_IDP)
  call blockj(sceq3,1,6,6,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

  call block0(vzt_eq,6,6,0,0,1,-1.0_IDP)

! Load 2nd Omega-d terms in fast ion density equation:
  if (myPE == 0) write(*,'("nfeq = ",1pe12.4)') nfeq(1)

  if (trapped_on .eq. 0) then
     do l=1,leqmax
        sceq1(:,l)=nfeq(:)*omdr(:,l)
        sceq2(:,l)=nfeq(:)*omdt(:,l)
        sceq3(:,l)=nfeq(:)*omdz(:,l)
     end do
  else
     do l=1,leqmax
        sceq1(:,l)=nfeq(:)*omdrprp(:,l)/vfova
        sceq2(:,l)=nfeq(:)*omdtprp(:,l)/vfova
        sceq3(:,l)=nfeq(:)*omdzprp(:,l)/vfova
     end do
  end if

  call blockj(sceq1,-1,5,2,0,1,0,-1.0_IDP)
  call blockj(sceq2,1,5,2,1,0,0,-1.0_IDP)
  call blockj(sceq3,1,5,2,0,0,1,-1.0_IDP)

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call block0_dlsq(sd1,5,5,stdifnf)
  else
     call block0_dlsqnr(fctr_dif,5,5,stdifnfn,dfctr_difdr)
  end if

! Load remaining terms in fast ion density equation

! Parallel gradient term

  do l=1,leqmax
     sceq1(:,l)=nfeq*bmod(:,l)/(feq-qqinv*cureq)
  end do
  
  call om(sceq1,1,5,6,0,0,0,-1.0_IDP)

! Omega* term

  sd1=dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
  sd2=dnfeqdr*feq/(feq-qqinv*cureq)
  
  call block0(sd1,5,2,0,0,1,1.0_IDP)
  call block0(sd2,5,2,1,0,0,-1.0_IDP)

! EP FLR effects

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then
     sd1=-dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
     sd2=-dnfeqdr*feq/(feq-qqinv*cureq)
     call block0(sd1,5,iw,0,0,1,1.0_IDP)
     call block0(sd2,5,iw,1,0,0,-1.0_IDP)
     ! sd1=epsq*omcyd*omegar*nfeq/vfova2
     ! call block0(sd1,5,iw,0,0,0,1.0_IDP)
  end if

! Load remaining terms in fast ion parallel velocity equation

! Landau closure term

  do l=1,leqmax
     sceq1(:,l)=1.414213*LcA1*vfova*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call blockj_landau_grad_parallel(sceq1,1,6,6,0,0,0,1.0_IDP)

! Parallel gradient terms

  do l=1,leqmax
     sceq1(:,l)=2*LcA0*vfova2*bmod(:,l)/(nfeq*(feq-qqinv*cureq))
  end do
  call om(sceq1,1,6,5,0,0,0,-1.0_IDP)

  do l=1,leqmax
     sceq1(:,l)=2*LcA0*epsq*omcyd*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call om(sceq1,1,6,2,0,0,0,-1.0_IDP)

! Omega* term

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then

!  EP FLR effects

     sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
     sd2=vfova2*dnfeqdr*rinv*feq/(nfeq*(feq-qqinv*cureq))
     call block0(sd1,6,ix1,0,0,0,1.0_IDP)
     call block0(sd2,6,ix2,0,0,0,-1.0_IDP)

!  EP FLR effects auxiliary equations

     sd1=1.0_IDP
     call block0(sd1,iw,iw,0,0,0,1.0_IDP)
     if (epflr_on == 1) then
        sd2=-r_epflr*r_epflr
     else if (epflr_on == 2) then
        sd2=-vfova2/(epsq*omcyd*omcyd)
     end if
     call block_dlsq(iw,iw,sd2)
     call block_dlsq(iw,2,-sd2)

     call block0(sd1,ix1,ix1,0,0,0,1.0_IDP)
     call block_dlsq(ix1,ix1,sd2)
     call block0(sd1,ix1,1,0,0,1,-1.0_IDP)

     call block0(sd1,ix2,ix2,0,0,0,1.0_IDP)
     call block_dlsq(ix2,ix2,sd2)
     call block0(r,ix2,1,1,0,0,-1.0_IDP)

  else 

     sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
     sd2=vfova2*dnfeqdr*feq/(nfeq*(feq-qqinv*cureq))
     call block0(sd1,6,1,0,0,1,1.0_IDP)
     call block0(sd2,6,1,1,0,0,-1.0_IDP)
   
  end if 

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call block0_dlsq(sd1,6,6,stdifvf)
  else
     call block0_dlsqnr(fctr_dif,6,6,stdifvfn,dfctr_difdr)
  end if

! End of fast ion moment equations

! Load terms of the thermal moment of energetic particles

! Pressure gradient term

  do l=1,leqmax
     sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
  end do
  call om(sceq1,1,7,3,0,0,0,-1.0_IDP)
  
! Magnetic field perturbation term  
  
  do l=1,leqmax
     sceq2(:,l)=sceq1(:,l)*dpreqdr
  end do  
  call blockj(sceq2,1,7,1,1,0,0,1.0_IDP)  
  
! Shared equilibrium toroidal flow velocity for the thermal moment of energetic particles equation

  call block0(vzt_eq,7,7,0,0,1,-1.0_IDP)

! diffusion term added
  sd1=1.0_IDP
  if (difnr_on == 0) then
     call block0_dlsq(sd1,7,7,stdifv)
  else
     call block0_dlsqnr(fctr_dif,7,7,stdifvn,dfctr_difdr)
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Alpha particle effects   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (alpha_on == 1) then

! Load 1st Omega-d terms in fast ion density equation:

     do l=1,leqmax
        sceq1(:,l)=valphaova2*omdr(:,l)/(epsq*omcydalp)
        sceq2(:,l)=valphaova2*omdt(:,l)/(epsq*omcydalp)
        sceq3(:,l)=valphaova2*omdz(:,l)/(epsq*omcydalp)
     end do
  
     call blockj(sceq1,-1,8,8,0,1,0,-1.0_IDP)
     call blockj(sceq2,1,8,8,1,0,0,-1.0_IDP)
     call blockj(sceq3,1,8,8,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion density equation

     call block0(vzt_eq,8,8,0,0,1,-1.0_IDP)

! Load 1st Omega-d terms in fast ion parallel velocity equation:

     call blockj(sceq1,-1,9,9,0,1,0,-1.0_IDP)
     call blockj(sceq2,1,9,9,1,0,0,-1.0_IDP)
     call blockj(sceq3,1,9,9,0,0,1,-1.0_IDP)

! Shared equilibrium toroidal flow velocity for fast ion parallel velocity equation

     call block0(vzt_eq,9,9,0,0,1,-1.0_IDP)

! Load 2nd Omega-d terms in fast ion density equation:

     do l=1,leqmax
        sceq1(:,l)=nalpeq(:)*omdr(:,l)
        sceq2(:,l)=nalpeq(:)*omdt(:,l)
        sceq3(:,l)=nalpeq(:)*omdz(:,l)
     end do

     call blockj(sceq1,-1,8,2,0,1,0,-1.0_IDP)
     call blockj(sceq2,1,8,2,1,0,0,-1.0_IDP)
     call blockj(sceq3,1,8,2,0,0,1,-1.0_IDP)

! diffusion term added
     sd1=1.0_IDP
     if (difnr_on == 0) then
        call block0_dlsq(sd1,8,8,stdifnalp)
     else
        call block0_dlsqnr(fctr_dif,8,8,stdifnalpn,dfctr_difdr)
     end if

! Load remaining terms in fast ion density equation

! Parallel gradient term

     do l=1,leqmax
        sceq1(:,l)=nalpeq*bmod(:,l)/(feq-qqinv*cureq)
     end do
  
     call om(sceq1,1,8,9,0,0,0,-1.0_IDP)

! Omega* term

     sd1=dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
     sd2=dnalpeqdr*feq/(feq-qqinv*cureq)
  
     call block0(sd1,8,2,0,0,1,1.0_IDP)
     call block0(sd2,8,2,1,0,0,-1.0_IDP)
   

! EP FLR effects

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then
        sd1=-dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
        sd2=-dnalpeqdr*feq/(feq-qqinv*cureq)
        call block0(sd1,8,iwa,0,0,1,1.0_IDP)
        call block0(sd2,8,iwa,1,0,0,-1.0_IDP)
        ! sd1=epsq*omcydalp*omegar*nalpeq/valphaova2
        ! call block0(sd1,8,iwa,0,0,0,1.0_IDP)
     end if

! Load remaining terms in fast ion parallel velocity equation

! Landau closure term

     do l=1,leqmax
        sceq1(:,l)=1.414213*LcA1alp*valphaova*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call blockj_landau_grad_parallel(sceq1,1,9,9,0,0,0,1.0_IDP)

! Parallel gradient terms

     do l=1,leqmax
        sceq1(:,l)=2*LcA0alp*valphaova2*bmod(:,l)/(nalpeq*(feq-qqinv*cureq))
     end do
     call om(sceq1,1,9,8,0,0,0,-1.0_IDP)

     do l=1,leqmax
        sceq1(:,l)=2*LcA0alp*epsq*omcydalp*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call om(sceq1,1,9,2,0,0,0,-1.0_IDP)

! Omega* term

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then

!  EP FLR effects

        sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
        sd2=valphaova2*dnalpeqdr*rinv*feq/(nalpeq*(feq-qqinv*cureq))
        call block0(sd1,9,ix1a,0,0,0,1.0_IDP)
        call block0(sd2,9,ix2a,0,0,0,-1.0_IDP)

!  EP FLR effects auxiliary equations

        sd1=1.0_IDP
        call block0(sd1,iwa,iwa,0,0,0,1.0_IDP)
        if (epflr_on == 1) then
           sd2=-r_epflralp*r_epflralp
        else if (epflr_on == 2) then
           sd2=-valphaova2/(epsq*omcydalp*omcydalp)
        end if
        call block_dlsq(iwa,iwa,sd2)
        call block_dlsq(iwa,2,-sd2)

        call block0(sd1,ix1a,ix1a,0,0,0,1.0_IDP)
        call block_dlsq(ix1a,ix1a,sd2)
        call block0(sd1,ix1a,1,0,0,1,-1.0_IDP)

        call block0(sd1,ix2a,ix2a,0,0,0,1.0_IDP)
        call block_dlsq(ix2a,ix2a,sd2)
        call block0(r,ix2a,1,1,0,0,-1.0_IDP)

     else
 
        sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
        sd2=valphaova2*dnalpeqdr*feq/(nalpeq*(feq-qqinv*cureq))
        call block0(sd1,9,1,0,0,1,1.0_IDP)
        call block0(sd2,9,1,1,0,0,-1.0_IDP)
   
     end if 

! diffusion term added
     sd1=1.0_IDP
     if (difnr_on == 0) then
        call block0_dlsq(sd1,9,9,stdifvalp)
     else
        call block0_dlsqnr(fctr_dif,9,9,stdifvalpn,dfctr_difdr)
     end if

! End of alpha particles moment equations

  end if

  amat=-dtd2*amat
  bmat=-dtd2*bmat
  cmat=-dtd2*cmat

! put in l.h.s. of equation

! psi equation

  sd1=1.0_IDP
  call block0(sd1,1,1,0,0,0,1.0_IDP)

! u-zeta equation

  call block0(sd1,2,4,0,0,0,1.0_IDP)

! Ion FLR effects

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then
     if (iflr_on == 1) then
        sd2=-iflr*iflr
     else
        sd2=-vthi*vthi*tieq/(epsq*omcyd*omcyd)
     end if
     call block_dlsq(2,4,sd2)
  end if

! p equation

  call block0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density moment equation

  call block0(sd1,5,5,0,0,0,1.0_IDP)

! EP FLR effects

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then
     if (epflr_on == 1) then
        sd2=nfeq/(r_epflr*r_epflr*omcyd)
     else if (epflr_on == 2) then
        sd2=epsq*omcyd*nfeq/vfova2
     end if
     call block0(sd2,5,iw,0,0,0,1.0_IDP)
  end if

! fast ion parallel velocity moment equation

  call block0(sd1,6,6,0,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-epsq*omcyd*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call blockj(sceq1,1,6,1,0,0,0,1.0_IDP)

! thermal moment equation

  call block0(sd1,7,7,0,0,0,1.0_IDP)  

  if (alpha_on == 1) then

! fast ion density moment equation

     call block0(sd1,8,8,0,0,0,1.0_IDP)

! EP FLR effects

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then
        if (epflr_on == 1) then
           sd2=nalpeq/(r_epflralp*r_epflralp*omcydalp)
        else if (epflr_on == 2) then
           sd2=epsq*omcydalp*nalpeq/valphaova2
        end if
        call block0(sd2,8,iwa,0,0,0,1.0_IDP)
     end if

! fast ion parallel velocity moment equation

     call block0(sd1,9,9,0,0,0,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=-epsq*omcydalp*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call blockj(sceq1,1,9,1,0,0,0,1.0_IDP)

  end if

  yt=0.0_IDP
  do i=n_start,n_end
     mnum3=noeqn*mnumn(i)
     do l1t=1,mnum3
        do l2t=1,mnum3
           do j=1,mjm1
              imt1=l1t+mnum3*(j-1)+nskpxn(i)
              imt2=l2t+mnum3*(j-1)+nskpxn(i)
              imat=l1t+mnum3*(l2t-1+mnum3*(j-1))+nskpn(i)
              if (j > 1) yt(imt1)=yt(imt1)+cmat(imat)*xt(imt2-mnum3)
              yt(imt1)=yt(imt1)+amat(imat)*xt(imt2)
              yt(imt1)=yt(imt1)+bmat(imat)*xt(imt2+mnum3)
           end do
        end do
     end do
  end do

! u-dlpersq(phi)=0

  do i=n_start,n_end
     mnum3=noeqn*mnumn(i)
     lskp=3*mnumn(i)
     do l1=1,mnumn(i)
        do j=1,mjm1
           l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
           yt(l1t)=0.
        end do
     end do
  end do

!  Ion FLR effects auxiliary equation

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then
     do i=n_start,n_end
        mnum3=noeqn*mnumn(i)
        lskp=(iq-1)*mnumn(i)
        do l1=1,mnumn(i)
           do j=1,mjm1
              l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
              yt(l1t)=0.
           end do
        end do
     end do
  end if

!  EP FLR effects auxiliary equations

  if (epflr_on > 0) then
     if (r_epflr > 0.0 .or. epflr_on == 2) then
        do i=n_start,n_end
           mnum3=noeqn*mnumn(i)
           lskp=(iw-1)*mnumn(i)
           do l1=1,3*mnumn(i)
              do j=1,mjm1
                 l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
                 yt(l1t)=0.
              end do
           end do
        end do
     end if
     if (alpha_on == 1 .and. (r_epflralp > 0.0 .or. epflr_on == 2)) then
        do i=n_start,n_end
           mnum3=noeqn*mnumn(i)
           lskp=(iwa-1)*mnumn(i)
           do l1=1,3*mnumn(i)
              do j=1,mjm1
                 l1t=l1+lskp+mnum3*(j-1)+nskpxn(i)
                 yt(l1t)=0.
              end do
           end do
        end do
     end if
  end if

  ! do i=n_start,n_end
  !    mnum6=noeqn*noeqn*mnumn(i)*mnumn(i)
  !    do l1=1,mnum6
  !       l1t=l1+nskpn(i)
  !       l2t=l1+mnum6*(mj-2)+nskpn(i)
  !       bmat(l2t)=0.
  !       cmat(l1t)=0.
  !    end do
  ! end do

  do it=n_start,n_end
     loca=nskpn(it)+1
     loci=nskpin(it)+1
     mnum3=noeqn*mnumn(it)
     call decbt(mnum3,mjm1,amat(loca:),bmat(loca:),cmat(loca:),ipc(loci:),ier)
     if (ier /= 0) then
        write(6,'("ier=",i5," ntor=",i5)') ier,it
        stop 19
     end if
  end do

  ! call cpu_time(time_em)
  ! write(6,'(/"Time spent in writing matrices:",1pe13.6," seconds")') time_em-time_sm
  ! time_c=0.0
  ! time_m=0.0

end subroutine linstart
