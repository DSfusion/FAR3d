subroutine solve

  use param
  use var_para
  use cotrol
  use domain
  use equil
  use dynamo
  use dbyd
  use block_mod
  use mult_mod
  use matrix
  use scratch

  implicit none

  integer :: l,i,mnum3,l1,l1t,it,loca,loci,locx,lskp,mnum2
  real(IDP) :: epsq,omcyd,omcydalp,scmn,scmx,uztmx,uztedg

  interface
     subroutine cnvt(idir)
        implicit none
        integer :: idir
     end subroutine cnvt
     subroutine delstar(ss,ff,itypf,wk1,wk2,wkeq1,c1,c2)
        use param
        use var_para
        implicit none
        integer :: itypf
        real(IDP) :: c1,c2
        real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2
        real(IDP), dimension(0:,0:) :: wkeq1
     end subroutine delstar
     subroutine bigf(f,g,itypeg,h,itypeh,sx1,sx2,c1,c2)
        use param
        use var_para
        implicit none
        integer :: itypeg,itypeh
        real(IDP) :: c1,c2
        real(IDP), dimension(mj_start:,0:) :: f,g,h,sx1,sx2
     end subroutine bigf
     subroutine divv(f,c1,c2)
        use param
        use var_para
        implicit none
        real(IDP) :: c1,c2
        real(IDP), dimension(mj_start:,0:) :: f
     end subroutine divv
  end interface

  epsq=eps*eps
  omcyd = omcy
  omcydalp = omcyalp

!  first half-step

!  From previous step we have yt=(L-R*Dt/2)Y(t), and we want to build yt=(L+R*Dt/2)Y(t)
!  So we change sign and add 2*L*Y(t).

  yt=-yt

! put in l.h.s. of equation

! psi equation

  sd1=2.0_IDP
  call b2lx0(sd1,1,1,0,0,0,1.0_IDP)

! u-zeta equation

  call b2lx0(sd1,2,4,0,0,0,1.0_IDP)

! Ion FLR effects

  if ((iflr_on == 1 .and. iflr > 0.0) .or. iflr_on == 2) then
     if (iflr_on == 1) then
        sd2=-2.*iflr*iflr
     else
        sd2=-2.*vthi*vthi*tieq/(epsq*omcyd*omcyd)
     end if
     call b2lx_dlsq(2,4,sd2)
  end if

! p equation

  call b2lx0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density equation

  call b2lx0(sd1,5,5,0,0,0,1.0_IDP)

! EP FLR effects

  if ((epflr_on == 1 .and. r_epflr > 0.0) .or. epflr_on == 2) then
     if (epflr_on == 1) then
        sd2=2.0*nfeq/(r_epflr*r_epflr*omcyd)
     else if (epflr_on == 2) then
        sd2=2.0*epsq*omcyd*nfeq/vfova2
     end if
     call b2lx0(sd2,5,iw,0,0,0,1.0_IDP)
  end if

! fast ion v-parallel equation

  call b2lx0(sd1,6,6,0,0,0,1.0_IDP)
  do l=1,leqmax
     sceq1(:,l)=-2.0*epsq*omcyd*bmod(:,l)/(feq-qqinv*cureq)
  end do
  call b2lx(sceq1,1,6,1,0,0,0,1.0_IDP)
  
! thermal moment equation

  call b2lx0(sd1,7,7,0,0,0,1.0_IDP)

  if (alpha_on == 1) then

! alpha density equation

     call b2lx0(sd1,8,8,0,0,0,1.0_IDP)

! EP FLR effects

     if ((epflr_on == 1 .and. r_epflralp > 0.0) .or. epflr_on == 2) then
        if (epflr_on == 1) then
           sd2=2.0*nalpeq/(r_epflralp*r_epflralp*omcydalp)
        else if (epflr_on == 2) then
           sd2=2.0*epsq*omcydalp*nalpeq/valphaova2
        end if
        call b2lx0(sd2,8,iwa,0,0,0,1.0_IDP)
     end if

! alpha v-parallel equation

     call b2lx0(sd1,9,9,0,0,0,1.0_IDP)
     do l=1,leqmax
        sceq1(:,l)=-2.0*epsq*omcydalp*bmod(:,l)/(feq-qqinv*cureq)
     end do
     call b2lx(sceq1,1,9,1,0,0,0,1.0_IDP)

  end if

  if (nonlin /= 0) then
  
! transfer (L+R*Dt/2)Y(t) to secondary arrays

     call cnvt(3)

! store values of secondary arrays

     sc5=uztp
     sc6=prp
     sc7=psip
     sc8=nfpp
     sc9=vprlfp
     sc10=vthprlfp
     if (alpha_on == 1) then
        sc11=nalpp
        sc12=vprlalpp
     end if
     
! psi equation

     call bigf(psi_nl(:,:,1),phi,-1,psi,1,sc3,sc4,0.0_IDP,dt)

! u-zeta equation

     call bigf(uzt_nl(:,:,1),phi,-1,uzt,-1,sc3,sc4,0.0_IDP,dt)
     call delstar(sc1,psi,1,sc2,sc3,sceq1,0.0_IDP,1.0_IDP)
     call bigf(uzt_nl(:,:,1),sc1,1,psi,1,sc2,sc3,1.0_IDP,dt)

! pr equation

     call bigf(pr_nl(:,:,1),phi,-1,pr,1,sc3,sc4,0.0_IDP,dt)
     if (gamma > 0.0_IDP) call divv(pr_nl(:,:,1),1.0_IDP,gamma*dt)

     call bigf(nf_nl(:,:,1),phi,-1,nf,1,sc3,sc4,0.0_IDP,dt)
 
     call bigf(vprlf_nl(:,:,1),phi,-1,vprlf,-1,sc3,sc4,0.0_IDP,dt)
 
     call bigf(vthprlf_nl(:,:,1),phi,-1,vthprlf,-1,sc3,sc4,0.0_IDP,dt)
     call bigf(sc1,psi,1,pr,1,sc2,sc3,0.0_IDP,dt)
     do l=1,leqmax
        sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
     end do
     call multed(vthprlf_nl(:,:,1),sceq1,1,sc1,-1,1.0_IDP,1.0_IDP)

     if (alpha_on == 1) then

        call bigf(nalp_nl(:,:,1),phi,-1,nalp,1,sc3,sc4,0.0_IDP,dt)
 
        call bigf(vprlalp_nl(:,:,1),phi,-1,vprlalp,-1,sc3,sc4,0.0_IDP,dt)
     end if

! parallel velocity terms
 
     do l=1,leqmax
        sceq1(:,l)=bmod(:,l)/(feq-qqinv*cureq)
     end do

! u-zeta equation

     call grdpar(sc1,uzt,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(uzt_nl(:,:,1),sc3,1,vthprlf,-1,1.0_IDP,-dt)

! pr equation

     ! call dbydth_par(sc2,psi,1,0.0_IDP,-1.0_IDP,0)
     ! do l=1,lmax
     !    sc1(:,l)=dpreqdr*sc2(:,l)
     ! end do
     ! call grdpar(sc1,pr,1,1.0_IDP,1.0_IDP)
     call grdpar(sc1,pr,1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
     call mult(pr_nl(:,:,1),sc3,-1,vthprlf,-1,1.0_IDP,-dt)

     ! do l=1,lmax
     !    sc1(:,l)=dnfeqdr*sc2(:,l)
     ! end do
     ! call grdpar(sc1,nf,1,1.0_IDP,1.0_IDP)
     call grdpar(sc1,nf,1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
     call mult(nf_nl(:,:,1),sc3,-1,vprlf,-1,1.0_IDP,-dt)

     call grdpar(sc1,vprlf,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(vprlf_nl(:,:,1),sc3,1,vprlf,-1,1.0_IDP,-dt)

     call grdpar(sc1,vthprlf,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(vthprlf_nl(:,:,1),sc3,1,vthprlf,-1,1.0_IDP,-dt)

     if (alpha_on == 1) then
        ! do l=1,lmax
        !    sc1(:,l)=dnalpeqdr*sc2(:,l)
        ! end do
        ! call grdpar(sc1,nalp,1,1.0_IDP,1.0_IDP)
        call grdpar(sc1,nalp,1,0.0_IDP,1.0_IDP)
        call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
        call mult(nalp_nl(:,:,1),sc3,-1,vprlalp,-1,1.0_IDP,-dt)

        call grdpar(sc1,vprlalp,-1,0.0_IDP,1.0_IDP)
        call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
        call mult(vprlalp_nl(:,:,1),sc3,1,vprlalp,-1,1.0_IDP,-dt)
     end if

! At this time, we have varp=(L+Dt*R/2)Y(t)+(Dt/12)*(23*NL[Y(t)]-16*NL[Y(t-Dt)]+5*NL[Y(t-2*Dt)])
 
     if (nopsievol_on == 1) psi_nl(:,l0,1)=0.0
     if (noprevol_on == 1) pr_nl(:,l0,1)=0.0
     if (nonfevol_on == 1) nf_nl(:,l0,1)=0.0
     if (alpha_on == 1 .and. nonalpevol_on == 1) nalp_nl(:,l0,1)=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source and sinks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (src_sink_th_on == 1) pr_nl(:,l0,1) = pr_nl(:,l0,1)+dt*src(mj_start:mj_end)
     if (src_sink_EP1_on == 1) nf_nl(:,l0,1) = nf_nl(:,l0,1)+dt*src_nf(mj_start:mj_end)
     if (alpha_on == 1 .and. src_sink_EP2_on == 1) nalp_nl(:,l0,1) = nalp_nl(:,l0,1)+dt*src_nalpha(mj_start:mj_end)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (nstep <= nstep1+2) then
        uztp=uztp+uzt_nl(:,:,1)
        prp=prp+pr_nl(:,:,1)
        psip=psip+psi_nl(:,:,1)
        nfpp=nfpp+nf_nl(:,:,1)
        vprlfp=vprlfp+vprlf_nl(:,:,1)
        vthprlfp=vthprlfp+vthprlf_nl(:,:,1)
        if (alpha_on == 1) then
           nalpp=nalpp+nalp_nl(:,:,1)
           vprlalpp=vprlalpp+vprlalp_nl(:,:,1)
        end if
     else
        uztp=uztp+(23.*uzt_nl(:,:,1)-16.*uzt_nl(:,:,2)+5.*uzt_nl(:,:,3))/12.
        prp=prp+(23.*pr_nl(:,:,1)-16.*pr_nl(:,:,2)+5.*pr_nl(:,:,3))/12.
        psip=psip+(23.*psi_nl(:,:,1)-16.*psi_nl(:,:,2)+5.*psi_nl(:,:,3))/12.
        nfpp=nfpp+(23.*nf_nl(:,:,1)-16.*nf_nl(:,:,2)+5.*nf_nl(:,:,3))/12.
        vprlfp=vprlfp+(23.*vprlf_nl(:,:,1)-16.*vprlf_nl(:,:,2)+5.*vprlf_nl(:,:,3))/12.
        vthprlfp=vthprlfp+(23.*vthprlf_nl(:,:,1)-16.*vthprlf_nl(:,:,2)+5.*vthprlf_nl(:,:,3))/12.
        if (alpha_on == 1) then
           nalpp=nalpp+(23.*nalp_nl(:,:,1)-16.*nalp_nl(:,:,2)+5.*nalp_nl(:,:,3))/12.
           vprlalpp=vprlalpp+(23.*vprlalp_nl(:,:,1)-16.*vprlalp_nl(:,:,2)+5.*vprlalp_nl(:,:,3))/12.
        end if
     end if 

     call cnvt(1)

  end if
  
! find predicted values

  xt=yt

  do i=n_start,n_end
     mnum3=noeqn*mnumn(i)
     do l1=1,mnum3
        l1t=l1+mnum3*mjm1+nskpxn(i)
        xt(l1t)=0.
     end do
  end do

  ! call cpu_time(time_sm)

  do it=n_start,n_end
     loca=nskpn(it)+1
     loci=nskpin(it)+1
     locx=nskpxn(it)+1
     mnum3=noeqn*mnumn(it)
     call solbt(mnum3,mjm1,amat(loca:),bmat(loca:),cmat(loca:),xt(locx:),ipc(loci:))
  end do

  ! call cpu_time(time_em)
  ! time_m=time_m+time_em-time_sm

! second half-step


  if (nonlin /= 0) then

! recover (L+R*Dt/2)Y(t) to secondary arrays

     uztp=sc5
     prp=sc6
     psip=sc7
     nfpp=sc8
     vprlfp=sc9
     vthprlfp=sc10
     if (alpha_on == 1) then
        nalpp=sc11
        vprlalpp=sc12
     end if

! transfer predicted values to original arrays

     call cnvt(2)

! psi equation

     call bigf(psi_nl(:,:,3),phi,-1,psi,1,sc1,sc2,0.0_IDP,dt)

! u-zeta equation

     call bigf(uzt_nl(:,:,3),phi,-1,uzt,-1,sc1,sc2,0.0_IDP,dt)
     call delstar(sc1,psi,1,sc2,sc3,sceq1,0.0_IDP,1.0_IDP)
     call bigf(uzt_nl(:,:,3),sc1,1,psi,1,sc2,sc3,1.0_IDP,dt)

! pr equation

     call bigf(pr_nl(:,:,3),phi,-1,pr,1,sc1,sc2,0.0_IDP,dt)
     if (gamma > 0.0_IDP) call divv(pr_nl(:,:,3),1.0_IDP,gamma*dt)

     call bigf(nf_nl(:,:,3),phi,-1,nf,1,sc1,sc2,0.0_IDP,dt)
 
     call bigf(vprlf_nl(:,:,3),phi,-1,vprlf,-1,sc1,sc2,0.0_IDP,dt)
 
     call bigf(vthprlf_nl(:,:,3),phi,-1,vthprlf,-1,sc1,sc2,0.0_IDP,dt)
     call bigf(sc1,psi,1,pr,1,sc2,sc3,0.0_IDP,dt)
     do l=1,leqmax
        sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
     end do
     call multed(vthprlf_nl(:,:,3),sceq1,1,sc1,-1,1.0_IDP,1.0_IDP)

     if (alpha_on == 1) then

        call bigf(nalp_nl(:,:,3),phi,-1,nalp,1,sc3,sc4,0.0_IDP,dt)
 
        call bigf(vprlalp_nl(:,:,3),phi,-1,vprlalp,-1,sc3,sc4,0.0_IDP,dt)
     end if

! parallel velocity terms
 
     do l=1,leqmax
        sceq1(:,l)=bmod(:,l)/(feq-qqinv*cureq)
     end do

! u-zeta equation

     call grdpar(sc1,uzt,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(uzt_nl(:,:,3),sc3,1,vthprlf,-1,1.0_IDP,-dt)

! pr equation

     ! call dbydth_par(sc2,psi,1,0.0_IDP,-1.0_IDP,0)
     ! do l=1,lmax
     !    sc1(:,l)=dpreqdr*sc2(:,l)
     ! end do
     ! call grdpar(sc1,pr,1,1.0_IDP,1.0_IDP)
     call grdpar(sc1,pr,1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
     call mult(pr_nl(:,:,3),sc3,-1,vthprlf,-1,1.0_IDP,-dt)

     ! do l=1,lmax
     !    sc1(:,l)=dnfeqdr*sc2(:,l)
     ! end do
     ! call grdpar(sc1,nf,1,1.0_IDP,1.0_IDP)
     call grdpar(sc1,nf,1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
     call mult(nf_nl(:,:,3),sc3,-1,vprlf,-1,1.0_IDP,-dt)

     call grdpar(sc1,vprlf,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(vprlf_nl(:,:,3),sc3,1,vprlf,-1,1.0_IDP,-dt)

     call grdpar(sc1,vthprlf,-1,0.0_IDP,1.0_IDP)
     call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
     call mult(vthprlf_nl(:,:,3),sc3,1,vthprlf,-1,1.0_IDP,-dt)

     if (alpha_on == 1) then
        ! do l=1,lmax
        !    sc1(:,l)=dnalpeqdr*sc2(:,l)
        ! end do
        ! call grdpar(sc1,nalp,1,1.0_IDP,1.0_IDP)
        call grdpar(sc1,nalp,1,0.0_IDP,1.0_IDP)
        call multed(sc3,sceq1,1,sc1,-1,0.0_IDP,1.0_IDP)
        call mult(nalp_nl(:,:,3),sc3,-1,vprlalp,-1,1.0_IDP,-dt)

        call grdpar(sc1,vprlalp,-1,0.0_IDP,1.0_IDP)
        call multed(sc3,sceq1,1,sc1,1,0.0_IDP,1.0_IDP)
        call mult(vprlalp_nl(:,:,3),sc3,1,vprlalp,-1,1.0_IDP,-dt)
     end if

     if (nopsievol_on == 1) psi_nl(:,l0,3)=0.0
     if (noprevol_on == 1) pr_nl(:,l0,3)=0.0
     if (nonfevol_on == 1) nf_nl(:,l0,3)=0.0
     if (alpha_on == 1 .and. nonalpevol_on == 1) nalp_nl(:,l0,3)=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Source and sinks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (src_sink_th_on == 1) pr_nl(:,l0,3) = pr_nl(:,l0,3)+dt*src(mj_start:mj_end)
     if (src_sink_EP1_on == 1) nf_nl(:,l0,3) = nf_nl(:,l0,3)+dt*src_nf(mj_start:mj_end)
     if (alpha_on == 1 .and. src_sink_EP2_on == 1) nalp_nl(:,l0,3) = nalp_nl(:,l0,3)+dt*src_nalpha(mj_start:mj_end)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (nstep <= nstep1+2) then
        uztp=uztp+0.5*(uzt_nl(:,:,1)+uzt_nl(:,:,3))
        prp=prp+0.5*(pr_nl(:,:,1)+pr_nl(:,:,3))
        psip=psip+0.5*(psi_nl(:,:,1)+psi_nl(:,:,3))
        nfpp=nfpp+0.5*(nf_nl(:,:,1)+nf_nl(:,:,3))
        vprlfp=vprlfp+0.5*(vprlf_nl(:,:,1)+vprlf_nl(:,:,3))
        vthprlfp=vthprlfp+0.5*(vthprlf_nl(:,:,1)+vthprlf_nl(:,:,3))
        if (alpha_on == 1) then
           nalpp=nalpp+0.5*(nalp_nl(:,:,1)+nalp_nl(:,:,3))
           vprlalpp=vprlalpp+0.5*(vprlalp_nl(:,:,1)+vprlalp_nl(:,:,3))
        end if
     else
        uztp=uztp+(8.*uzt_nl(:,:,1)-uzt_nl(:,:,2)+5.*uzt_nl(:,:,3))/12.
        prp=prp+(8.*pr_nl(:,:,1)-pr_nl(:,:,2)+5.*pr_nl(:,:,3))/12.
        psip=psip+(8.*psi_nl(:,:,1)-psi_nl(:,:,2)+5.*psi_nl(:,:,3))/12.
        nfpp=nfpp+(8.*nf_nl(:,:,1)-nf_nl(:,:,2)+5.*nf_nl(:,:,3))/12.
        vprlfp=vprlfp+(8.*vprlf_nl(:,:,1)-vprlf_nl(:,:,2)+5.*vprlf_nl(:,:,3))/12.
        vthprlfp=vthprlfp+(8.*vthprlf_nl(:,:,1)-vthprlf_nl(:,:,2)+5.*vthprlf_nl(:,:,3))/12.
        if (alpha_on == 1) then
           nalpp=nalpp+(8.*nalp_nl(:,:,1)-nalp_nl(:,:,2)+5.*nalp_nl(:,:,3))/12.
           vprlalpp=vprlalpp+(8.*vprlalp_nl(:,:,1)-vprlalp_nl(:,:,2)+5.*vprlalp_nl(:,:,3))/12.
        end if
     end if

! At this time, we have varp=(L+Dt*R/2)Y(t)+(Dt/12)*(5*NL[Y(p)]+8*NL[Y(t)]-NL[Y(t-Dt)])

     call cnvt(1)

! find time advanced values

     xt=yt
     do i=n_start,n_end
        mnum3=noeqn*mnumn(i)
        do l1=1,mnum3
           l1t=l1+mnum3*mjm1+nskpxn(i)
           xt(l1t)=0.
        end do
     end do

     ! call cpu_time(time_sm)

     do it=n_start,n_end
        loca=nskpn(it)+1
        loci=nskpin(it)+1
        locx=nskpxn(it)+1
        mnum3=noeqn*mnumn(it)
        call solbt(mnum3,mjm1,amat(loca:),bmat(loca:),cmat(loca:),xt(locx:),ipc(loci:))
     end do

     ! call cpu_time(time_em)
     ! time_m=time_m+time_em-time_sm

! save current values for the next step

     uzt_nl(:,:,3)=uzt_nl(:,:,2)
     pr_nl(:,:,3)=pr_nl(:,:,2)
     psi_nl(:,:,3)=psi_nl(:,:,2)
     nf_nl(:,:,3)=nf_nl(:,:,2)
     vprlf_nl(:,:,3)=vprlf_nl(:,:,2)
     vthprlf_nl(:,:,3)=vthprlf_nl(:,:,2)
     uzt_nl(:,:,2)=uzt_nl(:,:,1)
     pr_nl(:,:,2)=pr_nl(:,:,1)
     psi_nl(:,:,2)=psi_nl(:,:,1)
     nf_nl(:,:,2)=nf_nl(:,:,1)
     vprlf_nl(:,:,2)=vprlf_nl(:,:,1)
     vthprlf_nl(:,:,2)=vthprlf_nl(:,:,1)
     if (alpha_on == 1) then
        nalp_nl(:,:,3)=nalp_nl(:,:,2)
        vprlalp_nl(:,:,3)=vprlalp_nl(:,:,2)
        nalp_nl(:,:,2)=nalp_nl(:,:,1)
        vprlalp_nl(:,:,2)=vprlalp_nl(:,:,1)
     end if

  end if

! transfer time advanced values to original arrays

  call cnvt(2)

end subroutine solve
