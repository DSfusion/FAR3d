MODULE diagnostics

  USE param
  USE processor
  USE var_para
  USE domain
  USE cotrol
  USE equil
  USE dynamo
  USE dbyd
  USE mult_mod
  USE transfer
  USE tools
  USE scratch

  IMPLICIT NONE

CONTAINS

  subroutine energy(i)

    implicit none

    integer, save :: iflg = 0
    integer :: i,j,l,lp,l1,l2,mjp
    real(IDP) :: denom,scnorm,emetot,eketot,eprtot,ealptot,emetot0,eketot0,eprtot0,ealptot0
    character(len=1) :: t
    character(len=16) :: confil
    character(len=32) :: format1='("time",1200(a1,i4,"/",i4))'
    character(len=60) :: formatt='("time",a1,"total",a1,"total (w/o 0/0)",1200(a1,i4,"/",i4))'
    character(len=32) :: formatv='(1pe13.6,1202(a1,1pe15.8))'
    integer, dimension(lmax) :: lmap
    real(IDP), dimension(lmax) :: gampsi,gamalpha,gamphi,gampr
    real(IDP), dimension(mj) :: rint,xint,bb,cc,dd
    real(IDP), dimension(0:mj) :: xinte,bbe,cce,dde

    epsi(:,i)=0.0_IDP 
    ephi(:,i)=0.0_IDP 
    epr(:,i)=0.0_IDP 
    eprnc(:,i)=0.0_IDP 
    ekenc(:,i)=0.0_IDP 
    eke(:,i)=0.0_IDP 
    emenc(:,i)=0.0_IDP 
    eme(:,i)=0.0_IDP 
    ealp(:,i)=0.0_IDP 
    ealpnc(:,i)=0.0_IDP 

    !  store uzt-values into sc4

    do l=1,lmax
       sc4(:,l)=uzt(mj_start:mj_end,l)
    end do

    !  vr up
    call dbydth_par(sc2,phi,-1,0.0_IDP,-1.0_IDP,0)
    !  vth up
    call dbydr_par(sc3,phi,0.0_IDP,1.0_IDP,0)

    mjp=mj+1
    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=jbgrr(mj_start:mj_end,leq0)*sc2(:,lln(l))**2+jbgtt(mj_start:mj_end,leq0)*sc3(:,lln(l))**2
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*denseq*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ekenc(l,i))
          if (l == l0) ekenc(l,i)=2.*ekenc(l,i)
       end do
    end if

    call multed(sc5,jbgrr,1,sc2,1,0.0_IDP,1.0_IDP)
    call multed(sc5,jbgrt,-1,sc3,-1,1.0_IDP,1.0_IDP)
    call multed(sc6,jbgrt,-1,sc2,1,0.0_IDP,1.0_IDP)
    call multed(sc6,jbgtt,1,sc3,-1,1.0_IDP,1.0_IDP)
    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=sc2(:,lln(l))*sc5(:,lln(l))+sc3(:,lln(l))*sc6(:,lln(l))
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*denseq*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eke(l,i))
          if (l == l0) eke(l,i)=2.*eke(l,i)
       end do
    end if

    !  br up
    call dbydth_par(sc2,psi,1,0.0_IDP,-1.0_IDP,0)
    !  bth up
    call dbydr_par(sc3,psi,0.0_IDP,1.0_IDP,0)

    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=grroj(mj_start:mj_end,leq0)*sc2(:,lln(l))**2+gttoj(mj_start:mj_end,leq0)*sc3(:,lln(l))**2
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,emenc(l,i))
          if (l == l0) emenc(l,i)=2.*emenc(l,i)
       end do
    end if

    call multed(sc5,grroj,1,sc2,-1,0.0_IDP,1.0_IDP)
    call multed(sc5,grtoj,-1,sc3,1,1.0_IDP,1.0_IDP)
    call multed(sc6,grtoj,-1,sc2,-1,0.0_IDP,1.0_IDP)
    call multed(sc6,gttoj,1,sc3,1,1.0_IDP,1.0_IDP)
    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=sc2(:,lln(l))*sc5(:,lln(l))+sc3(:,lln(l))*sc6(:,lln(l))
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eme(l,i))
          if (l == l0) eme(l,i)=2.*eme(l,i)
       end do
    end if

    call multed(sc2,sqg,1,vprlf,-1,0.0_IDP,1.0_IDP)
    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=nfeq(mj_start:mj_end)*sqg(mj_start:mj_end,leq0)*vprlf(mj_start:mj_end,lln(l))**2
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eprnc(l,i))
          if (l == l0) eprnc(l,i)=2.*eprnc(l,i)
       end do
    end if
    do l=1,lmaxn
       uzt(mj_start:mj_end,l)=nfeq(mj_start:mj_end)*sc2(:,lln(l))*vprlf(mj_start:mj_end,lln(l))
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do l=1,lmaxn
          xinte=r*uzt(:,l)
          call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,epr(l,i))
          if (l == l0) epr(l,i)=2.*epr(l,i)
       end do
    end if

    if (alpha_on == 1) then
       call multed(sc2,sqg,1,vprlalp,-1,0.0_IDP,1.0_IDP)
       do l=1,lmaxn
          uzt(mj_start:mj_end,l)=nalpeq(mj_start:mj_end)*sqg(mj_start:mj_end,leq0)*vprlalp(mj_start:mj_end,lln(l))**2
       end do
       call trnsfr0(uzt,1)
       if (myPE == 0) then
          do l=1,lmaxn
             xinte=r*uzt(:,l)
             call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealpnc(l,i))
             if (l == l0) ealpnc(l,i)=2.*ealpnc(l,i)
          end do
       end if
       do l=1,lmaxn
          uzt(mj_start:mj_end,l)=nalpeq(mj_start:mj_end)*sc2(:,lln(l))*vprlalp(mj_start:mj_end,lln(l))
       end do
       call trnsfr0(uzt,1)
       if (myPE == 0) then
          do l=1,lmaxn
             xinte=r*uzt(:,l)
             call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealp(l,i))
             if (l == l0) ealp(l,i)=2.*ealp(l,i)
          end do
       end if
    end if

    !  restore uzt-values

    do l=1,lmax
       uzt(mj_start:mj_end,l)=sc4(:,l)
    end do

    if (i == 1) return

    if (myPE == 0) then

       write(6,'(/"energy:numrun=",2a2,a1,",numruno=",2a2,a1,",nstep=",i10,",time=",1pe12.5,",dt=",1pe12.5)') numrun,numruno, &
            nstep,time,dt

       lp=0
       emetot=0.0
       eketot=0.0
       eprtot=0.0
       ealptot=0.0
       emetot0=0.0
       eketot0=0.0
       eprtot0=0.0
       ealptot0=0.0
       do l=1,lmaxn
          l1=lln(l)
          if (signl(l1) < 0) cycle
          if (signl(l1) > 0 .and. lo(l) > 0) then
             l2=lo(l)
             ekenc(l,1)=ekenc(l,1)+ekenc(l2,1)
             eke(l,1)=eke(l,1)+eke(l2,1)
             emenc(l,1)=emenc(l,1)+emenc(l2,1)
             eme(l,1)=eme(l,1)+eme(l2,1)
             eprnc(l,1)=eprnc(l,1)+eprnc(l2,1)
             epr(l,1)=epr(l,1)+epr(l2,1)
             ealpnc(l,1)=ealpnc(l,1)+ealpnc(l2,1)
             ealp(l,1)=ealp(l,1)+ealp(l2,1)
             ekenc(l,2)=ekenc(l,2)+ekenc(l2,2)
             eke(l,2)=eke(l,2)+eke(l2,2)
             emenc(l,2)=emenc(l,2)+emenc(l2,2)
             eme(l,2)=eme(l,2)+eme(l2,2)
             eprnc(l,2)=eprnc(l,2)+eprnc(l2,2)
             epr(l,2)=epr(l,2)+epr(l2,2)
             ealpnc(l,2)=ealpnc(l,2)+ealpnc(l2,2)
             ealp(l,2)=ealp(l,2)+ealp(l2,2)
          end if
          gampsi(l)=0.0_IDP
          denom=dt*(emenc(l,1)+emenc(l,2))
          if (denom /= 0.0_IDP) gampsi(l)=(emenc(l,2)-emenc(l,1))/denom
          gamphi(l)=0.0_IDP
          denom=dt*(ekenc(l,1)+ekenc(l,2))
          if (denom /= 0.0_IDP) gamphi(l)=(ekenc(l,2)-ekenc(l,1))/denom
          gampr(l)=0.0_IDP
          denom=dt*(eprnc(l,1)+eprnc(l,2))
          if (denom /= 0.0_IDP) gampr(l)=(eprnc(l,2)-eprnc(l,1))/denom
          emetot=emetot+eme(l,2)
          eketot=eketot+eke(l,2)
          eprtot=eprtot+epr(l,2)
          ealptot=ealptot+ealp(l,2)
          if (lln(l) /= l0) then
             emetot0=emetot0+eme(l,2)
             eketot0=eketot0+eke(l,2)
             eprtot0=eprtot0+epr(l,2)
             ealptot0=ealptot0+ealp(l,2)
          end if
          lp=lp+1
          lmap(lp)=l
       end do
       t=char(9)
       if (iflg == 0) then
          write(confil,'("eme_",2a2)') numrun(1),numrun(2)
          open(unit=83,file=confil,recl=19384)
          if (nstres == 0) write(83,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          write(confil,'("eke_",2a2)') numrun(1),numrun(2)
          open(unit=84,file=confil,recl=19384)
          if (nstres == 0) write(84,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          write(confil,'("evprlf_",2a2)') numrun(1),numrun(2)
          open(unit=85,file=confil,recl=19384)
          if (nstres == 0) write(85,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          write(confil,'("emenc_",2a2)') numrun(1),numrun(2)
          open(unit=93,file=confil,recl=19384)
          if (nstres == 0) write(93,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          write(confil,'("ekenc_",2a2)') numrun(1),numrun(2)
          open(unit=94,file=confil,recl=19384)
          if (nstres == 0) write(94,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          write(confil,'("evprlfnc_",2a2)') numrun(1),numrun(2)
          open(unit=95,file=confil,recl=19384)
          if (nstres == 0) write(95,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
       end if
       write(83,formatv) time,t,eps*eps*emetot,t,eps*eps*emetot0,(t,eps*eps*eme(lmap(l),2),l=1,lp) 
       write(84,formatv) time,t,eps*eps*eketot,t,eps*eps*eketot0,(t,eps*eps*eke(lmap(l),2),l=1,lp) 
       write(85,formatv) time,t,bet0_f*eprtot/2.0,t,bet0_f*eprtot0/2.0,(t,bet0_f*epr(lmap(l),2)/2.0,l=1,lp) 
       write(93,formatv) time,(t,eps*eps*emenc(lmap(l),2),l=1,lp) 
       write(94,formatv) time,(t,eps*eps*ekenc(lmap(l),2),l=1,lp) 
       write(95,formatv) time,(t,bet0_f*eprnc(lmap(l),2)/2.0,l=1,lp) 

       if (alpha_on == 1) then
          do l=1,lmaxn
             gamalpha(l)=0.0_IDP
             denom=dt*(ealpnc(l,1)+ealpnc(l,2))
             if (denom /= 0.0_IDP) gamalpha(l)=(ealpnc(l,2)-ealpnc(l,1))/denom
          end do
          write(6,'(/"   l   m/  n        ke        me     vprlf   vprlalp     gamke     gamme    gamvpr    gamalp")')
          do l=1,lmaxn
             l1=lln(l)
             if (signl(l1) < 0) cycle
             write(6,'(3i4,1p8e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l,2),eps*eps*emenc(l,2),bet0_f*eprnc(l,2)/2.0, &
                  bet0_alp*ealpnc(l,2)/2.0,gamphi(l),gampsi(l),gampr(l),gamalpha(l)
          end do
          if (iflg == 0) then
             write(confil,'("evalp_",2a2)') numrun(1),numrun(2)
             open(unit=86,file=confil,recl=19384)
             if (nstres == 0) write(86,formatt) t,t,(t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
             write(confil,'("evalpnc_",2a2)') numrun(1),numrun(2)
             open(unit=96,file=confil,recl=19384)
             if (nstres == 0) write(96,format1) (t,mm(lln(lmap(l))),nn(lln(lmap(l))),l=1,lp)
          end if
          write(86,formatv) time,t,bet0_alp*ealptot/2.0,t,bet0_alp*ealptot0/2.0,(t,bet0_alp*ealp(lmap(l),2)/2.0,l=1,lp) 
          write(96,formatv) time,(t,bet0_alp*ealpnc(lmap(l),2)/2.0,l=1,lp) 
       else
          write(6,'(/"   l   m/  n        ke        me     vprlf     gamke     gamme    gamvpr")')
          do l=1,lmaxn
             l1=lln(l)
             if (signl(l1) < 0) cycle
             write(6,'(3i4,1p6e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l,2),eps*eps*emenc(l,2),bet0_f*eprnc(l,2)/2.0, &
                  gamphi(l),gampsi(l),gampr(l)
          end do
       end if

       iflg=1

       if (nonlin == 0 .and. nstep == nstep1+maxstp) then
          scnorm=max(maxval(ekenc),maxval(eke),maxval(emenc),maxval(eme))
          if (scnorm /= 0.0_IDP) then
             ekenc(1:lmaxn,2)=ekenc(1:lmaxn,2)/scnorm
             eke(1:lmaxn,2)=eke(1:lmaxn,2)/scnorm
             emenc(1:lmaxn,2)=emenc(1:lmaxn,2)/scnorm
             eme(1:lmaxn,2)=eme(1:lmaxn,2)/scnorm
             eprnc(1:lmaxn,2)=bet0_f*eprnc(1:lmaxn,2)/(2.0*eps*eps*scnorm)
             epr(1:lmaxn,2)=bet0_f*epr(1:lmaxn,2)/(2.0*eps*eps*scnorm)
             ealpnc(1:lmaxn,2)=bet0_alp*ealpnc(1:lmaxn,2)/(2.0*eps*eps*scnorm)
             ealp(1:lmaxn,2)=bet0_alp*ealp(1:lmaxn,2)/(2.0*eps*eps*scnorm)
          end if
          write(confil,'("spctr_",2a2)') numrun(1),numrun(2)
          open(unit=92,file=confil)
          if (alpha_on == 1) then
             write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf", &
                  a1,"vprlalpnc",a1,"vprlalp")') (t,l=1,10)
             do l=1,lmaxn
                if (signl(lln(l)) < 0) cycle
                write(92,'(i4,2(a1,i4),8(a1,1pe13.6))') lln(l),t,mm(lln(l)),t,nn(lln(l)),t,ekenc(l,2), &
                     t,eke(l,2),t,emenc(l,2),t,eme(l,2),t,eprnc(l,2),t,epr(l,2),t,ealpnc(l,2),t,ealp(l,2)
             end do
          else
             write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf")') (t,l=1,8)
             do l=1,lmaxn
                if (signl(lln(l)) < 0) cycle
                write(92,'(i4,2(a1,i4),6(a1,1pe13.6))') lln(l),t,mm(lln(l)),t,nn(lln(l)),t,ekenc(l,2), &
                     t,eke(l,2),t,emenc(l,2),t,eme(l,2),t,eprnc(l,2),t,epr(l,2)
             end do
          end if
          close(92)
       end if

    end if

  end subroutine energy

!
!     Routine for high frequency (every time step) output of a few quantities
!
  subroutine hifreq(last)

    implicit none

    real(IDP) :: bthta, btht1, btht2, btht3, btht4, phi00, &
         psi00,pr00,vprl00,epsq,omcyd
    real(IDP) :: flux_nf1,flux_nf2,flux_nf3,flux_nf4,flux_nf5,   &
         flux_nalf1,flux_nalf2,flux_nalf3,flux_nalf4,    &
         flux_nalf5
    real(IDP) :: phi_rd2,phi_rd4,phi_rd6,phi_rd8,phi_rda
    integer :: rd2, rd4, rd6, rd8, rda, ic
    integer ::  i, j, l, mjp, lsign, ltype
    character(len=24) :: formath='(1pe13.6,14(a1,1pe15.8))'
    character(len=24) :: formatj='(1pe13.6,19(a1,1pe15.8))'
    character(len=23) :: formatz='(1pe15.8,5(a1,1pe15.8))'
    character(len=23) :: formatz1='(1pe15.8,6(a1,1pe15.8))'
    character(len=32) :: confil
    character(len=6) :: i_diag_str
    character(len=1) :: t
    logical :: last

    epsq=eps*eps
    omcyd = omcy
    vfova2 = vfova*vfova
    !
    !   Calculate delta-B_pol/B at 5 radial locations for output at each time step
    !
    call dbydr_par(sc1,psi,0.0_IDP,1.0_IDP,0)
    call multed(sc3,bmod,1,sc1,1,0.0_IDP,1.0_IDP)
    do l=1,lmax
       sc2(:,l)=uzt(mj_start:mj_end,l)
       uzt(mj_start:mj_end,l)=sc3(:,l)
    end do
    call trnsfr0e(uzt,1)
    if (myPE == 0) then
       t=char(9)
       bthta = 0.0_IDP; btht1 = 0.0_IDP; btht2 = 0.0_IDP
       btht3 = 0.0_IDP; btht4 = 0.0_IDP
       rda=mj;rd2=mj/5;rd4=2*mj/5;rd6=3*mj/5;rd8=4*mj/5
       ltype = +1
       do l=1,lmax
          lsign=signl(l)
          if(ltype*lsign .gt. 0) bthta=bthta+uzt(rda,l)*eps/(feq(rda)-qqinv(rda)*cureq(rda))
          if(ltype*lsign .gt. 0) btht1=btht1+uzt(rd2,l)*eps/(feq(rd2)-qqinv(rd2)*cureq(rd2))
          if(ltype*lsign .gt. 0) btht2=btht2+uzt(rd4,l)*eps/(feq(rd4)-qqinv(rd4)*cureq(rd4))
          if(ltype*lsign .gt. 0) btht3=btht3+uzt(rd6,l)*eps/(feq(rd6)-qqinv(rd6)*cureq(rd6))
          if(ltype*lsign .gt. 0) btht4=btht4+uzt(rd8,l)*eps/(feq(rd8)-qqinv(rd8)*cureq(rd8))
       end do
       !
       !   Calculate integrals of m,n = 0,0 amplitudes for output at each time step
       !
       do j=0,mj
          sd1(j) = phi(j,l0)*r(j)
          sd2(j) = psi(j,l0)*r(j)
          sd3(j) = nf(j,l0)*r(j)
          sd4(j) = vprlf(j,l0)*r(j)
       end do
       mjp = mj + 1
       call quadq(mjp,r,sd1,1._IDP,sd5,sd6,sd7,phi00)
       call quadq(mjp,r,sd2,1._IDP,sd5,sd6,sd7,psi00)
       call quadq(mjp,r,sd3,1._IDP,sd5,sd6,sd7,pr00)
       call quadq(mjp,r,sd4,1._IDP,sd5,sd6,sd7,vprl00)
    end if
    !
    !    Calculate macro transport fluxes, taking into account both potential and magnetic fluctuations.
    !    Three types of data gathering are done: (1) flux surface averaged fluxes are generated on 5
    !    flux surfaces for output at every time step, (2) Radially resolved flux surface averaged flows
    !    are generated for output only at the final time step, (3) 2D transport flows are formed for output
    !    only at the final time step (eventually this should be extended to 3D transport)
    !
    sc5(:,:) = 0.0_IDP; sc6(:,:) = 0.0_IDP
    sc7(:,:) = 0.0_IDP; sc8(:,:) = 0.0_IDP
    sc9(:,:) = 0.0_IDP; sc10(:,:) = 0.0_IDP
    sc11(:,:) = 0.0_IDP; sc12(:,:) = 0.0_IDP

    !  v_ExB_r * nf  (goes into sc6,sc10)

    call dbydth_par(sc5,phi,-1,0.0_IDP,-1.0_IDP,0)
    call mult(sc6,sc5,+1,nf,+1,0.0_IDP,1.0_IDP)
    if(alpha_on == 1) then
       call mult(sc10,sc5,+1,nalp,+1,0.0_IDP,1.0_IDP)
    end if

    !  vf0 * (deltaBr/B) * nf  (goes into sc6,sc12) ---- Need to check the dimensional factors on sc8

    call dbydth_par(sc7,psi,+1,0.0_IDP,-1.0_IDP,0)
    call mult(sc8,sc7,-1,nf,+1,0.0_IDP,1.0_IDP)
    if(alpha_on == 1) then
       call mult(sc12,sc7,-1,nalp,+1,0.0_IDP,1.0_IDP)
    end if

    !  Sum (goes into sc1, sc2)
    sc1(:,:) = 0.0_IDP
    sc2(:,:) = 0.0_IDP
    do l = 1,lmax
       sc8(:,l) = vfova(:)*LcA3*sc8(:,l)/(feq(:)-qqinv(:)*cureq(:))
       sc1(:,l) = sc6(:,l) + sc8(:,l)
    end do

    if(alpha_on == 1)then
       do l = 1,lmax
          sc12(:,l) = vfova(:)*LcA3*sc12(:,l)/(feq(:)-qqinv(:)*cureq(:))
          sc2(:,l) = sc10(:,l) + sc12(:,l)
       end do
    endif

    do l=1,lmax
       uzt(mj_start:mj_end,l)=sc1(:,l)
    end do
    call trnsfr0e(uzt,1)
    rda=mj/2;rd2=mj/10;rd4=2*mj/10;rd6=3*mj/10;rd8=4*mj/10
    flux_nf1 = uzt(rd2,l0)
    flux_nf2 = uzt(rd4,l0)
    flux_nf3 = uzt(rd6,l0)
    flux_nf4 = uzt(rd8,l0)
    flux_nf5 = uzt(rda,l0)

    if(alpha_on == 1) then
       do l=1,lmax
          uzt(mj_start:mj_end,l)=sc2(:,l)
       end do
       call trnsfr0e(uzt,1)
       flux_nalf1 = uzt(rd2,l0)
       flux_nalf2 = uzt(rd4,l0)
       flux_nalf3 = uzt(rd6,l0)
       flux_nalf4 = uzt(rd8,l0)
       flux_nalf5 = uzt(rda,l0)        
    endif

    if(myPE ==0) then

       if(alpha_on == 0) then
          write(82,formath) time,t,btht1,t,btht2,t,btht3,t,btht4,  &
               t,bthta,t,phi00,t,psi00,t,pr00,t,vprl00,      &
               t,flux_nf1,t,flux_nf2,t,flux_nf3,             &
               t,flux_nf4,t,flux_nf5
       else if(alpha_on == 1) then
          write(82,formatj) time,t,btht1,t,btht2,t,btht3,t,btht4,     &
               t,bthta,t,phi00,t,psi00,t,pr00,t,vprl00,         &
               t,flux_nf1,t,flux_nf2,t,flux_nf3,                &
               t,flux_nf4,t,flux_nf5,t,flux_nalf1,t,flux_nalf2, &
               t,flux_nalf3,t,flux_nalf4,t,flux_nalf5
       endif

       if (last .or. mod(nstep,ndiag) == 0) then
          write(i_diag_str, '("_",I5.5)') idiag
          confil = "final_n0_data_"//numrun(1)//numrun(2)//i_diag_str
          open(unit=83,file=confil,status="unknown")
          if(alpha_on == 0) then
             write(83,'("r",a1,"flux",a1,"phi_00",a1,"psi_00",a1,"nf_00",a1,"vprlr_00")') t,t,t,t,t
             do j=1,mj
                write(83,formatz) r(j),t,uzt(j,l0),t,phi(j,l0),t,psi(j,l0),t,nf(j,l0),t,vprlf(j,l0)
             end do
          else if(alpha_on == 1) then
             write(83,'("r",a1,"flux",a1,"phi_00",a1,"psi_00",a1,  &
                  "nf_00",a1,"nalp_00",a1,"vprlr_00")') t,t,t,t,t,t
             do j=1,mj
                write(83,formatz1) r(j),t,uzt(j,l0),t,phi(j,l0),t,psi(j,l0),t,nf(j,l0),t,nalp(j,l0),t,vprlf(j,l0)
             end do
          end if
          close(unit=83)

          !confil = "final_2D_flux_ExB_"//numrun(1)//numrun(2)//i_diag_str
          !open(unit=84,file=confil,status="unknown")
          !confil = "final_2D_flux_deltaB_"//numrun(1)//numrun(2)//i_diag_str
          !open(unit=85,file=confil,status="unknown")
          
          !ic = 0
          !do l=1,lmax
          !   if(nn(l) .eq. 0) ic = ic + 1
          !end do
          !write(84,*) mj,ic
          !do l=1,lmax
          !   do j=1,mj
          !      if(nn(l) .eq. 0) write(84,'(i5,2x,i4,2x,i3,2(2x,e15.8))') l,mm(l),signl(l),r(j),sc6(j,l)
          !      if(nn(l) .eq. 0) write(85,'(i5,2x,i4,2x,i3,2(2x,e15.8))') l,mm(l),signl(l),r(j),sc8(j,l)
          !   end do
          !end do
          !close(unit=84)
          !close(unit=85)
          idiag=idiag+1

       end if  !if(last)
    end if !myPE == 0

    do l=1,lmax
       uzt(mj_start:mj_end,l)= phi(:,l)
    end do
    call trnsfr0e(uzt,-1)

    if(myPE == 0) then
       phi_rd2 = 0; phi_rd4 = 0; phi_rd6 = 0; phi_rd8 = 0; phi_rda = 0
       !    Here we are adding up phi at theta = 0, zeta = 0
       do l=1,lmax
          phi_rd2 = phi_rd2 + uzt(rd2,l)
          phi_rd4 = phi_rd4 + uzt(rd4,l)
          phi_rd6 = phi_rd6 + uzt(rd6,l)
          phi_rd8 = phi_rd8 + uzt(rd8,l)
          phi_rda = phi_rda + uzt(rda,l)
       end do

       write(77,'(1pe13.6,10(a1,e15.8))')  &
            time,t,phi_rd2,t,phi_rd4,t,phi_rd6,t,phi_rd8,t,phi_rda,  &
            t,uzt(rd2,l0),t,uzt(rd4,l0),t,uzt(rd6,l0),t,uzt(rd8,l0), &
            t,uzt(rda,l0)
    end if

    do l=1,lmax
       uzt(mj_start:mj_end,l)=sc2(:,l)
    end do

    sc1(:,:) = 0.0_IDP; sc2(:,:) = 0.0_IDP; sc3(:,:) = 0.0_IDP
    sc4(:,:) = 0.0_IDP; sc5(:,:) = 0.0_IDP; sc6(:,:) = 0.0_IDP 

  end subroutine hifreq

END MODULE diagnostics
