MODULE output_mod

  use param
  use cotrol
  use domain
  use equil
  use dynamo
  use dbyd
  use scratch

  IMPLICIT NONE

CONTAINS

  subroutine output

    implicit none

    integer :: i,nd,l,m,n,lh,lheq,le,j

    write(6,'(//52("*"),2a8,51("*")/)') ndevice
    if (nocpl == 0) then
       write(6,'(50("*")," toroidal geometry ",50("*")/)')
    else
       write(6,'(49("*")," cylindrical geometry ",48("*")/)')
    end if
    if (nonlin == 0) write(6,'(56("*")," linear ",55("*")/)')
    if (nonlin /= 0) write(6,'(54("*")," nonlinear ",54("*")//)')

    write(6,'(/"control:"/"numrun  =     ",2a2,a1," numruno =     ",2a2,a1," numruns =     ",2a2,a1, &  
         " nstep   =",i10," maxstp  =",i10," nonlin  =",i10/"ndump   =",i10," nprint  =",i10," ndiag   =",i10, &
         " lplots  =",i10," itime   =",i10," dt0     =",1pe10.3/)') numrun,numruno,numruns,nstep,maxstp,nonlin, &
         ndump,nprint,ndiag,lplots,itime,dt0

    write(6,'("numhist =",12("   ",a6))') (numhist(i),i=1,ihist)

    if (Auto_grid_on == 0) then
       write(6,'(/"domain:"/"mj      =",i10," ni      =",i10," nis     =",i10," ne      =",i10," rc      =",0pf10.7, &
            " delta   =",0pf10.7/"lmax    =",i10," leqmax  =",i10," mmin    =",i10," mmax    =",i10, &
            " nmin    =",i10," nmax    =",i10/"mmineq  =",i10," mmaxeq  =",i10," nmineq  =",i10," nmaxeq  =",i10, &
            " l0      =",i10," leq0    =",i10)') &
            mj,ni,nis,ne,rc,delta,lmax,leqmax,mmin,mmax,nmin,nmax,mmineq,mmaxeq,nmineq,nmaxeq,l0,leq0
    else
       write(6,'(/"domain:"/"mj      =",i10,"  grid with evenly spaced points" &
            /"lmax    =",i10," leqmax  =",i10," mmin    =",i10," mmax    =",i10, &
            " nmin    =",i10," nmax    =",i10/"mmineq  =",i10," mmaxeq  =",i10," nmineq  =",i10," nmaxeq  =",i10, &
            " l0      =",i10," leq0    =",i10)') &
            mj,lmax,leqmax,mmin,mmax,nmin,nmax,mmineq,mmaxeq,nmineq,nmaxeq,l0,leq0
    end if

    write(6,'(/"equil:")')
    if (ext_prof == 0) then
       write(6,'(a)') "eq name:  "//trim(eq_name)
    else
       write(6,'(a)') "eq name:  "//trim(eq_name)//"   ext prof name:  "//trim(ext_prof_name)
    end if
    if (alpha_on == 0) then
       write(6,'("ext_prof        =",i2," alpha_on        =",i2," trapped_on      =",i2," EP_dens_on      =",i2)') &
            ext_prof,alpha_on,trapped_on,EP_dens_on
       write(6,'("R0      =",0pf10.6," B0      =",0pf10.6," eps     =",0pf10.6," bet0    =",1pe10.3," bet0_f  =",1pe10.3, &
            " omcy    =",1pe10.3)') bigrn,bmodn,eps,bet0,bet0_f,omcy
       if (trapped_on == 1) write(6,'("omcyb   =",1pe10.3," rbound  =",1pe10.3)') omcyb,rbound
       if (ext_prof == 0) then
          if (EP_dens_on == 1) then
             write(6,'("Adens   =",1pe10.3," Bdens   =",1pe10.3)') Adens,Bdens
          else
             if (cnfp(0) > 0.0_IDP) write(6,'("cnfp            =",1p6e17.9)') (cnfp(i),i=0,10)
          end if
          if (cnep(0) > 0.0_IDP) write(6,'("cnep            =",1p6e17.9)') (cnep(i),i=0,10)
          if (ctep(0) > 0.0_IDP) write(6,'("ctep            =",1p6e17.9)') (ctep(i),i=0,10)
          if (cvep(0) > 0.0_IDP) write(6,'("cvep            =",1p6e17.9)') (cvep(i),i=0,10)
          if (cvfp(0) > 0.0_IDP) write(6,'("cvfp            =",1p6e17.9)') (cvfp(i),i=0,10)
          if (eqvt(0) > 0.0_IDP) write(6,'("eqvt            =",1p6e17.9)') (eqvt(i),i=0,10)
          if (eqvp(0) > 0.0_IDP) write(6,'("eqvp            =",1p6e17.9)') (eqvp(i),i=0,10)
       else
          write(6,'("DIIID_u         =",i2," q_prof_on       =",i2," EP_vel_on       =",i2," Alpha_vel_on    =",i2, &
               " Eq_presseq_on   =",i2," Eq_presstot_on  =",i2/"Eq_vel_on       =",i2," Eq_velp_on      =",i2, &
               " spe1            =",i2," spe2            =",i2)') DIIID_u,q_prof_on,EP_vel_on,Alpha_vel_on, &
               Eq_presseq_on,Eq_presstot_on,Eq_vel_on,Eq_velp_on,spe1,spe2
       end if
    else
       write(6,'("ext_prof        =",i2," alpha_on        =",i2," trapped_on      =",i2," EP_dens_on      =",i2, &
            " Alpha_dens_on   =",i2)') ext_prof,alpha_on,trapped_on,EP_dens_on,Alpha_dens_on
       if (trapped_on == 0) then
          write(6,'("R0      =",0pf10.6," B0      =",0pf10.6," eps     =",0pf10.6," bet0    =",1pe10.3," bet0_f  =", &
               1pe10.3," bet0_alp=",1pe10.3/"omcy    =",1pe10.3," omcyalp =",1pe10.3)') bigrn,bmodn,eps,bet0, &
               bet0_f,bet0_alp,omcy,omcyalp
       else
          write(6,'("R0      =",0pf10.6," B0      =",0pf10.6," eps     =",0pf10.6," bet0    =",1pe10.3," bet0_f  =", &
               1pe10.3," bet0_alp=",1pe10.3/"omcy    =",1pe10.3," omcyalp =",1pe10.3," omcyb   =",1pe10.3, &
               " rbound  =",1pe10.3)') bigrn,bmodn,eps,bet0,bet0_f,bet0_alp,omcy,omcyalp,omcyb,rbound
       end if
       if (ext_prof == 0) then
          if (EP_dens_on == 1) then
             if (Alpha_dens_on == 1) then
                write(6,'("Adens   =",1pe10.3," Bdens   =",1pe10.3," Adensalp=",1pe10.3," Bdensalp=",1pe10.3)') &
                     Adens,Bdens,Adensalp,Bdensalp
             else
                write(6,'("Adens   =",1pe10.3," Bdens   =",1pe10.3)') Adens,Bdens
                if (cnfpalp(0) > 0.0_IDP) write(6,'("cnfpalp         =",1p6e17.9)') (cnfpalp(i),i=0,10)
             end if
          else
             if (Alpha_dens_on == 1) then
                write(6,'("Adensalp=",1pe10.3," Bdensalp=",1pe10.3)') Adensalp,Bdensalp
                if (cnfp(0) > 0.0_IDP) write(6,'("cnfp            =",1p6e17.9)') (cnfp(i),i=0,10)
             else
                if (cnfp(0) > 0.0_IDP) write(6,'("cnfp            =",1p6e17.9)') (cnfp(i),i=0,10)
                if (cnfpalp(0) > 0.0_IDP) write(6,'("cnfpalp         =",1p6e17.9)') (cnfpalp(i),i=0,10)
             end if
          end if
          if (cnep(0) > 0.0_IDP) write(6,'("cnep            =",1p6e17.9)') (cnep(i),i=0,10)
          if (ctep(0) > 0.0_IDP) write(6,'("ctep            =",1p6e17.9)') (ctep(i),i=0,10)
          if (cvep(0) > 0.0_IDP) write(6,'("cvep            =",1p6e17.9)') (cvep(i),i=0,10)
          if (cvfp(0) > 0.0_IDP) write(6,'("cvfp            =",1p6e17.9)') (cvfp(i),i=0,10)
          if (cvfpalp(0) > 0.0_IDP) write(6,'("cvfpalp         =",1p6e17.9)') (cvfpalp(i),i=0,10)
          if (eqvt(0) > 0.0_IDP) write(6,'("eqvt            =",1p6e17.9)') (eqvt(i),i=0,10)
          if (eqvp(0) > 0.0_IDP) write(6,'("eqvp            =",1p6e17.9)') (eqvp(i),i=0,10)
       else
          write(6,'("DIIID_u         =",i2," q_prof_on       =",i2," EP_vel_on       =",i2, &
               " Eq_presseq_on   =",i2," Eq_presstot_on  =",i2/"Eq_vel_on       =",i2," Eq_velp_on      =",i2, &
               " spe1            =",i2," spe2            =",i2)') DIIID_u,q_prof_on,EP_vel_on, &
               Eq_presseq_on,Eq_presstot_on,Eq_vel_on,Eq_velp_on,spe1,spe2
       end if
    end if
    write(6,'("Edge_on =",i10," edge_p  =",i10)') Edge_on,edge_p

    write(6,'(/"dynamo:"/"s       =",1pe10.3," ietaeq  =",i10," etascl  =",1pe10.3," eta0    =",1pe10.3, &
         " reta    =",0pf10.7," etalmb  =",1pe10.3/"gamma   =",1pe10.3," LcA0    =",1pe10.3, &
         " LcA1    =",1pe10.3," LcA2    =",1pe10.3," LcA3    =",1pe10.3," omegar  =",1pe10.3/"epflr_on=",i10, &
         " r_epflr =",1pe10.3," iflr_on =",i10," iflr    =",1pe10.3," twofl_on=",i10," dpres   =",1pe10.3/ &
         "B_par_on=",i10)') s,ietaeq,etascl,eta0,reta,etalmb,gamma,LcA0,LcA1,LcA2,LcA3,omegar,epflr_on,r_epflr, &
         iflr_on,iflr,twofl_on,dpres,B_par_on
    if (alpha_on == 1) write(6,'("LcA0alp =",1pe10.3," LcA1alp =",1pe10.3," LcA2alp =",1pe10.3," LcA3alp =",1pe10.3, &
         " r_flralp=",1pe10.3)') LcA0alp,LcA1alp,LcA2alp,LcA3alp,r_epflralp
    if (ext_prof == 1) write(6,'("ieldamp_on =",i7," uion    =",1pe10.3," xnuelc0 =",1pe10.3)') ieldamp_on,uion,xnuelc0
    write(6,'("m0dy    =",i10," ipert   =",i10," pertscl =",1pe10.3," xle     =",1pe10.3," dt      =",1pe10.3, &
         " time    =",1pe10.3/"difnr_on=",i10," stdifp  =",1pe10.3," stdifu  =",1pe10.3," stdifv  =",1pe10.3, &
         " stdifnf =",1pe10.3," stdifvf =",1pe10.3)') m0dy,ipert,pertscl,xle,dt,time,difnr_on,stdifp,stdifu, &
         stdifv,stdifnf,stdifvf
    if (alpha_on == 1) write(6,'("difnalp =",1pe10.3," difvalp =",1pe10.3)') stdifnalp,stdifvalp  
    write(6,'("nopsievol_on    =",i2," noprevol_on     =",i2," nonfevol_on     =",i2," nonalpevol_on   =",i2/ &
         "src_sink_th_on  =",i2," src_sink_EP1_on =",i2," src_sink_EP2_on =",i2," src_sink_DIIID  =",i2, &
         " src_sink_ITER   =",i2)') nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on,src_sink_th_on, &
         src_sink_EP1_on,src_sink_EP2_on,src_sink_DIIID_on,src_sink_ITER_on
    write(6,'("asrc    =",1pe10.3," rsrc    =",1pe10.3," wsrc    =",1pe10.3," asrc_EP1=",1pe10.3," rsrc_EP1=",1pe10.3, &
         " wsrc_EP1=",1pe10.3)') asrc,rsrc,wsrc,asrc_EP1,rsrc_EP1,wsrc_EP1
    if (alpha_on == 1) write(6,'("asrc_EP2=",1pe10.3," rsrc_EP2=",1pe10.3," wsrc_EP2=",1pe10.3)') asrc_EP2,rsrc_EP2, &
         wsrc_EP2
    if (srcsinkth(0) > 0.0_IDP) write(6,'("srcsinkth       =",1p6e17.9)') (srcsinkth(i),i=0,10)
    if (srcsinkEP1(0) > 0.0_IDP) write(6,'("srcsinkEP1      =",1p6e17.9)') (srcsinkEP1(i),i=0,10)
    if (srcsinkEP2(0) > 0.0_IDP) write(6,'("srcsinkEP2      =",1p6e17.9)') (srcsinkEP2(i),i=0,10)

    if (betath_factor /= 1.0_IDP) write(6,'(" WARNING THERMAL BETA FACTOR ACTIVE = ",1pe13.6)') betath_factor 
    if (deltaq /= 0) write(6,'(" WARNING SAFETY FACTOR PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaq
    if (deltaiota /= 0) write(6,'(" WARNING IOTA PROFILE DISPLACEMET ACTIVE = ",1pe13.6)') deltaiota

    write(6,'(/"ll(m,n)")')
    do m=mmax,mmin,-1
       write(6,'(i4,"  ",31i4)') m,(ll(m,n),n=nmin,nmax)
    end do
    write(6,'("   m/n",31i4)') (n,n=nmin,nmax)

    write(6,'(/"helicities"/"   lh   mh   nh")')
    write(6,'(3i5)') (lh,mh(lh),nh(lh),lh=1,lhmax)

    write(6,'(/"    l   mm   nn          signl             rs         widthi         gammai")')
    write(6,'(3i5,"     ",i10,"     ",0pf10.5,1p2e15.3)') (l,mm(l),nn(l),signl(l),rs(l),widthi(l),gammai(l),l=1,lmax)

    write(6,'(/"lleq(m,n)")')
    do m=mmaxeq,mmineq,-1
       write(6,'(i3,"    ",37i3)') m,(lleq(m,n),n=nmineq,nmaxeq)
    end do
    write(6,'("meq/neq",37i3)') (n,n=nmineq,nmaxeq)

    write(6,'(/"eq helicities"/" lheq mheq nheq")')
    write(6,'(3i5)') (lheq,mheq(lheq),nheq(lheq),lheq=1,lheqmx)

    write(6,'(/"  leq mmeq nneq         sgnleq")')
    write(6,'(3i5,"     ",i10)') (le,mmeq(le),nneq(le),sgnleq(le),le=1,leqmax)

    write(6,'(/"  linear matrix mode entries: lmaxn,lln-",i5/("                                        ",12i5))') lmaxn, &
         (lln(i),i=1,lmaxn)

    if (alpha_on == 0) then
       write(6,'(/"           r            q            J            I         preq       denseq         teeq", &
            "         nfeq        vfova          eta")')
       do j=0,mj
          write(6,'(0pf12.8,1p9e13.5)') r(j),qq(j),feq(j),cureq(j),preq(j),denseq(j),teeq(j),nfeq(j),vfova(j),eta(j)
       end do
    else
       write(6,'(/"           r            q            J            I         preq       denseq         teeq", &
            "         nfeq        vfova       nalpeq    valphaova          eta")')
       do j=0,mj
          write(6,'(0pf12.8,1p11e13.5)') r(j),qq(j),feq(j),cureq(j),preq(j),denseq(j),teeq(j),nfeq(j),vfova(j), &
               nalpeq(j),valphaova(j),eta(j)
       end do
    end if

  end subroutine output

  subroutine endrun

    use processor
    use var_para
    use dbyd
    use mult_mod
    use transfer

    implicit none

    integer :: i,j,l,lwrt,l2t,mnum
    real(IDP) :: scmn,scmx,sclmn,sclmx
    real(IDP), dimension(:), allocatable :: scnorm
    character(len=1) :: tb
    character(len=12) :: confil
    character(len=32) :: formatt='("r",1200(a1,i4,"/",i4))'
    character(len=32) :: formatv='(1pe13.6,1200(a1,1pe15.8))'

    interface
       subroutine delstar(ss,ff,itypf,wk1,wk2,wkeq1,c1,c2)
         use param
         use var_para
         implicit none
         integer :: itypf
         real(IDP) :: c1,c2
         real(IDP), dimension(mj_start:,0:) :: ss,ff,wk1,wk2
         real(IDP), dimension(0:,0:) :: wkeq1
       end subroutine delstar
    end interface

    tb=char(9)

    !  store uzt-values into sc4

    do l=1,lmax
       sc4(:,l)=uzt(mj_start:mj_end,l)
    end do

    allocate (scnorm(nnum))
    scnorm=1.0_IDP

    if (myPE == 0) then

       !  vr up
       call dbydth(uzt,phi,-1,0.0_IDP,-eps,0)

       do i=1,nnum
          mnum=mnumn(i)
          if (nonlin == 0) then
             scmn=0.
             scmx=0.
             do l2t=1,mnum
                l=l2t+lnumn(i-1)
                sclmn=minval(uzt(:,lln(l)))
                sclmx=maxval(uzt(:,lln(l)))
                scmn=min(scmn,sclmn)
                scmx=max(scmx,sclmx)
             end do
             if (scmx > abs(scmn)) then
                scnorm(i)=scmx
             else
                scnorm(i)=scmn
             endif
          end if
          if(scnorm(i) /= 0.0_IDP) then
             do l2t=1,mnum
                l=l2t+lnumn(i-1)
                psi(:,lln(l))=psi(:,lln(l))/scnorm(i)
                where (abs(psi(:,lln(l))) < 1.e-50_IDP) psi(:,lln(l))=0
                phi(:,lln(l))=phi(:,lln(l))/scnorm(i)
                where (abs(phi(:,lln(l))) < 1.e-50_IDP) phi(:,lln(l))=0
                pr(:,lln(l))=pr(:,lln(l))/scnorm(i)
                where (abs(pr(:,lln(l))) < 1.e-50_IDP) pr(:,lln(l))=0
                nf(:,lln(l))=nf(:,lln(l))/scnorm(i)
                where (abs(nf(:,lln(l))) < 1.e-50_IDP) nf(:,lln(l))=0
                vprlf(:,lln(l))=vprlf(:,lln(l))/scnorm(i)
                where (abs(vprlf(:,lln(l))) < 1.e-50_IDP) vprlf(:,lln(l))=0
                vthprlf(:,lln(l))=vthprlf(:,lln(l))/scnorm(i)
                where (abs(vthprlf(:,lln(l))) < 1.e-50_IDP) vthprlf(:,lln(l))=0
                uzt(:,lln(l))=uzt(:,lln(l))/scnorm(i)
                where (abs(uzt(:,lln(l))) < 1.e-50_IDP) uzt(:,lln(l))=0
             end do
             if (alpha_on == 1) then
                do l2t=1,mnum
                   l=l2t+lnumn(i-1)
                   nalp(:,lln(l))=nalp(:,lln(l))/scnorm(i)
                   where (abs(nalp(:,lln(l))) < 1.e-50_IDP) nalp(:,lln(l))=0
                   vprlalp(:,lln(l))=vprlalp(:,lln(l))/scnorm(i)
                   where (abs(vprlalp(:,lln(l))) < 1.e-50_IDP) vprlalp(:,lln(l))=0
                end do
             end if
          end if
       end do
       lwrt=abs(lplots)
       if (lwrt > lmax) lwrt=lmax
       write(confil,'("phi_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,phi(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("psi_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,psi(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("pr_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,pr(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("vr_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("nf_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,nf(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("vprlf_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,vprlf(j,l),l=1,lwrt)
       end do
       close(92)
       write(confil,'("vthprlf_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,vthprlf(j,l),l=1,lwrt)
       end do
       close(92)
       if (alpha_on == 1) then
          write(confil,'("nalp_",2a2)') numrun(1),numrun(2)
          open(unit=92,file=confil,recl=19384)
          write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
          do j=0,mj
             write(92,formatv) r(j),(tb,nalp(j,l),l=1,lwrt)
          end do
          close(92)
          write(confil,'("vprlalp_",2a2)') numrun(1),numrun(2)
          open(unit=92,file=confil,recl=19384)
          write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
          do j=0,mj
             write(92,formatv) r(j),(tb,vprlalp(j,l),l=1,lwrt)
          end do
          close(92)
       end if

       !  vth up
       call dbydr(uzt,phi,0.0_IDP,eps,0)
       write(confil,'("vth_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
    end if

    !  br up
    call dbydth_par(sc3,psi,1,0.0_IDP,-eps,0)
    !  bth up
    call dbydr_par(sc5,psi,0.0_IDP,eps,0)

    call multed(uzt,sqgi,1,sc3,-1,0.0_IDP,1.0_IDP)
    call trnsfr0(uzt,-1)
    if (myPE == 0) then
       write(confil,'("br_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
    end if

    call multed(uzt,sqgi,1,sc5,1,0.0_IDP,1.0_IDP)
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       write(confil,'("bth_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
    end if

    call delstar(uzt,psi,1,sc2,sc3,sceq1,0.0_IDP,1.0_IDP)
    uzt=eps*uzt
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       write(confil,'("curzt_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
    end if

    !  restore uzt-values

    do l=1,lmax
       uzt(mj_start:mj_end,l)=sc4(:,l)
    end do
    call trnsfr0(uzt,1)
    if (myPE == 0) then
       do i=1,nnum
          mnum=mnumn(i)
          if (nonlin == 0) then
             scmn=0.
             scmx=0.
             do l2t=1,mnum
                l=l2t+lnumn(i-1)
                sclmn=minval(uzt(:,lln(l)))
                sclmx=maxval(uzt(:,lln(l)))
                scmn=min(scmn,sclmn)
                scmx=max(scmx,sclmx)
             end do
             if (scmx > abs(scmn)) then
                scnorm(i)=scmx
             else
                scnorm(i)=scmn
             endif
          end if
          if(scnorm(i) /= 0.0_IDP) then
             do l2t=1,mnum
                l=l2t+lnumn(i-1)
                uzt(:,lln(l))=uzt(:,lln(l))/scnorm(i)
                where (abs(uzt(:,lln(l))) < 1.e-50_IDP) uzt(:,lln(l))=0
             end do
          end if
       end do
       write(confil,'("uzt_",2a2)') numrun(1),numrun(2)
       open(unit=92,file=confil,recl=19384)
       write(92,formatt) (tb,mm(l),nn(l),l=1,lwrt)
       do j=0,mj
          write(92,formatv) r(j),(tb,uzt(j,l),l=1,lwrt)
       end do
       close(92)
    end if

  end subroutine endrun

  subroutine rddump

    implicit none

    integer :: i,l,j,ihisto,idum
    real(IDP) :: pertsclo,dum
    character(len=2), dimension(3) :: numrunp
    character(len=5) :: numvac

    if (old_rd) then
       read(8) ihisto,numruno,numrunp,numruns,idum,numvac,(numhist(i),i=1,ihisto),maxstp,nstep, &
            ndump,nprint,lplots,idum,idum,idum,idum,itime,dt0,nstep1,nonlin,m0dy,nocpl,noprevol_on,idum, &
            Adens,Bdens,LcA0,LcA1,LcA2,LcA3,omcy,dpres,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,twofl_on,ieldamp_on, &
            omegar,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp
    else
       read(8) ihisto,numruno,numrunp,numruns,(numhist(i),i=1,ihisto),maxstp,nstep, &
            ndump,nprint,lplots,itime,dt0,nstep1,nonlin,m0dy,nocpl,noprevol_on, &
            Adens,Bdens,LcA0,LcA1,LcA2,LcA3,omcy,dpres,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,twofl_on,ieldamp_on, &
            omegar,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp,ext_prof_name, &
            trapped_on,omcyb,rbound,B_par_on
    end if
    read(8) mj,lmaxo,leqmax,mjm1
    if (lmax < lmaxo) then
       read(8) (r(j),j=0,mj),(mm(l),l=1,lmax),(idum,l=lmax+1,lmaxo),(nn(l),l=1,lmax),(idum,l=lmax+1,lmaxo), &
            (mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &  
            (dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
    else
       read(8) (r(j),j=0,mj),(mm(l),l=1,lmaxo),(nn(l),l=1,lmaxo),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &  
            (dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
    end if

    allocate (sqgi(0:mj,0:leqmax),sqg(0:mj,0:leqmax),bst(0:mj,0:leqmax),grr(0:mj,0:leqmax),grt(0:mj,0:leqmax), &
         gtt(0:mj,0:leqmax),grroj(0:mj,0:leqmax),grtoj(0:mj,0:leqmax),gttoj(0:mj,0:leqmax),jbgrr(0:mj,0:leqmax), &
         jbgrt(0:mj,0:leqmax),jbgtt(0:mj,0:leqmax),lplrr(0:mj,0:leqmax),lplrt(0:mj,0:leqmax),lplrz(0:mj,0:leqmax), &
         lpltt(0:mj,0:leqmax),lpltz(0:mj,0:leqmax),lplzz(0:mj,0:leqmax),lplr(0:mj,0:leqmax),lplt(0:mj,0:leqmax), &
         lplz(0:mj,0:leqmax),djroj(0:mj,0:leqmax),djtoj(0:mj,0:leqmax),djzoj(0:mj,0:leqmax),omdr(0:mj,0:leqmax), &
         omdt(0:mj,0:leqmax),omdz(0:mj,0:leqmax),bmod(0:mj,0:leqmax),dbsjtoj(0:mj,0:leqmax),dbsjzoj(0:mj,0:leqmax))

    read(8) ni,nis,ne,delta,rc,fti,fte
    if (old_rd) then
       read(8) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(dum,j=0,mj),(dum,j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
            (cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=0,mj),(vfova(j),j=0,mj),(vzt_eq(j),j=0,mj), &
            (tieq(j),j=0,mj)
    else
       read(8) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
            (cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=0,mj),(vfova(j),j=0,mj),(vzt_eq(j),j=0,mj), &
            (tieq(j),j=0,mj),(vth_eq(j),j=0,mj)
    end if
    if (alpha_on == 1) read(8) (nalpeq(j),j=0,mj),(valphaova(j),j=0,mj)
    if (ext_prof == 1) then
       allocate (vAlfven(0:mj),vtherm_elecP(0:mj))
       read(8) (vAlfven(j),j=0,mj),(vtherm_elecP(j),j=0,mj)
    else
       read(8) (vtherm_elc(j),j=0,mj)
    end if
    if (old_rd) then
       if (lmax < lmaxo) then
          read(8) idum,ndevice(1),ndevice(2),idum,eps,bet0,bet0_f,idum,dum,dum,dum,(rs(l),l=1,lmax)
       else
          read(8) idum,ndevice(1),ndevice(2),idum,eps,bet0,bet0_f,idum,dum,dum,dum,(rs(l),l=1,lmaxo)
       end if
    else
       if (lmax < lmaxo) then
          read(8) ndevice(1),ndevice(2),eps,bet0,bet0_f,(rs(l),l=1,lmax)
       else
          read(8) ndevice(1),ndevice(2),eps,bet0,bet0_f,(rs(l),l=1,lmaxo)
       end if
    end if
    read(8) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
         ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
         ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
         ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
         ((lplrr(j,l),j=0,mj),l=1,leqmax),((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
         ((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax),((lplzz(j,l),j=0,mj),l=1,leqmax), &
         ((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &
         ((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &
         ((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &
         ((bmod(j,l),j=0,mj),l=1,leqmax),((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax)
    if (ieldamp_on == 1) then
       allocate (eildr(0:mj,0:leqmax),eildt(0:mj,0:leqmax),eildz(0:mj,0:leqmax),eildrr(0:mj,0:leqmax), &
            eildrt(0:mj,0:leqmax),eildrz(0:mj,0:leqmax),eildtt(0:mj,0:leqmax),eildtz(0:mj,0:leqmax), &
            eildzz(0:mj,0:leqmax))
       read(8) ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax), &
            ((eildz(j,l),j=0,mj),l=1,leqmax),((eildrr(j,l),j=0,mj),l=1,leqmax), &
            ((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &
            ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax), &
            ((eildzz(j,l),j=0,mj),l=1,leqmax)
    end if
    if (trapped_on == 1) then
       allocate (omdrprp(0:mj,0:leqmax),omdtprp(0:mj,0:leqmax),omdzprp(0:mj,0:leqmax))
       read(8) ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax), &
            ((omdzprp(j,l),j=0,mj),l=1,leqmax)
    end if

    read(8) (eta(j),j=0,mj),dt,time

  end subroutine rddump

  subroutine wrdump(i_diag, is_restart)

    implicit none

    integer, intent(in) :: i_diag
    logical, intent(in) :: is_restart
    integer :: i,l,j
    character(len=12) :: confil
    character(len=6) :: i_diag_str

    if (is_restart) then
       confil="fs"//numrun(1)//numrun(2)//numrun(3)
       open(unit=7,file=confil,status='new',convert='big_endian',form='unformatted')
    else
       write(i_diag_str, '("_",I5.5)') i_diag
       confil="fs"//numrun(1)//numrun(2)//i_diag_str
       open(unit=7,file=confil,status='new',convert='big_endian',form='unformatted')   
    end if

    write(7) ihist,numrun,numruno,numruns,(numhist(i),i=1,ihist),maxstp,nstep, &
         ndump,nprint,lplots,itime,dt0,nstep1,nonlin,m0dy,nocpl,noprevol_on, &
         Adens,Bdens,LcA0,LcA1,LcA2,LcA3,omcy,dpres,ext_prof,epflr_on,r_epflr,alpha_on,iflr_on,iflr,twofl_on,ieldamp_on, &
         omegar,bet0_alp,Adensalp,Bdensalp,LcA0alp,LcA1alp,LcA2alp,LcA3alp,omcyalp,r_epflralp,ext_prof_name, &
         trapped_on,omcyb,rbound,B_par_on
    write(7) mj,lmax,leqmax,mjm1
    write(7) (r(j),j=0,mj),(mm(l),l=1,lmax),(nn(l),l=1,lmax),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &  
         (dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
    write(7) ni,nis,ne,delta,rc,fti,fte
    write(7) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
         (cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=0,mj),(vfova(j),j=0,mj),(vzt_eq(j),j=0,mj), &
         (tieq(j),j=0,mj),(vth_eq(j),j=0,mj)
    if (alpha_on == 1) write(7) (nalpeq(j),j=0,mj),(valphaova(j),j=0,mj)
    if (ext_prof == 1) then
       write(7) (vAlfven(j),j=0,mj),(vtherm_elecP(j),j=0,mj)
    else
       write(7) (vtherm_elc(j),j=0,mj)
    end if
    write(7) ndevice(1),ndevice(2),eps,bet0,bet0_f,(rs(l),l=1,lmax)
    write(7) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
         ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
         ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
         ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
         ((lplrr(j,l),j=0,mj),l=1,leqmax),((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
         ((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax),((lplzz(j,l),j=0,mj),l=1,leqmax), &
         ((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &
         ((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &
         ((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &
         ((bmod(j,l),j=0,mj),l=1,leqmax),((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax)
    if (ieldamp_on == 1) write(7) ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax), &
         ((eildz(j,l),j=0,mj),l=1,leqmax),((eildrr(j,l),j=0,mj),l=1,leqmax), &
         ((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &
         ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax), &
         ((eildzz(j,l),j=0,mj),l=1,leqmax)
    if (trapped_on == 1) write(7) ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax), &
         ((omdzprp(j,l),j=0,mj),l=1,leqmax)

    write(7) (eta(j),j=0,mj),dt,time

    write(7) ((psi(j,l),j=0,mj),l=1,lmax)
    write(7) ((phi(j,l),j=0,mj),l=1,lmax)
    write(7) ((pr(j,l),j=0,mj),l=1,lmax)
    write(7) ((nf(j,l),j=0,mj),l=1,lmax)
    write(7) ((vprlf(j,l),j=0,mj),l=1,lmax)
    write(7) ((vthprlf(j,l),j=0,mj),l=1,lmax)
    if (alpha_on == 1) then
       write(7) ((nalp(j,l),j=0,mj),l=1,lmax)
       write(7) ((vprlalp(j,l),j=0,mj),l=1,lmax)
    end if

    write(7) etascl,reta,eta0,etalmb,ietaeq,stdifp,stdifu,stdifnf,stdifv,stdifvf, &
         stdifnalp,stdifvalp,s,gamma,xle,ipert,(widthi(l),l=1,lmax),(gammai(l),l=1,lmax),pertscl
    write(7) ext_prof_name,difnr_on,nopsievol_on,noprevol_on,nonfevol_on,nonalpevol_on, &
         src_sink_th_on,src_sink_EP1_on,src_sink_EP2_on,src_sink_DIIID_on,src_sink_ITER_on,rsrc,wsrc,asrc, &
         rsrc_EP1,wsrc_EP1,asrc_EP1,rsrc_EP2,wsrc_EP2,asrc_EP2,AWfctr,Nfctr,AWfctr_dif,Rfctr,Wfctr
    write(7) (srcsinkth(i),i=0,10),(srcsinkEP1(i),i=0,10),(srcsinkEP2(i),i=0,10)

    rewind(7)
    close(unit=7)
    write (6,'(/"  wrdump: have written fs",2a2,a1)') (numrun(i),i=1,3)

  end subroutine wrdump

END MODULE output_mod
