	subroutine vmec

		use param
		use cotrol
		use domain
		use equil
		use scratch
		implicit none
 
		integer :: j,l,l1,lq,mjeqp,mjeqm1,mpol,ntor,lb0,m,m2,mm1,mm2,mm12,mm22,mm2p1,mm2p12,mp2,mp22,isg,idummy,mband,n,nnn,mmm,lbm2,ir,jp,lp,nbst
		real(IDP) :: bigrn,pror,dummy,bet0in,twopi,twophip,one,gc,iota,iotap,p,pprimel,pprime,qmin,qmax,xl
		integer, dimension(:), allocatable :: llc,lls,lbst	
		real(IDP), dimension(:), allocatable :: rbinv,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,rsb,Abst,Wbst
		real(IDP), dimension(:,:), allocatable :: rmnb,sqgib,sqgb,grrb,grtb,gttb,bmodb,grrojb,grtojb,gttojb,jbgrrb,jbgrtb,jbgttb
		real(IDP), dimension(:,:), allocatable :: sqgieq,sqgeq,grreq,grteq,gtteq,grzeq,gtzeq,gzzeq,bmodeq,bsteq,jbgrreq,jbgrteq,jbgtteq,jbgrzeq,jbgtzeq, &
							  grrojeq,grtojeq,gttojeq,grzojeq,gtzojeq,gzzojeq,omdreq,omdteq,omdzeq, &
                                                          djrojeq,djtojeq,djzojeq,omdrprpeq,omdtprpeq,omdzprpeq, &
                                                          grrupeq,grtupeq,grzupeq,gttupeq,gtzupeq,gzzupeq,sb1,sb2,sb3,sb4,sb5,sb6,sb7, &
                                                          dbsjtojeq,dbsjzojeq,dbsjtbjeq,dgttreq,dgrrteq,dgrtteq,dgttteq,dgrrzeq, &
							  dgrtzeq,dgttzeq,dgrtpeq,dgttpeq,jsqeq,bsgrteq,bsgtteq,bsqeq,bsqgtteq, &
							  lplrreq,lplrteq,lplrzeq,lpltteq,lpltzeq,lplzzeq,lplreq,lplteq,lplzeq, &
							  eildreq,eildteq,eildzeq,eildrreq,eildrteq,eildrzeq,eildtteq,eildtzeq,eildzzeq, &
                                                          sqgdrojeq,sqgdthojeq,sqgdztojeq,sqgdthojbsteq,sqgdztojbsteq,sqgibmoditheq,sqgibmodizteq, &
                                                          testeq,testreq,testteq,testrreq,testtteq,testrteq, &
                                                          bmodieq,bdceq,brtdceq,bbbdceq,bbdceq 
		LOGICAL :: lasym
        complex(IDP) :: zetai,zetae,zetai2,zetai3,zetai4,zi,y0i,y1i, &
                   y2i,ddii,zetae2,zetae3,zetae4,ze,y0e,y1e,y2e,ddee,rei,sei, &
                   stfe,stfi,reii,sei1,sei2,sei3,zetaiinv,zetaeinv,cmplx1
        real(IDP) :: xnuion0, xkprl, abkprl, sgkprl, xsii, xsie, xsi2, xsi3, xsi4, &
                     ztr, zti, zir, zii, zer, zei, &
					 vther, vthir, tauie, xnuelc, xnuion	

		integer :: lwrt
		character(len=32) :: formatt='("r",127(a1,i4,"/",i4))'
		character(len=32) :: formatv='(1pe13.6,127(a1,1pe15.8))'
		character(len=1) :: t
		character(len=12) :: confil
					 	  
		interface
			subroutine mmblims
			end subroutine mmblims
			subroutine multb(f,g,typeg,h,typeh,c1,c2)
				use param
				implicit none
				integer :: typeg,typeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine multb
			subroutine grparb(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine grparb
			subroutine dbydthb(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydthb
			subroutine dbydrb(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydrb
			subroutine dbydztb(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydztb
			subroutine dbydtheq(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydtheq
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydzteq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzteq
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine eqsplns(yeq,yeqr,yfar,m,m2,choice,smoo)
				use param
				implicit none
				integer :: m,m2
				character(len=6) :: choice
				real(IDP) :: smoo
				real(IDP), dimension(0:) :: yeq,yeqr
				real(IDP), dimension(0:) :: yfar
			end subroutine eqsplns
			subroutine fitter(yeq,yeqr,yfar,ir,m)
				use param
				implicit none
				integer :: ir,m
				real(IDP), dimension(0:) :: yeq,yeqr,yfar
			end subroutine fitter	
 	        subroutine zzdisp(x,y,zzr,zzi)
 	            use param
 	            implicit none
                real(IDP) :: x,y,zzr,zzi
 	        end subroutine zzdisp			
		end interface

		write(0,'(" ====> Equilibria set up ... ")')			

!  		reads in high beta equilibrium data from file eq_name
!  		on unit 25.  then sets up equilibrium quantities
!  		required by FAR3d.
		
		bet0in=bet0

		allocate (sqg(0:mj,0:leqmax),sqgi(0:mj,0:leqmax),bmod(0:mj,0:leqmax),bst(0:mj,0:leqmax), &
			  grr(0:mj,0:leqmax),grt(0:mj,0:leqmax),gtt(0:mj,0:leqmax), &
                          grz(0:mj,0:leqmax),gtz(0:mj,0:leqmax),gzz(0:mj,0:leqmax), &
			  grroj(0:mj,0:leqmax),grtoj(0:mj,0:leqmax),gttoj(0:mj,0:leqmax), &
                          grzoj(0:mj,0:leqmax),gtzoj(0:mj,0:leqmax),gzzoj(0:mj,0:leqmax), &
                          grrup(0:mj,0:leqmax),grtup(0:mj,0:leqmax),grzup(0:mj,0:leqmax), &
                          gttup(0:mj,0:leqmax),gtzup(0:mj,0:leqmax), gzzup(0:mj,0:leqmax), &
			  jbgrr(0:mj,0:leqmax),jbgrt(0:mj,0:leqmax),jbgtt(0:mj,0:leqmax),jbgrz(0:mj,0:leqmax),jbgtz(0:mj,0:leqmax), &
			  djroj(0:mj,0:leqmax),djtoj(0:mj,0:leqmax),djzoj(0:mj,0:leqmax), &
			  jsq(0:mj,0:leqmax),omdr(0:mj,0:leqmax),omdt(0:mj,0:leqmax),omdz(0:mj,0:leqmax), &
			  dbsjtoj(0:mj,0:leqmax),dbsjzoj(0:mj,0:leqmax),dbsjtbj(0:mj,0:leqmax),dgttr(0:mj,0:leqmax), &
			  dgrrt(0:mj,0:leqmax),dgrtt(0:mj,0:leqmax),dgttt(0:mj,0:leqmax),dgrrz(0:mj,0:leqmax), &
			  dgrtz(0:mj,0:leqmax),dgttz(0:mj,0:leqmax),dgrtp(0:mj,0:leqmax),dgttp(0:mj,0:leqmax), &
			  bsgrt(0:mj,0:leqmax),bsgtt(0:mj,0:leqmax),bsq(0:mj,0:leqmax),bsqgtt(0:mj,0:leqmax), &
			  lplrr(0:mj,0:leqmax),lplrt(0:mj,0:leqmax),lplrz(0:mj,0:leqmax),lplr(0:mj,0:leqmax),lpltt(0:mj,0:leqmax), &
			  lpltz(0:mj,0:leqmax),lplt(0:mj,0:leqmax),lplz(0:mj,0:leqmax),lplzz(0:mj,0:leqmax), &
			  sqgdroj(0:mj,0:leqmax),sqgdthoj(0:mj,0:leqmax),sqgdztoj(0:mj,0:leqmax),sqgdthojbst(0:mj,0:leqmax), &
			  sqgdztojbst(0:mj,0:leqmax),sqgibmodith(0:mj,0:leqmax),sqgibmodizt(0:mj,0:leqmax),test(0:mj,0:leqmax), &
                          testr(0:mj,0:leqmax),testt(0:mj,0:leqmax),testrr(0:mj,0:leqmax),testtt(0:mj,0:leqmax),testrt(0:mj,0:leqmax))
		if (ieldamp_on == 1) allocate (eildr(0:mj,0:leqmax),eildt(0:mj,0:leqmax),eildz(0:mj,0:leqmax),eildrr(0:mj,0:leqmax), &
					       eildrt(0:mj,0:leqmax),eildrz(0:mj,0:leqmax),eildtt(0:mj,0:leqmax),eildtz(0:mj,0:leqmax), &
					       eildzz(0:mj,0:leqmax))	
                if(Trapped_on == 1) allocate (omdrprp(0:mj,0:leqmax),omdtprp(0:mj,0:leqmax),omdzprp(0:mj,0:leqmax))

!		Here we open the equilibria file and read the equilibria parameters		
 		open(unit=25,file=eq_name,status='old',convert='big_endian',form='unformatted')	
		
		read(25) nfp,lbmax,mjeqp,lasym
		if (nfp == 1) then
			ndevice(1)="    toka"
			ndevice(2)="mak     "
		else
			ndevice(1)="  stella"
			ndevice(2)="rator   "
		end if

		mjeq=mjeqp-1
		mjeqm1=mjeq-1
		lbm2=2*lbmax-1

		
		if (lasym) then
			allocate (mmb(lbm2),nnb(lbm2),llc(lbmax),lls(lbmax))
		else
			allocate (mmb(lbmax),nnb(lbmax))
		end if		

		allocate (rfar(0:mjeq),rbinv(0:mjeq),qfar(0:mjeq),pfar(0:mjeq),phip(0:mjeq),curfar(0:mjeq),ffar(0:mjeq), &
			  sfar1(0:mjeq),sfar2(0:mjeq),sfar3(0:mjeq),sfar4(0:mjeq))
		if (lasym) then
			allocate (rmnb(0:mjeq,0:lbm2),sqgb(0:mjeq,0:lbm2),sqgib(0:mjeq,0:lbm2),bmodb(0:mjeq,0:lbm2), &
				  grrb(0:mjeq,0:lbm2),grtb(0:mjeq,0:lbm2),gttb(0:mjeq,0:lbm2), &
				  grrojb(0:mjeq,0:lbm2),grtojb(0:mjeq,0:lbm2),gttojb(0:mjeq,0:lbm2), &
				  jbgrrb(0:mjeq,0:lbm2),jbgrtb(0:mjeq,0:lbm2),jbgttb(0:mjeq,0:lbm2))
		else
			allocate (rmnb(0:mjeq,0:lbmax),sqgb(0:mjeq,0:lbmax),sqgib(0:mjeq,0:lbmax),bmodb(0:mjeq,0:lbmax), &
				  grrb(0:mjeq,0:lbmax),grtb(0:mjeq,0:lbmax),gttb(0:mjeq,0:lbmax), &
				  grrojb(0:mjeq,0:lbmax),grtojb(0:mjeq,0:lbmax),gttojb(0:mjeq,0:lbmax), &
				  jbgrrb(0:mjeq,0:lbmax),jbgrtb(0:mjeq,0:lbmax),jbgttb(0:mjeq,0:lbmax))
		end if

		rmnb(:,0)=0.0_IDP
		sqgib(:,0)=0.0_IDP
		sqgb(:,0)=0.0_IDP
		grrb(:,0)=0.0_IDP
		grtb(:,0)=0.0_IDP
		gttb(:,0)=0.0_IDP
		grrojb(:,0)=0.0_IDP
		grtojb(:,0)=0.0_IDP
		gttojb(:,0)=0.0_IDP
 	
		read(25) (mmb(l),nnb(l),l=1,lbmax)
		read(25) (phip(j),qfar(j),curfar(j),ffar(j),pfar(j),j=1,mjeq)
		do j=1,mjeq
			read(25) (rmnb(j,l),dummy,dummy,bmodb(j,l),sqgb(j,l),sqgib(j,l),l=1,lbmax)
		end do
		do j=1,mjeq
			read(25) (grrb(j,l),grtb(j,l),gttb(j,l),grrojb(j,l),grtojb(j,l),gttojb(j,l),jbgrrb(j,l),jbgrtb(j,l),jbgttb(j,l), &
				  l=1,lbmax)
		end do

		lb0=0
		do l=1,lbmax
			if (mmb(l) == 0 .and. nnb(l) == 0) exit
		end do
		lb0=l
		if (lb0 == 0 .or. lb0 > lbmax) then
			write (6,'("  vmec: lb0=0")')
			stop
		end if

		if (lasym) then
			l1=lbmax
			do l=1,lbmax
				if (l == lb0) cycle
				l1=l1+1
				llc(l)=l1
				lls(l)=l1
			end do
			llc(lb0)=lb0
			lls(lb0)=0

			do j=1,mjeq
				read(25) (rmnb(j,lls(l)),dummy,dummy,bmodb(j,lls(l)),sqgb(j,lls(l)),sqgib(j,lls(l)),l=1,lbmax)
			end do
			do j=1,mjeq
				read(25) (grrb(j,lls(l)),grtb(j,llc(l)),gttb(j,lls(l)),grrojb(j,lls(l)),grtojb(j,llc(l)),gttojb(j,lls(l)), &
					  jbgrrb(j,lls(l)),jbgrtb(j,llc(l)),jbgttb(j,lls(l)),l=1,lbmax)
			end do
		end if

		do l=1,lbmax
			if (nnb(l) == 0 .and. mmb(l) < 0) then
				mmb(l)=-mmb(l)
				grtb(:,l)=-grtb(:,l)
				grtojb(:,l)=-grtojb(:,l)
				jbgrtb(:,l)=-jbgrtb(:,l)
				if (lasym) then
					rmnb(:,lls(l))=-rmnb(:,lls(l))
					bmodb(:,lls(l))=-bmodb(:,lls(l))
					sqgb(:,lls(l))=-sqgb(:,lls(l))
					sqgib(:,lls(l))=-sqgib(:,lls(l))
					grrb(:,lls(l))=-grrb(:,lls(l))
					gttb(:,lls(l))=-gttb(:,lls(l))
					grrojb(:,lls(l))=-grrojb(:,lls(l))
					gttojb(:,lls(l))=-gttojb(:,lls(l))
					jbgrrb(:,lls(l))=-jbgrrb(:,lls(l))
					jbgttb(:,lls(l))=-jbgttb(:,lls(l))
				end if
			end if
			if (nnb(l) < 0) then
				mmb(l)=-mmb(l)
				nnb(l)=-nnb(l)
				grtb(:,l)=-grtb(:,l)
				grtojb(:,l)=-grtojb(:,l)
				jbgrtb(:,l)=-jbgrtb(:,l)
				if (lasym) then
					rmnb(:,lls(l))=-rmnb(:,lls(l))
					bmodb(:,lls(l))=-bmodb(:,lls(l))
					sqgb(:,lls(l))=-sqgb(:,lls(l))
					sqgib(:,lls(l))=-sqgib(:,lls(l))
					grrb(:,lls(l))=-grrb(:,lls(l))
					gttb(:,lls(l))=-gttb(:,lls(l))
					grrojb(:,lls(l))=-grrojb(:,lls(l))
					gttojb(:,lls(l))=-gttojb(:,lls(l))
					jbgrrb(:,lls(l))=-jbgrrb(:,lls(l))
					jbgttb(:,lls(l))=-jbgttb(:,lls(l))
				end if
			end if
		end do

		if (lasym) then
			l1=lbmax
			do l=1,lbmax
				if (l == lb0) cycle
				l1=l1+1
				mmb(l1)=-mmb(l)
				nnb(l1)=-nnb(l)
			end do
			lbmax=lbm2
		end if
	
!	normalizations

		twopi=8.0_IDP*atan(1.0_IDP)
		ffar=-ffar
		pfar=2.e-7*twopi*pfar

		bigrn=1.5*rmnb(1,lb0)-0.5*rmnb(2,lb0)
		bmodn=1.5*bmodb(1,lb0)-0.5*bmodb(2,lb0)

		twophip=2.0*phip(mjeq)
		curfar=sign(1.0_IDP,twophip)*curfar/(bigrn*bmodn)
		ffar=sign(1.0_IDP,twophip)*ffar/(bigrn*bmodn)
		bmodb=bmodb/bmodn
		sqgb=-bmodn*sqgb/bigrn
		sqgib=-bigrn*sqgib/bmodn
		grrb=bmodn*grrb/abs(twophip)
		grtb=bmodn*grtb/abs(twophip)
		gttb=bmodn*gttb/abs(twophip)
		grrojb=bigrn*grrojb/abs(twophip)
		grtojb=bigrn*grtojb/abs(twophip)
		gttojb=bigrn*gttojb/abs(twophip)
		jbgrrb=bmodn*bmodn*jbgrrb/(bigrn*abs(twophip))
		jbgrtb=bmodn*bmodn*jbgrtb/(bigrn*abs(twophip))
		jbgttb=bmodn*bmodn*jbgttb/(bigrn*abs(twophip))

		eps=sqrt(abs(twophip))/(bigrn*sqrt(bmodn))

		write(6,'(/"R_0 =",1pe13.6," B0 =",1pe13.6," eps=",1pe13.6)') bigrn,bmodn,eps

		do j=0,mjeq
			rfar(j)=sqrt(1.0_IDP*j/mjeq)
		end do

		rbinv(0)=0.0_IDP
		do j=1,mjeq
			rbinv(j)=1.0_IDP/rfar(j)
		end do

		nnb=nnb/nfp

		mmaxb=maxval(mmb(1:lbmax))
		mminb=minval(mmb(1:lbmax))
		nmaxb=maxval(nnb(1:lbmax))
		nminb=minval(nnb(1:lbmax))
		mmaxxb=max(mmaxb,abs(mminb))

		allocate (llb(mminb:mmaxb,nminb:nmaxb))

		llb=0
		do l=1,lbmax
			llb(mmb(l),nnb(l))=l
		end do

		mband=0
		do m=0,mmaxb
			if (llb(m,0) /= 0) mband=mband+1
		end do
		mxmbandb=2*mband-1
		do n=1,nmaxb
			mband=0
			do m=mminb,mmaxb
				if (llb(m,n) /= 0) mband=mband+1
			end do
			if (mband > mxmbandb) mxmbandb=mband
		end do
		call mmblims

		nnb=nnb*nfp

		sqgib(0,:)=0.0_IDP
		sqgb(0,:)=0.0_IDP
		grrb(0,:)=0.0_IDP
		grtb(0,:)=0.0_IDP
		gttb(0,:)=0.0_IDP
		bmodb(0,:)=0.0_IDP
		grrojb(0,:)=0.0_IDP
		grtojb(0,:)=0.0_IDP
		gttojb(0,:)=0.0_IDP
		jbgrrb(0,:)=0.0_IDP
		jbgrtb(0,:)=0.0_IDP
		jbgttb(0,:)=0.0_IDP
		qfar(0)=1.5*qfar(1)-0.5*qfar(2)
		pfar(0)=1.5*pfar(1)-0.5*pfar(2)
		ffar(0)=1.5*ffar(1)-0.5*ffar(2)
		do j=1,mjeq
			curfar(j)=curfar(j)*mjeq/(j-0.5_IDP)
		end do
		curfar(0)=1.5*curfar(1)-0.5*curfar(2)

		do l=1,lbmax
			if (mmb(l) == 0) then
				sqgib(0,l)=2.*sqgib(1,l)-sqgib(2,l)
				sqgb(0,l)=1.5*sqgb(1,l)-0.5*sqgb(2,l)
				bmodb(0,l)=1.5*bmodb(1,l)-0.5*bmodb(2,l)
			endif
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				grrb(0,l)=2.*grrb(1,l)-grrb(2,l)
				grtb(0,l)=2.*grtb(1,l)-grtb(2,l)
				gttb(0,l)=2.*gttb(1,l)-gttb(2,l)
				grrojb(0,l)=2.*grrojb(1,l)-grrojb(2,l)
				grtojb(0,l)=2.*grtojb(1,l)-grtojb(2,l)
				gttojb(0,l)=2.*gttojb(1,l)-gttojb(2,l)
				jbgrrb(0,l)=2.*jbgrrb(1,l)-jbgrrb(2,l)
				jbgrtb(0,l)=2.*jbgrtb(1,l)-jbgrtb(2,l)
				jbgttb(0,l)=2.*jbgttb(1,l)-jbgttb(2,l)
			endif
		end do

		bet0=2.*pfar(0)/bmodn**2
		if (bet0 /= 0.0_IDP) then
			pror=pfar(0)
			pfar=pfar/pror
		end if

		if (bet0in /= 1.0_IDP) bet0=bet0in

!		Change the thermal beta

		bet0=bet0*betath_factor

!		Here we use splines to transfer the equilibria parameters to code parameters
		
		dummy=1.5*qfar(mjeq)-0.5*qfar(mjeqm1)
		do j=1,mjeqm1
			qfar(j)=0.5*(qfar(j+1)+qfar(j))
		end do
		qfar(mjeq)=dummy
		call eqsplns(qqinv,sd1,qfar,0,0,"spline",1.e-07_IDP)
		dummy=1.5*pfar(mjeq)-0.5*pfar(mjeqm1)
		do j=1,mjeqm1
			pfar(j)=0.5*(pfar(j+1)+pfar(j))
		end do
		pfar(mjeq)=dummy
		call eqsplns(preq,sd1,pfar,0,0,"spline",1.e-07_IDP)
		dummy=1.5*ffar(mjeq)-0.5*ffar(mjeqm1)
		do j=1,mjeqm1
			ffar(j)=0.5*(ffar(j+1)+ffar(j))
		end do
		ffar(mjeq)=dummy
		call eqsplns(feq,sd1,ffar,0,0,"spline",1.e-07_IDP)
		dummy=1.5*curfar(mjeq)-0.5*curfar(mjeqm1)
		do j=1,mjeqm1
			curfar(j)=0.5*(curfar(j+1)+curfar(j))
		end do
		curfar(mjeq)=dummy
		call eqsplns(cureq,sd1,curfar,0,0,"spline",1.e-07_IDP)
		do j=1,mj
			cureq(j)=cureq(j)*r(j)**2
		end do
		cureq(0)=0.0_IDP

!		Displacment included to the safety factor / iota

		qq=(1.0_IDP/qqinv) + deltaq
		qqinv = (1.0_IDP/qq) + deltaiota	

!		Introduce the external profiles	

	           if(ext_prof .eq. 1) then	
                         call ae_profiles
		end if	

		allocate (sqgeq(0:mj,0:lbmax),sqgieq(0:mj,0:lbmax),bmodeq(0:mj,0:lbmax),bsteq(0:mj,0:lbmax), &
			  grreq(0:mj,0:lbmax),grteq(0:mj,0:lbmax),gtteq(0:mj,0:lbmax), &
                          grzeq(0:mj,0:lbmax),gtzeq(0:mj,0:lbmax),gzzeq(0:mj,0:lbmax), &
			  grrupeq(0:mj,0:lbmax),grtupeq(0:mj,0:lbmax),gttupeq(0:mj,0:lbmax), &
                          grzupeq(0:mj,0:lbmax),gtzupeq(0:mj,0:lbmax),gzzupeq(0:mj,0:lbmax), &
			  grrojeq(0:mj,0:lbmax),grtojeq(0:mj,0:lbmax),gttojeq(0:mj,0:lbmax), &
			  grzojeq(0:mj,0:lbmax),gtzojeq(0:mj,0:lbmax),gzzojeq(0:mj,0:lbmax), &
			  jbgrreq(0:mj,0:lbmax),jbgrteq(0:mj,0:lbmax),jbgtteq(0:mj,0:lbmax),jbgrzeq(0:mj,0:lbmax),jbgtzeq(0:mj,0:lbmax), &
			  djrojeq(0:mj,0:lbmax),djtojeq(0:mj,0:lbmax),djzojeq(0:mj,0:lbmax), &
			  jsqeq(0:mj,0:lbmax),omdreq(0:mj,0:lbmax),omdteq(0:mj,0:lbmax),omdzeq(0:mj,0:lbmax), &
			  dbsjtojeq(0:mj,0:lbmax),dbsjzojeq(0:mj,0:lbmax),dbsjtbjeq(0:mj,0:lbmax),dgttreq(0:mj,0:lbmax), &
			  dgrrteq(0:mj,0:lbmax),dgrtteq(0:mj,0:lbmax),dgttteq(0:mj,0:lbmax),dgrrzeq(0:mj,0:lbmax), &
			  dgrtzeq(0:mj,0:lbmax),dgttzeq(0:mj,0:lbmax),dgrtpeq(0:mj,0:lbmax),dgttpeq(0:mj,0:lbmax), &
			  bsgrteq(0:mj,0:lbmax),bsgtteq(0:mj,0:lbmax),bsqeq(0:mj,0:lbmax),bsqgtteq(0:mj,0:lbmax), &
			  lplrreq(0:mj,0:lbmax),lplrteq(0:mj,0:lbmax),lplrzeq(0:mj,0:lbmax),lplreq(0:mj,0:lbmax),lpltteq(0:mj,0:lbmax), &
			  lpltzeq(0:mj,0:lbmax),lplteq(0:mj,0:lbmax),lplzeq(0:mj,0:lbmax),lplzzeq(0:mj,0:lbmax), &
			  sb1(0:mj,0:lbmax),sb2(0:mj,0:lbmax),sb3(0:mj,0:lbmax),sb4(0:mj,0:lbmax),sb5(0:mj,0:lbmax),sb6(0:mj,0:lbmax),sb7(0:mj,0:lbmax), &
			  sqgdrojeq(0:mj,0:lbmax),sqgdthojeq(0:mj,0:lbmax),sqgdztojeq(0:mj,0:lbmax),sqgdthojbsteq(0:mj,0:lbmax), &	
			  sqgdztojbsteq(0:mj,0:lbmax),sqgibmoditheq(0:mj,0:lbmax),sqgibmodizteq(0:mj,0:lbmax), &
                          testeq(0:mj,0:lbmax),testreq(0:mj,0:lbmax),testteq(0:mj,0:lbmax),testrreq(0:mj,0:lbmax), &
                          testtteq(0:mj,0:lbmax),testrteq(0:mj,0:lbmax), &
                          bmodieq(0:mj,0:lbmax),bdceq(0:mj,0:lbmax),brtdceq(0:mj,0:lbmax),bbbdceq(0:mj,0:lbmax),bbdceq(0:mj,0:lbmax)) 
                        if(Trapped_on == 1) allocate (omdrprpeq(0:mj,0:lbmax),omdtprpeq(0:mj,0:lbmax),omdzprpeq(0:mj,0:lbmax))

		sqgieq=0.0_IDP
		sqgeq=0.0_IDP
		bsteq=0.0_IDP
		grreq=0.0_IDP
		grteq=0.0_IDP
		grzeq=0.0_IDP
		gtteq=0.0_IDP
		gtzeq=0.0_IDP
		gzzeq=0.0_IDP
		bmodeq=0.0_IDP
		grrojeq=0.0_IDP
		grtojeq=0.0_IDP
		grzojeq=0.0_IDP
		gttojeq=0.0_IDP
		gtzojeq=0.0_IDP
		gzzojeq=0.0_IDP
		jbgrreq=0.0_IDP
		jbgrteq=0.0_IDP
		jbgtteq=0.0_IDP	
		jbgrzeq=0.0_IDP
		jbgtzeq=0.0_IDP		

                bmodieq=0.0_IDP
                bdceq=0.0_IDP
                brtdceq=0.0_IDP
                bbbdceq=0.0_IDP
                bbdceq=0.0_IDP

                lplrreq=0.0_IDP
                lplrteq=0.0_IDP
                lplrzeq=0.0_IDP
                lpltzeq=0.0_IDP
                lpltteq=0.0_IDP
                lplzzeq=0.0_IDP
                lplreq=0.0_IDP
                lplteq=0.0_IDP
                lplzeq=0.0_IDP
  
		ir=max(15,mjeq/10)
		do l=1,lbmax
			m=min(abs(mmb(l)),30)
			m2=m
			call eqsplns(sqgieq(:,l),sd1,sqgib(:,l),m,m2,"icsscu",1.e-07_IDP)
			dummy=1.5*sqgb(mjeq,l)-0.5*sqgb(mjeqm1,l)
			do j=1,mjeqm1
				sqgb(j,l)=0.5*(sqgb(j+1,l)+sqgb(j,l))
			end do
			sqgb(mjeq,l)=dummy
			call eqsplns(sqgeq(:,l),sd1,sqgb(:,l),m,m2,"icsscu",1.e-07_IDP)
			dummy=1.5*bmodb(mjeq,l)-0.5*bmodb(mjeqm1,l)
			do j=1,mjeqm1
				bmodb(j,l)=0.5*(bmodb(j+1,l)+bmodb(j,l))
			end do
			bmodb(mjeq,l)=dummy
			call eqsplns(bmodeq(:,l),sd1,bmodb(:,l),m,m2,"icsscu",1.e-07_IDP)
			mm1=abs(m-1)
			mm2=abs(mm1-1)
			mm12=mm1
			mm22=mm2
			mm2p1=mm2+1
			mm2p12=mm2p1
			mp2=m+2
			mp22=mp2
			isg=sign(1,mmb(l))

			sfar1=grrb(:,l)+isg*grtb(:,l)
			sfar2=grrb(:,l)-isg*grtb(:,l)
			sfar3=gttb(:,l)-isg*grtb(:,l)
			sfar4=sfar3-sfar1
			call fitter(sd1,sd6,sfar1,ir,m)
			call fitter(sd2,sd6,sfar2,ir,mm2)
			call fitter(sd4,sd6,sfar4,ir,mp2)
			grreq(:,l)=(sd1+sd2)/2.
			grteq(:,l)=isg*(sd1-sd2)/2.
			gtteq(:,l)=(2.*sd4+3.*sd1-sd2)/2.

			sfar1=grrojb(:,l)+isg*grtojb(:,l)
			sfar2=grrojb(:,l)-isg*grtojb(:,l)
			sfar3=gttojb(:,l)-isg*grtojb(:,l)
			sfar4=sfar3-sfar1
			call fitter(sd1,sd6,sfar1,ir,m)
			call fitter(sd2,sd6,sfar2,ir,mm2)
			call fitter(sd4,sd6,sfar4,ir,mp2)
			grrojeq(:,l)=(sd1+sd2)/2.
			grtojeq(:,l)=isg*(sd1-sd2)/2.
			gttojeq(:,l)=(2.*sd4+3.*sd1-sd2)/2.

			sfar1=jbgrrb(:,l)+isg*jbgrtb(:,l)
			sfar2=jbgrrb(:,l)-isg*jbgrtb(:,l)
			sfar3=jbgttb(:,l)-isg*jbgrtb(:,l)
			sfar4=sfar3-sfar1
			call fitter(sd1,sd6,sfar1,ir,m)
			call fitter(sd2,sd6,sfar2,ir,mm2)
			call fitter(sd4,sd6,sfar4,ir,mp2)
			jbgrreq(:,l)=(sd1+sd2)/2.
			jbgrteq(:,l)=isg*(sd1-sd2)/2.
			jbgtteq(:,l)=(2.*sd4+3.*sd1-sd2)/2.
			
		end do

!       Edge extrapolation

                if(Edge_on .eq. 1) then      
		  do j=edge_p,mj

		      sqgieq(j,:)=sqgieq(j-2,:) + (r(j)-r(j-2))*(sqgieq(j-1,:)-sqgieq(j-2,:))/(r(j-1)-r(j-2))
		      sqgeq(j,:)=sqgeq(j-2,:) + (r(j)-r(j-2))*(sqgeq(j-1,:)-sqgeq(j-2,:))/(r(j-1)-r(j-2))
		      grreq(j,:)=grreq(j-2,:) + (r(j)-r(j-2))*(grreq(j-1,:)-grreq(j-2,:))/(r(j-1)-r(j-2))
		      grteq(j,:)=grteq(j-2,:) + (r(j)-r(j-2))*(grteq(j-1,:)-grteq(j-2,:))/(r(j-1)-r(j-2))
		      gtteq(j,:)=gtteq(j-2,:) + (r(j)-r(j-2))*(gtteq(j-1,:)-gtteq(j-2,:))/(r(j-1)-r(j-2))
		      bmodeq(j,:)=bmodeq(j-2,:) + (r(j)-r(j-2))*(bmodeq(j-1,:)-bmodeq(j-2,:))/(r(j-1)-r(j-2))
		      grrojeq(j,:)=grrojeq(j-2,:) + (r(j)-r(j-2))*(grrojeq(j-1,:)-grrojeq(j-2,:))/(r(j-1)-r(j-2))
		      grtojeq(j,:)=grtojeq(j-2,:) + (r(j)-r(j-2))*(grtojeq(j-1,:)-grtojeq(j-2,:))/(r(j-1)-r(j-2))
		      gttojeq(j,:)=gttojeq(j-2,:) + (r(j)-r(j-2))*(gttojeq(j-1,:)-gttojeq(j-2,:))/(r(j-1)-r(j-2))
		      jbgrreq(j,:)=jbgrreq(j-2,:) + (r(j)-r(j-2))*(jbgrreq(j-1,:)-jbgrreq(j-2,:))/(r(j-1)-r(j-2))
		      jbgrteq(j,:)=jbgrteq(j-2,:) + (r(j)-r(j-2))*(jbgrteq(j-1,:)-jbgrteq(j-2,:))/(r(j-1)-r(j-2))
		      jbgtteq(j,:)=jbgtteq(j-2,:) + (r(j)-r(j-2))*(jbgtteq(j-1,:)-jbgtteq(j-2,:))/(r(j-1)-r(j-2))

		      preq(j)=preq(j-2) + (r(j)-r(j-2))*(preq(j-1)-preq(j-2))/(r(j-1)-r(j-2))
		      feq(j)=feq(j-2) + (r(j)-r(j-2))*(feq(j-1)-feq(j-2))/(r(j-1)-r(j-2))
		      cureq(j)=cureq(j-2) + (r(j)-r(j-2))*(cureq(j-1)-cureq(j-2))/(r(j-1)-r(j-2))

                  end do

		  t=char(9)
                  lwrt=abs(lplots)
		  write(confil,'("sqgieq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
		  	write(92,formatv) r(j),(t,sqgieq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("sqgeq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,sqgeq(j,l),l=1,lwrt)
		  end do
		  close(92)

		  write(confil,'("grreq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,grreq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("grteq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,grteq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("gtteq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,gtteq(j,l),l=1,lwrt)
		  end do
		  close(92)

		  write(confil,'("bmodeq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,bmodeq(j,l),l=1,lwrt)
		  end do
		  close(92)

		  write(confil,'("grrojeq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,grrojeq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("grtojeq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,grtojeq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("gttojeq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,gttojeq(j,l),l=1,lwrt)
		  end do
		  close(92)

		  write(confil,'("jbgrreq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,jbgrreq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("jbgrteq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,jbgrteq(j,l),l=1,lwrt)
		  end do
		  close(92)
		  write(confil,'("jbgtteq_",2a2)') numrun(1),numrun(2)
		  open(unit=92,file=confil,recl=2048)
		  write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		  do j=0,mj
			write(92,formatv) r(j),(t,jbgtteq(j,l),l=1,lwrt)
		  end do
		  close(92)

                end if


!	Calculation of bst		
!	bst is modified for resonant modes

		! allocate (lbst(lbmax))
		! allocate (rsb(lbmax), Abst(lbmax), Wbst(lbmax))

		! bsteq=0
		! one = 1.0_IDP
		! lbst=0
		! Abst=0
		! qmin=minval(qqinv(0:mj))
		! qmax=maxval(qqinv(0:mj))
		! call dbydr0(sd1,preq,0.0_IDP,-1.0_IDP,0)
		! call dbydr0(sd3,qqinv,0.0_IDP,1.0_IDP,0)
		! lp=0
		! do l=1,lbmax
			! if (mmb(l) == 0 .or. nnb(l) == 0) cycle
			! p=one*nnb(l)/mmb(l)
			! if ((p-qmin)*(p-qmax) >= 0.0_IDP) cycle
			! lp=lp+1
			! lbst(lp)=l
			! do j=1,mj
				! if ((qqinv(j-1)-p)*(qqinv(j)-p) <= 0.) exit
			! end do
			! jp=j
			! rsb(lp)=r(jp-1)+(r(jp)-r(jp-1))*(qqinv(jp-1)-p)/(qqinv(jp-1)-qqinv(jp))
			! if (jp == 1) then
				! Abst(lp) = (preq(0)-preq(1))/r(1)
				! iotap = (qqinv(1)-qqinv(0))/r(1)
			! else if (jp < mj) then
				! Abst(lp) = ((r(jp)-rsb(lp))*sd1(jp-1)+(rsb(lp)-r(jp-1))*sd1(jp))/(r(jp)-r(jp-1))
				! iotap = ((r(jp)-rsb(lp))*sd3(jp-1)+(rsb(lp)-r(jp-1))*sd3(jp))/(r(jp)-r(jp-1))
			! else
				! Abst(lp) = (preq(mjm1)-preq(mj))/(r(mj)-r(mjm1))
				! iotap = (qqinv(mj)-qqinv(mjm1))/(r(mj)-r(mjm1))
			! end if
			! Wbst(lp) = abs(p/(mmb(l)*iotap))
		! end do
		! nbst=lp
		! open (unit=6,file="farprt",status="old",POSITION="APPEND")
		! write(6,'(/"bst is modified for",i5," resonant modes")') nbst
		! do l=1,nbst
			! write(6,'("m =",i5," n =",i5," rs =",1pe13.6," A =",1pe13.6," W =",1pe13.6)') mmb(lbst(l)), &
				! nnb(lbst(l)),rsb(l),Abst(l),Wbst(l)
		! end do
		! close(6)

		! sd2=rinv*sd1/2.0_IDP
		! sd2(0)=(preq(0)-preq(1))/r(1)**2
		! do l=1,nbst
			! do j=0,mj
				! pprime=sd2(j)
				! iota=qqinv(j)
				! gc=sqgeq(j,lbst(l))
				! xl=r(j)-rsb(l)
				! pprimel=pprime-0.5*rinv(j)*Abst(l)*(1-(xl/Wbst(l))**2)*exp(-0.5*(xl/Wbst(l))**2)
				! bsteq(j,lbst(l))=bet0*pprimel*gc/(mmb(lbst(l))*iota-nnb(lbst(l)))
			! end do
		! end do

		! do l=1,lbmax
			! if(mmb(l) == 0 .and. nnb(l) == 0) cycle
			! p=0.0
			! if (mmb(l) /= 0) p=one*nnb(l)/mmb(l)
			! if ((p-qmin)*(p-qmax) < 0.0_IDP) cycle
			! bsteq(:,l)=bet0*sd2*sqgeq(:,l)/(mmb(l)*qqinv-nnb(l))
		! end do
		
		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
		sd2=-rinv*sd1/2.0_IDP
		sd2(0)=(preq(0)-preq(1))/r(1)**2
		do l=1,lbmax
			if(mmb(l) == 0 .and. nnb(l) == 0) cycle
			bsteq(:,l)=bet0*sd2*sqgeq(:,l)/(mmb(l)*qqinv-nnb(l))
		end do	

!       Definition of contravariant metric tensors	

		call multb(sb1,sqgeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			grzeq(:,l)=eps*r*qqinv*grteq(:,l)+r*sb1(:,l)/eps
			gtzeq(:,l)=eps*r*qqinv*gtteq(:,l)+cureq*sqgeq(:,l)*rinv/eps			
			gzzeq(:,l)=(feq+qqinv*cureq)*sqgeq(:,l)+ (eps*r*qqinv)*(eps*r*qqinv)*gtteq(:,l)			
			grrupeq(:,l)=(feq+qqinv*cureq)*gttojeq(:,l)-(cureq*rinv/eps)*(cureq*rinv/eps)
		end do

		grzeq(0,:)=0.0_IDP
		gtzeq(0,:)=0.0_IDP
		gzzeq(0,:)=0.0_IDP
		grrupeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				grzeq(0,l)=2.*grzeq(1,l)-grzeq(2,l)
				gtzeq(0,l)=2.*gtzeq(1,l)-gtzeq(2,l)
				gzzeq(0,l)=2.*gzzeq(1,l)-gzzeq(2,l)
				grrupeq(0,l)=2.*grrupeq(1,l)-grrupeq(2,l)
			endif
		end do

		call multb(grzojeq,grzeq,-1,sqgieq,1,0.0_IDP,1.0_IDP)	
		call multb(gtzojeq,gtzeq,1,sqgieq,1,0.0_IDP,1.0_IDP)

		grzojeq(0,:)=0.0_IDP
		gtzojeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				grzojeq(0,l)=2.*grzojeq(1,l)-grzojeq(2,l)
				gtzojeq(0,l)=2.*gtzojeq(1,l)-gtzojeq(2,l)
			endif
		end do
			
		call multb(sb1,gttojeq,1,bsteq,-1,0.0_IDP,1.0_IDP)		
		do l=1,lbmax		
			grtupeq(:,l)=-feq*grtojeq(:,l)+r*r*qqinv*sb1(:,l)+cureq*bsteq(:,l)/(eps*eps)		
			grzupeq(:,l)=((cureq*grtojeq(:,l)*rinv)-r*sb1(:,l))/eps
		end do

		grtupeq(0,:)=0.0_IDP
		grzupeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				grtupeq(0,l)=2.*grtupeq(1,l)-grtupeq(2,l)
				grzupeq(0,l)=2.*grzupeq(1,l)-grzupeq(2,l)
			endif
		end do

		call multb(sb1,grtojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb2,bsteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)	
		call multb(sb3,sb2,1,gttojeq,1,0.0_IDP,1.0_IDP)	
		call multb(sb4,grzojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)	
		do l=1,lbmax							
			gzzupeq(:,l)=(1/(feq-qqinv*cureq))*(sqgieq(:,l)+(cureq*cureq*grrojeq(:,l)*rinv*rinv/(eps*eps))-2*cureq*sb1(:,l)/(eps*eps)+r*r*sb3(:,l)/(eps*eps))
		end do	

		gzzupeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				gzzupeq(0,l)=2.*gzzupeq(1,l)-gzzupeq(2,l)
			endif
		end do

		do l=1,lbmax							
                        gtzupeq(:,l)=-(cureq*grrojeq(:,l)*rinv -r*sb1(:,l))/eps -eps*r*qqinv*gzzupeq(:,l)
			gttupeq(:,l)=(feq+qqinv*cureq)*grrojeq(:,l)-2*r*r*qqinv*sb1(:,l)-(r*r/(eps*eps))*sb2(:,l)+(eps*r*qqinv)*(eps*r*qqinv)*gzzupeq(:,l)
		end do	

		gtzupeq(0,:)=0.0_IDP
		gttupeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				gtzupeq(0,l)=2.*gtzupeq(1,l)-gtzupeq(2,l)
				gttupeq(0,l)=2.*gttupeq(1,l)-gttupeq(2,l)
			endif
		end do

		call multb(jbgrzeq,grzeq,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(jbgtzeq,gtzeq,1,sqgeq,1,0.0_IDP,1.0_IDP)

		jbgrzeq(0,:)=0.0_IDP
		jbgtzeq(0,:)=0.0_IDP
		do l=1,lbmax
			if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
				jbgrzeq(0,l)=2.*jbgrzeq(1,l)-jbgrzeq(2,l)
				jbgtzeq(0,l)=2.*jbgtzeq(1,l)-jbgtzeq(2,l)
			endif
		end do


!               Boundary condition variables

		call multb(testeq,grreq,1,grtupeq,-1,0.0_IDP,1.0_IDP)

		call multb(sb1,sqgeq,1,gtteq,1,0.0_IDP,1.0_IDP)
		call dbydrb(sb3,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb3,1,gtteq,1,0.0_IDP,1.0_IDP)
		call dbydrb(sb4,gtteq,0.0_IDP,1.0_IDP,0)
		call multb(sb3,sb4,1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb5,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb4,sb5,-1,grteq,-1,0.0_IDP,1.0_IDP)
		call dbydthb(sb6,grteq,-1,0.0_IDP,1.0_IDP,0)
		call multb(sb5,sb6,1,sqgeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			testreq(:,l)=(sb1(:,l)-sb4(:,l)-sb5(:,l))*rinv + sb2(:,l) + sb3(:,l)
		end do

		call multb(sb1,sqgeq,1,grteq,-1,0.0_IDP,1.0_IDP)
		call dbydrb(sb3,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb3,1,grteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb3,sqgeq,1,grteq,-1,0.0_IDP,1.0_IDP)
		call dbydrb(sb5,grteq,0.0_IDP,1.0_IDP,0)
		call multb(sb4,sqgeq,1,sb5,-1,0.0_IDP,1.0_IDP)
		call dbydthb(sb6,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb5,sb6,-1,grreq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb7,grreq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb6,sb7,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			testteq(:,l)=(sb1(:,l)-sb3(:,l)-sb5(:,l)-sb6(:,l))*rinv + sb2(:,l) + sb4(:,l)
		end do

		call multb(testrreq,sqgeq,1,gtteq,1,0.0_IDP,1.0_IDP)

		call multb(testtteq,sqgeq,1,grreq,1,0.0_IDP,1.0_IDP)

		call multb(testrteq,sqgeq,1,grteq,-1,0.0_IDP,1.0_IDP)

!		Definition of new variables used in linstart
		
		sb1=sqgeq
		call multb(jsqeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
		
		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(djrojeq,sqgieq,1,sb1,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(djtojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(djzojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)

!		Calculation of Omega d operator components		

		sd1=feq-qqinv*cureq	
		call dbydztb(sb1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(omdteq,bsteq,-1,sb1,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb2(:,l)=rinv*cureq*sb1(:,l)
		end do
		call dbydthb(sb1,sqgieq,1,0.0_IDP,1.0_IDP,0)
		call multb(omdzeq,bsteq,-1,sb1,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb2(:,l)=sb2(:,l)-feq*sb1(:,l)
		end do
		call multb(omdreq,sqgeq,1,sb2,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			omdreq(:,l)=omdreq(:,l)/(2.0_IDP*sd1)
		end do

                if(Trapped_on .eq. 1) then

		    do l=1,lbmax
			    if (mm(l) .eq. 0)  then
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound)
			    else
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound*mm(l))
			    endif
		     end do
		     call multb(omdrprpeq,omdreq,-1,sb1,1,0.0_IDP,1.0_IDP)

                end if

		do l=1,lbmax
			sb1(:,l)=sd1*sqgieq(:,l)
		end do
		call dbydrb(sb2,sb1,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			sb1(:,l)=feq*sb2(:,l)-sd1*r*omdteq(:,l)
		end do
		call multb(omdteq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			omdteq(:,l)=omdteq(:,l)/(2.0_IDP*sd1*sd1)
		end do

                if(Trapped_on .eq. 1) then

		    do l=1,lbmax
			    if (mm(l) .eq. 0)  then
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound)
			    else
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound*mm(l))
			    endif
		    end do
		    call multb(omdtprpeq,omdteq,1,sb1,1,0.0_IDP,1.0_IDP)

                end if

		do l=1,lbmax
			sb1(:,l)=sd1*r*omdzeq(:,l)-rinv*cureq*sb2(:,l)
		end do
		call multb(omdzeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			omdzeq(:,l)=omdzeq(:,l)/(2.0_IDP*sd1*sd1)
		end do

                if(Trapped_on .eq. 1) then

		    do l=1,lbmax
			    if (mm(l) .eq. 0)  then
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound)
			    else
			     sb1(:,l)=-eps*eps*r*r*omcyb*sqgeq(:,l)/(2*rbound*mm(l))
			    endif
		    end do
		    call multb(omdzprpeq,omdzeq,1,sb1,1,0.0_IDP,1.0_IDP)

                end if

		call multb(sb1,bsteq,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb2,sb1,-1,0.0_IDP,1.0_IDP,0)
		call multb(dbsjtojeq,sb2,1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(dbsjtbjeq,sb2,1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydztb(sb2,sb1,-1,0.0_IDP,1.0_IDP)
		call multb(dbsjzojeq,sb2,1,sqgieq,1,0.0_IDP,1.0_IDP)

		call dbydrb(sb1,gtteq,0.0_IDP,1.0_IDP,2)
		call multb(dgttreq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,grreq,1,0.0_IDP,1.0_IDP,2)
		call multb(dgrrteq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,grteq,-1,0.0_IDP,1.0_IDP,2)
		call multb(dgrtteq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,gtteq,1,0.0_IDP,1.0_IDP,2)
		call multb(dgttteq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,grreq,1,0.0_IDP,1.0_IDP)
		call multb(dgrrzeq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,grteq,-1,0.0_IDP,1.0_IDP)
		call multb(dgrtzeq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,gtteq,1,0.0_IDP,1.0_IDP)
		call multb(dgttzeq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
		call grparb(sb1,grteq,-1,0.0_IDP,1.0_IDP)
		call multb(dgrtpeq,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)
		call grparb(sb1,gtteq,1,0.0_IDP,1.0_IDP)
		call multb(dgttpeq,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)

		do l=1,lbmax
			dbsjtojeq(:,l)=r*dbsjtojeq(:,l)/sd1
			dbsjzojeq(:,l)=dbsjzojeq(:,l)/sd1
			dbsjtbjeq(:,l)=r*dbsjtbjeq(:,l)/sd1
			dgttreq(:,l)=dgttreq(:,l)/sd1
			dgrrteq(:,l)=r*dgrrteq(:,l)/sd1
			dgrtteq(:,l)=r*dgrtteq(:,l)/sd1
			dgttteq(:,l)=r*dgttteq(:,l)/sd1
			dgrrzeq(:,l)=dgrrzeq(:,l)/sd1
			dgrtzeq(:,l)=dgrtzeq(:,l)/sd1
			dgttzeq(:,l)=dgttzeq(:,l)/sd1
			dgrtpeq(:,l)=dgrtpeq(:,l)/sd1
			dgttpeq(:,l)=dgttpeq(:,l)/sd1
		end do

		sb1=bsteq
		call multb(bsqeq,sb1,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bsgrteq,grtojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bsgtteq,gttojeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bsqgtteq,bsgtteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)

!               FLR components
                
		call multb(sb1,bmodeq,1,sqgeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			bmodieq(:,l)=sb1(:,l)/(feq-qqinv*cureq)
		end do
		do l=1,lbmax
			brtdceq(:,l)=bmodeq(:,l)*r*qqinv/(feq-qqinv*cureq)
		end do
		do l=1,lbmax
			bdceq(:,l)=bmodeq(:,l)/(feq-qqinv*cureq)
		end do
		call multb(sb1,bmodeq,1,bmodeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			bbdceq(:,l)=sb1(:,l)/((feq-qqinv*cureq)*(feq-qqinv*cureq))
		end do
		call multb(sb2,sb1,1,bmodeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			bbbdceq(:,l)=sb2(:,l)/((feq-qqinv*cureq)*(feq-qqinv*cureq))
		end do

!               FLR component d/dr
                
		call dbydrb(sb1,brtdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=-sb3(:,l)*cureq*rinv
		end do
		call dbydthb(sb1,brtdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)*cureq*rinv
		end do
		call dbydztb(sb1,brtdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)*cureq*rinv*eps
		end do
		call dbydrb(sb1,bdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb3(:,l)*feq
		end do
		call dbydthb(sb1,bdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb3(:,l)*feq
		end do
		call dbydztb(sb1,bdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)*feq*eps
		end do
		call dbydrb(sb1,bmodieq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)
		end do
		call dbydthb(sb1,bmodieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)
		end do
		call dbydztb(sb1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb3(:,l)*eps
		end do
		call dbydztb(sb1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)-sb2(:,l)*r/(feq-qqinv*cureq)
		end do
		call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)
		call dbydr0(sd2,feq,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+(grrupeq(:,l)*sd2/(feq-qqinv*cureq))-(grrupeq(:,l)*qqinv*sd1/(feq-qqinv*cureq))
		end do
		call dbydthb(sb1,bsteq,-1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb2(:,l)*r/(feq-qqinv*cureq)
		end do
		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grrupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb3(:,l)+grrupeq(:,l)*rinv
		end do
		call dbydrb(sb1,grrupeq,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb1(:,l)
		end do
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb3(:,l)
		end do
		call dbydthb(sb1,grtupeq,-1,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb1(:,l)
		end do
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb3(:,l)*eps
		end do
		call dbydztb(sb1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplreq(:,l)=lplreq(:,l)+sb1(:,l)*eps
		end do

!               FLR component d/dth

		call dbydthb(sb1,bmodieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bbbdceq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=-sb2(:,l)*2*eps*eps*r*r*qqinv*qqinv
		end do
		call dbydztb(sb1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bbbdceq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb2(:,l)*2*eps*eps*r*qqinv
		end do
		call dbydrb(sb1,brtdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)*cureq*rinv
		end do
		call dbydthb(sb1,brtdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gttupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)*cureq*rinv
		end do
		call dbydztb(sb1,brtdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)*cureq*rinv*eps
		end do
		call dbydrb(sb1,bdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)*feq
		end do
		call dbydthb(sb1,bdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gttupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)*feq
		end do
		call dbydztb(sb1,bdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)*feq*eps
		end do
		call dbydrb(sb1,bmodieq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)
		end do
		call dbydthb(sb1,bmodieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gttupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)
		end do
		call dbydztb(sb1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-sb3(:,l)*eps
		end do
		call dbydztb(sb1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)-(sb2(:,l)*r/(feq-qqinv*cureq))+(grtupeq(:,l)*sd2/(feq-qqinv*cureq))-(grtupeq(:,l)*qqinv*sd1/(feq-qqinv*cureq))
		end do
		call dbydthb(sb1,bsteq,-1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb2(:,l)*r*r*qqinv/(feq-qqinv*cureq)
		end do
		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grtupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)
		end do
		call dbydrb(sb1,grtupeq,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb1(:,l)
		end do
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gttupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)
		end do
		call dbydthb(sb1,gttupeq,1,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb1(:,l)
		end do
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb3(:,l)*eps
		end do
		call dbydztb(sb1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplteq(:,l)=lplteq(:,l)+sb1(:,l)*eps
		end do

!               FLR component d/dz

		call dbydthb(sb1,bmodieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bbbdceq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=sb2(:,l)*2*eps*eps*r*qqinv
		end do
		call dbydztb(sb1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bbbdceq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb2(:,l)*2*eps*eps
		end do
		call dbydrb(sb1,brtdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*cureq*rinv*eps
		end do
		call dbydthb(sb1,brtdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*cureq*rinv*eps
		end do
		call dbydztb(sb1,brtdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gzzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*cureq*rinv*eps*eps
		end do
		call dbydrb(sb1,bdceq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*feq*eps
		end do
		call dbydthb(sb1,bdceq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*feq*eps
		end do
		call dbydztb(sb1,bdceq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gzzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*feq*eps*eps
		end do
		call dbydrb(sb1,bmodieq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*eps
		end do
		call dbydthb(sb1,bmodieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*eps
		end do
		call dbydztb(sb1,bmodieq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,bmodeq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gzzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-sb3(:,l)*eps*eps
		end do
		call dbydztb(sb1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)-(sb2(:,l)*eps*r/(feq-qqinv*cureq))+(grzupeq(:,l)*sd2*eps/(feq-qqinv*cureq))-(grzupeq(:,l)*sd1*eps*qqinv/(feq-qqinv*cureq))
		end do
		call dbydthb(sb1,bsteq,-1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb2(:,l)*eps*r*r*qqinv/(feq-qqinv*cureq)
		end do
		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,1,grzupeq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*eps+grzupeq(:,l)*eps*rinv
		end do
		call dbydrb(sb1,grzupeq,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb1(:,l)*eps
		end do
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gtzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*eps
		end do
		call dbydthb(sb1,gtzupeq,1,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*eps
		end do
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sb3,sb2,-1,gzzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb3(:,l)*eps*eps
		end do
		call dbydztb(sb1,gzzupeq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			lplzeq(:,l)=lplzeq(:,l)+sb1(:,l)*eps*eps
		end do

!               FLR component d2/drdt

		do l=1,lbmax
			lplrteq(:,l)=2*grtupeq(:,l)*rinv
		end do

!               FLR component d2/drdz

		do l=1,lbmax
			lplrzeq(:,l)=2*grzupeq(:,l)*eps
		end do

!               FLR component d2/dtdz

		do l=1,lbmax
			lpltzeq(:,l)=2*bbdceq(:,l)*eps*eps*r*qqinv+2*gtzupeq(:,l)*eps
		end do

!               FLR component d2/drdr

		do l=1,lbmax
			lplrreq(:,l)=grrupeq(:,l)
		end do

!               FLR component d2/dtdt

		do l=1,lbmax
			lpltteq(:,l)=-bbdceq(:,l)*eps*eps*r*r*qqinv*qqinv+gttupeq(:,l)
		end do

!               FLR component d2/dzdz

		do l=1,lbmax
			lplzzeq(:,l)=-bbdceq(:,l)*eps*eps+gzzupeq(:,l)*eps*eps
		end do

                sqgdrojeq=0.0_IDP	
                sqgdthojeq=0.0_IDP	
                sqgdztojeq=0.0_IDP	
                sqgdthojbsteq=0.0_IDP	
                sqgdztojbsteq=0.0_IDP	
                sqgibmoditheq=0.0_IDP	
                sqgibmodizteq=0.0_IDP	

		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)	
		call multb(sb2,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)	
		do l=1,lbmax				
                    sqgdrojeq(:,l)=sb2(:,l)
		end do					
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)			
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax				
                    sqgdthojeq(:,l)=sb2(:,l)	
		end do 				
		call multb(sqgdthojbsteq,sqgdthojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)				
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)	
		call multb(sb2,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)	
		do l=1,lbmax			
                     sqgdztojeq(:,l)=sb2(:,l)	
		end do		 			
		call multb(sqgdztojbsteq,sqgdztojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,sqgieq,1,0.0_IDP,1.0_IDP,0)
		call multb(sb2,bmodeq,1,sqgeq,1,0.0_IDP,1.0_IDP)	
		call multb(sqgibmoditheq,sb1,-1,sb2,1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(sqgibmodizteq,sb1,-1,sb2,1,0.0_IDP,1.0_IDP)	
		
!! Calculate the Landau damping elements		
        if(ieldamp_on .eq. 1) then	    
	  
          cmplx1 = cmplx(0._IDP,1._IDP,kind=IDP)                                  !! complex number i 
          tauie = tinn(0)/tenn(0)
	  xnuelc0 = xnuelc0/omegar    
          xnuion0 = 0.0165*xnuelc0/(tauie*sqrt(uion*tauie))                       !! ion-ion collision FR axis (Te=Ti axis)

          sb5=0.0_IDP	  
          do l=1,lbmax
            if (l == lb0) cycle     
            m=mm(l)
            n=nn(l)                                                                
            do j=0,mj
              vther = vthe*sqrt(te(j))
              vthir = vthi*sqrt(ti(j))			  
              if(vthir .le. 0.) vthir=.001
              if(vther .le. 0.) vther=.01
              if(te(j) .le. 0.) te(j) = .01
              tauie = tinn(j)/tenn(j)			  
              xnuelc = xnuelc0*dni(j)/(te(j)*sqrt(te(j)))                          !! electron-ion collision FR
              xnuion = xnuion0*dni(j)/(ti(j)*sqrt(ti(j)))                      !! ion-ion collision FR
              xkprl = real(n,kind=IDP) - real(m,kind=IDP)*qqinv(j)
              abkprl = abs(xkprl)
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
                     zi = cmplx(zir,zii,kind=IDP)
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
                     ze = cmplx(zer,zei,kind=IDP)
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
              stfi = 1. + (y0i/ddii)/tauie
              reii = 1. + stfe + stfi
              rei = 1./reii
              sei1 = y2e + y2i*tauie
              stfe = (y1e*y1e)/ddee
              stfi = (y1i*y1i)/ddii
              sei2 = -cmplx1*(xnuion*stfi + xnuelc*stfe)
              sei3 = rei*(((y1e/sqrt(ddee)) - (y1i/sqrt(ddii)))**2)
              sei = sei1 + sei2 + sei3 
              sb5(j,l) = aimag(sei)	  
            end do
           end do			   
					
		        allocate (eildreq(0:mj,0:lbmax),eildteq(0:mj,0:lbmax),eildzeq(0:mj,0:lbmax),eildrreq(0:mj,0:lbmax), &
			          eildrteq(0:mj,0:lbmax),eildrzeq(0:mj,0:lbmax),eildtteq(0:mj,0:lbmax),eildtzeq(0:mj,0:lbmax), &
			          eildzzeq(0:mj,0:lbmax))

			call dbydrb(sb4,omdreq,0.0_IDP,1.0_IDP,3)
			call multb(sb1,omdreq,-1,sb4,-1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdreq,-1,0.0_IDP,1.0_IDP,3)
			call multb(sb1,omdteq,1,sb4,1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(sb1,omdzeq,1,sb4,1,1.0_IDP,1.0_IDP)

			call dbydrb(sb4,omdteq,0.0_IDP,1.0_IDP,3)
			do l=1,lbmax
				sb4(:,l)=sb4(:,l)-rinv*omdteq(:,l)
				if (mmb(l) == 0 .or. abs(mmb(l)) == 2) sb4(0,l)=sb4(0,l)-rinv(1)*omdteq(1,l)
			end do
			call multb(sb2,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdteq,1,0.0_IDP,1.0_IDP,3)
			call multb(sb2,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(sb2,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

			call dbydrb(sb4,omdzeq,0.0_IDP,1.0_IDP,0)
			call multb(sb3,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdzeq,1,0.0_IDP,1.0_IDP,0)
			call multb(sb3,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdzeq,1,0.0_IDP,1.0_IDP)
			call multb(sb3,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

			call multb(eildreq,sb5,1,sb1,1,0.0_IDP,1.0_IDP)
			call multb(eildteq,sb5,1,sb2,-1,0.0_IDP,1.0_IDP)
			call multb(eildzeq,sb5,1,sb3,-1,0.0_IDP,1.0_IDP)

			call multb(sb1,sb5,1,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(eildrreq,sb1,-1,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(sb1,sb5,1,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(eildrteq,sb1,-1,omdteq,1,0.0_IDP,2.0_IDP)
			call multb(sb1,sb5,1,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(eildrzeq,sb1,-1,omdzeq,1,0.0_IDP,2.0_IDP)
			call multb(sb1,sb5,1,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(eildtteq,sb1,1,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(sb1,sb5,1,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(eildtzeq,sb1,1,omdzeq,1,0.0_IDP,2.0_IDP)
			call multb(sb1,sb5,1,omdzeq,1,0.0_IDP,1.0_IDP)
			call multb(eildzzeq,sb1,1,omdzeq,1,0.0_IDP,1.0_IDP)

			eildr=0.0_IDP
			eildt=0.0_IDP
			eildz=0.0_IDP
			eildrr=0.0_IDP
			eildrt=0.0_IDP
			eildrz=0.0_IDP
			eildtt=0.0_IDP
			eildtz=0.0_IDP
			eildzz=0.0_IDP

			do lq=1,leqmax
				do l=1,lbmax
					if(mmb(l).eq.mmeq(lq).and.nnb(l).eq.nneq(lq)) exit
				end do
				if (l > lbmax) cycle
				eildr(:,lq)=sd1*preq*eildreq(:,l)
				eildt(:,lq)=sd1*preq*eildteq(:,l)
				eildz(:,lq)=sd1*preq*eildzeq(:,l)
				eildrr(:,lq)=sd1*preq*eildrreq(:,l)
				eildrt(:,lq)=sd1*preq*eildrteq(:,l)
				eildrz(:,lq)=sd1*preq*eildrzeq(:,l)
				eildtt(:,lq)=sd1*preq*eildtteq(:,l)
				eildtz(:,lq)=sd1*preq*eildtzeq(:,l)
				eildzz(:,lq)=sd1*preq*eildzzeq(:,l)
			end do

			deallocate (eildreq,eildteq,eildzeq,eildrreq,eildrteq,eildrzeq,eildtteq,eildtzeq,eildzzeq)
				
		end if
		
		sqgi=0.0_IDP
		sqg=0.0_IDP
		bst=0.0_IDP
		grr=0.0_IDP
		grt=0.0_IDP
		gtt=0.0_IDP
		bmod=0.0_IDP
		grroj=0.0_IDP
		grtoj=0.0_IDP
		grzoj=0.0_IDP
		gttoj=0.0_IDP
		gtzoj=0.0_IDP
                djroj=0.0_IDP
                djtojeq=0.0_IDP
                djzojeq=0.0_IDP
		jbgrr=0.0_IDP
		jbgrt=0.0_IDP
		jbgtt=0.0_IDP
		omdr=0.0_IDP
		omdt=0.0_IDP
		omdz=0.0_IDP
		dbsjtoj=0.0_IDP
		dbsjzoj=0.0_IDP
		dbsjtbj=0.0_IDP
 		dgttr=0.0_IDP
 		dgrrt=0.0_IDP
 		dgrtt=0.0_IDP
 		dgttt=0.0_IDP
 		dgrrz=0.0_IDP
 		dgrtz=0.0_IDP
 		dgttz=0.0_IDP
 		dgrtp=0.0_IDP
		dgttp=0.0_IDP
		jsq=0.0_IDP
		bsq=0.0_IDP
		bsgrt=0.0_IDP
		bsgtt=0.0_IDP
		bsqgtt=0.0_IDP
		lplrr=0.0_IDP
		lplrt=0.0_IDP
		lplrz=0.0_IDP
		lpltz=0.0_IDP
		lpltt=0.0_IDP
		lplzz=0.0_IDP
		lplr=0.0_IDP
		lplt=0.0_IDP
		lplz=0.0_IDP

                sqgdroj=0.0_IDP
		sqgdthoj=0.0_IDP
		sqgdztoj=0.0_IDP
  		sqgdthojbst=0.0_IDP
		sqgdztojbst=0.0_IDP
		sqgibmodith=0.0_IDP
		sqgibmodizt=0.0_IDP

		test=0.0_IDP
		testr=0.0_IDP
		testt=0.0_IDP
		testrr=0.0_IDP
		testtt=0.0_IDP
		testrt=0.0_IDP

                if(Trapped_on .eq. 1) then
		    omdrprp=0.0_IDP
		    omdtprp=0.0_IDP
		    omdzprp=0.0_IDP
		end if	
	
		do lq=1,leqmax
			do l=1,lbmax
				if(mmb(l).eq.mmeq(lq).and.nnb(l).eq.nneq(lq)) exit
			end do
			if (l > lbmax) cycle
			sqgi(:,lq)=sqgieq(:,l)
			sqg(:,lq)=sqgeq(:,l)
			bst(:,lq)=bsteq(:,l)
			grr(:,lq)=grreq(:,l)
			grt(:,lq)=grteq(:,l)
			gtt(:,lq)=gtteq(:,l)
			grz(:,lq)=grzeq(:,l)
			gtz(:,lq)=gtzeq(:,l)
			gzz(:,lq)=gzzeq(:,l)			
			bmod(:,lq)=bmodeq(:,l)
			grroj(:,lq)=grrojeq(:,l)
			grtoj(:,lq)=grtojeq(:,l)
			grzoj(:,lq)=grzojeq(:,l)
			gttoj(:,lq)=gttojeq(:,l)
			gtzoj(:,lq)=gtzojeq(:,l)
			gzzoj(:,lq)=gzzojeq(:,l)	
			jbgrr(:,lq)=jbgrreq(:,l)
			jbgrt(:,lq)=jbgrteq(:,l)
			jbgtt(:,lq)=jbgtteq(:,l)
			jbgrz(:,lq)=jbgrzeq(:,l)
			jbgtz(:,lq)=jbgtzeq(:,l)
			djroj(:,lq)=djrojeq(:,l)
			djtoj(:,lq)=djtojeq(:,l)
			djzoj(:,lq)=djzojeq(:,l)
			omdr(:,lq)=omdreq(:,l)
			omdt(:,lq)=omdteq(:,l)
			omdz(:,lq)=omdzeq(:,l)
			grrup(:,lq)=grrupeq(:,l)	
			grtup(:,lq)=grtupeq(:,l)			
			grzup(:,lq)=grzupeq(:,l)
			gttup(:,lq)=gttupeq(:,l)
			gtzup(:,lq)=gtzupeq(:,l)
			gzzup(:,lq)=gzzupeq(:,l)
			dbsjtoj(:,lq)=dbsjtojeq(:,l)
			dbsjzoj(:,lq)=dbsjzojeq(:,l)
			dbsjtbj(:,lq)=dbsjtbjeq(:,l)
 			dgttr(:,lq)=dgttreq(:,l)
 			dgrrt(:,lq)=dgrrteq(:,l)
 			dgrtt(:,lq)=dgrtteq(:,l)
 			dgttt(:,lq)=dgttteq(:,l)
 			dgrrz(:,lq)=dgrrzeq(:,l)
 			dgrtz(:,lq)=dgrtzeq(:,l)
 			dgttz(:,lq)=dgttzeq(:,l)
 			dgrtp(:,lq)=dgrtpeq(:,l)
			dgttp(:,lq)=dgttpeq(:,l)
			jsq(:,lq)=jsqeq(:,l)
			bsq(:,lq)=bsqeq(:,l)
			bsgrt(:,lq)=bsgrteq(:,l)
			bsgtt(:,lq)=bsgtteq(:,l)
			bsqgtt(:,lq)=bsqgtteq(:,l)
			lplrr(:,lq)=lplrreq(:,l)
			lplrt(:,lq)=lplrteq(:,l)
			lplrz(:,lq)=lplrzeq(:,l)
			lpltt(:,lq)=lpltteq(:,l)
			lpltz(:,lq)=lpltzeq(:,l)
			lplzz(:,lq)=lplzzeq(:,l)	
			lplr(:,lq)=lplreq(:,l)
			lplt(:,lq)=lplteq(:,l)
			lplz(:,lq)=lplzeq(:,l)

                        sqgdroj(:,lq)=sqgdrojeq(:,l)
                        sqgdthoj(:,lq)=sqgdthojeq(:,l)
                        sqgdztoj(:,lq)=sqgdztojeq(:,l)
                        sqgdthojbst(:,lq)=sqgdthojbsteq(:,l)
                        sqgdthojbst(:,lq)=sqgdztojbsteq(:,l)
                        sqgibmodith(:,lq)=sqgibmoditheq(:,l)
                        sqgibmodizt(:,lq)=sqgibmodizteq(:,l)

                        test(:,lq)=testeq(:,l)
                        testr(:,lq)=testreq(:,l)
                        testt(:,lq)=testteq(:,l)
                        testrr(:,lq)=testrreq(:,l)
                        testtt(:,lq)=testtteq(:,l)
                        testrt(:,lq)=testrteq(:,l)

                         if(Trapped_on .eq. 1) then

			    omdrprp(:,lq)=omdrprpeq(:,l)
			    omdtprp(:,lq)=omdtprpeq(:,l)
			    omdzprp(:,lq)=omdzprpeq(:,l)

		         end if	

		
		end do							

		write(0,'(" ====> Equilibria set up DONE !! ")')
		
		deallocate (sqgieq)
		deallocate (sqgeq)
		deallocate (grreq)
		deallocate (grteq)
		deallocate (gtteq)
		deallocate (grzeq)
		deallocate (gtzeq)
		deallocate (gzzeq)		
		deallocate (bmodeq)
		deallocate (grrojeq)
		deallocate (grtojeq)
		deallocate (grzojeq)
		deallocate (gttojeq)	
		deallocate (gtzojeq)
		deallocate (gzzojeq)
		deallocate (jbgrreq)
		deallocate (jbgrteq)
		deallocate (jbgtteq)
		deallocate (jbgrzeq)
		deallocate (jbgtzeq)
		deallocate (bsteq)
		deallocate (lplrreq)
		deallocate (lplrteq)
		deallocate (lplrzeq)
		deallocate (lpltzeq)
		deallocate (lpltteq)
		deallocate (lplzzeq)
		deallocate (lplreq)
		deallocate (lplteq)
		deallocate (lplzeq)
		deallocate (omdreq)
		deallocate (omdteq)
		deallocate (omdzeq)
		deallocate (sb1)
		deallocate (sb2)
		deallocate (sb3)		
		deallocate (sb4)
		deallocate (sb5)	
		deallocate (sb6)			
		deallocate (sfar1)
		deallocate (sfar2)
		deallocate (sfar3)
		deallocate (sfar4)

		deallocate (mmstartb)
		deallocate (mmendb)
		deallocate (llb)

		deallocate (mmb)
		deallocate (nnb)
		deallocate (rfar)
		deallocate (rbinv)
		deallocate (qfar)
		deallocate (pfar)
		deallocate (phip)
		deallocate (curfar)
		deallocate (ffar)
		deallocate (sqgib)
		deallocate (sqgb)
		deallocate (grrb)
		deallocate (grtb)
		deallocate (gttb)
		deallocate (bmodb)
		deallocate (grrojb)
		deallocate (grtojb)
		deallocate (gttojb)
		deallocate (jbgrrb)
		deallocate (jbgrtb)
		deallocate (jbgttb)
		
		deallocate (grrupeq)	
		deallocate (grtupeq)	
		deallocate (grzupeq)	
		deallocate (gttupeq)	
		deallocate (gtzupeq)	
		deallocate (gzzupeq)

		deallocate (sqgdrojeq)
		deallocate  (sqgdthojeq)
		deallocate (sqgdztojeq)
		deallocate (sqgdthojbsteq)
		deallocate (sqgdztojbsteq)
		deallocate (sqgibmoditheq)
		deallocate (sqgibmodizteq)

		deallocate (testeq)
		deallocate (testreq)
		deallocate (testteq)
		deallocate (testrreq)
		deallocate (testtteq)
		deallocate (testrteq)

                if(Trapped_on .eq. 1) then

		    deallocate (omdrprpeq)
		    deallocate (omdtprpeq)
		    deallocate (omdzprpeq)

                end if
		
		close(25)
		
	end subroutine vmec

	subroutine mmblims

	!   Calculate the minimum and maximum m values for each n value
	!   in the complex exponential representation.

		use param
		use domain
		implicit none

		integer :: n,l

		allocate (mmstartb(-nmaxb:nmaxb))
		allocate (mmendb(-nmaxb:nmaxb))

		do n = 0,nmaxb
			mmstartb(n) = 1000
			mmendb(n) = -1000
		end do

		do l = 1,lbmax

			n = nnb(l)
			if (n >= 0) then
				if (mmb(l) < mmstartb(n)) mmstartb(n) = mmb(l)
				if (mmb(l) > mmendb(n)) mmendb(n) = mmb(l)
			else
				if (-mmb(l) < mmstartb(-n)) mmstartb(-n) = -mmb(l)
				if (-mmb(l) > mmendb(-n)) mmendb(-n) = -mmb(l)
			end if

		end do

		!   in the complex exponential format, [fr,fi](-m,-n) = [fr,-fi](m,n),
		!   so limits for negative n rows are just negatives of corresponding
		!   positive n limits.

		do n = 0,nmaxb

			mmstartb(-n) = -mmendb(n)
			mmendb(-n) = -mmstartb(n)

		end do

	end subroutine mmblims
