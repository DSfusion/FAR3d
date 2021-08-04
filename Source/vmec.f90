	subroutine vmec

		use param
		use cotrol
		use domain
		use equil
		use scratch
		implicit none
 
		integer :: j,l,l1,lq,mjeqp,mjeqm1,mpol,ntor,lb0,m,m2,mm1,mm2,mm12,mm22,mm2p1,mm2p12,mp2,mp22,isg,idummy,mband,n,lbm2,ir,jp,lp,nbst
		real(IDP) :: bigrn,pror,dummy,bet0in,twopi,twophip,one,gc,iota,iotap,pprime,qmin,qmax,p,xl,pprimel
		integer, dimension(:), allocatable :: llc,lls,lbst
		real(IDP), dimension(:), allocatable :: rbinv,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,rsb,Abst,Wbst
		real(IDP), dimension(:,:), allocatable :: lpltr,lplzr,lplzt
		real(IDP), dimension(:,:), allocatable :: rmnb,sqgib,sqgb,grrb,grtb,gttb,bmodb,grrojb,grtojb,gttojb,jbgrrb,jbgrtb,jbgttb
		real(IDP), dimension(:,:), allocatable :: sqgieq,sqgeq,grreq,grteq,gtteq,bmodeq,bsteq,jbgrreq,jbgrteq,jbgtteq, &
							  grrojeq,grtojeq,gttojeq,lplreq,lplteq,djrojeq,djtojeq,djzojeq,omdreq,omdteq,&
							  omdzeq,dbsjtojeq,dbsjzojeq,dbsjtbjeq,dgttreq,dgrrteq,dgrtteq,dgttteq,dgrrzeq, &
							  dgrtzeq,dgttzeq,dgrtpeq,dgttpeq,jsqeq,bsgrteq,bsgtteq,bsqeq,bsqgtteq, &
							  lplrreq,lplrteq,lplrzeq,lpltreq,lpltteq,lpltzeq,lplzreq,lplzteq,lplzzeq, &
							  bsjgrteq,bsjgtteq,bstgeq,bsqgeq,eildreq,eildteq,eildzeq,eildrreq,eildrteq, &
							  eildrzeq,eildtteq,eildtzeq,eildzzeq,sb1,sb2,sb3,sb4,sb5,sb6
		LOGICAL :: lasym

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
		end interface

		bet0in=bet0

		allocate (sqg(0:mj,0:leqmax),sqgi(0:mj,0:leqmax),bmod(0:mj,0:leqmax),bst(0:mj,0:leqmax), &
			  grr(0:mj,0:leqmax),grt(0:mj,0:leqmax),gtt(0:mj,0:leqmax), &
			  grroj(0:mj,0:leqmax),grtoj(0:mj,0:leqmax),gttoj(0:mj,0:leqmax), &
			  jbgrr(0:mj,0:leqmax),jbgrt(0:mj,0:leqmax),jbgtt(0:mj,0:leqmax),lplr_r(0:mj,0:leqmax),lplt_r(0:mj,0:leqmax), &
			  lplr(0:mj,0:leqmax),lplt(0:mj,0:leqmax),lplz(0:mj,0:leqmax),djroj(0:mj,0:leqmax),djtoj(0:mj,0:leqmax), &
			  djzoj(0:mj,0:leqmax),jsq(0:mj,0:leqmax),omdr(0:mj,0:leqmax),omdt(0:mj,0:leqmax),omdz(0:mj,0:leqmax), &
			  dbsjtoj(0:mj,0:leqmax),dbsjzoj(0:mj,0:leqmax),dbsjtbj(0:mj,0:leqmax),dgttr(0:mj,0:leqmax), &
			  dgrrt(0:mj,0:leqmax),dgrtt(0:mj,0:leqmax),dgttt(0:mj,0:leqmax),dgrrz(0:mj,0:leqmax), &
			  dgrtz(0:mj,0:leqmax),dgttz(0:mj,0:leqmax),dgrtp(0:mj,0:leqmax),dgttp(0:mj,0:leqmax), &
			  bsgrt(0:mj,0:leqmax),bsgtt(0:mj,0:leqmax),bsq(0:mj,0:leqmax),bsqgtt(0:mj,0:leqmax), &
			  lplrr(0:mj,0:leqmax),lplrt(0:mj,0:leqmax),lplrz(0:mj,0:leqmax),lpltr(0:mj,0:leqmax),lpltt(0:mj,0:leqmax), &
			  lpltz(0:mj,0:leqmax),lplzr(0:mj,0:leqmax),lplzt(0:mj,0:leqmax),lplzz(0:mj,0:leqmax), &
			  bsjgrt(0:mj,0:leqmax),bsjgtt(0:mj,0:leqmax),bstg(0:mj,0:leqmax),bsqg(0:mj,0:leqmax))
		if (ieldamp_on == 1) allocate (eildr(0:mj,0:leqmax),eildt(0:mj,0:leqmax),eildz(0:mj,0:leqmax),eildrr(0:mj,0:leqmax), &
					       eildrt(0:mj,0:leqmax),eildrz(0:mj,0:leqmax),eildtt(0:mj,0:leqmax),eildtz(0:mj,0:leqmax), &
					       eildzz(0:mj,0:leqmax))

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

		qq=1.0_IDP/qqinv
	
		allocate (sqgeq(0:mj,0:lbmax),sqgieq(0:mj,0:lbmax),bmodeq(0:mj,0:lbmax),bsteq(0:mj,0:lbmax), &
			  grreq(0:mj,0:lbmax),grteq(0:mj,0:lbmax),gtteq(0:mj,0:lbmax), &
			  grrojeq(0:mj,0:lbmax),grtojeq(0:mj,0:lbmax),gttojeq(0:mj,0:lbmax), &
			  jbgrreq(0:mj,0:lbmax),jbgrteq(0:mj,0:lbmax),jbgtteq(0:mj,0:lbmax), &
			  lplreq(0:mj,0:lbmax),lplteq(0:mj,0:lbmax),djrojeq(0:mj,0:lbmax),djtojeq(0:mj,0:lbmax),djzojeq(0:mj,0:lbmax), &
			  jsqeq(0:mj,0:lbmax),omdreq(0:mj,0:lbmax),omdteq(0:mj,0:lbmax),omdzeq(0:mj,0:lbmax), &
			  dbsjtojeq(0:mj,0:lbmax),dbsjzojeq(0:mj,0:lbmax),dbsjtbjeq(0:mj,0:lbmax),dgttreq(0:mj,0:lbmax), &
			  dgrrteq(0:mj,0:lbmax),dgrtteq(0:mj,0:lbmax),dgttteq(0:mj,0:lbmax),dgrrzeq(0:mj,0:lbmax), &
			  dgrtzeq(0:mj,0:lbmax),dgttzeq(0:mj,0:lbmax),dgrtpeq(0:mj,0:lbmax),dgttpeq(0:mj,0:lbmax), &
			  bsgrteq(0:mj,0:lbmax),bsgtteq(0:mj,0:lbmax),bsqeq(0:mj,0:lbmax),bsqgtteq(0:mj,0:lbmax), &
			  lplrreq(0:mj,0:lbmax),lplrteq(0:mj,0:lbmax),lplrzeq(0:mj,0:lbmax),lpltreq(0:mj,0:lbmax),lpltteq(0:mj,0:lbmax), &
			  lpltzeq(0:mj,0:lbmax),lplzreq(0:mj,0:lbmax),lplzteq(0:mj,0:lbmax),lplzzeq(0:mj,0:lbmax), &
			  bsjgrteq(0:mj,0:lbmax),bsjgtteq(0:mj,0:lbmax),bstgeq(0:mj,0:lbmax),bsqgeq(0:mj,0:lbmax), &
			  sb1(0:mj,0:lbmax),sb2(0:mj,0:lbmax),sb3(0:mj,0:lbmax),sb4(0:mj,0:lbmax),sb5(0:mj,0:lbmax),sb6(0:mj,0:lbmax))

		sqgieq=0.0_IDP
		sqgeq=0.0_IDP
		bsteq=0.0_IDP
		grreq=0.0_IDP
		grteq=0.0_IDP
		gtteq=0.0_IDP
		bmodeq=0.0_IDP
		grrojeq=0.0_IDP
		grtojeq=0.0_IDP
		gttojeq=0.0_IDP
		jbgrreq=0.0_IDP
		jbgrteq=0.0_IDP
		jbgtteq=0.0_IDP

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

		if (ext_prof == 1) call ae_profiles

!	bst is modified for resonant modes

		allocate (lbst(lbmax))
		allocate (rsb(lbmax), Abst(lbmax), Wbst(lbmax))

		bsteq=0
		one = 1.0_IDP
		lbst=0
		Abst=0
		qmin=minval(qqinv(0:mj))
		qmax=maxval(qqinv(0:mj))
		call dbydr0(sd1,preq,0.0_IDP,-1.0_IDP,0)
		call dbydr0(sd3,qqinv,0.0_IDP,1.0_IDP,0)
		lp=0
		do l=1,lbmax
			if (mmb(l) == 0 .or. nnb(l) == 0) cycle
			p=one*nnb(l)/mmb(l)
			if ((p-qmin)*(p-qmax) >= 0.0_IDP) cycle
			lp=lp+1
			lbst(lp)=l
			do j=1,mj
				if ((qqinv(j-1)-p)*(qqinv(j)-p) <= 0.) exit
			end do
			jp=j
			rsb(lp)=r(jp-1)+(r(jp)-r(jp-1))*(qqinv(jp-1)-p)/(qqinv(jp-1)-qqinv(jp))
!			write(6,'("m =",i5," n =",i5," rs =",1pe13.6," r(",i3,") =",1pe13.6," r(",i3,") =",1pe13.6)') mmb(l), &
!				nnb(l),rsb(lp),jp-1,r(jp-1),jp,r(jp)
			if (jp == 1) then
				Abst(lp) = (preq(0)-preq(1))/r(1)
				iotap = (qqinv(1)-qqinv(0))/r(1)
			else if (jp < mj) then
				Abst(lp) = ((r(jp)-rsb(lp))*sd1(jp-1)+(rsb(lp)-r(jp-1))*sd1(jp))/(r(jp)-r(jp-1))
				iotap = ((r(jp)-rsb(lp))*sd3(jp-1)+(rsb(lp)-r(jp-1))*sd3(jp))/(r(jp)-r(jp-1))
			else
				Abst(lp) = (preq(mjm1)-preq(mj))/(r(mj)-r(mjm1))
				iotap = (qqinv(mj)-qqinv(mjm1))/(r(mj)-r(mjm1))
			end if
			Wbst(lp) = abs(p/(mmb(l)*iotap))
		end do
		nbst=lp
		write(6,'(/"bst is modified for",i5," resonant modes")') nbst
		do l=1,nbst
			write(6,'("m =",i5," n =",i5," rs =",1pe13.6," A =",1pe13.6," W =",1pe13.6)') mmb(lbst(l)), &
				nnb(lbst(l)),rsb(l),Abst(l),Wbst(l)
		end do

		sd2=rinv*sd1/2.0_IDP
		sd2(0)=(preq(0)-preq(1))/r(1)**2
		do l=1,nbst
			do j=0,mj
				pprime=sd2(j)
				iota=qqinv(j)
				gc=sqgeq(j,lbst(l))
				xl=r(j)-rsb(l)
				pprimel=pprime-0.5*rinv(j)*Abst(l)*(1-(xl/Wbst(l))**2)*exp(-0.5*(xl/Wbst(l))**2)
				bsteq(j,lbst(l))=bet0*pprimel*gc/(mmb(lbst(l))*iota-nnb(lbst(l)))
			end do
		end do

		do l=1,lbmax
			if(mmb(l) == 0 .and. nnb(l) == 0) cycle
			p=0.0
			if (mmb(l) /= 0) p=one*nnb(l)/mmb(l)
			if ((p-qmin)*(p-qmax) < 0.0_IDP) cycle
			bsteq(:,l)=bet0*sd2*sqgeq(:,l)/(mmb(l)*qqinv-nnb(l))
		end do
		
		sb1=sqgeq
		call multb(jsqeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
		call dbydrb(sb1,sqgieq,0.0_IDP,1.0_IDP,0)
		call multb(lplreq,gtteq,1,sb1,1,0.0_IDP,-1.0_IDP)
		call dbydthb(sb2,sqgieq,1,0.0_IDP,1.0_IDP,0)
		call multb(lplreq,grteq,-1,sb2,-1,1.0_IDP,1.0_IDP)
		call multb(lplteq,grteq,-1,sb1,1,0.0_IDP,1.0_IDP)
		call multb(lplteq,grreq,1,sb2,-1,1.0_IDP,-1.0_IDP)
		
		call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
		call multb(djrojeq,sqgieq,1,sb1,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
		call multb(djtojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)
		call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
		call multb(djzojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)

!	WARNING: sd1 is feq-qqinv*cureq

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
		do l=1,lbmax
			sb1(:,l)=sd1*r*omdzeq(:,l)-rinv*cureq*sb2(:,l)
		end do
		call multb(omdzeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			omdzeq(:,l)=omdzeq(:,l)/(2.0_IDP*sd1*sd1)
		end do

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
		call dbydrb(sb1,sqgieq,0.0_IDP,1.0_IDP,0)
		do l=1,lbmax
			sb3(:,l)=sd1*gtteq(:,l)-rinv*rinv*cureq*cureq*sqgeq(:,l)/(eps*eps)
		end do
		call multb(lplrreq,sb3,1,sb1,1,0.0_IDP,-1.0_IDP)
		call multb(sb4,gtteq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb5,sqgeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=feq*grteq(:,l)-r*r*qqinv*sb4(:,l)-cureq*sb5(:,l)/(eps*eps)
		end do
		call multb(lplrteq,sb3,-1,sb1,1,0.0_IDP,1.0_IDP)
		call dbydthb(sb2,sqgieq,1,0.0_IDP,1.0_IDP,0)
		call multb(lpltreq,sb3,-1,sb2,-1,0.0_IDP,1.0_IDP)
		call multb(sb6,sb5,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(sb5,grteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=feq*(feq*grreq(:,l)-2*r*r*qqinv*sb5(:,l))/sd1-r*r*sb6(:,l)/(eps*eps)
		end do
		call multb(sb6,sb4,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=sb3(:,l)+r*r*r*r*qqinv*qqinv*sb6(:,l)/sd1
		end do
		call multb(lpltteq,sb3,1,sb2,-1,0.0_IDP,-1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=rinv*cureq*grteq(:,l)-r*sb4(:,l)
		end do
		call multb(lplrzeq,sb3,-1,sb1,1,0.0_IDP,-1.0_IDP)
		call dbydztb(sb1,sqgieq,1,0.0_IDP,1.0_IDP)
		call multb(lplzreq,sb3,-1,sb1,-1,0.0_IDP,-1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=(feq*(rinv*cureq*grreq(:,l)-r*sb5(:,l))-r*qqinv*(cureq*sb5(:,l)-r*r*sb6(:,l)))/sd1
		end do
		call multb(lpltzeq,sb3,1,sb2,-1,0.0_IDP,1.0_IDP)
		call multb(lplzteq,sb3,1,sb1,-1,0.0_IDP,1.0_IDP)
		do l=1,lbmax
			sb3(:,l)=(rinv*cureq*(rinv*cureq*grreq(:,l)-r*sb5(:,l))-(cureq*sb5(:,l)-r*r*sb6(:,l)))/sd1
		end do
		call multb(lplzzeq,sb3,1,sb1,-1,0.0_IDP,-1.0_IDP)

		call multb(bsjgrteq,jbgrteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bsjgtteq,jbgtteq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bstgeq,jsqeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
		call multb(bsqgeq,bstgeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)

		sqgi=0.0_IDP
		sqg=0.0_IDP
		bst=0.0_IDP
		grr=0.0_IDP
		grt=0.0_IDP
		gtt=0.0_IDP
		bmod=0.0_IDP
		grroj=0.0_IDP
		grtoj=0.0_IDP
		gttoj=0.0_IDP
		jbgrr=0.0_IDP
		jbgrt=0.0_IDP
		jbgtt=0.0_IDP
		lplr=0.0_IDP
		lplt=0.0_IDP
		djroj=0.0_IDP
		djtoj=0.0_IDP
		djzoj=0.0_IDP
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
		lpltr=0.0_IDP
		lpltt=0.0_IDP
		lpltz=0.0_IDP
		lplzr=0.0_IDP
		lplzt=0.0_IDP
		lplzz=0.0_IDP
		bsjgrt=0.0_IDP
		bsjgtt=0.0_IDP
		bstg=0.0_IDP
		bsqg=0.0_IDP

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
			bmod(:,lq)=bmodeq(:,l)
			grroj(:,lq)=grrojeq(:,l)
			grtoj(:,lq)=grtojeq(:,l)
			gttoj(:,lq)=gttojeq(:,l)
			jbgrr(:,lq)=jbgrreq(:,l)
			jbgrt(:,lq)=jbgrteq(:,l)
			jbgtt(:,lq)=jbgtteq(:,l)
			lplr(:,lq)=lplreq(:,l)
			lplt(:,lq)=lplteq(:,l)
			djroj(:,lq)=djrojeq(:,l)
			djtoj(:,lq)=djtojeq(:,l)
			djzoj(:,lq)=djzojeq(:,l)
			omdr(:,lq)=omdreq(:,l)
			omdt(:,lq)=omdteq(:,l)
			omdz(:,lq)=omdzeq(:,l)
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
			lpltr(:,lq)=lpltreq(:,l)
			lpltt(:,lq)=lpltteq(:,l)
			lpltz(:,lq)=lpltzeq(:,l)
			lplzr(:,lq)=lplzreq(:,l)
			lplzt(:,lq)=lplzteq(:,l)
			lplzz(:,lq)=lplzzeq(:,l)
			bsjgrt(:,lq)=bsjgrteq(:,l)
			bsjgtt(:,lq)=bsjgtteq(:,l)
			bstg(:,lq)=bstgeq(:,l)
			bsqg(:,lq)=bsqgeq(:,l)
		end do

		if (ieldamp_on == 1) then

			allocate (eildreq(0:mj,0:lbmax),eildteq(0:mj,0:lbmax),eildzeq(0:mj,0:lbmax),eildrreq(0:mj,0:lbmax), &
				  eildrteq(0:mj,0:lbmax),eildrzeq(0:mj,0:lbmax),eildtteq(0:mj,0:lbmax),eildtzeq(0:mj,0:lbmax), &
				  eildzzeq(0:mj,0:lbmax))

			call dbydrb(sb4,omdreq,0.0_IDP,1.0_IDP,3)
			call multb(eildreq,omdreq,-1,sb4,-1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdreq,-1,0.0_IDP,1.0_IDP,3)
			call multb(eildreq,omdteq,1,sb4,1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(eildreq,omdzeq,1,sb4,1,1.0_IDP,1.0_IDP)

			call dbydrb(sb4,omdteq,0.0_IDP,1.0_IDP,3)
			do l=1,lbmax
				sb4(:,l)=sb4(:,l)-rinv*omdteq(:,l)
				if (mmb(l) == 0 .or. abs(mmb(l)) == 2) sb4(0,l)=sb4(0,l)-rinv(1)*omdteq(1,l)
			end do
			call multb(eildteq,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdteq,1,0.0_IDP,1.0_IDP,3)
			call multb(eildteq,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(eildteq,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

			call dbydrb(sb4,omdzeq,0.0_IDP,1.0_IDP,0)
			call multb(eildzeq,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
			call dbydthb(sb4,omdzeq,1,0.0_IDP,1.0_IDP,0)
			call multb(eildzeq,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
			call dbydztb(sb4,omdzeq,1,0.0_IDP,1.0_IDP)
			call multb(eildzeq,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

			sb1=omdreq
			call multb(eildrreq,sb1,-1,omdreq,-1,0.0_IDP,1.0_IDP)
			call multb(eildrteq,sb1,-1,omdteq,1,0.0_IDP,2.0_IDP)
			call multb(eildrzeq,sb1,-1,omdzeq,1,0.0_IDP,2.0_IDP)
			sb1=omdteq
			call multb(eildtteq,sb1,1,omdteq,1,0.0_IDP,1.0_IDP)
			call multb(eildtzeq,sb1,1,omdzeq,1,0.0_IDP,2.0_IDP)
			sb1=omdzeq
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

!	equilibrium arrays for delperpsq

		lplr_r=lplr
		lplt_r=lplt
		call dbydreq(lplr_r,gttoj,1.0_IDP,1.0_IDP,2)
		call dbydtheq(lplr_r,grtoj,-1,1.0_IDP,-1.0_IDP,2)
		do l=1,leqmax
			lplr_r(:,l)=lplr_r(:,l)+rinv*gttoj(:,l)
		end do
		call dbydtheq(lplt_r,grroj,1,1.0_IDP,1.0_IDP,2)
		call dbydreq(lplt_r,grtoj,1.0_IDP,-1.0_IDP,2)

		lplr=lplrr+lpltr+lplzr
		lplt=lplrt+lpltt+lplzt
		lplz=lplrz+lpltz+lplzz

		do l=1,leqmax
			lplrr(:,l)=sd1*gttoj(:,l)
		end do
		lplrr(:,leq0)=lplrr(:,leq0)-rinv*rinv*cureq*cureq/(eps*eps)
		call dbydreq(lplr,lplrr,1.0_IDP,1.0_IDP,2)
		do l=1,leqmax
			lplr(:,l)=lplr(:,l)+rinv*lplrr(:,l)
		end do

		do l=1,leqmax
			lplrt(:,l)=-feq*grtoj(:,l)+r*r*qqinv*bsgtt(:,l)+cureq*bst(:,l)/(eps*eps)
		end do
		call dbydtheq(lplr,lplrt,-1,1.0_IDP,1.0_IDP,2)
		call dbydreq(lplt,lplrt,1.0_IDP,1.0_IDP,2)
		lplrt=2.0*lplrt

		do l=1,leqmax
			lplrz(:,l)=rinv*cureq*grtoj(:,l)-r*bsgtt(:,l)
		end do
		call dbydzteq(lplr,lplrz,-1,1.0_IDP,1.0_IDP)
		call dbydreq(lplz,lplrz,1.0_IDP,1.0_IDP,3)
		do l=1,leqmax
			lplz(:,l)=lplz(:,l)+rinv*lplrz(:,l)
		end do
		lplrz=2.0*lplrz

		do l=1,leqmax
			lpltt(:,l)=(feq*feq*grroj(:,l)-r*r*qqinv*(2*feq*bsgrt(:,l)-r*r*qqinv*bsqgtt(:,l)))/sd1-r*r*bsq(:,l)/(eps*eps)
		end do
		call dbydtheq(lplt,lpltt,1,1.0_IDP,1.0_IDP,2)

		do l=1,leqmax
			lpltz(:,l)=-(rinv*feq*cureq*grroj(:,l)-r*((feq+qqinv*cureq)*bsgrt(:,l)-r*r*qqinv*bsqgtt(:,l)))/sd1
		end do
		call dbydzteq(lplt,lpltz,1,1.0_IDP,1.0_IDP)
		call dbydtheq(lplz,lpltz,1,1.0_IDP,1.0_IDP,3)
		lpltz=2.0*lpltz

		do l=1,leqmax
			lplzz(:,l)=(cureq*(rinv*rinv*cureq*grroj(:,l)-2*bsgrt(:,l))+r*r*bsqgtt(:,l))/sd1
		end do
		call dbydzteq(lplz,lplzz,1,1.0_IDP,1.0_IDP)

!	END WARNING

		write(0,'(" ====> Equilibria set up DONE !! ")')

		deallocate (llb,mmb,nnb,rfar,rbinv,qfar,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,lbst,rsb,Abst,Wbst)
		if (lasym) deallocate (llc,lls)
		deallocate (rmnb,sqgb,sqgib,bmodb,grrb,grtb,gttb,grrojb,grtojb,gttojb,jbgrrb,jbgrtb,jbgttb)
		deallocate (sqgeq,sqgieq,bmodeq,bsteq,grreq,grteq,gtteq,grrojeq,grtojeq,gttojeq,jbgrreq,jbgrteq,jbgtteq, &
			    lplreq,lplteq,djrojeq,djtojeq,djzojeq,jsqeq,omdreq,omdteq,omdzeq,dbsjtojeq,dbsjzojeq,dbsjtbjeq,dgttreq, &
			    dgrrteq,dgrtteq,dgttteq,dgrrzeq,dgrtzeq,dgttzeq,dgrtpeq,dgttpeq,bsgrteq,bsgtteq,bsqeq,bsqgtteq, &
			    lplrreq,lplrteq,lplrzeq,lpltreq,lpltteq,lpltzeq,lplzreq,lplzteq,lplzzeq,bsjgrteq,bsjgtteq, &
			    bstgeq,bsqgeq,sb1,sb2,sb3,sb4,sb5,sb6)

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
