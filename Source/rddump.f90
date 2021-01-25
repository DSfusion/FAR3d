	subroutine rddump

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

		integer :: i,l,j,lmaxo,m,n,ihisto,l1
		integer, dimension(ldim) :: mmo,nno
		integer, dimension(ldim) :: lmap
		real(IDP) :: pertsclo
		real(IDP), dimension(ldim) :: widthio,gammaio
		character(len=2), dimension(3) :: numrunp
		character(len=8) :: confil		

!		The data of the previous run is read to continue the simulation		

		confil="fs"//numruno(1)//numruno(2)//numruno(3)		
		open(unit=8,file=confil,status='unknown',form='unformatted',POSITION="REWIND")

		read(8) ihisto,numruno,numrunp,numruns,(numhist(i),i=1,ihisto),nstep,maxstp, &
				ndump,nprint,lplots,itime,dt0,nstep1,nonlin,m0dy,nocpl
		read(8) mj,lmaxo,leqmax,mjm1,mjm2,leqdim,ldim,jdim
		read(8) (r(j),j=0,mj),(mmo(l),l=1,lmaxo),(nno(l),l=1,lmaxo),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &
				(dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
				
		if (lmax == 0) lmax=lmaxo
		do l=1,lmaxo
			if (mm(l) /= 0 .or. nn(l) /= 0) exit
			mm(l)=mmo(l)
			nn(l)=nno(l)
		end do

		lmap=0
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			do l1=1,lmaxo
				if (mmo(l1) == m .and. nno(l1) == n) exit
			end do
			if (l1 <= lmaxo) lmap(l1)=l
		end do

		allocate (sqgi(0:mj,0:leqmax))
		allocate (sqg(0:mj,0:leqmax))
		allocate (bst(0:mj,0:leqmax))
		allocate (grr(0:mj,0:leqmax))
		allocate (grt(0:mj,0:leqmax))
		allocate (gtt(0:mj,0:leqmax))
		allocate (grroj(0:mj,0:leqmax))
		allocate (grtoj(0:mj,0:leqmax))
		allocate (gttoj(0:mj,0:leqmax))		
		allocate (bmod(0:mj,0:leqmax))				
		allocate (jbgrr(0:mj,0:leqmax))
		allocate (jbgrt(0:mj,0:leqmax))
		allocate (jbgtt(0:mj,0:leqmax))
		allocate (lplr(0:mj,0:leqmax))
		allocate (lplt(0:mj,0:leqmax))		
		allocate (lplz(0:mj,0:leqmax))	
		allocate (lplr_r(0:mj,0:leqmax))
		allocate (lplt_r(0:mj,0:leqmax))		
		allocate (djroj(0:mj,0:leqmax))	
		allocate (djtoj(0:mj,0:leqmax))	
		allocate (djzoj(0:mj,0:leqmax))	
		allocate (omdr(0:mj,0:leqmax))
		allocate (omdt(0:mj,0:leqmax))
		allocate (omdz(0:mj,0:leqmax))
		
		allocate (dbsjtoj(0:mj,0:leqmax))	
		allocate (dbsjzoj(0:mj,0:leqmax))
		allocate (dbsjtbj(0:mj,0:leqmax))
		allocate (dgttr(0:mj,0:leqmax))
		allocate (dgrrt(0:mj,0:leqmax))
		allocate (dgrtt(0:mj,0:leqmax))
		allocate (dgttt(0:mj,0:leqmax))
		allocate (dgrrz(0:mj,0:leqmax))
		allocate (dgrtz(0:mj,0:leqmax))
		allocate (dgttz(0:mj,0:leqmax))
		allocate (dgrtp(0:mj,0:leqmax))
		allocate (dgttp(0:mj,0:leqmax))
		allocate (jsq(0:mj,0:leqmax))
		allocate (bsgrt(0:mj,0:leqmax))
		allocate (bsgtt(0:mj,0:leqmax))
		allocate (bsq(0:mj,0:leqmax))
		allocate (bsqgtt(0:mj,0:leqmax))
		allocate (lplrr(0:mj,0:leqmax))
		allocate (lplrt(0:mj,0:leqmax))
		allocate (lplrz(0:mj,0:leqmax))
		allocate (lpltt(0:mj,0:leqmax))
		allocate (lpltz(0:mj,0:leqmax))
		allocate (lplzz(0:mj,0:leqmax))
					
        if(ieldamp_on .eq. 1 .and. ext_prof .eq. 1) then		
		  allocate (eildrr(0:mj,0:leqmax))		
		  allocate (eildrt(0:mj,0:leqmax))		
		  allocate (eildrz(0:mj,0:leqmax))			
		  allocate (eildtt(0:mj,0:leqmax))		
		  allocate (eildtz(0:mj,0:leqmax))		
		  allocate (eildzz(0:mj,0:leqmax))			
		  allocate (eildr(0:mj,0:leqmax))		
		  allocate (eildt(0:mj,0:leqmax))		
		  allocate (eildz(0:mj,0:leqmax))			  			  
		end if	

        if(Trapped_on .eq. 1) then	 
		allocate (omdrprp(0:mj,0:leqmax))
		allocate (omdtprp(0:mj,0:leqmax))
		allocate (omdzprp(0:mj,0:leqmax))	
		end if	
		
		sqgi=0.0_IDP
		sqg=0.0_IDP
		bst=0.0_IDP
		grr=0.0_IDP
		grt=0.0_IDP
		gtt=0.0_IDP
		grroj=0.0_IDP
		grtoj=0.0_IDP
		gttoj=0.0_IDP
		bmod=0.0_IDP
		jbgrr=0.0_IDP
		jbgrt=0.0_IDP
		jbgtt=0.0_IDP
		lplr=0.0_IDP
		lplt=0.0_IDP
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
		bsgrt=0.0_IDP
		bsgtt=0.0_IDP
		bsq=0.0_IDP
		bsqgtt=0.0_IDP
		lplrr=0.0_IDP
		lplrt=0.0_IDP
		lplrz=0.0_IDP
		lpltt=0.0_IDP
		lpltz=0.0_IDP
		lplzz=0.0_IDP

        if(Trapped_on .eq. 1) then
		omdrprp=0.0_IDP
		omdtprp=0.0_IDP
		omdzprp=0.0_IDP
		end if	

        if(ieldamp_on .eq. 1 .and. ext_prof .eq. 1) then		
		  eildrr=0.0_IDP
		  eildrt=0.0_IDP	
		  eildrz=0.0_IDP	
		  eildtt=0.0_IDP	
		  eildtz=0.0_IDP	
		  eildzz=0.0_IDP	
		  eildr=0.0_IDP	
		  eildt=0.0_IDP	
		  eildz=0.0_IDP		  
		end if	 				

		read(8) ni,nis,ne,delta,rc,fti,fte
		read(8) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(psieq(j),j=0,mj),(chieq(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
				(cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=1,mj),(vfova(j),j=1,mj), &
                                (vtherm_ion(j),j=1,mj),(vzt_eq(j),j=1,mj)

                if(alpha_on .eq. 1) then
		  read(8) (nalpeq(j),j=1,mj),(valphaova(j),j=1,mj)
		end if	

		read(8) ngeneq,ndevice(1),ndevice(2),eps,bet0,(rs(lmap(l)),l=1,lmaxo),leq,mjeq,nfp  

		read(8) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
				 ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
				 ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
		 
				 ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
				 ((lplr(j,l),j=0,mj),l=1,leqmax),((lplt(j,l),j=0,mj),l=1,leqmax),((lplz(j,l),j=0,mj),l=1,leqmax), &	
				 ((djroj(j,l),j=0,mj),l=1,leqmax),((djtoj(j,l),j=0,mj),l=1,leqmax),((djzoj(j,l),j=0,mj),l=1,leqmax), &				 
				 ((omdr(j,l),j=0,mj),l=1,leqmax),((omdt(j,l),j=0,mj),l=1,leqmax),((omdz(j,l),j=0,mj),l=1,leqmax), &	
				 ((dbsjtoj(j,l),j=0,mj),l=1,leqmax),((dbsjzoj(j,l),j=0,mj),l=1,leqmax),((dbsjtbj(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttr(j,l),j=0,mj),l=1,leqmax),((dgrrt(j,l),j=0,mj),l=1,leqmax),((dgrtt(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttt(j,l),j=0,mj),l=1,leqmax),((dgrrz(j,l),j=0,mj),l=1,leqmax),((dgrtz(j,l),j=0,mj),l=1,leqmax), &
				 ((dgttz(j,l),j=0,mj),l=1,leqmax),((dgrtp(j,l),j=0,mj),l=1,leqmax),((dgttp(j,l),j=0,mj),l=1,leqmax), &
				 ((jsq(j,l),j=0,mj),l=1,leqmax),((bsgrt(j,l),j=0,mj),l=1,leqmax),((bsgtt(j,l),j=0,mj),l=1,leqmax), &
				 ((bsq(j,l),j=0,mj),l=1,leqmax),((bsqgtt(j,l),j=0,mj),l=1,leqmax),((lplrr(j,l),j=0,mj),l=1,leqmax), &
				 ((lplrt(j,l),j=0,mj),l=1,leqmax),((lplrz(j,l),j=0,mj),l=1,leqmax), &
				 ((lpltt(j,l),j=0,mj),l=1,leqmax),((lpltz(j,l),j=0,mj),l=1,leqmax), &
				 ((lplzz(j,l),j=0,mj),l=1,leqmax),((bmod(j,l),j=0,mj),l=1,leqmax), &
				 ((lplr_r(j,l),j=0,mj),l=1,leqmax),((lplt_r(j,l),j=0,mj),l=1,leqmax)	
				 
        if(ieldamp_on .eq. 1) then
		  read(8)  ((eildrr(j,l),j=0,mj),l=1,leqmax),((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &					 
				   ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax),((eildzz(j,l),j=0,mj),l=1,leqmax), &		
				   ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax),((eildz(j,l),j=0,mj),l=1,leqmax)		
		end if	

        if(Trapped_on .eq. 1) then
		  read(8)  ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax),((omdzprp(j,l),j=0,mj),l=1,leqmax)
		end if	
				 
		read(8) (eta(j),j=0,mj),dt,time

		read(8) ((psi(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((phi(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((pr(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((nf(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((vprlf(j,lmap(l)),j=0,mj),l=1,lmaxo)
		read(8) ((vthprlf(j,lmap(l)),j=0,mj),l=1,lmaxo)		
        if(alpha_on .eq. 1) then 
		  read(8) ((nalp(j,lmap(l)),j=0,mj),l=1,lmaxo)	
		  read(8) ((vprlalp(j,lmap(l)),j=0,mj),l=1,lmaxo)
        end if		 
 
		read(8) etascl,reta,eta0,etalmb,ietaeq,stdifp,stdifu,stdifv,stdifnf,stdifvf,stdifnalp,stdifvalp, &
		        s,gamma,ipert,(widthio(l),l=1,lmaxo),(gammaio(l),l=1,lmaxo),pertsclo,norm_eildump

		rewind(8)
		close(unit=8)
		
!		open (unit=6,file="farprt",status="old",POSITION="APPEND")		
		write(6,'("  rddump: have read fs",2a2,a1)') (numruno(i),i=1,3)
!		close(6)

	end subroutine rddump
