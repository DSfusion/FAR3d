	subroutine wrdump

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

		integer :: i,l,j
		character(len=8) :: confil

!		The data of the run is stored in case there is a continuation run	
		
		confil="fs"//numrun(1)//numrun(2)//numrun(3)
                if (nonlin == 1) then
                  open(unit=7,file=confil,status='new',form='unformatted')
                else
		  open(unit=7,file=confil,status='unknown',form='unformatted')
                end if

		write(7) ihist,numrun,numruno,numruns,(numhist(i),i=1,ihist),nstep,maxstp, &
				 ndump,nprint,lplots,itime,dt0,nstep1,nonlin,m0dy,nocpl	

		write(7) mj,lmax,leqmax,mjm1,mjm2,leqdim,ldim,jdim
		write(7) (r(j),j=0,mj),(mm(l),l=1,lmax),(nn(l),l=1,lmax),(mmeq(l),l=1,leqmax),(nneq(l),l=1,leqmax),(rinv(j),j=0,mj), &		
				 (dc1m(j),j=1,mj),(dc1p(j),j=1,mj),(dc2m(j),j=1,mj),(dc2p(j),j=1,mj),(del2cm(j),j=1,mj),(del2cp(j),j=1,mj)
	
		write(7) ni,nis,ne,delta,rc,fti,fte
		write(7) (qq(j),j=0,mj),(qqinv(j),j=0,mj),(psieq(j),j=0,mj),(chieq(j),j=0,mj),(preq(j),j=0,mj),(feq(j),j=0,mj), &
				 (cureq(j),j=0,mj),(denseq(j),j=0,mj),(teeq(j),j=0,mj),(nfeq(j),j=1,mj),(vfova(j),j=1,mj), &
                                 (vtherm_ion(j),j=1,mj),(vzt_eq(j),j=1,mj)
				 
                if(alpha_on .eq. 1) then
		  write(7) (nalpeq(j),j=1,mj),(valphaova(j),j=1,mj)
		end if	

		write(7) ngeneq,ndevice(1),ndevice(2),eps,bet0,(rs(l),l=1,lmax),leq,mjeq,nfp

		write(7) ((sqgi(j,l),j=0,mj),l=1,leqmax),((sqg(j,l),j=0,mj),l=1,leqmax),((bst(j,l),j=0,mj),l=1,leqmax), &
				 ((grr(j,l),j=0,mj),l=1,leqmax),((grt(j,l),j=0,mj),l=1,leqmax),((gtt(j,l),j=0,mj),l=1,leqmax), &
				 ((grz(j,l),j=0,mj),l=1,leqmax),((gtz(j,l),j=0,mj),l=1,leqmax),((gzz(j,l),j=0,mj),l=1,leqmax), &
				 ((grroj(j,l),j=0,mj),l=1,leqmax),((grtoj(j,l),j=0,mj),l=1,leqmax),((gttoj(j,l),j=0,mj),l=1,leqmax), &
				 ((grzoj(j,l),j=0,mj),l=1,leqmax),((gtzoj(j,l),j=0,mj),l=1,leqmax),((gzzoj(j,l),j=0,mj),l=1,leqmax), &
				 ((grrup(j,l),j=0,mj),l=1,leqmax),((grtup(j,l),j=0,mj),l=1,leqmax),((grzup(j,l),j=0,mj),l=1,leqmax), &
				 ((gttup(j,l),j=0,mj),l=1,leqmax),((gtzup(j,l),j=0,mj),l=1,leqmax),((gzzup(j,l),j=0,mj),l=1,leqmax), &				 
				 ((jbgrr(j,l),j=0,mj),l=1,leqmax),((jbgrt(j,l),j=0,mj),l=1,leqmax),((jbgtt(j,l),j=0,mj),l=1,leqmax), &
                                 ((jbgrz(j,l),j=0,mj),l=1,leqmax),((jbgtz(j,l),j=0,mj),l=1,leqmax), &
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
				 ((sqgdroj(j,l),j=0,mj),l=1,leqmax),((sqgdthoj(j,l),j=0,mj),l=1,leqmax),((sqgdztoj(j,l),j=0,mj),l=1,leqmax), &
				 ((sqgdthojbst(j,l),j=0,mj),l=1,leqmax),((sqgdztojbst(j,l),j=0,mj),l=1,leqmax), &
				 ((sqgibmodith(j,l),j=0,mj),l=1,leqmax),((sqgibmodizt(j,l),j=0,mj),l=1,leqmax) 				 
 	                
        if(ieldamp_on .eq. 1) then
		  write(7) ((eildrr(j,l),j=0,mj),l=1,leqmax),((eildrt(j,l),j=0,mj),l=1,leqmax),((eildrz(j,l),j=0,mj),l=1,leqmax), &					 
				   ((eildtt(j,l),j=0,mj),l=1,leqmax),((eildtz(j,l),j=0,mj),l=1,leqmax),((eildzz(j,l),j=0,mj),l=1,leqmax), &		
				   ((eildr(j,l),j=0,mj),l=1,leqmax),((eildt(j,l),j=0,mj),l=1,leqmax),((eildz(j,l),j=0,mj),l=1,leqmax)	
		end if		

        if(Trapped_on .eq. 1) then
		  write(7)  ((omdrprp(j,l),j=0,mj),l=1,leqmax),((omdtprp(j,l),j=0,mj),l=1,leqmax),((omdzprp(j,l),j=0,mj),l=1,leqmax)
		end if	
		
		write(7) (eta(j),j=0,mj),dt,time

		write(7) ((psi(j,l),j=0,mj),l=1,lmax)
		write(7) ((phi(j,l),j=0,mj),l=1,lmax)
		write(7) ((pr(j,l),j=0,mj),l=1,lmax)
		write(7) ((nf(j,l),j=0,mj),l=1,lmax)
		write(7) ((vprlf(j,l),j=0,mj),l=1,lmax)
		write(7) ((vthprlf(j,l),j=0,mj),l=1,lmax)	
        if(alpha_on .eq. 1) then 
		  write(7) ((nalp(j,l),j=0,mj),l=1,lmax)
		  write(7) ((vprlalp(j,l),j=0,mj),l=1,lmax)
        end if		
 
		write(7) etascl,reta,eta0,etalmb,ietaeq,stdifp,stdifu,stdifv,stdifnf,stdifvf,stdifnalp,stdifvalp, &
		         s,gamma,ipert,(widthi(l),l=1,lmax),(gammai(l),l=1,lmax),pertscl,norm_eildump

		rewind(7)
		close(unit=7)
		
		open (unit=6,file="farprt",status="old",POSITION="APPEND")	
		write (6,'("  wrdump: have written fs",2a2,a1)') (numrun(i),i=1,3)
		close(6)

	end subroutine wrdump