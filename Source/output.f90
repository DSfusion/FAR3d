	subroutine far3d_output

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: i,nd,l,m,n,lh,lheq,le,j

		interface
			subroutine eqtodyn(adyn,aeq,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: adyn,aeq
			end subroutine eqtodyn
			subroutine dyntoeq(aeq,adyn,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: aeq,adyn
			end subroutine dyntoeq
			subroutine dbydr(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydr
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydrl(d,a,c1,c2,k,l)
				use param
				implicit none
				integer :: k,l
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydrl
		end interface	

!		Summary of run paramer in the output file "farprt"	
		
!		open (unit=6,file="farprt",status="old",POSITION="APPEND")
		
		write(6,'(//52("*"),2a8,51("*")/)') ndevice
		if (ngeneq == 0) write(6,'(50("*"),"cylindrical geometry",50("*")/)')
		if (ngeneq /= 0) write(6,'(51("*"),"toroidal geometry",51("*")/)')
		if (nonlin == 0) write(6,'(57("*"),"linear",56("*")/)')
		if (nonlin /= 0) write(6,'(55("*"),"nonlinear",55("*")//)')

		write(6,'(/"control:"/"numrun  =     ",2a2,a1," numruno =     ",2a2,a1," numruns =     ",2a2,a1," numvac  =     ",a5, &  
				     " maxstp  =",i10," nstep   =",i10/"ndump   =",i10," nprint  =",i10," lplots  =",i10, &
				     " itime   =",i10," dt0     =",1pe10.3," nonlin  =",i10/)') numrun,numruno,numruns,numvac,maxstp,nstep, &
				ndump,nprint,lplots,itime,dt0,nonlin

		write(6,'("numhist =",11("    ",a6))') (numhist(i),i=1,ihist)

		write(6,'(/"domain:"/"mj      =",i10," lmax    =",i10," leqmax  =",i10," mmin    =",i10," mmax    =",i10, &
				    " nmin    =",i10/"nmax    =",i10," mmineq  =",i10," mmaxeq  =",i10," nmineq  =",i10," nmaxeq  =",i10, &
				    " l0      =",i10/"leq0    =",i10," lhmax   =",i10," lheqmx  =",i10," ni      =",i10," nis     =",i10, &
				    " ne      =",i10/"delta   =",1pe10.3," rc      =",0pf10.7," fti     =",0pf10.7," fte     =",0pf10.7)') &
				mj,lmax,leqmax,mmin,mmax,nmin,nmax,mmineq,mmaxeq,nmineq,nmaxeq,l0,leq0,lhmax,lheqmx,ni,nis,ne,delta,rc,fti,fte

		write(6,'(/"equil:"/"ngeneq  =",i10," ndevice =",2a8/"eps     =",0pf10.7," bet0    =",1pe10.3," bet0_f  =",1pe10.3, &
				   " bet0_al =",1pe10.3," omcy    =",1pe10.3," omcyalp =",1pe10.3/"Adens   =",1pe10.3," Bdens   =",1pe10.3, &
				   " Adensal =",1pe10.3," Bdensal =",1pe10.3)') &				 
				ngeneq,ndevice,eps,bet0,bet0_f,bet0_alp,omcy,omcyalp,Adens,Bdens,Adensalp,Bdensalp

		write(6,'(/"dynamo:"/"etascl  =",1pe10.3," reta    =",0pf10.7," eta0    =",1pe10.3," etalmb  =",1pe10.3," ietaeq  =",i10 &
				    /"s       =",1pe10.3," gamma   =",1pe10.3," stdifp  =",1pe10.3," stdifu  =",1pe10.3," stdifnf =",1pe10.3, &
				    " stdifvf =",1pe10.3/"stdifv  =",1pe10.3," stdifnal=",1pe10.3," stdifval=",1pe10.3/"LcA0    =",1pe10.3, &
				    " LcA1    =",1pe10.3," LcA2    =",1pe10.3," LcA3    =",1pe10.3/"LcA0alp =",1pe10.3, &
				    " LcA1alp =",1pe10.3," LcA2alp =",1pe10.3," LcA3alp =",1pe10.3/"iflr    =",i10," epflr   =",i10, &
				    " ieldamp =",i10," twofl   =",i10," alpha   =",i10/"omegar  =",1pe10.3," r_iflr  =",1pe10.3, &
				    " r_epflr =",1pe10.3," r_epalp =",1pe10.3," dpres   =",1pe10.3/"ipert   =",i10, &
				    " pertscl =",1pe10.3," dt      =",1pe10.3," time    =",1pe10.3/)') etascl,reta,eta0, &
				etalmb,ietaeq,s,gamma,stdifp,stdifu,stdifnf,stdifvf,stdifv,stdifnalp,stdifvalp,LcA0,LcA1,LcA2,LcA3,LcA0alp, &
				LcA1alp,LcA2alp,LcA3alp,iflr_on,epflr_on,ieldamp_on,twofl_on,alpha_on,omegar,iflr,r_epflr,r_epflralp,dpres, &
				ipert,pertscl,dt,time

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

		write(6,'(/"  linear matrix bldg. blocks: lmax0,ll0-",i5/("                                        ",12i5))') lmax0, &
			  (ll0(i),i=1,lmax0)
		write(6,'(/"  linear matrix mode entries: lmaxn,lln-",i5/("                                        ",12i5))') lmaxn, &
			  (lln(i),i=1,lmaxn)

!		close(6)

	end subroutine far3d_output
