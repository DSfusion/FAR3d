	subroutine energy(i)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: i,j,l,lmaxe,mjp
		real(IDP) :: denom,scnorm
		character(len=1) :: t
		character(len=12) :: confil
		real(IDP), dimension(ldim,2) :: epsi,ephi,epr,eprnc,ekenc,eke,emenc,eme
		real(IDP), dimension(ldim) :: gampsi,gamchi,gamlamb,gamalpha,gamphi,gampr
		real(IDP), dimension(jdim) :: rint,xint,bb,cc,dd
		real(IDP), dimension(0:jdim) :: xinte,bbe,cce,dde
		
		interface
			subroutine dbydth(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydth
			subroutine dbydr(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydr
			subroutine dbydzt(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzt
			subroutine eqtodyn(adyn,aeq,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: adyn,aeq
			end subroutine eqtodyn
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine mult
			subroutine quadq(n,x,y,x0,b,c,d,result)
				use param
				implicit none
				integer :: n
				real(IDP) :: x0,result
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine quadq
		end interface

!		Caluclation of the system energy and growth rates 		
		
		epsi(:,i)=0.0_IDP 
		ephi(:,i)=0.0_IDP 
		epr(:,i)=0.0_IDP 
		eprnc(:,i)=0.0_IDP 
		ekenc(:,i)=0.0_IDP 
		eke(:,i)=0.0_IDP 
		emenc(:,i)=0.0_IDP 
		eme(:,i)=0.0_IDP 
		
		lmaxe=lmax
		if (nonlin == 0) lmaxe=lmaxn

		do l=1,lmaxe
			if (mm(l) == 0) then
				psi(0,l)=r(2)**2*psi(1,l)/(r(2)**2-r(1)**2)-r(1)**2*psi(2,l)/(r(2)**2-r(1)**2)
				phi(0,l)=r(2)**2*phi(1,l)/(r(2)**2-r(1)**2)-r(1)**2*phi(2,l)/(r(2)**2-r(1)**2)
			end if
			do j=1,mj
				rint(j)=(r(j)+r(j-1))*.5
				xint(j)=rint(j)*(((psi(j,l)-psi(j-1,l))/(r(j)-r(j-1)))**2+(mm(l)/rint(j)*.5*(psi(j,l)+psi(j-1,l)))**2)
			end do
			call quadq(mj,rint,xint,rint(mj),bb,cc,dd,epsi(l,i))
			do j=1,mj
				xint(j)=rint(j)*(((phi(j,l)-phi(j-1,l))/(r(j)-r(j-1)))**2+(mm(l)/rint(j)*.5*(phi(j,l)+phi(j-1,l)))**2)
			end do
			call quadq(mj,rint,xint,rint(mj),bb,cc,dd,ephi(l,i))
		end do

!  vr up
		call dbydth(sc2,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr(sc3,phi,0.0_IDP,1.0_IDP,0)
		
		mjp=mj+1
		do l=1,lmaxe
			sc1(:,l)=jbgrr(:,leq0)*sc2(:,l)**2+jbgtt(:,leq0)*sc3(:,l)**2
			xinte=r*denseq*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ekenc(l,i))
		end do

		call eqtodyn(sc8,jbgrr,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,jbgrt,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,-1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,jbgtt,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,-1,1.0_IDP,1.0_IDP)
		sc8=0.0_IDP
		do l=1,lmaxe
			sc1(:,l)=sc2(:,l)*sc5(:,l)+sc3(:,l)*sc6(:,l)
			xinte=r*denseq*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eke(l,i))
		end do

!  br up
		call dbydth(sc2,psi,1,0.0_IDP,-s,0)
!  bth up
		call dbydr(sc3,psi,0.0_IDP,s,0)

		do l=1,lmaxe
			sc1(:,l)=grroj(:,leq0)*sc2(:,l)**2+gttoj(:,leq0)*sc3(:,l)**2
			xinte=r*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,emenc(l,i))
		end do

		call eqtodyn(sc8,grroj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,grtoj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,gttoj,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,1,1.0_IDP,1.0_IDP)
		sc8=0.0_IDP
		do l=1,lmaxe
			sc1(:,l)=sc2(:,l)*sc5(:,l)+sc3(:,l)*sc6(:,l)
			xinte=r*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eme(l,i))
		end do

		call eqtodyn(sc8,sqg,0.0_IDP,1.0_IDP)
		call mult(sc2,sc8,1,pr,1,0.0_IDP,1.0_IDP)
		do l=1,lmaxe
			xinte=r*sqg(:,leq0)*pr(:,l)**2
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eprnc(l,i))
			xinte=r*sc2(:,l)*pr(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,epr(l,i))
		end do

		if (i == 1) return

		open (unit=6,file="farprt",status="old",POSITION="APPEND")		
		write(6,'(/"energy:numrun=",2a2,a1,",numruno=",2a2,a1,",nstep=",i10,",time=",1pe12.5,",dt=",1pe12.5)') numrun,numruno, &
				nstep,time,dt
		close(6)

		do l=1,lmaxe
			gampsi(l)=0.0_IDP
			denom=dt*(emenc(l,1)+emenc(l,2))
			if (denom /= 0.0_IDP) gampsi(l)=(emenc(l,2)-emenc(l,1))/denom
			gamchi(l)=0.0_IDP
			denom=dt*(eme(l,1)+eme(l,2))
			if (denom /= 0.0_IDP) gamchi(l)=(eme(l,2)-eme(l,1))/denom
			gamphi(l)=0.0_IDP
			denom=dt*(ekenc(l,1)+ekenc(l,2))
			if (denom /= 0.0_IDP) gamphi(l)=(ekenc(l,2)-ekenc(l,1))/denom
			gamlamb(l)=0.0_IDP
			denom=dt*(eke(l,1)+eke(l,2))
			if (denom /= 0.0_IDP) gamlamb(l)=(eke(l,2)-eke(l,1))/denom
			gamalpha(l)=0.0_IDP
			denom=dt*(eprnc(l,1)+eprnc(l,2))
			if (denom /= 0.0_IDP) gamalpha(l)=(eprnc(l,2)-eprnc(l,1))/denom
			gampr(l)=0.0_IDP
			denom=dt*(epr(l,1)+epr(l,2))
			if (denom /= 0.0_IDP) gampr(l)=(epr(l,2)-epr(l,1))/denom
		end do

		open (unit=6,file="farprt",status="old",POSITION="APPEND")		
		write(6,'(/"  l  m/ n      kenc        ke      menc        me      prnc        pr   gamkenc     gamke   gammenc     ", &
				  "gamme   gamprnc     gampr")')
		do l=1,lmaxe
		   write(6,'(3i3,1p12e10.3)') l,mm(l),nn(l),ekenc(l,2),eke(l,2),emenc(l,2),eme(l,2),eprnc(l,2),epr(l,2), &
					gamphi(l)/s,gamlamb(l)/s,gampsi(l)/s,gamchi(l)/s,gamalpha(l)/s,gampr(l)/s
		end do
		close(6)

		if (nstep == nstep1+maxstp) then
			t=char(9)
			scnorm=max(maxval(ekenc),maxval(eke),maxval(emenc),maxval(eme))
			if (scnorm /= 0.0_IDP) then
				ekenc(1:lmaxe,2)=ekenc(1:lmaxe,2)/scnorm
				eke(1:lmaxe,2)=eke(1:lmaxe,2)/scnorm
				emenc(1:lmaxe,2)=emenc(1:lmaxe,2)/scnorm
				eme(1:lmaxe,2)=eme(1:lmaxe,2)/scnorm
				eprnc(1:lmaxe,2)=eprnc(1:lmaxe,2)/scnorm
				epr(1:lmaxe,2)=epr(1:lmaxe,2)/scnorm
			end if
			write(confil,'("spctr_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil)
			write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"prnc",a1,"pr")') (t,l=1,6)
			do l=1,lmaxe
				write(92,'(i4,2(a1,i4),6(a1,1pe13.6))') l,t,mm(l),t,nn(l),t,ekenc(l,2),t,eke(l,2),t,emenc(l,2),t,eme(l,2), &
					t,eprnc(l,2),t,epr(l,2)
			end do
			close(92)
		end if

	end subroutine energy