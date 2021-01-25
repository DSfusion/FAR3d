	subroutine energy(i)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: i,j,l,l1,l2,mjp
		real(IDP) :: denom,scmn,scmx,sclmn,sclmx,scnorm
		character(len=1) :: t
		character(len=12) :: confil
		real(IDP), dimension(ldim) :: gampsi,gamalpha,gamphi,gampr
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

!  vr up
		call dbydth(sc2,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr(sc3,phi,0.0_IDP,1.0_IDP,0)
		
		mjp=mj+1
		do l=1,lmaxn
			sc1(:,l)=jbgrr(:,leq0)*sc2(:,l)**2+jbgtt(:,leq0)*sc3(:,l)**2
			xinte=r*denseq*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ekenc(l,i))
			if (l == l0) ekenc(l,i)=2.*ekenc(l,i)
		end do

		call eqtodyn(sc8,jbgrr,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,jbgrt,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,-1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,jbgtt,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,-1,1.0_IDP,1.0_IDP)
		do l=1,lmaxn
			sc1(:,l)=sc2(:,l)*sc5(:,l)+sc3(:,l)*sc6(:,l)
			xinte=r*denseq*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eke(l,i))
			if (l == l0) eke(l,i)=2.*eke(l,i)
		end do

!  br up
		call dbydth(sc2,psi,1,0.0_IDP,-1.0_IDP,0)
!  bth up
		call dbydr(sc3,psi,0.0_IDP,1.0_IDP,0)

		do l=1,lmaxn
			sc1(:,l)=grroj(:,leq0)*sc2(:,l)**2+gttoj(:,leq0)*sc3(:,l)**2
			xinte=r*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,emenc(l,i))
			if (l == l0) emenc(l,i)=2.*emenc(l,i)
		end do

		call eqtodyn(sc8,grroj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,grtoj,0.0_IDP,1.0_IDP)
		call mult(sc5,sc8,-1,sc3,1,1.0_IDP,1.0_IDP)
		call mult(sc6,sc8,-1,sc2,-1,0.0_IDP,1.0_IDP)
		call eqtodyn(sc8,gttoj,0.0_IDP,1.0_IDP)
		call mult(sc6,sc8,1,sc3,1,1.0_IDP,1.0_IDP)
		do l=1,lmaxn
			sc1(:,l)=sc2(:,l)*sc5(:,l)+sc3(:,l)*sc6(:,l)
			xinte=r*sc1(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eme(l,i))
			if (l == l0) eme(l,i)=2.*eme(l,i)
		end do

		call eqtodyn(sc8,sqg,0.0_IDP,1.0_IDP)
		call mult(sc2,sc8,1,vprlf,-1,0.0_IDP,1.0_IDP)
		do l=1,lmaxn
			xinte=r*nfeq*sqg(:,leq0)*vprlf(:,l)**2
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,eprnc(l,i))
			xinte=r*nfeq*sc2(:,l)*vprlf(:,l)
			call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,epr(l,i))
		end do

		if (alpha_on == 1) then
			call mult(sc2,sc8,1,vprlalp,-1,0.0_IDP,1.0_IDP)
			do l=1,lmaxn
				xinte=r*nalpeq*sqg(:,leq0)*vprlalp(:,l)**2
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealpnc(l,i))
				xinte=r*nalpeq*sc2(:,l)*vprlalp(:,l)
				call quadq(mjp,r,xinte,r(mj),bbe,cce,dde,ealp(l,i))
			end do
		end if

		if (i == 1) return

		write(6,'(/"energy:numrun=",2a2,a1,",numruno=",2a2,a1,",nstep=",i10,",time=",1pe12.5,",dt=",1pe12.5)') numrun,numruno, &
				nstep,time,dt
		do l=1,lmaxn
			l1=lln(l)
			if (signl(l1) < 0) cycle
			l2=0
			if (signl(l1) > 0 .and. lo(l) > 0) l2=lln(lo(l))
			gampsi(l1)=0.0_IDP
			if (l2 == 0) then
				denom=dt*(emenc(l1,1)+emenc(l1,2))
				if (denom /= 0.0_IDP) gampsi(l1)=(emenc(l1,2)-emenc(l1,1))/denom
			else
				denom=dt*(emenc(l1,1)+emenc(l2,1)+emenc(l1,2)+emenc(l2,2))
				if (denom /= 0.0_IDP) gampsi(l1)=(emenc(l1,2)+emenc(l2,2)-emenc(l1,1)-emenc(l2,1))/denom
			end if
			gamphi(l1)=0.0_IDP
			if (l2 == 0) then
				denom=dt*(ekenc(l1,1)+ekenc(l1,2))
				if (denom /= 0.0_IDP) gamphi(l1)=(ekenc(l1,2)-ekenc(l1,1))/denom
			else
				denom=dt*(ekenc(l1,1)+ekenc(l2,1)+ekenc(l1,2)+ekenc(l2,2))
				if (denom /= 0.0_IDP) gamphi(l1)=(ekenc(l1,2)+ekenc(l2,2)-ekenc(l1,1)-ekenc(l2,1))/denom
			end if
			gampr(l1)=0.0_IDP
			if (l2 == 0) then
				denom=dt*(eprnc(l1,1)+eprnc(l1,2))
				if (denom /= 0.0_IDP) gampr(l1)=(eprnc(l1,2)-eprnc(l1,1))/denom
			else
				denom=dt*(eprnc(l1,1)+eprnc(l2,1)+eprnc(l1,2)+eprnc(l2,2))
				if (denom /= 0.0_IDP) gampr(l1)=(eprnc(l1,2)+eprnc(l2,2)-eprnc(l1,1)-eprnc(l2,1))/denom
			end if
		end do

		if (alpha_on == 1) then
			do l=1,lmaxn
				l1=lln(l)
				if (signl(l1) < 0) cycle
				l2=0
				if (signl(l1) > 0 .and. lo(l) > 0) l2=lln(lo(l))
				gamalpha(l1)=0.0_IDP
				if (l2 == 0) then
					denom=dt*(ealpnc(l1,1)+ealpnc(l1,2))
					if (denom /= 0.0_IDP) gamalpha(l1)=(ealpnc(l1,2)-ealpnc(l1,1))/denom
				else
					denom=dt*(ealpnc(l1,1)+ealpnc(l2,1)+ealpnc(l1,2)+ealpnc(l2,2))
					if (denom /= 0.0_IDP) gamalpha(l1)=(ealpnc(l1,2)+ealpnc(l2,2)-ealpnc(l1,1)-ealpnc(l2,1))/denom
				end if
			end do
			write(6,'(/"   l   m/  n        ke        me     vprlf   vprlalp     gamke     gamme    gamvpr    gamalp")')
			do l=1,lmaxn
				l1=lln(l)
				if (signl(l1) < 0) cycle
				write(6,'(3i4,1p6e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l1,2),eps*eps*emenc(l1,2),bet0_f*eprnc(l1,2)/bet0, &
							  bet0_alp*ealpnc(l1,2)/bet0,gamphi(l1),gampsi(l1),gampr(l1),gamalpha(l1)
			end do
		else
			write(6,'(/"   l   m/  n        ke        me     vprlf     gamke     gamme    gamvpr")')
			do l=1,lmaxn
				l1=lln(l)
				if (signl(l1) < 0) cycle
				write(6,'(3i4,1p6e10.3)') l1,mm(l1),nn(l1),eps*eps*ekenc(l1,2),eps*eps*emenc(l1,2),bet0_f*eprnc(l1,2)/bet0, &
							  gamphi(l1),gampsi(l1),gampr(l1)
			end do
		end if

		scmn=0.0_IDP
		scmx=0.0_IDP
		do l=1,lmaxn
			sclmn=minval(phi(:,l))
			sclmx=maxval(phi(:,l))
			scmn=min(scmn,sclmn)
			scmx=max(scmx,sclmx)
		end do
		if(scmx > abs(scmn)) then
			scnorm=scmx
		else
			scnorm=scmn
		end if
		if (scnorm /= 0.0_IDP) then
			xt=xt/scnorm
			yt=yt/scnorm
		end if

		if (nstep == nstep1+maxstp) then
			do l=1,lmaxn
				l1=lln(l)
				if (signl(l1) < 0) cycle
				if (lo(l) /= 0) then
					l2=lln(lo(l))
					ekenc(l1,2)=ekenc(l1,2)+ekenc(l2,2)
					eke(l1,2)=eke(l1,2)+eke(l2,2)
					emenc(l1,2)=emenc(l1,2)+emenc(l2,2)
					eme(l1,2)=eme(l1,2)+eme(l2,2)
					eprnc(l1,2)=eprnc(l1,2)+eprnc(l2,2)
					epr(l1,2)=epr(l1,2)+epr(l2,2)
					ealpnc(l1,2)=ealpnc(l1,2)+ealpnc(l2,2)
					ealp(l1,2)=ealp(l1,2)+ealp(l2,2)
				end if
			end do
			t=char(9)
			scnorm=max(maxval(ekenc),maxval(eke),maxval(emenc),maxval(eme))
			if (scnorm /= 0.0_IDP) then
				ekenc(1:lmaxn,2)=ekenc(1:lmaxn,2)/scnorm
				eke(1:lmaxn,2)=eke(1:lmaxn,2)/scnorm
				emenc(1:lmaxn,2)=emenc(1:lmaxn,2)/scnorm
				eme(1:lmaxn,2)=eme(1:lmaxn,2)/scnorm
				eprnc(1:lmaxn,2)=bet0_f*eprnc(1:lmaxn,2)/(bet0*eps*eps*scnorm)
				epr(1:lmaxn,2)=bet0_f*epr(1:lmaxn,2)/(bet0*eps*eps*scnorm)
				ealpnc(1:lmaxn,2)=bet0_alp*ealpnc(1:lmaxn,2)/(bet0*eps*eps*scnorm)
				ealp(1:lmaxn,2)=bet0_alp*ealp(1:lmaxn,2)/(bet0*eps*eps*scnorm)
			end if
			write(confil,'("spctr_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil)
			if (alpha_on == 1) then
				write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf",a1,"vprlalpnc",a1,"vprlalp")') &
					 (t,l=1,10)
				do l=1,lmaxn
					if (signl(l) < 0) cycle
					write(92,'(i4,2(a1,i4),8(a1,1pe13.6))') l,t,mm(l),t,nn(l),t,ekenc(l,2),t,eke(l,2),t,emenc(l,2),t,eme(l,2), &
						t,eprnc(l,2),t,epr(l,2),t,ealpnc(l,2),t,ealp(l,2)
				end do
			else
				write(92,'("l",a1,"m",a1,"n",a1,"kenc",a1,"ke",a1,"menc",a1,"me",a1,"vprlfnc",a1,"vprlf")') (t,l=1,8)
				do l=1,lmaxn
					if (signl(l) < 0) cycle
					write(92,'(i4,2(a1,i4),6(a1,1pe13.6))') l,t,mm(l),t,nn(l),t,ekenc(l,2),t,eke(l,2),t,emenc(l,2),t,eme(l,2), &
						t,eprnc(l,2),t,epr(l,2)
				end do
			end if
			close(92)
		end if

	end subroutine energy
