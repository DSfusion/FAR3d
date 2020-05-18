	subroutine endrun

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: j,l,lwrt
		real(IDP) :: scmn,scmx,sclmn,sclmx,scnorm
		character(len=1) :: t
		character(len=12) :: confil
		character(len=32) :: formatt='("r",127(a1,i4,"/",i4))'
		character(len=32) :: formatv='(1pe13.6,127(a1,1pe15.8))'

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
			subroutine blockp(tx,itx,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockp
			subroutine b2lxp(tx,itx,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lxp
			subroutine decbt(m,n,a,b,c,ip,ier)
				use param
				implicit none
				integer :: m,n,ier
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
			end subroutine decbt
			subroutine solbt(m,n,a,b,c,y,ip)
				use param
				implicit none
				integer :: m,n
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
				real(IDP), dimension(m,n) :: y
			end subroutine solbt
		end interface

!		Creates the radial and poloidal components of the velocity and magnetic fields		
		
		t=char(9)

!  vr up
		call dbydth(sc1,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr(sc2,phi,0.0_IDP,1.0_IDP,0)

		scmn=0.0_IDP
		scmx=0.0_IDP
		do l=1,lmaxn
			sclmn=minval(sc1(:,l))
			sclmx=maxval(sc1(:,l))
			scmn=min(scmn,sclmn)
			scmx=max(scmx,sclmx)
		end do
		if(scmx > abs(scmn)) then
			scnorm=scmx
		else
			scnorm=scmn
		end if
		if (scnorm /= 0.0_IDP) then
			do l=1,lmaxn
				psi(:,l)=psi(:,l)/scnorm
				where (abs(psi(:,l)) < 1.e-50_IDP) psi(:,l)=0
				phi(:,l)=phi(:,l)/scnorm
				where (abs(phi(:,l)) < 1.e-50_IDP) phi(:,l)=0
				pr(:,l)=pr(:,l)/scnorm
				where (abs(pr(:,l)) < 1.e-50_IDP) pr(:,l)=0
				uzt(:,l)=uzt(:,l)/scnorm
				where (abs(uzt(:,l)) < 1.e-50_IDP) uzt(:,l)=0
				sc1(:,l)=sc1(:,l)/scnorm
				where (abs(sc1(:,l)) < 1.e-50_IDP) sc1(:,l)=0
				sc2(:,l)=sc2(:,l)/scnorm
				where (abs(sc2(:,l)) < 1.e-50_IDP) sc2(:,l)=0
				nf(:,l)=nf(:,l)/scnorm
				where (abs(nf(:,l)) < 1.e-50_IDP) nf(:,l)=0
				vprlf(:,l)=vprlf(:,l)/scnorm
				where (abs(vprlf(:,l)) < 1.e-50_IDP) vprlf(:,l)=0
				vthprlf(:,l)=vthprlf(:,l)/scnorm
				where (abs(vthprlf(:,l)) < 1.e-50_IDP) vthprlf(:,l)=0
				if(alpha_on .eq. 1) then
				  nalp(:,l)=nalp(:,l)/scnorm
				  where (abs(nalp(:,l)) < 1.e-50_IDP) nalp(:,l)=0
				  vprlalp(:,l)=vprlalp(:,l)/scnorm
				  where (abs(vprlalp(:,l)) < 1.e-50_IDP) vprlalp(:,l)=0
				end if
			end do
			xt=xt/scnorm
		end if
		lwrt=abs(lplots)
		write(confil,'("phi_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,phi(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("psi_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,psi(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("pr_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,pr(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("uzt_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,uzt(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("nf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,nf(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("vprlf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,vprlf(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("vthprlf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,vthprlf(j,l),l=1,lwrt)
		end do
		close(92)

		write(confil,'("vr_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,sc1(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("vth_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,sc2(j,l),l=1,lwrt)
		end do
		close(92)

!  br up
		call dbydth(sc3,psi,1,0.0_IDP,-s,0)
!  bth up
		call dbydr(sc4,psi,0.0_IDP,s,0)

		call eqtodyn(sc8,sqgi,0.0_IDP,1.0_IDP)

		call mult(sc1,sc8,1,sc3,-1,0.0_IDP,1.0_IDP)
		call mult(sc2,sc8,1,sc4,1,0.0_IDP,1.0_IDP)

		write(confil,'("br_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,sc1(j,l),l=1,lwrt)
		end do
		close(92)
		write(confil,'("bth_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (t,mm(l),nn(l),l=1,lwrt)
		do j=0,mj
			write(92,formatv) r(j),(t,sc2(j,l),l=1,lwrt)
		end do
		close(92)

	end subroutine endrun