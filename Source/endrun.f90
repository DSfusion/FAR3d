	subroutine endrun

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: j,l,l1,l2,lvrmax,jvrmax,lwrt,lwrto2,iwrt
		integer, dimension(1) :: lvrmx,jvrmx
		integer, dimension(:), allocatable :: mwrt,nwrt
		real(IDP), dimension(:), allocatable :: vrsqmx
		real(IDP) :: theta,csth,snth,scnorm
		character(len=1) :: tb
		character(len=12) :: confil
		character(len=32) :: formatt='("r",127(a1,a1,i3,"/",i2))'
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
		end interface

!		Calculates the radial and poloidal components of the velocity and magnetic field

		tb=char(9)

!  vr up
		call dbydth(sc1,phi,-1,0.0_IDP,-1.0_IDP,0)
!  vth up
		call dbydr(sc2,phi,0.0_IDP,1.0_IDP,0)

		allocate (vrsqmx(lmaxn))
		vrsqmx=0.0
		sc3=0.0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			do j=0,mj
				sc3(j,l)=sc1(j,l1)*sc1(j,l1)+sc1(j,l2)*sc1(j,l2)
			end do
			vrsqmx(l)=maxval(sc3(:,l))
		end do
		lvrmx=maxloc(vrsqmx)
		lvrmax=lvrmx(1)
		jvrmx=maxloc(sc3(:,lvrmax))
		jvrmax=jvrmx(1)

		theta=atan2(sc1(jvrmax,lln(lo(lvrmax))),sc1(jvrmax,lln(lvrmax)))
		csth=cos(theta)
		snth=sin(theta)
		scnorm=sc1(jvrmax,lln(lvrmax))*csth+sc1(jvrmax,lln(lo(lvrmax)))*snth

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
				nf(:,l)=nf(:,l)/scnorm
				where (abs(nf(:,l)) < 1.e-50_IDP) nf(:,l)=0
				vprlf(:,l)=vprlf(:,l)/scnorm
				where (abs(vprlf(:,l)) < 1.e-50_IDP) vprlf(:,l)=0
				vthprlf(:,l)=vthprlf(:,l)/scnorm
				where (abs(vthprlf(:,l)) < 1.e-50_IDP) vthprlf(:,l)=0
				sc1(:,l)=eps*sc1(:,l)/scnorm
				where (abs(sc1(:,l)) < 1.e-50_IDP) sc1(:,l)=0
				sc2(:,l)=eps*sc2(:,l)/scnorm
				where (abs(sc2(:,l)) < 1.e-50_IDP) sc2(:,l)=0
			end do
			if (alpha_on == 1) then
				do l=1,lmaxn
					nalp(:,l)=nalp(:,l)/scnorm
					where (abs(nalp(:,l)) < 1.e-50_IDP) nalp(:,l)=0
					vprlalp(:,l)=vprlalp(:,l)/scnorm
					where (abs(vprlalp(:,l)) < 1.e-50_IDP) vprlalp(:,l)=0
				end do
			end if
			xt=xt/scnorm
		end if
		lwrt=min(abs(lplots),lmaxn)
		lwrto2=lwrt/2
		allocate (mwrt(lwrto2),nwrt(lwrto2))
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			mwrt(iwrt)=mm(l1)
			nwrt(iwrt)=nn(l1)
			sc3(:,iwrt)=csth*phi(:,l2)+snth*phi(:,l1)
			sc3(:,iwrt+lwrto2)=snth*phi(:,l2)-csth*phi(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("phi_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*psi(:,l1)+snth*psi(:,l2)
			sc3(:,iwrt+lwrto2)=snth*psi(:,l1)-csth*psi(:,l2)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("psi_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*pr(:,l1)+snth*pr(:,l2)
			sc3(:,iwrt+lwrto2)=snth*pr(:,l1)-csth*pr(:,l2)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("pr_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*uzt(:,l2)+snth*uzt(:,l1)
			sc3(:,iwrt+lwrto2)=snth*uzt(:,l2)-csth*uzt(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("uzt_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*nf(:,l1)+snth*nf(:,l2)
			sc3(:,iwrt+lwrto2)=snth*nf(:,l1)-csth*nf(:,l2)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("nf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*vprlf(:,l2)+snth*vprlf(:,l1)
			sc3(:,iwrt+lwrto2)=snth*vprlf(:,l2)-csth*vprlf(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("vprlf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		if (alpha_on == 1) then
			iwrt=0
			do l=1,lmaxn
				l1=lln(l)
				if (nn(l1) < 0) cycle
				l2=lln(lo(l))
				iwrt=iwrt+1
				sc3(:,iwrt)=csth*nalp(:,l1)+snth*nalp(:,l2)
				sc3(:,iwrt+lwrto2)=snth*nalp(:,l1)-csth*nalp(:,l2)
				if (iwrt == lwrto2) exit
			end do
			write(confil,'("nalp_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=2048)
			write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
			do j=0,mj
				write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
			end do
			close(92)
			iwrt=0
			do l=1,lmaxn
				l1=lln(l)
				if (nn(l1) < 0) cycle
				l2=lln(lo(l))
				iwrt=iwrt+1
				sc3(:,iwrt)=csth*vprlalp(:,l2)+snth*vprlalp(:,l1)
				sc3(:,iwrt+lwrto2)=snth*vprlalp(:,l2)-csth*vprlalp(:,l1)
				if (iwrt == lwrto2) exit
			end do
			write(confil,'("vprlalp_",2a2)') numrun(1),numrun(2)
			open(unit=92,file=confil,recl=2048)
			write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
			do j=0,mj
				write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
			end do
			close(92)
		end if
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*vthprlf(:,l2)+snth*vthprlf(:,l1)
			sc3(:,iwrt+lwrto2)=snth*vthprlf(:,l2)-csth*vthprlf(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("vthprlf_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*sc1(:,l1)+snth*sc1(:,l2)
			sc3(:,iwrt+lwrto2)=snth*sc1(:,l1)-csth*sc1(:,l2)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("vr_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*sc2(:,l2)+snth*sc2(:,l1)
			sc3(:,iwrt+lwrto2)=snth*sc2(:,l2)-csth*sc2(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("vth_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)

!  br up
		call dbydth(sc3,psi,1,0.0_IDP,-eps,0)
!  bth up
		call dbydr(sc4,psi,0.0_IDP,eps,0)

		call eqtodyn(sc8,sqgi,0.0_IDP,1.0_IDP)

		call mult(sc1,sc8,1,sc3,-1,0.0_IDP,1.0_IDP)
		call mult(sc2,sc8,1,sc4,1,0.0_IDP,1.0_IDP)

		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*sc1(:,l2)+snth*sc1(:,l1)
			sc3(:,iwrt+lwrto2)=snth*sc1(:,l2)-csth*sc1(:,l1)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("br_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)
		iwrt=0
		do l=1,lmaxn
			l1=lln(l)
			if (nn(l1) < 0) cycle
			l2=lln(lo(l))
			iwrt=iwrt+1
			sc3(:,iwrt)=csth*sc2(:,l1)+snth*sc2(:,l2)
			sc3(:,iwrt+lwrto2)=snth*sc2(:,l1)-csth*sc2(:,l2)
			if (iwrt == lwrto2) exit
		end do
		write(confil,'("bth_",2a2)') numrun(1),numrun(2)
		open(unit=92,file=confil,recl=2048)
		write(92,formatt) (tb,"R",mwrt(l),nwrt(l),l=1,lwrto2),(tb,"I",mwrt(l),nwrt(l),l=1,lwrto2)
		do j=0,mj
			write(92,formatv) r(j),(tb,sc3(j,l),l=1,lwrt)
		end do
		close(92)

	end subroutine endrun
