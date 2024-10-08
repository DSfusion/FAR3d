	subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt
		real(IDP) :: coef,x
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(jdim) :: xa
		integer :: ity,iph,ich,l,ln,m,n,l1,l2,l3,l3n,lp,lpn,lpp,ibnd,j,il

!		Operator in lincheck subroutine	
!		Add the model terms to Y=R*X	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		tx = multiplicative factor (with radial and angular dependencies)
!		itz = multiplicative factor sign
!		ieqn = equation index where the term is included
!		ivar = perturbation variable index
!		ith = number of derivation in theta
!		ir = number of devivation in rho
!		izt = number of derivations in zeta
!		coef = coefficient that multiplies tx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		ity=1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
		ity=ity*(-1)**(ith+izt)
		ich=itx*ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do l=1,lmaxn
			ln=l
!			ln=lln(l)
			l2=l+(ivar-1)*lmaxn
			m=mm(lln(l))
			n=nn(lln(l))
			ibnd=1
			if (m == 0) ibnd=2
			do lp=1,lmaxn
				lpn=lp
!				lpn=lln(lp)
				lpp=lp
				if (ich == -1) lpp=lo(lp)
				if (lpp == 0) cycle
				l1=lpp+(ieqn-1)*lmaxn
				do l3=1,leqmax
					l3n=l3
					if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n,ln)
					if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n,ln)
					if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n,ln)
					if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n,ln)
					do il=1,ith
						x=x*m
					end do
					do il=1,izt
						x=x*n
					end do
					if (x == 0.0_IDP) cycle
					x=x*iph
					do j=1,mjm1
						xa(j)=x*(rinv(j))**ith
					end do
					if (ir == 0) then
						do j=1,mjm1
							yt(l1,j)=yt(l1,j)+xa(j)*tx(j,l3n)*xt(l2,j)
						end do
					else if (ir == 1) then
						yt(l1,1)=yt(l1,1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
						yt(l1,1)=yt(l1,1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm2
							yt(l1,j)=yt(l1,j)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
!						if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						else if (ivar == 4) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
							yt(l1,mjm1)=yt(l1,mjm1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
						else
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						end if
					else
						yt(l1,1)=yt(l1,1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
						yt(l1,1)=yt(l1,1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm2
							yt(l1,j)=yt(l1,j)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
!						if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						else if (ivar == 4) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
							yt(l1,mjm1)=yt(l1,mjm1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
						else
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						end if
					end if
				end do
			end do
		end do

	end subroutine b2lx

	subroutine b2lxl(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt
		real(IDP) :: x
		real(IDP), dimension(0:,0:) :: tx
		real(IDP), dimension(:,:) :: coef
		real(IDP), dimension(jdim) :: xa
		integer :: ity,iph,ich,l,ln,m,n,l1,l2,l3,l3n,lp,lpn,lpp,ibnd,j,il

!		Operator in lincheck subroutine	
!		Add the model terms to Y=R*X	
!		Similar to blockj but coef depends on the l.h.s. l-value and r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		tx = multiplicative factor (with radial and angular dependencies)
!		itz = multiplicative factor sign
!		ieqn = equation index where the term is included
!		ivar = perturbation variable index
!		ith = number of derivation in theta
!		ir = number of devivation in rho
!		izt = number of derivations in zeta
!		coef = coefficient (with radial and mode dependecies) that multiplies tx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			

		ity=1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
		ity=ity*(-1)**(ith+izt)
		ich=itx*ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do l=1,lmaxn
			ln=l
!			ln=lln(l)
			l2=l+(ivar-1)*lmaxn
			m=mm(lln(l))
			n=nn(lln(l))
			ibnd=1
			if (m == 0) ibnd=2
			do lp=1,lmaxn
				lpn=lp
!				lpn=lln(lp)
				lpp=lp
				if (ich == -1) lpp=lo(lp)
				if (lpp == 0) cycle
				l1=lpp+(ieqn-1)*lmaxn
				do l3=1,leqmax
					l3n=l3
					if (ity ==  1 .and. itx ==  1) x=cmapp(lpn,l3n,ln)
					if (ity ==  1 .and. itx == -1) x=cmamp(lpn,l3n,ln)
					if (ity == -1 .and. itx == -1) x=cmamm(lpn,l3n,ln)
					if (ity == -1 .and. itx ==  1) x=cmapm(lpn,l3n,ln)
					do il=1,ith
						x=x*m
					end do
					do il=1,izt
						x=x*n
					end do
					if (x == 0.0_IDP) cycle
					x=x*iph
					do j=1,mjm1
						xa(j)=x*coef(j,lp)*(rinv(j))**ith
					end do
					if (ir == 0) then
						do j=1,mjm1
							yt(l1,j)=yt(l1,j)+xa(j)*tx(j,l3n)*xt(l2,j)
						end do
					else if (ir == 1) then
						yt(l1,1)=yt(l1,1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
						yt(l1,1)=yt(l1,1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm2
							yt(l1,j)=yt(l1,j)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
!						if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						else if (ivar == 4) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
							yt(l1,mjm1)=yt(l1,mjm1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
						else
							yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						end if
					else
						yt(l1,1)=yt(l1,1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
						yt(l1,1)=yt(l1,1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm2
							yt(l1,j)=yt(l1,j)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
!						if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
						if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						else if (ivar == 4) then
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
							yt(l1,mjm1)=yt(l1,mjm1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3n)* &
										(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
						else
							yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm1)
							yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3n)*xt(l2,mjm2)
						end if
					end if
				end do
			end do
		end do

	end subroutine b2lxl

	subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: ieqn,ivar,ith,ir,izt
		real(IDP) :: coef,x
		real(IDP), dimension(0:) :: tx
		real(IDP), dimension(jdim) :: xa
		integer :: ity,iph,ich,l,m,n,l1,l2,lp,ibnd,j,il

!		Operator in lincheck subroutine	
!		Add the model terms to Y=R*X	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		tx = multiplicative factor (with only radial dependency)
!		itz = multiplicative factor sign
!		ieqn = equation index where the term is included
!		ivar = perturbation variable index
!		ith = number of derivation in theta
!		ir = number of devivation in rho
!		izt = number of derivations in zeta
!		coef = coefficient that multiplies tx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

		ity=1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
		ity=ity*(-1)**(ith+izt)
		ich=ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do l=1,lmaxn
			l2=l+(ivar-1)*lmaxn
			m=mm(lln(l))
			n=nn(lln(l))
			ibnd=1
			if (m == 0) ibnd=2
			lp=l
			if (ich == -1) lp=lo(l)
			if (lp == 0) cycle
			l1=lp+(ieqn-1)*lmaxn
			x=coef
			do il=1,ith
				x=x*m
			end do
			do il=1,izt
				x=x*n
			end do
			if (x == 0.0_IDP) cycle
			x=x*iph
			do j=1,mjm1
				xa(j)=x*(rinv(j))**ith
			end do
			if (ir == 0) then
				do j=1,mjm1
					yt(l1,j)=yt(l1,j)+xa(j)*tx(j)*xt(l2,j)
				end do
			else if (ir == 1)then
				yt(l1,1)=yt(l1,1)+wt10(1,ibnd)*xa(1)*tx(1)*xt(l2,1)
				yt(l1,1)=yt(l1,1)+wt1p(1,ibnd)*xa(1)*tx(1)*xt(l2,2)
				do j=2,mjm2
					yt(l1,j)=yt(l1,j)+wt10(j,ibnd)*xa(j)*tx(j)*xt(l2,j)
					yt(l1,j)=yt(l1,j)+wt1m(j,ibnd)*xa(j)*tx(j)*xt(l2,j-1)
					yt(l1,j)=yt(l1,j)+wt1p(j,ibnd)*xa(j)*tx(j)*xt(l2,j+1)
				end do
!				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
				if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
					yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
				else if (ivar == 4) then
					yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
					yt(l1,mjm1)=yt(l1,mjm1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)* &
								(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)* &
								(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
				else
					yt(l1,mjm1)=yt(l1,mjm1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
				end if
			else
				yt(l1,1)=yt(l1,1)+wt20(1,ibnd)*xa(1)*tx(1)*xt(l2,1)
				yt(l1,1)=yt(l1,1)+wt2p(1,ibnd)*xa(1)*tx(1)*xt(l2,2)
				do j=2,mjm2
					yt(l1,j)=yt(l1,j)+wt20(j,ibnd)*xa(j)*tx(j)*xt(l2,j)
					yt(l1,j)=yt(l1,j)+wt2m(j,ibnd)*xa(j)*tx(j)*xt(l2,j-1)
					yt(l1,j)=yt(l1,j)+wt2p(j,ibnd)*xa(j)*tx(j)*xt(l2,j+1)
				end do
!				if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
				if ((ivar > 5 .and. ivar < 8) .or. (alpha_on == 1 .and. ivar == 9)) then
					yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
				else if (ivar == 4) then
					yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
					yt(l1,mjm1)=yt(l1,mjm1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)* &
								(r(mj)-r(mjm2))/(r(mjm1)-r(mjm2))*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)* &
								(r(mj)-r(mjm1))/(r(mjm1)-r(mjm2))*xt(l2,mjm2)
				else
					yt(l1,mjm1)=yt(l1,mjm1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm1)
					yt(l1,mjm1)=yt(l1,mjm1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(l2,mjm2)
				end if
			end if
		end do

	end subroutine b2lx0

	subroutine b2lx_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itx,ieqn,ivar,ith,ir,izt,ith_zt
		real(kind=IDP) :: coef,x
		real(kind=IDP), dimension(0:,0:) :: tx
		real(kind=IDP), dimension(jdim) :: xa
		integer :: ity,iph,ich,l,ln,m,n,l1,l2,l3,l3n,lp,lpn,lpp,ibnd,j,il

!		Operator in lincheck subroutine	
!		Add the model terms to Y=R*X	
!		Similar to blockj but including the factor abs(dble(n)-dble(m)*qqinv(j))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		tx = multiplicative factor (with radial and angular dependencies)
!		itz = multiplicative factor sign
!		ieqn = equation index where the term is included
!		ivar = perturbation variable index
!		ith = number of derivation in theta
!		ir = number of devivation in rho
!		izt = number of derivations in zeta
!		coef = coefficient that multiplies tx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			

		ity=1
		if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
		
		iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
		ity=ity*(-1)**(ith+izt)

		ich=itx*ity
		if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
		do l=1,lmaxn
			ln=l
!			ln=lln(l)
			l2=l+(ivar-1)*lmaxn
			m=mm(lln(l))
			n=nn(lln(l))
			ibnd=1
			if (m == 0) ibnd=2
			do lp=1,lmaxn
				lpn=lp
!				lpn=lln(lp)
				lpp=lp
				if (ich == -1) lpp=lo(lp)
				if (lpp == 0) cycle
				l1=lpp+(ieqn-1)*lmaxn
				do l3=1,leqmax
					l3n=l3
					if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3n,ln)
					if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3n,ln)
					if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3n,ln)
					if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3n,ln)
					do il=1,ith
						x=x*m
					end do
					do il=1,izt
						x=x*n
					end do
					if (x == 0.0_IDP) cycle
					x=x*iph
					do j=1,mjm1
						xa(j)=x*(rinv(j))**ith*abs(dble(n)-dble(m)*qqinv(j))
					end do
					if (ir == 0) then
						do j=1,mjm1
							yt(l1,j)=yt(l1,j)+xa(j)*tx(j,l3n)*xt(l2,j)
						end do
					else if (ir == 1) then
							yt(l1,1)=yt(l1,1)+wt10(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
							yt(l1,1)=yt(l1,1)+wt1p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm1
							yt(l1,j)=yt(l1,j)+wt10(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
					else
							yt(l1,1)=yt(l1,1)+wt20(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,1)
							yt(l1,1)=yt(l1,1)+wt2p(1,ibnd)*xa(1)*tx(1,l3n)*xt(l2,2)
						do j=2,mjm1
							yt(l1,j)=yt(l1,j)+wt20(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j)
							yt(l1,j)=yt(l1,j)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j-1)
							yt(l1,j)=yt(l1,j)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)*xt(l2,j+1)
						end do
					end if
				end do
			end do
		end do

	end subroutine b2lx_landau_grad_parallel

	subroutine b2lx_dlsq(ieqn,ivar,coef)

		use param
		use equil
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef

		interface
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		call b2lx(lplrr,1,ieqn,ivar,0,2,0,coef)
		call b2lx(lplrt,-1,ieqn,ivar,1,1,0,coef)
		call b2lx(lplrz,-1,ieqn,ivar,0,1,1,coef)
		call b2lx(lpltt,1,ieqn,ivar,2,0,0,coef)
		call b2lx(lpltz,1,ieqn,ivar,1,0,1,coef)
		call b2lx(lplzz,1,ieqn,ivar,0,0,2,coef)
		call b2lx(lplr,1,ieqn,ivar,0,1,0,coef)
		call b2lx(lplt,-1,ieqn,ivar,1,0,0,coef)
		call b2lx(lplz,-1,ieqn,ivar,0,0,1,coef)

	end subroutine b2lx_dlsq

	subroutine b2lx_dlsq_r(ieqn,ivar,coef)

		use param
		use equil
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef

		interface
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		call b2lx(gttoj,1,ieqn,ivar,0,2,0,coef)
		call b2lx(grtoj,-1,ieqn,ivar,1,1,0,-2.*coef)
		call b2lx(grroj,1,ieqn,ivar,2,0,0,coef)
		call b2lx(lplr_r,1,ieqn,ivar,0,1,0,coef)
		call b2lx(lplt_r,-1,ieqn,ivar,1,0,0,coef)

	end subroutine b2lx_dlsq_r
