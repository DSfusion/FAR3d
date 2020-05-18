	subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)

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

!		Operator in linstart subroutine	
!		Add the model terms to the tridiagonal matrix	
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
							amat(l1,l2,j)=amat(l1,l2,j)+xa(j)*tx(j,l3n)
						end do
					else if(ir == 1) then
						do j=1,mjm1
							amat(l1,l2,j)=amat(l1,l2,j)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
							cmat(l1,l2,j)=cmat(l1,l2,j)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
							bmat(l1,l2,j)=bmat(l1,l2,j)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
						end do
					else
						do j=1,mjm1
							amat(l1,l2,j)=amat(l1,l2,j)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
							cmat(l1,l2,j)=cmat(l1,l2,j)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
							bmat(l1,l2,j)=bmat(l1,l2,j)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
						end do
					end if
				end do
			end do
		end do

	end subroutine blockj

	subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)

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

!		Operator in linstart subroutine	
!		Add the model terms to the tridiagonal matrix	
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
					amat(l1,l2,j)=amat(l1,l2,j)+xa(j)*tx(j)
				end do
			else if (ir == 1) then
				do j=1,mjm1
					amat(l1,l2,j)=amat(l1,l2,j)+wt10(j,ibnd)*xa(j)*tx(j)
					cmat(l1,l2,j)=cmat(l1,l2,j)+wt1m(j,ibnd)*xa(j)*tx(j)
					bmat(l1,l2,j)=bmat(l1,l2,j)+wt1p(j,ibnd)*xa(j)*tx(j)
				end do
			else
				do j=1,mjm1
					amat(l1,l2,j)=amat(l1,l2,j)+wt20(j,ibnd)*xa(j)*tx(j)
					cmat(l1,l2,j)=cmat(l1,l2,j)+wt2m(j,ibnd)*xa(j)*tx(j)
					bmat(l1,l2,j)=bmat(l1,l2,j)+wt2p(j,ibnd)*xa(j)*tx(j)
				end do
			end if
		end do

	end subroutine block0

	subroutine blockj_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

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

!		Operator in linstart subroutine	
!		Add the model terms to the tridiagonal matrix	
!		Similar than blockj but including the factor abs(dble(n)-dble(m)*qqinv(j))
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
							amat(l1,l2,j)=amat(l1,l2,j)+xa(j)*tx(j,l3n)
						end do
					else if(ir == 1) then
						do j=1,mjm1
							amat(l1,l2,j)=amat(l1,l2,j)+wt10(j,ibnd)*xa(j)*tx(j,l3n)
							cmat(l1,l2,j)=cmat(l1,l2,j)+wt1m(j,ibnd)*xa(j)*tx(j,l3n)
							bmat(l1,l2,j)=bmat(l1,l2,j)+wt1p(j,ibnd)*xa(j)*tx(j,l3n)
						end do
					else
						do j=1,mjm1
							amat(l1,l2,j)=amat(l1,l2,j)+wt20(j,ibnd)*xa(j)*tx(j,l3n)
							cmat(l1,l2,j)=cmat(l1,l2,j)+wt2m(j,ibnd)*xa(j)*tx(j,l3n)
							bmat(l1,l2,j)=bmat(l1,l2,j)+wt2p(j,ibnd)*xa(j)*tx(j,l3n)
						end do
					end if
				end do
			end do
		end do

	end subroutine blockj_landau_grad_parallel

	subroutine block_dlsq(ieqn,ivar,coef)

		use param
		use equil
		implicit none

		integer :: ieqn,ivar
		real(IDP) :: coef

		interface
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
		end interface

		call blockj(lplrr,1,ieqn,ivar,0,2,0,coef)
		call blockj(lplrt,-1,ieqn,ivar,1,1,0,coef)
		call blockj(lplrz,-1,ieqn,ivar,0,1,1,coef)
		call blockj(lpltt,1,ieqn,ivar,2,0,0,coef)
		call blockj(lpltz,1,ieqn,ivar,1,0,1,coef)
		call blockj(lplzz,1,ieqn,ivar,0,0,2,coef)
		call blockj(lplr,1,ieqn,ivar,0,1,0,coef)
		call blockj(lplt,-1,ieqn,ivar,1,0,0,coef)
		call blockj(lplz,-1,ieqn,ivar,0,0,1,coef)

	end subroutine block_dlsq
