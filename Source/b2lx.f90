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

!		Operator in linscheck subroutine	
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

	end subroutine b2lx

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

!		Operator in linscheck subroutine	
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
					yt(l1,j)=yt(l1,j)+xa(j)*tx(j)*xt(l2,j)
				end do
			else if (ir == 1)then
					yt(l1,1)=yt(l1,1)+wt10(1,ibnd)*xa(1)*tx(1)*xt(l2,1)
					yt(l1,1)=yt(l1,1)+wt1p(1,ibnd)*xa(1)*tx(1)*xt(l2,2)
				do j=2,mjm1
					yt(l1,j)=yt(l1,j)+wt10(j,ibnd)*xa(j)*tx(j)*xt(l2,j)
					yt(l1,j)=yt(l1,j)+wt1m(j,ibnd)*xa(j)*tx(j)*xt(l2,j-1)
					yt(l1,j)=yt(l1,j)+wt1p(j,ibnd)*xa(j)*tx(j)*xt(l2,j+1)
				end do
			else
					yt(l1,1)=yt(l1,1)+wt20(1,ibnd)*xa(1)*tx(1)*xt(l2,1)
					yt(l1,1)=yt(l1,1)+wt2p(1,ibnd)*xa(1)*tx(1)*xt(l2,2)
				do j=2,mjm1
					yt(l1,j)=yt(l1,j)+wt20(j,ibnd)*xa(j)*tx(j)*xt(l2,j)
					yt(l1,j)=yt(l1,j)+wt2m(j,ibnd)*xa(j)*tx(j)*xt(l2,j-1)
					yt(l1,j)=yt(l1,j)+wt2p(j,ibnd)*xa(j)*tx(j)*xt(l2,j+1)
				end do
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

!		Operator in linscheck subroutine	
!		Similar than b2lx but including the factor abs(dble(n)-dble(m)*qqinv(j))
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

!		ith_zt=1
!		iph=ity**(ith_zt)*(-1)**((ith_zt+1)/2)  !change on 5/26/2015 - grad-parallel-theta and zeta derivatives
!		ity=ity*(-1)**(ith_zt)                  !change on 5/26/2015 - grad-parallel-theta and zeta derivatives
		
!		ich=-itx*ity
		ich=itx*ity                             !change on 5/26/2015 - grad-parallel-theta and zeta derivatives
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

	subroutine b2lx_dlsq_r(ieqn,ivar,wk1,wk2,coef)

		use param
		use domain
		use equil
		implicit none

		integer :: ieqn,ivar,itypf,l,j
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: wk1,wk2

		interface
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine mult
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
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		wk1=0. 
		wk2=0. 

		do l=1,lmax
			do j=1,mjm1
				wk1(j,l)=gttoj(j,l)*rinv(j)
			end do
		end do
		call b2lx(wk1,1,ieqn,ivar,0,1,0,coef)
		call dbydr(wk1,gttoj,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk2,wk1,1,sqgi,1,0.0_IDP,1.0_IDP)
		call b2lx(wk2,1,ieqn,ivar,0,1,0,coef)
		call dbydth(wk1,grtoj,-1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk2,wk1,1,sqgi,1,0.0_IDP,1.0_IDP)
		call b2lx(wk2,1,ieqn,ivar,0,1,0,-coef)

		call dbydth(wk1,grroj,1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk2,wk1,-1,sqgi,1,0.0_IDP,1.0_IDP)
		call b2lx(wk2,-1,ieqn,ivar,1,0,0,coef)
		call dbydr(wk1,grtoj,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk2,wk1,-1,sqgi,1,0.0_IDP,1.0_IDP)
		call b2lx(wk2,-1,ieqn,ivar,1,0,0,-coef)

		call b2lx(grtoj,-1,ieqn,ivar,1,1,0,-2*coef)

		call b2lx(gttoj,1,ieqn,ivar,0,2,0,coef)

		call b2lx(grroj,1,ieqn,ivar,2,0,0,coef)

		wk1(:,0)=0. 
		wk2(:,0)=0. 

	end subroutine b2lx_dlsq_r

	subroutine b2lx_dlsq_x(ieqn,ivar,wk1,wk2,coef)

		use param
		use domain
		use equil
		implicit none

		integer :: ieqn,ivar,itypf,l,j
		real(IDP) :: coef
		real(IDP), dimension(0:,0:) :: wk1,wk2

		interface
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine mult
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
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
		end interface

		wk1=0. 

		do l=1,lmax
			wk1(:,l)=gtt(:,l)*rinv*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,1,ieqn,ivar,0,1,0,coef)
		call dbydr(wk1,gtt,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		do l=1,lmax
		        wk1(:,l)=wk1(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,1,ieqn,ivar,0,1,0,coef)
		call dbydth(wk1,grt,-1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		do l=1,lmax
		        wk1(:,l)=wk1(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,1,ieqn,ivar,0,1,0,-coef)

		call dbydth(wk1,grr,1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		do l=1,lmax
		        wk1(:,l)=wk1(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,-1,ieqn,ivar,1,0,0,coef)
		call dbydr(wk1,grt,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		do l=1,lmax
		        wk1(:,l)=wk1(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,-1,ieqn,ivar,1,0,0,-coef)

		do l=1,lmax
		        wk1(:,l)=grt(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,-1,ieqn,ivar,1,1,0,-2*coef)

		do l=1,lmax
		        wk1(:,l)=gtt(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,1,ieqn,ivar,0,2,0,coef)

		do l=1,lmax
		        wk1(:,l)=grr(:,l)*ti/(feq-qqinv*cureq)
		end do
		call b2lx(wk1,1,ieqn,ivar,2,0,0,coef)

		wk1(:,0)=0.

	end subroutine b2lx_dlsq_x