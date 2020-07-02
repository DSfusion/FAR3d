	subroutine dlsq_r(ss,ff,itypf,wk1,wk2,wk3,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: itypf,l,j,le
		real(IDP) :: c1,c2,xm,xn,coefs
		real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3

		interface
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
		end interface

		ff(:,0)=0. 
		ss(:,0)=0. 
		wk1=0. 
		wk2=0. 

		call dbydr(wk2,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		do l=1,lmax
			do j=1,mjm1
				wk2(j,l)=wk2(j,l)*rinv(j)
			end do
		end do
		call mult(ss,wk2,itypf,gttoj,1,1.0_IDP,c2)
		call dbydr(wk1,gttoj,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk3,wk1,1,sqgi,1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,1,1.0_IDP,c2)
		call dbydth(wk1,grtoj,-1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk3,wk1,1,sqgi,1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,1,1.0_IDP,-c2)

		call dbydth(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call dbydth(wk1,grroj,1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk3,wk1,-1,sqgi,1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,-itypf,wk3,-1,1.0_IDP,c2)
		call dbydr(wk1,grtoj,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(wk3,wk1,-1,sqgi,1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,-itypf,wk3,-1,1.0_IDP,-c2)

		call dbydr(wk1,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call dbydth(wk2,wk1,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,-itypf,grtoj,-1,1.0_IDP,-2*c2)

		call dbydr(wk2,wk1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(ss,wk2,itypf,gttoj,1,1.0_IDP,c2)

		call dbydth(wk1,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call dbydth(wk2,wk1,-itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,itypf,grroj,1,1.0_IDP,c2)

		wk1(:,0)=0. 
		wk2(:,0)=0. 
		wk3(:,0)=0. 
		ss(:,0)=0. 
		ff(:,0)=0. 

	end subroutine dlsq_r