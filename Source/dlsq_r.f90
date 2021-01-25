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
		end interface

		ff(:,0)=0. 
		ss(:,0)=0. 
		wk1=0. 
		wk2=0. 

		do l=1,lmax
			do j=1,mjm1
				wk2(j,l)=dc2p(j)*(ff(j+1,l)-ff(j,l))+dc2m(j)*(ff(j-1,l)-ff(j,l))
			end do
		end do
		call eqtodyn(wk3,gttoj,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,2,c1,c2)
		call eqtodyn(wk1,lplr_r,0.0_IDP,1.0_IDP)
		call dbydr(wk2,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call eqtodyn(wk3,grtoj,0.0_IDP,1.0_IDP)
		call dbydth(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,-2.*c2)
		call eqtodyn(wk1,lplt_r,0.0_IDP,1.0_IDP)
		call dbydth(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)
		call eqtodyn(wk3,grroj,0.0_IDP,1.0_IDP)
		do l=1,lmax
			xm=mm(l)
			wk2(1:mj,l)=-(rinv(1:mj)*xm)**2*ff(1:mj,l)
		end do
		call mult(ss,wk2,itypf,wk3,2,1.0_IDP,c2)
		wk1(:,0)=0. 
		wk2(:,0)=0. 
		wk3(:,0)=0. 
		ss(:,0)=0. 
		ff(:,0)=0. 

	end subroutine dlsq_r
