	subroutine dlstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: itypf,l,j
		real(IDP) :: c1,c2,xm
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
				wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
			end do
		end do
		call eqtodyn(wk3,jbgtt,0.0_IDP,1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq*wk3(:,l)
		end do
		call mult(ss,wk2,itypf,wk3,2,c1,c2)
		call dbydr(wk2,wk3,0.0_IDP,1.0_IDP,2)
		call eqtodyn(wk3,jbgrt,0.0_IDP,-1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq*wk3(:,l)
		end do
		call dbydth(wk1,wk3,-1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
		call dbydr(wk2,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call dbydth(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(ss,wk1,-itypf,wk3,-2,1.0_IDP,2.*c2)
		call dbydr(wk1,wk3,0.0_IDP,1.0_IDP,0)
		call eqtodyn(wk3,jbgrr,0.0_IDP,1.0_IDP)
		do l=1,lmax
			wk3(:,l)=denseq*wk3(:,l)
		end do
		call dbydth(wk2,wk3,1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk1(:,l)=wk1(:,l)+wk2(:,l)
		end do
		call dbydth(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)
		do l=1,lmax
			xm=mm(l)
			wk2(1:mj,l)=-(rinv(1:mj)*xm)**2*ff(1:mj,l)
		end do
		call mult(ss,wk3,2,wk2,itypf,1.0_IDP,c2)
		do l=1,lmax
			ss(mj,l)=(ss(mjm1,l)*(r(mj)-r(mjm2))-ss(mjm2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mjm2))
		end do
		do l=1,lmax
			if (mm(l) == 0) ss(0,l)=(r(2)**2*ss(1,l)-r(1)**2*ss(2,l))/(r(2)**2-r(1)**2)
		end do
		wk1(:,0)=0. 
		wk2(:,0)=0. 
		wk3(:,0)=0. 
		ss(:,0)=0. 
		ff(:,0)=0. 

	end subroutine dlstar

	subroutine dlstar_ext(ss,ff,itypf,wk1,wk2,wk3,wkeq1,wkeq2,wkeq3,tx,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: itypf,l,j
		real(IDP) :: c1,c2,xm,epsq
		real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3
		real(IDP), dimension(0:,0:) :: wkeq1,wkeq2,wkeq3
		real(IDP), dimension(0:) ::  tx

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
		end interface

		ff(:,0)=0. 
		ss(:,0)=0. 
		wk1=0. 
		wk2=0. 

		epsq=eps*eps
		tx=feq-qqinv*cureq
		do l=1,leqmax
			wkeq1(:,l)=denseq*(jbgtt(:,l)-rinv*rinv*cureq*cureq*jsq(:,l)/(epsq*tx))
		end do
		do l=1,lmax
			do j=1,mjm1
				wk2(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
			end do
		end do
		call eqtodyn(wk1,wkeq1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk1,2,c1,c2)
		do l=1,leqmax
			wkeq2(:,l)=-denseq*(jbgrt(:,l)-cureq*bstg(:,l)/(epsq*tx))
		end do
		call dbydreq(wkeq3,wkeq1,0.0_IDP,1.0_IDP,2)
		call dbydtheq(wkeq3,wkeq2,-1,1.0_IDP,1.0_IDP,2)
		call eqtodyn(wk1,wkeq3,0.0_IDP,1.0_IDP)
		call dbydr(wk2,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,itypf,wk1,2,1.0_IDP,c2)
		call dbydth(wk1,wk2,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		do l=1,leqmax
			wkeq1(:,l)=-denseq*(feq*jbgrt(:,l)-r*r*qqinv*bsjgtt(:,l)-cureq*bstg(:,l)/epsq)/tx
		end do
		wkeq3=wkeq1+wkeq2
		call eqtodyn(wk2,wkeq3,0.0_IDP,1.0_IDP)
		call mult(ss,wk1,-itypf,wk2,-2,1.0_IDP,c2)
		do l=1,leqmax
			wkeq2(:,l)=denseq*(feq*jbgrr(:,l)-r*r*(qqinv*bsjgrt(:,l)+bsqg(:,l)/epsq))/tx
		end do
		call eqtodyn(wk1,wkeq2,0.0_IDP,1.0_IDP)
		do l=1,lmax
			xm=mm(l)
			wk2(1:mj,l)=-(rinv(1:mj)*xm)**2*ff(1:mj,l)
		end do
		call mult(ss,wk1,2,wk2,itypf,1.0_IDP,c2)
		call dbydreq(wkeq3,wkeq1,0.0_IDP,1.0_IDP,2)
		call dbydtheq(wkeq3,wkeq2,1,1.0_IDP,1.0_IDP,2)
		call eqtodyn(wk1,wkeq3,0.0_IDP,1.0_IDP)
		call dbydth(wk2,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		call mult(ss,wk2,-itypf,wk1,-2,1.0_IDP,c2)
		call dbydzt(wk1,ff,itypf,0.0_IDP,1.0_IDP)
		call dbydr(wk2,wk1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		do l=1,leqmax
			wkeq1(:,l)=denseq*(rinv*cureq*jbgrt(:,l)-r*bsjgtt(:,l))/tx
		end do
		call eqtodyn(wk3,wkeq1,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,-itypf,wk3,-2,1.0_IDP,c2)
		call dbydth(wk2,wk1,-itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk2(0,l)=0.
		end do
		do l=1,leqmax
			wkeq2(:,l)=-denseq*(rinv*cureq*jbgrr(:,l)-r*bsjgrt(:,l))/tx
		end do
		call eqtodyn(wk3,wkeq2,0.0_IDP,1.0_IDP)
		call mult(ss,wk2,itypf,wk3,2,1.0_IDP,c2)
		call dbydreq(wkeq3,wkeq1,0.0_IDP,1.0_IDP,3)
		do l=1,leqmax
			wkeq3(:,l)=wkeq3(:,l)+rinv*wkeq1(:,l)
		end do
		call dbydtheq(wkeq3,wkeq2,1,1.0_IDP,1.0_IDP,3)
		call eqtodyn(wk2,wkeq3,0.0_IDP,1.0_IDP)
		do l=1,lmax
			wk1(0,l)=0.
		end do
		call mult(ss,wk1,-itypf,wk2,-2,1.0_IDP,c2)
		do l=1,lmax
			ss(mj,l)=(ss(mjm1,l)*(r(mj)-r(mjm2))-ss(mjm2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mjm2))
		end do
		do l=1,lmax
			if (mm(l) == 0) ss(0,l)=(r(2)**2*ss(1,l)-r(1)**2*ss(2,l))/(r(2)**2-r(1)**2)
		end do
		wkeq1(:,0)=0. 
		wkeq2(:,0)=0. 
		wkeq3(:,0)=0. 
		wk1(:,0)=0. 
		wk2(:,0)=0. 
		wk3(:,0)=0. 
		ss(:,0)=0. 
		ff(:,0)=0. 

	end subroutine dlstar_ext
