	subroutine dlstar(ss,ff,itypf,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,c1,c2)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itypf,l,j
		real(IDP) :: c1,c2,xm
		real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11

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

!		Move equilibrium variables to dynamic

		call eqtodyn(wk1,jbgrr,0.0_IDP,1.0_IDP)
		call eqtodyn(wk2,jbgrt,0.0_IDP,1.0_IDP)
		call eqtodyn(wk3,jbgtt,0.0_IDP,1.0_IDP)
		call eqtodyn(wk4,jbgrz,0.0_IDP,1.0_IDP)
		call eqtodyn(wk5,jbgtz,0.0_IDP,1.0_IDP)
		call eqtodyn(wk6,bst,0.0_IDP,1.0_IDP)

!		d2Phi / drho2 components

		do l=1,lmax
			do j=1,mjm1
				wk7(j,l)=del2cp(j)*(ff(j+1,l)-ff(j,l))+del2cm(j)*(ff(j-1,l)-ff(j,l))
			end do
		end do
		do l=1,lmax
			wk8(:,l)=-denseq*(feq*wk3(:,l)+(rinv*cureq*wk5(:,l)/eps))/(feq-qqinv*cureq)
		end do
		call mult(ss,wk7,itypf,wk8,1,c1,c2)

!		d2Phi / dtheta2 components

		do l=1,lmax
			xm=mm(l)
			wk7(1:mj,l)=-denseq(1:mj)*feq(1:mj)*(rinv(1:mj)*xm)**2*ff(1:mj,l)/(feq(1:mj)-qqinv(1:mj)*cureq(1:mj))
		end do
		do l=1,lmax
			wk7(0,l)=0.
		end do
		call mult(ss,wk1,1,wk7,itypf,1.0_IDP,c2)
		do l=1,lmax
			xm=mm(l)
			wk8(1:mj,l)=-denseq(1:mj)*r(1:mj)*(rinv(1:mj)*xm)**2*ff(1:mj,l)/(eps*(feq(1:mj)-qqinv(1:mj)*cureq(1:mj)))
		end do
		do l=1,lmax
			wk8(0,l)=0.
		end do
		call mult(wk9,wk6,-1,wk8,itypf,0.0_IDP,c2)
		call mult(ss,wk4,-1,wk9,-itypf,1.0_IDP,c2)

! !		d2Phi / drhodzeta components

		call dbydr(wk7,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(0,l)=0.
		end do
		call dbydth(wk8,wk7,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk8(0,l)=0.
		end do
		call mult(wk9,wk3,1,wk6,-1,0.0_IDP,c2)
		do l=1,lmax
			wk10(:,l)=denseq*(cureq*rinv*wk2(:,l)-r*wk9(:,l))/(feq-qqinv*cureq)
		end do
		call mult(ss,wk8,-itypf,wk10,-1,1.0_IDP,c2)

! !		d2Phi / dthetadzeta components

        call dbydzt(wk7,ff,itypf,0.0_IDP,1.0_IDP)
		call dbydth(wk8,wk7,-itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk8(0,l)=0.
		end do
		call mult(wk9,wk2,-1,wk6,-1,0.0_IDP,c2)
		do l=1,lmax
			wk10(:,l)=denseq*(r*wk9(:,l)-rinv*cureq*wk1(:,l))/(feq-qqinv*cureq)
		end do
		call mult(ss,wk8,itypf,wk10,1,1.0_IDP,c2)

! !		d2Phi / drhodtheta components

		call dbydth(wk7,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(0,l)=0.
		end do
		call dbydr(wk8,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk8(0,l)=0.
		end do
                call mult(wk7,wk5,1,wk6,-1,0.0_IDP,c2)
		do l=1,lmax
			wk10(:,l)=denseq*((cureq*rinv*wk4(:,l)/eps)+2*feq*wk2(:,l)+(r*wk7(:,l)/eps))/(feq-qqinv*cureq)
		end do
		call mult(ss,wk8,-itypf,wk10,-1,1.0_IDP,c2)

! !		dPhi / drho components

		do l=1,lmax
			wk7(:,l)=denseq*r*feq*wk3(:,l)/(feq-qqinv*cureq)
		end do
		call dbydr(wk8,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=denseq*cureq*wk5(:,l)/(eps*(feq-qqinv*cureq))
		end do
		call dbydr(wk9,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=denseq*feq*wk2(:,l)/(feq-qqinv*cureq)
		end do
		call dbydth(wk10,wk7,-1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk7(:,l)=denseq*rinv*cureq*wk4(:,l)/(eps*(feq-qqinv*cureq))
		end do
		call dbydth(wk11,wk7,-1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk7(:,l)=-rinv*(wk8(:,l)+wk9(:,l)) + wk10(:,l) + wk11(:,l)
		end do
		call dbydr(wk8,ff,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk8(0,l)=0.
		end do
		call mult(ss,wk8,itypf,wk7,1,1.0_IDP,c2)

! !		dPhi / dtheta components

		call mult(wk9,wk6,-1,wk5,-1,0.0_IDP,c2)
		do l=1,lmax
			wk7(:,l)=denseq*r*wk9(:,l)/(eps*(feq-qqinv*cureq))
		end do
		call dbydr(wk8,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=denseq*feq*wk2(:,l)/(feq-qqinv*cureq)
		end do
		call dbydr(wk9,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=denseq*feq*wk1(:,l)/(feq-qqinv*cureq)
		end do
		call dbydth(wk10,wk7,1,0.0_IDP,1.0_IDP,2)
		call mult(wk11,wk6,-1,wk4,-1,0.0_IDP,c2)
		do l=1,lmax
			wk7(:,l)=denseq*r*wk11(:,l)/(eps*(feq-qqinv*cureq))
		end do
		call dbydth(wk11,wk7,1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk7(:,l)=wk8(:,l) + wk9(:,l) - wk10(:,l) - wk11(:,l)
		end do
		call dbydth(wk8,ff,itypf,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk8(0,l)=0.
		end do
		call mult(ss,wk8,-itypf,wk7,-1,1.0_IDP,c2)

! !		dPhi / dzeta components

		do l=1,lmax
			wk7(:,l)=denseq*cureq*wk2(:,l)/(feq-qqinv*cureq)
		end do
		call dbydr(wk8,wk7,0.0_IDP,1.0_IDP,0)
		call mult(wk9,wk6,-1,wk3,1,0.0_IDP,c2)
		do l=1,lmax
			wk7(:,l)=denseq*r*r*wk9(:,l)/(feq-qqinv*cureq)
		end do
		call dbydr(wk9,wk7,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=rinv*denseq*cureq*wk1(:,l)/(feq-qqinv*cureq)
		end do
		call dbydth(wk10,wk7,1,0.0_IDP,1.0_IDP,2)
		call mult(wk11,wk6,-1,wk2,-1,0.0_IDP,c2)
		do l=1,lmax
			wk7(:,l)=r*denseq*wk11(:,l)/(feq-qqinv*cureq)
		end do
		call dbydth(wk11,wk7,1,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk7(:,l)=rinv*(wk8(:,l) - wk9(:,l)) - wk10(:,l) + wk11(:,l)
		end do
                call dbydzt(wk8,ff,itypf,0.0_IDP,1.0_IDP)
		call mult(ss,wk8,-itypf,wk7,-1,1.0_IDP,c2)

! !     Extra terms of the vorticity from the parallel thermal plasma velocity

		do l=1,lmax
			wk7(:,l)=r*denseq*wk5(:,l)/eps
		end do
		call dbydr(wk8,wk7,0.0_IDP,1.0_IDP,2)
		do l=1,lmax
			wk7(:,l)=denseq*wk4(:,l)/eps
		end do
		call dbydth(wk9,wk7,-1,0.0_IDP,1.0_IDP,0)
		do l=1,lmax
			wk7(:,l)=rinv*wk8(:,l)-wk9(:,l)
		end do
		call mult(ss,wk7,1,vthprlf,-1,1.0_IDP,1.0_IDP)

		do l=1,lmax
			wk7(:,l)=denseq*wk5(:,l)/eps
		end do
		call dbydr(wk8,vthprlf,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk8,-1,wk7,1,1.0_IDP,1.0_IDP)
		do l=1,lmax
			wk7(:,l)=-denseq*wk4(:,l)/eps
		end do
		call dbydth(wk8,vthprlf,-1,0.0_IDP,1.0_IDP,0)
		call mult(ss,wk8,1,wk7,-1,1.0_IDP,1.0_IDP)

		do l=1,lmax
			ss(mj,l)=(ss(mjm1,l)*(r(mj)-r(mj-2))-ss(mj-2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mj-2))
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