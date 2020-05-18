	subroutine dbydzt(d,a,ltype,c1,c2)

		use param
		use domain
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d

!		Operator: derivative in zeta 	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dzt
	
		a(:,0)=0. 
		do l=1,lmax
			d(:,l)=c1*d(:,l)-ltype*nn(l)*c2*a(:,l)
		end do
		d(:,0)=0. 

	end subroutine dbydzt

	subroutine dbydth(d,a,ltype,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,ltype,l,m
		real(IDP) :: c1,c2,temp
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in theta 	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*rinv*d(a)/dth		
		
		a(:,0)=0.0_IDP
		do l=1,lmax
			dold=d(:,l) 
			temp=-mm(l)*ltype
			d(:,l)=temp*rinv*a(:,l) 
			d(0,l)=0.0_IDP
			m=abs(mm(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
			d(:,l)=c1*dold+c2*d(:,l) 
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydth

	subroutine grdpar(d,a,ltype,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d

!		Operator: d/dzt - qqinv*rinv*d/dth  	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*(d(a)/dzt - qqinv*rinv*d(a)/dth)		
		
		a(:,0)=0.0_IDP
		do l=1,lmax
			d(:,l)=c1*d(:,l)-ltype*(nn(l)-mm(l)*qqinv)*c2*a(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine grdpar

	subroutine dbydr(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in rho 	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr		
		
		a(:,0)=0.0_IDP
		do l=1,lmax
			dold=d(:,l) 
			do j=1,mjm1
				d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
			end do
			d(0,l)=0.0_IDP
			m=abs(mm(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=rinv(1)*a(1,l)
			d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
			d(:,l)=c1*dold+c2*d(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydr

	subroutine dbydr0(d,a,c1,c2,k)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: k,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in rho (variables with only radial dependency)	
!		d = term name
!		a = element to be derivated
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr				
		
		dold=d
		do j=1,mjm1
			d(j)=dc1m(j)*(a(j-1)-a(j))+dc1p(j)*(a(j+1)-a(j))
		end do
		d(0)=0.
		m=0
		if (k > 1) m=abs(abs(m-1)-1)
		if (k == 1 .or. k == 3) m=m+1
		if (m == 1) d(0)=rinv(1)*a(1)
		d(mj)=(a(mj)-a(mjm1))/(r(mj)-r(mjm1))
		d=c1*dold+c2*d

	end subroutine dbydr0

	subroutine d2bydr20(d,a,c1,c2,k)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: k,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: second order derivative in rho (variables with only radial dependency)	
!		d = term name
!		a = element to be derivated
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr				
		
		dold=d
		do j=1,mjm1
			d(j)=dc2m(j)*(a(j-1)-a(j))+dc2p(j)*(a(j+1)-a(j))
		end do
		d(0)=0.
		m=0
		if (k > 1) m=abs(abs(m-1)-1)
		if (k == 1 .or. k == 3) m=m+1
		if (m == 0 .or. m == 2) d(0)=d(1)
		d(mj)=d(mjm1)
		d=c1*dold+c2*d

	end subroutine d2bydr20

	subroutine dbydtheq(d,a,ltype,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,ltype,l,m
		real(IDP) :: c1,c2,temp
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in theta (equilibrium variables)	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*rinv*d(a)/dth				
		
		a(:,0)=0.0_IDP
		do l=1,leqmax
			dold=d(:,l) 
			temp=-mmeq(l)*ltype
			d(:,l)=temp*rinv*a(:,l) 
			d(0,l)=0.0_IDP
			m=abs(mmeq(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (k == 4) m=m+2
			if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
			d(:,l)=c1*dold+c2*d(:,l) 
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydtheq

	subroutine grpareq(d,a,ltype,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d

!		Operator: d/dzt - qqinv*rinv*d/dth (equilibrium variables)	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*(d(a)/dzt - qqinv*rinv*d(a)/dth)			
		
		a(:,0)=0.0_IDP
		do l=1,leqmax
			d(:,l)=c1*d(:,l)-ltype*(nneq(l)-mmeq(l)*qqinv)*c2*a(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine grpareq

	subroutine dbydreq(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Operator: derivative in rho (equilibrium variables)	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dr				
		
		a(:,0)=0.0_IDP
		do l=1,leqmax
			dold=d(:,l) 
			do j=1,mjm1
				d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
			end do
			d(0,l)=0.0_IDP
			m=abs(mmeq(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (k == 4) m=m+2
			if (m == 1) d(0,l)=rinv(1)*a(1,l)
			d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
			d(:,l)=c1*dold+c2*d(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydreq

	subroutine dbydzteq(d,a,ltype,c1,c2)

		use param
		use domain
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d

!		Operator: derivative in zeta (equilibrium variables)	
!		d = term name
!		a = element to be derivated
!		ltype = sign of the element
!		c1 = multiplicative coefficient
!		c2 = includes in d a previous d term component
!		d = c2*d + c1*ltype*d(a)/dzt		
		
		a(:,0)=0.0_IDP 
		do l=1,leqmax
			d(:,l)=c1*d(:,l)-ltype*nneq(l)*c2*a(:,l)
		end do
		d(:,0)=0.0_IDP

	end subroutine dbydzteq

	subroutine dbydrl(d,a,c1,c2,k,l)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Special operator to calculate the magnetic shear	
!		used in the output subroutine in qqinv	
		
		dold=d 
		do j=1,mjm1
			d(j)=dc1m(j)*(a(j-1)-a(j))+dc1p(j)*(a(j+1)-a(j))
		end do
		d(0)=0.
		m=abs(mm(l))
		if (k > 1) m=abs(abs(m-1)-1)
		if (k == 1 .or. k == 3) m=m+1
		if (m == 1) d(0)=rinv(1)*a(1)
		d(mj)=(a(mj)-a(mjm1))/(r(mj)-r(mjm1))
		d=c1*dold+c2*d 

	end subroutine dbydrl


	subroutine dbydztb(d,a,ltype,c1,c2)

		use param
		use domain
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		
!		Derivative in zeta for equilibrium variables		

		a(:,0)=0.0_IDP 
		do l=1,lbmax
			d(:,l)=c1*d(:,l)-ltype*nnb(l)*c2*a(:,l)
		end do
		d(:,0)=0.0_IDP 

	end subroutine dbydztb

	subroutine dbydthb(d,a,ltype,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,ltype,l,m
		real(IDP) :: c1,c2,temp
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold
		
!		Derivative in theta for equilibrium variables			

		a(:,0)=0.0_IDP 
		do l=1,lbmax
			dold=d(:,l) 
			temp=-mmb(l)*ltype
			d(:,l)=temp*rinv*a(:,l) 
			d(0,l)=0.0_IDP
			m=abs(mmb(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
			d(:,l)=c1*dold+c2*d(:,l) 
		end do
		d(:,0)=0.0_IDP 

	end subroutine dbydthb

	subroutine grparb(d,a,ltype,c1,c2)

		use param
		use domain
		use equil
		implicit none

		integer :: ltype,l
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		
!		Operator d/dzt - qqinv*d/dth for equilibrium variables

		a(:,0)=0.0_IDP 
		do l=1,lbmax
			d(:,l)=c1*d(:,l)-ltype*(nnb(l)-mmb(l)*qqinv)*c2*a(:,l)
		end do
		d(:,0)=0.0_IDP 

	end subroutine grparb

	subroutine dbydrb(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold
		
!		Derivative in rho for equilibrium variables			

		a(:,0)=0.0_IDP 
		do l=1,lbmax
			dold=d(:,l) 
			do j=1,mjm1
				d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
			end do
			d(0,l)=0.0_IDP
			m=abs(mmb(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=rinv(1)*a(1,l)
			d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
			d(:,l)=c1*dold+c2*d(:,l)
		end do
		d(:,0)=0.0_IDP 

	end subroutine dbydrb

	subroutine dbydrrb(d,a,c1,c2,k)

		use param
		use domain
		implicit none

		integer :: k,l,j,m
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: a,d
		real(IDP), dimension(0:mj) :: dold

!		Second order derivative in rho for equilibrium variables			
		
		a(:,0)=0.0_IDP 
		do l=1,lbmax
			dold=d(:,l) 
			do j=1,mjm1
				d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))+rinv(j)*a(j,l)
			end do
			d(0,l)=0.0_IDP
			m=abs(mmb(l))
			if (k > 1) m=abs(abs(m-1)-1)
			if (k == 1 .or. k == 3) m=m+1
			if (m == 1) d(0,l)=2.*rinv(1)*a(1,l)
			d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))+rinv(mj)*a(mj,l)
			d(:,l)=c1*dold+c2*d(:,l)
		end do
		d(:,0)=0.0_IDP 

	end subroutine dbydrrb
