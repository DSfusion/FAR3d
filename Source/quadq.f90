	subroutine quadq(n,x,y,x0,b,c,d,result)

		use param
		implicit none

		integer :: n
		real(IDP) :: x0,result
		real(IDP), dimension(:) :: x,y,b,c,d

!  sub. quadq

!  	a fortran sub. for computing an approximation
!  	to the integral of a fcn. y(x), given values of
!  	y(x) at the n distinct points  x(1).lt.x(2)...lt.x(n).

!  	n.b. this routine has been modified (from quad2) so that
!  	y(x) satisfies the constraint d**2 y/d x**2 at x=x(1)=0.

!  	the sub. uses a cubic spline fcn. to
!  	interpolate y(x). the resultant spline is then
!  	integrated, from x=x(1) to x=x0, to obtain
!  	the required approximation.
!  	the coefficients of the interpolating spline
!  	are obtained using an in-line modified version
!  	of the sub. "spline" which is published in
!  	the book 'computer methods for mathematical
!  	computations' by forsythe,malcolm, and moler.

!  input
!  	n.........the number of data points
!  	x.........the array of data points
!  	y.........the array of fcn. values
!  	x0........the upper limit of integration
!  	b,c,d.....work arrays


!  output
!  	result....the approximation to the integral of y(x) from
!  				 x = x(1)  to  x = x0.


!  p.w.gaffney   24th april 1978.


		real(IDP) :: sum,t,dj,alpha,a1,a3,a4,a20,b1,di
		integer :: nm1,ib,i,iup,j

!  	check that the upper limit of integration is sensible

		if (x0 < x(1)) return

!  	compute the coefficients of the interpolating spline

		if (n > 2) then

			nm1 = n-1

!  	set up tridiagonal system

!  	b = diagonal, d = offdiagonal, c = right hand side.

			d(1) = x(2) - x(1)
			c(2) = (y(2) - y(1))/d(1)
			do i = 2, nm1
				d(i) = x(i+1) - x(i)
				b(i) = 2.*(d(i-1) + d(i))
				c(i+1) = (y(i+1) - y(i))/d(i)
				c(i) = c(i+1) - c(i)
			end do

!  	end conditions.  third derivatives at  x(1)  and  x(n)
!  	obtained from divided differences

			b(1) = 1.
			d(1) = 0.
			b(n) = -d(nm1)
			c(1) = 0.
			c(n) = 0.
			if (n > 3) then
				c(n) = c(nm1)/(x(n)-x(n-2)) - c(n-2)/(x(nm1)-x(n-3))
				c(n) = -c(n)*d(nm1)**2/(x(n)-x(n-3))
			end if

!  	forward elimination

			do i = 2, n
				t = d(i-1)/b(i-1)
				b(i) = b(i) - t*d(i-1)
				c(i) = c(i) - t*c(i-1)
			end do

!  	back substitution

			c(n) = c(n)/b(n)
			do ib = 1, nm1
				i = n-ib
				c(i) = (c(i) - d(i)*c(i+1))/b(i)
			end do

!  	c(i) is now the sigma(i) of the text

!  	compute second derivative array c

			do i = 1,n
				c(i) = 3.*c(i)
			end do

!  	compute the integral of y(x) from x=x(1) to x=x0

			result = 0.0
			iup = nm1
			if (x0 < x(n)) then

!  	compute j such that x(j).le.x.lt.x(j+1)

				do j=1,nm1
					if (x0 < x(j+1)) exit
				end do

!  	first compute the integral from x=x(j) to x=x0

				iup = j-1
				dj = x(j+1)-x(j)
				alpha = (x0 - x(j))/dj
				a1 = 0.5*alpha*alpha
				a20=alpha-a1
				a3=a1*(a1-1.0)
				a4=alpha*alpha*(alpha-1.-0.5*a1)
				result=dj*(a1*y(j+1)+a20*y(j))+dj**3*(a3*c(j+1)+a4*c(j))/3.
			end if

!  	then compute the integral from x=x(1) to x=x(j)

			sum = 0.0
			do i=1,iup
				di = x(i+1)-x(i)
				sum = sum + 0.5*di*(y(i)+y(i+1)) - (c(i)+c(i+1))*di**3/12.
			end do
			result = sum + result

		else if (n == 2) then

			if (x0 < x(n)) then
				b1 = (y(2)-y(1))/(x(2)-x(1))
				result = (x0-x(1))*(y(1)+0.5*b1*(x0-x(1)))
			else
				result = 0.5*(x(2)-x(1))*(y(2)+y(1))
			end if
		
		end if

	end subroutine quadq