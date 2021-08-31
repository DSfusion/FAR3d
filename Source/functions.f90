module ffunctions
  implicit none
  contains
	function zeroin(ax,bx,f,tol)

		use param
		implicit none
		
		interface
			function f(x)
				use param
				implicit none
				real(IDP) :: f,x
			end function f
		end interface
		
		real(IDP) :: ax,bx,tol,zeroin

!  	a zero of the function  f(x)  is computed in the interval ax,bx .

!  input..

!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!  		 the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!  		 final result ( .ge. 0.0)


!  output..

!  zeroin abcissa approximating a zero of  f  in the interval ax,bx


!  	it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!  	this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).


		real(IDP) :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s

!  compute eps, the relative machine precision

		eps = 1.0
		tol1 = 2.0
		do while (tol1 > 1.0)
			eps = eps/2.0
			tol1 = 1.0 + eps
		end do

!  initialization

		a = ax
		b = bx
		fa = f(a)
		fb = f(b)

!  begin step

		c = a
		fc = fa
		d = b - a
		e = d
		do while (.true.)
			if (abs(fc) < abs(fb)) then
      		a = b
      		b = c
      		c = a
      		fa = fb
      		fb = fc
      		fc = fa
			end if

!  convergence test

			tol1 = 2.0*eps*abs(b) + 0.5*tol
			xm = .5*(c - b)
			if (abs(xm) <= tol1) exit
			if (fb == 0.0) exit

!  is bisection necessary

      	if (abs(e) < tol1 .or. abs(fa) <= abs(fb)) then

!  bisection

				d = xm
				e = d
				
			else

!  is quadratic interpolation possible

				if (a == c) then

!  linear interpolation

					s = fb/fa
					p = 2.0*xm*s
					q = 1.0 - s

				else

!  inverse quadratic interpolation

					q = fa/fc
					r = fb/fc
					s = fb/fa
					p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
					q = (q - 1.0)*(r - 1.0)*(s - 1.0)

				end if

!  adjust signs

				if (p > 0.0) q = -q
				p = abs(p)

!  is interpolation acceptable

      		if ((2.0*p) >= (3.0*xm*q - abs(tol1*q)) .or. p >= abs(0.5*e*q)) then
					d = xm
					e = d
				else
					e = d
      			d = p/q
				end if
				
			end if

!  complete step

			a = b
			fa = fb
			if (abs(d) > tol1) then
				b = b + d
			else
				b = b + sign(tol1, xm)
			end if
			fb = f(b)
			if ((fb*(fc/abs(fc))) > 0.0) then
				c = a
				fc = fa
				d = b - a
				e = d
			end if
		end do

!  done

		zeroin = b

	end function zeroin

	function seval(n, u, x, y, b, c, d)

		use param
		implicit none

		integer :: n
		real(IDP) :: seval,u
		real(IDP), dimension(:) :: x, y, b, c, d

!  	this subroutine evaluates the cubic spline function

!  	seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3

!  	where  x(i) .lt. u .lt. x(i+1), using horner's rule

!  	if  u .lt. x(1) then  i = 1  is used.
!  	if  u .ge. x(n) then  i = n  is used.

!  	input..

!  	n = the number of data points
!  	u = the abscissa at which the spline is to be evaluated
!  	x,y = the arrays of data abscissas and ordinates
!  	b,c,d = arrays of spline coefficients computed by spline

!  	if  u  is not in the same interval as the previous call, then a
!  	binary search is performed to determine the proper interval.

		integer :: j, k
		real(IDP) :: dx
		integer, save :: i = 1

		if ( i >= n ) i = 1
		if ( u < x(i) .or. u > x(i+1) ) then

!  	binary search

			i = 1
			j = n+1
			do while ( j > i+1 )
				k = (i+j)/2
				if ( u < x(k) ) then
					j = k
				else
					i = k
				end if
			end do

		end if

!  	evaluate spline

		dx = u - x(i)
		seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

	end function seval
	function erf(x)

		use param
		implicit none

		real(IDP) :: erf,x,sign_x,xh,one,four,pi,dn,du,ai,u,fac,arg,aih
		integer :: n,np,i

!  	this routine calculates the standard error fun..

		if (x == 0.0) then
			erf=0.0
		else
			sign_x=sign(1.0_IDP,x)
			xh=x
			if (abs(x) > 5.0_IDP) then
				erf=sign_x
			else
				if (x < 0.0) x=-x
				one=1.0
				four=4.0
				pi=atan(one)*four
				n=1000
				dn=n
				du=x/dn
				if (x > 3.0_IDP) du=3.0/dn
				np=n+1
				ai=0.0
				u=0.0

				do i=1,np
					fac=1.0
					if (i == 1 .or. i == np) fac=0.5
					arg=exp(-u*u)
					ai=ai+arg*fac
					u=u+du
				end do
				aih=ai
				erf=ai*du*2.0/sqrt(pi)*sign_x
				if (x > 3.0_IDP) then

					du=(x-3.0)/dn
					u=3.0
					ai=0.0
					do i=1,np
						fac=1.0
						if (i == 1 .or. i == np) fac=0.5
						arg=exp(-u*u)
						ai=ai+arg*fac
						u=u+du
					end do
					erf=ai*du*2.0/sqrt(pi)*sign_x+erf
				end if
				x=xh

			end if
		end if

	end function erf
end module ffunctions
