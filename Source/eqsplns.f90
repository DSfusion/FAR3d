	subroutine eqsplns(yeq,yeqr,yfar,m,m2,choice,smoo)

		use param
		use domain
		implicit none

		integer :: m,m2,mjeqp,mjeqm,ier,j,jj
		character(len=6) :: choice
		real(IDP) :: smoo,fac1,fac2,dspl,rjm,rjmp,fac,facp
		real(IDP), dimension(0:) :: yeq,yeqr
		real(IDP), dimension(0:) :: yfar
		real(IDP), dimension(:), allocatable :: xspl,dfspl,fnspl,yspl
		real(IDP), dimension(:,:), allocatable :: cspl
		real(IDP), dimension(:,:), allocatable :: wkspl

		interface
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
			subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
				use param
				implicit none
				integer :: nx,ic,ier
				real(IDP) :: sm
				real(IDP), dimension(:) :: x,f,df,y
				real(IDP), dimension(:,:) :: c
				real(IDP), dimension(:,:) :: wk
			end subroutine icsscu
		end interface

!		Calculates splines for the equilibrium variables in vmec		
		
		allocate (xspl(0:mjeq))
		allocate (dfspl(0:mjeq))
		allocate (fnspl(0:mjeq))
		allocate (yspl(0:mjeq))
		allocate (cspl(0:mjeq,3))
		allocate (wkspl(7,mjeq+3))

		mjeqp=mjeq+1
		mjeqm=mjeq-1
		xspl=rfar**2
		fac1=xspl(2)/(xspl(2)-xspl(1))
		fac2=-xspl(1)/(xspl(2)-xspl(1))
		dfspl(1:mjeq)=1./rfar(1:mjeq)**m2
		fnspl(1:mjeq)=yfar(1:mjeq)/rfar(1:mjeq)**m
		dfspl(0)=dfspl(1)
		fnspl(0)=fac1*fnspl(1)+fac2*fnspl(2)
		if (m == 0) then
			dfspl(0)=1.0_IDP
			fnspl(0)=yfar(0)
		end if
		yspl=fnspl
		ier=0
		if (choice == "spline") call spline(mjeqp,xspl,yspl,cspl(:,1),cspl(:,2),cspl(:,3))
		if (choice == "icsscu") call icsscu(xspl,fnspl,dfspl,mjeqp,smoo,yspl,cspl,mjeqp,wkspl,ier)
!		open (unit=6,file="farprt",status="old",POSITION="APPEND")
		if (ier == 129 .or. ier == 130 .or. ier == 131) write(6,'("  eqsplns: ier=",i5)') ier
!		close(6)
		do j=0,mj
			do jj=1,mjeqm
				if (rfar(jj) > r(j)) exit
			end do
			jj=jj-1
			dspl=r(j)**2-xspl(jj)
			rjm=1.0_IDP
			if (m /= 0) rjm=0.
			if (r(j) /= 0.) rjm=r(j)**m
			rjmp=0.
			if (m == 1) rjmp=1.0_IDP
			if (m /= 0 .and. r(j) /= 0.) rjmp=m*r(j)**(m-1)
			fac=((cspl(jj,3)*dspl+cspl(jj,2))*dspl+cspl(jj,1))*dspl+yspl(jj)
			facp=r(j)*((6.*cspl(jj,3)*dspl+4.*cspl(jj,2))*dspl+2.*cspl(jj,1))
			yeq(j)=rjm*fac
			yeqr(j)=rjmp*fac+rjm*facp
		end do

		deallocate (wkspl)
		deallocate (cspl)
		deallocate (yspl)
		deallocate (fnspl)
		deallocate (dfspl)
		deallocate (xspl)

	end subroutine eqsplns

	subroutine spline (n, x, y, b, c, d)

		use param
		implicit none

		integer :: n
		real(IDP), dimension(:) :: x, y, b, c, d

!  	the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  	for a cubic interpolating spline

!  	s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

!  	for  x(i) .le. x .le. x(i+1)

!  	input..

!  	n = the number of data points or knots (n.ge.2)
!  	x = the abscissas of the knots in strictly increasing order
!  	y = the ordinates of the knots

!  	output..

!  	b, c, d  = arrays of spline coefficients as defined above.

!  	using  p  to denote differentiation,

!  	y(i) = s(x(i))
!  	b(i) = sp(x(i))
!  	c(i) = spp(x(i))/2
!  	d(i) = sppp(x(i))/6  (derivative from the right)

!  	the accompanying function subprogram  seval  can be used
!  	to evaluate the spline.


		integer :: nm1, ib, i
		real(IDP) :: t

		nm1 = n-1

		if ( n > 2 ) then

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

			b(1) = -d(1)
			b(n) = -d(n-1)
			c(1) = 0.
			c(n) = 0.
			if ( n > 3 ) then
				c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
				c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
				c(1) = c(1)*d(1)**2/(x(4)-x(1))
				c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
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

!  	compute polynomial coefficients

			b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
			do i = 1, nm1
				b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
				d(i) = (c(i+1) - c(i))/d(i)
				c(i) = 3.*c(i)
			end do
			c(n) = 3.*c(n)
			d(n) = d(n-1)

		else if (n == 2) then

			b(1) = (y(2)-y(1))/(x(2)-x(1))
			c(1) = 0.
			d(1) = 0.
			b(2) = b(1)
			c(2) = 0.
			d(2) = 0.

		end if

	end subroutine spline

	subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)

!  icsscu-------s-------library 2---------------------------------------

!  function            - cubic spline data smoothing
!  usage               - call icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
!  parameters   x      - vector of length nx containing the abscissae
!  							 of the nx data points (x(i),f(i)) i=1,...,
!  							 nx (input). x must be ordered so that
!  							 x(i) .lt. x(i+1).
!  				 f      - vector of length nx containing the ordinates
!  							 (or function values) of the nx data points
!  							 (input).
!  				 df     - vector of length nx (input).
!  							 df(i) is the relative weight of data
!  							 point i (see parameter sm below).
!  				 nx     - number of elements in x, f, df, and y (input).
!  							 nx must be .ge. 2.
!  				 sm     - a non-negative number which controls the
!  							 extent of smoothing (input). the spline
!  							 function s is determined such that the
!  							 sum from 1 to nx of
!  							 ((s(x(i))-f(i))/df(i))**2 .le. sm,
!  							 where equality holds unless s describes
!  							 a straight line.
!  				 y,c    - spline coefficients (output). y is a vector
!  							 of length nx. c is an nx-1 by 3 matrix.
!  							 the value of the spline approximation
!  							 at t is
!  							 s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!  							 where x(i) .le. t .lt. x(i+1) and
!  							 d = t-x(i).
!  				 ic     - row dimension of matrix c in the calling
!  							 program (input). ic must be .ge. nx-1.
!  				 wk     - work area vector of length greater than or
!  							 equal to 7*nx+14.
!  				 ier    - error parameter
!  							 terminal error
!  							 	ier = 129. ic is less than nx-1.
!  							 	ier = 130, nx is less than 2.
!  							 	ier = 131, input abscissae are not ordered
!  							 				  so that x(1) .lt. x(2) ... .lt. x(nx).
!  precision           - single
!  req'd imsl routines - uertst
!  language            - fortran

!  latest revision     - february 16, 1976
!  							 dec

		use param
		implicit none

		integer :: nx,ic,ier,m2,np1,i,np3,j
		real(IDP) :: sm,p,h,f2,ff,g,onedh,e,hmg
		real(IDP), dimension(:) :: x,f,df,y
		real(IDP), dimension(:,:) :: c
		real(IDP), dimension(:,:) :: wk

		interface
			subroutine uertst(ier,name)
				implicit none
				integer :: ier
				character(len=6) :: name
			end subroutine uertst
		end interface

!  	check error conditions
		ier = 0
		if (ic < nx-1) then
			ier = 129
			call uertst(ier,"icsscu")
      else if (nx < 2) then
			ier = 130
			call uertst(ier,"icsscu")
		else
!  	set up working areas
      	m2 = nx+2
      	np1 = nx+1
      	wk(1,1) = 0.0
      	wk(1,2) = 0.0
      	wk(2,np1) = 0.0
      	wk(3,m2) = 0.0
      	wk(3,np1) = 0.0
      	wk(6,1) = 0.0
      	wk(6,2) = 0.0
      	wk(6,m2) = 0.0
      	wk(6,np1) = 0.0
      	p = 0.0
      	h = x(2)-x(1)
      	if (h <= 0.0) then
				ier = 131
				call uertst(ier,"icsscu")
			else
				f2 = -sm
				ff = (f(2)-f(1))/h
				if (nx > 2) then
					do i=3,nx
						g = h
						h = x(i)-x(i-1)
      				if (h <= 0.0) then
							ier = 131
							call uertst(ier,"icsscu")
							return
						end if
						onedh = 1.0/h
						e = ff
						ff = (f(i)-f(i-1))*onedh
						y(i) = ff-e
						wk(4,i) = 2.*(g+h)/3.
						wk(5,i) = h/3.0
						wk(3,i) = df(i-2)/g
						wk(1,i) = df(i)*onedh
						wk(2,i) = -df(i-1)/g-df(i-1)*onedh
					end do
					do i=3,nx
						c(i-1,1) = wk(1,i)*wk(1,i)+wk(2,i)*wk(2,i)+wk(3,i)*wk(3,i)
						c(i-1,2) = wk(1,i)*wk(2,i+1)+wk(2,i)*wk(3,i+1)
						c(i-1,3) = wk(1,i)*wk(3,i+2)
					end do
				end if
!  	next iteration
				do while (.true.)
					if (nx > 2) then
						do i=3,nx
							wk(2,i-1) = ff*wk(1,i-1)
							wk(3,i-2) = g*wk(1,i-2)
							wk(1,i) = 1.0/(p*c(i-1,1)+wk(4,i)-ff*wk(2,i-1)-g*wk(3,i-2))
							wk(6,i) = y(i)-wk(2,i-1)*wk(6,i-1)-wk(3,i-2)*wk(6,i-2)
							ff = p*c(i-1,2)+wk(5,i)-h*wk(2,i-1)
							g = h
							h = c(i-1,3)*p
						end do
						np3 = nx+3
						do i=3,nx
							j = np3-i
							wk(6,j) = wk(1,j)*wk(6,j)-wk(2,j)*wk(6,j+1)-wk(3,j)*wk(6,j+2)
						end do
					end if
					e = 0.0
					h = 0.0
!  	compute u and accumulate e
					do i=2,nx
						g = h
						h = (wk(6,i+1)-wk(6,i))/(x(i)-x(i-1))
						hmg = h-g
						wk(7,i) = hmg*df(i-1)*df(i-1)
						e = e+wk(7,i)*hmg
					end do
					g = -h*df(nx)*df(nx)
					wk(7,np1) = g
					e = e-g*h
					g = f2
					f2 = e*p*p
					if (f2 >= sm .or. f2 <= g) exit
					ff = 0.0
					h = (wk(7,3)-wk(7,2))/(x(2)-x(1))
					if (nx > 2) then
						do i=3,nx
							g = h
							h = (wk(7,i+1)-wk(7,i))/(x(i)-x(i-1))
							g = h-g-wk(2,i-1)*wk(1,i-1)-wk(3,i-2)*wk(1,i-2)
							ff = ff+g*wk(1,i)*g
							wk(1,i) = g
						end do
					end if
					h = e-p*ff
					if (h <= 0.0) exit
!  	update the lagrange multiplier p for the next iteration
					p = p+(sm-f2)/((sqrt(sm/e)+p)*h)
				end do
!  	if e less than or equal to s, compute the coefficients and return.
				np1 = nx-1
				do i=1,np1
					y(i) = f(i)-p*wk(7,i+1)
					c(i,2) = wk(6,i+1)
					wk(1,i) = y(i)
				end do
				wk(1,nx) = f(nx)-p*wk(7,nx+1)
				y(nx) = wk(1,nx)
				do i=2,nx
					h = x(i)-x(i-1)
					c(i-1,3) = (wk(6,i+1)-c(i-1,2))/(h+h+h)
					c(i-1,1) = (wk(1,i)-y(i-1))/h-(h*c(i-1,3)+c(i-1,2))*h
				end do
			end if
		end if

	end subroutine icsscu

	subroutine uertst(ier,name)

!  uertst---------------library 2---------------------------------------

!  function            - error message generation
!  usage               - call uertst(ier,name)
!  parameters   ier    - error parameter. type + n  where
!  								 type= 128 implies terminal error
!  										  64 implies warning with fix
!  										  32 implies warning
!  								 n   = error code relevant to calling routine
!  				 name   - input scalar (double precision on dec)
!  							 containing the name of the calling routine
!  							 as a 6-character literal string.
!  language            - fortran

!  latest revision     - october 1,1975
!  							 dec

		implicit none
		
		integer :: ier,ier1,ier2,i
		character(len=6) :: name
		integer, dimension(4) :: ibit=(/32,64,128,0/)
		character(len=5), dimension(4,4) :: ityp
!  	character(len=5), dimension(4,4) :: ityp= (/"warni","ng   ","     ","     ","warni","ng(wi","th fi","x)   ", &
!  															  "termi","nal  ","     ","     ","non-d","efine","d    ","     "/)

		ityp(:,1)=(/"warni","ng   ","     ","     "/)
		ityp(:,2)=(/"warni","ng(wi","th fi","x)   "/)
		ityp(:,3)=(/"termi","nal  ","     ","     "/)
		ityp(:,4)=(/"non-d","efine","d    ","     "/)

 		ier2=ier
		if (ier2 < ibit(1)) then
!  	non-defined
			ier1=4
		else if (ier2 < ibit(2)) then
!  	warning
			ier1=1
		else if (ier2 < ibit(3)) then
!  	warning(with fix)
			ier1=2
		else
!  	terminal
			ier1=3
		end if
!  	extract 'n'
		ier2=ier2-ibit(ier1)
!  	print error message
      write (6,'(" *** i m s l(uertst) ***  ",4a5,2x,a6,2x,i2," (ier = ",i3,")")') (ityp(i,ier1),i=1,4),name,ier2,ier

	end subroutine uertst
