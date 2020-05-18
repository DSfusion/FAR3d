	subroutine fitter(yeq,yeqr,yfar,ir,m)

		use param
		use domain
		implicit none

		integer :: ir,m
		real(IDP), dimension(0:) :: yeq,yeqr,yfar
		integer :: i,i1,j,jj,js,je,l,nr,irp1,kplus1,nrows,mf,lyf,np1,lwrk,liwrk,ifail
		real(IDP) :: pmin,pmax,ave,x,aux,dspl
		integer, dimension(1) :: ip
		integer, dimension(30) :: iwrk
		real(IDP), dimension(ir+1) :: flux2,w,temp
		real(IDP), dimension(mjeq-ir+3) :: xspl,yspl
		real(IDP), dimension(mjeq-ir+3,3) :: cspl
		real(IDP), dimension(1) :: xf
		real(IDP), dimension(3) :: yf
		real(IDP), dimension(5,5) :: a
		real(IDP), dimension(5) :: s
		real(IDP), dimension(0:4) :: xhfit
		real(IDP), dimension(600) :: wrk

		interface
			subroutine e02age(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, kplus1, liwrk, lwrk, lyf, m, mf, np1, nrows
				real(IDP), dimension(nrows,kplus1) :: a
				real(IDP), dimension(kplus1) :: s
				real(IDP), dimension(m) :: w, x, y
				real(IDP), dimension(lwrk) :: wrk
				real(IDP), dimension(mf) :: xf
				real(IDP), dimension(lyf) :: yf
				integer, dimension(mf) :: ip
				integer, dimension(liwrk) :: iwrk
			end subroutine e02age
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
		end interface

		irp1=ir+1
		flux2(1)=0.0_IDP
		do i=1,ir
			flux2(i+1)=rfar(i)*rfar(i)
		end do
!		write(*,'("irp1=",i3)') irp1
!		write(*,'("flux2=",1p10e13.6)') (flux2(i),i=1,irp1)

		do i=1,ir
			temp(i+1)=yfar(i)/rfar(i)**m
		end do
		temp(1)=5.*temp(5)-4.*temp(6)
!		write(*,'("temp=",1p10e13.6)') (temp(i),i=1,irp1)
		ave = 0.0_IDP
	 	do i=1,irp1
			ave = ave + abs(temp(i))
		end do
		ave = ave/irp1
		do i=1,irp1
			w(i) = 10.0_IDP
			if (abs(temp(i)) > 10.*ave) w(i) = ave/abs(temp(i))
		end do
		do i=1,4
			w(i)=1.e-10_IDP
		end do
!		write(*,'("w=",1p10e13.6)') (w(i),i=1,irp1)
		pmin   = flux2(1)
		pmax   = flux2(irp1)
		xf(1)  = flux2(irp1)
		yf(1)  = temp(irp1)
		aux    = yfar(irp1)/rfar(irp1)**m
		yf(2)  = (aux-temp(ir))/(2.*flux2(2))
		yf(3)  = (aux-2.*temp(irp1)+temp(ir))/(flux2(2)*flux2(2))
!		write(*,'("yf=",1p10e13.6)') (yf(i),i=1,3)
		w(irp1)= 0.0_IDP
		kplus1 = 5
		nrows  = 5
		lyf    = 3
		mf     = 1
		ip(1)  = 2
		np1    = 4
		lwrk   = 600
		liwrk  = 30
		ifail  = 0
		call e02age(ir,kplus1,nrows,pmin,pmax,flux2,temp,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

		xhfit(0) = 0.5*a(5,1) - a(5,3) + a(5,5)
		xhfit(1) = a(5,2) - 3.*a(5,4)
		xhfit(2) = 2.*a(5,3) - 8.*a(5,5)
		xhfit(3) = 4.*a(5,4)
		xhfit(4) = 8.*a(5,5)

		do j=0,mj
			if (r(j) > rfar(ir)) exit
			x = ((r(j)*r(j)-pmin)-(pmax-r(j)*r(j)))/(pmax-pmin)
			yeq(j) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
			yeqr(j) = ((4.*xhfit(4)*x+3.*xhfit(3))*x+2.*xhfit(2))*x+xhfit(1)
			yeqr(j) = 2.0*yeqr(j)/(pmax-pmin)
		end do
		je=j-1
!		write(*,'("r=",1p10e13.6)') (r(j),j=0,je)
!		write(*,'("yeq=",1p10e13.6)') (yeq(j),j=0,je)
		xspl(1)=rfar(ir-2)
		x = ((xspl(1)*xspl(1)-pmin)-(pmax-xspl(1)*xspl(1)))/(pmax-pmin)
		yspl(1) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
		xspl(2)=rfar(ir-1)
		x = ((xspl(2)*xspl(2)-pmin)-(pmax-xspl(2)*xspl(2)))/(pmax-pmin)
		yspl(2) = (((xhfit(4)*x+xhfit(3))*x+xhfit(2))*x+xhfit(1))*x+xhfit(0)
		if (m == 0) then
			do j=0,je
				yeqr(j)=2.*r(j)*yeqr(j)
			end do
		else if (m == 1) then
			do j=0,je
				yeqr(j)=yeq(j)+2.*r(j)*r(j)*yeqr(j)
				yeq(j)=r(j)*yeq(j)
			end do
			yspl(1)=xspl(1)*yspl(1)
			yspl(2)=xspl(2)*yspl(2)
		else
			do j=0,je
				yeqr(j)=m*yeq(j)*r(j)**(m-1)+2.*yeqr(j)*r(j)**(m+1)
				yeq(j)=yeq(j)*r(j)**m
			end do
			yspl(1)=yspl(1)*xspl(1)**m
			yspl(2)=yspl(2)*xspl(2)**m
		end if

	 	do i=ir,mjeq
			xspl(i-ir+3)=rfar(i)
			yspl(i-ir+3)=yfar(i)
		end do
		nr=mjeq-ir+3
!		write(*,'("xspl=",1p10e13.6)') (xspl(i),i=1,nr)
!		write(*,'("yspl=",1p10e13.6)') (yspl(i),i=1,nr)
		call spline(nr,xspl,yspl,cspl(:,1),cspl(:,2),cspl(:,3))

		js=3
		do j=je+1,mj
			do jj=js,nr
				if (xspl(jj) > r(j)) exit
			end do
			js=jj
			jj=jj-1
			dspl=r(j)-xspl(jj)
			yeq(j)=((cspl(jj,3)*dspl+cspl(jj,2))*dspl+cspl(jj,1))*dspl+yspl(jj)
			yeqr(j)=(3.*cspl(jj,3)*dspl+2.*cspl(jj,2))*dspl+cspl(jj,1)
		end do
!		write(*,'("r=",1p10e13.6)') (r(j),j=je+1,mj)
!		write(*,'("yeq=",1p10e13.6)') (yeq(j),j=je+1,mj)

	end subroutine fitter	

	subroutine e02ade(m,kplus1,nrows,x,y,w,work1,work2,a,s,ifail)

		use param
		implicit none

!  	nag library subroutine  e02ade

!  	e02ade  computes weighted least-squares polynomial
!  	approximations to an arbitrary set of data points.

!  	forsythe-clenshaw method with modifications due to
!  	reinsch and gentleman.

!  	uses nag library routine  p01aae.
!  	uses basic external function  sqrt.

!  	started - 1973.
!  	completed - 1976.
!  	author - mgc and jgh.

!  	work1  and  work2  are workspace areas.
!  	work1(1, r)  contains the value of the  r th  weighted
!  	residual for the current degree  i.
!  	work1(2, r)  contains the value of  x(r)  transformed
!  	to the range  -1  to  +1.
!  	work1(3, r)  contains the weighted value of the current
!  	orthogonal polynomial (of degree  i)  at the  r th
!  	data point.
!  	work2(1, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the current orthogonal polynomial
!  	(of degree  i).
!  	work2(2, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the previous orthogonal polynomial
!  	(of degree  i - 1).

!  	nag copyright 1975
!  	mark 5 release
!  	mark 6 revised  ier-84
!  	mark 11.5(f77) revised. (sept 1985.)

!  	check that the values of  m  and  kplus1  are reasonable


!  	.. parameters ..
		character(len=6), parameter :: srname="e02ade"
!  	.. scalar arguments ..
		integer :: ifail, kplus1, m, nrows
!  	.. array arguments ..
		real(IDP), dimension(nrows,kplus1) :: a
		real(IDP), dimension(kplus1) :: s
		real(IDP), dimension(m) :: w, x, y
		real(IDP), dimension(3,m) :: work1
		real(IDP), dimension(2,kplus1) :: work2
!  	.. local scalars ..
		real(IDP) :: alpip1, betai, bj, bjp1, bjp2, ci, d, df, di, dim1, dj, epsr, factor, pij, sigmai, wrpr, wrprsq, x1, xcapr, xm
		integer :: i, ierror, iplus1, iplus2, j, jplus1, jplus2, jrev, k, mdist, r
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic sqrt
!  	.. executable statements ..
		if (kplus1 < 1 .or. m < kplus1) then
			ierror = 4
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		k = kplus1 - 1

!  	test the validity of the data.

!  	check that the weights are strictly positive.

		do r = 1, m
			if (w(r) <= 0.0_IDP) then
				ierror = 1
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
		end do

!  	check that the values of  x(r)  are non-decreasing and
!  	determine
!  	the number  (mdist)  of distinct values of  x(r).

		mdist = 1
		do r = 2, m
			if (x(r) < x(r-1)) then
				ierror = 2
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			if (x(r) == x(r-1)) cycle
			mdist = mdist + 1
		end do

!  	if the  x(r)  all have the same value, i.e.  mdist = 1,
!  	the normalization of the independent variable is not
!  	possible.

		if (mdist == 1) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

!  	if the number of distinct values of  x(r)  fails to exceed
!  	the maximum degree  k  there is no unique polynomial
!  	approximation of that degree.

		if (mdist <= k) then
			ierror = 4
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

!  	check that  nrows  has been set sufficiently large.

		if (nrows < kplus1) then
			ierror = 5
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if

		x1 = x(1)
		xm = x(m)
		d = xm - x1

!  	the initial values  work1(1, r)  (r = 1, 2, ..., m)  of the
!  	weighted residuals and the values  work1(2, r)  (r = 1, 2,
!  	..., m)
!  	of the normalized independent variable are computed.  note
!  	that
!  	work1(2, r)  is computed from the expression below rather
!  	than the
!  	more natural form

!  	(2.0*x(r) - x1 - xm)/d,

!  	since the former guarantees the computed value to differ from
!  	the true value by at most  4.0*machine accuracy,  whereas the
!  	latter has no such guarantee.

		do r = 1, m
			work1(1,r) = w(r)*y(r)
			work1(2,r) = ((x(r)-x1)-(xm-x(r)))/d
		end do
		i = 1
		betai = 0.0_IDP
		do iplus1 = 1, kplus1

!  		set starting values for degree  i.

			iplus2 = iplus1 + 1
			if (iplus1 < kplus1) then
				do jplus1 = iplus2, kplus1
					a(iplus1,jplus1) = 0.0_IDP
				end do
				work2(1,iplus2) = 0.0_IDP
				work2(2,iplus2) = 0.0_IDP
			end if
			alpip1 = 0.0_IDP
			ci = 0.0_IDP
			di = 0.0_IDP
			a(i,iplus1) = 0.0_IDP
			work2(1,iplus1) = 1.0_IDP
			if (kplus1 > 1) work2(2,1) = work2(1,2)
			do r = 1, m
				xcapr = work1(2,r)

!  			the weighted value  work1(3, r)  of the orthogonal
!  			polynomial of
!  			degree  i  at  x = x(r)  is computed by recurrence from its
!  			chebyshev-series representation.

				if (iplus1 == 1) then
					wrpr = w(r)*0.5*work2(1,1)
					work1(3,r) = wrpr
				else
					j = iplus2
					if (xcapr > 0.5) then

!  			reinsch*s modified recurrence.

						factor = 2.0*(1.0-xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							dj = work2(1,j) + dj - factor*bj
							bj = bj + dj
						end do
						wrpr = w(r)*(0.5*work2(1,1)+dj-0.5*factor*bj)
						work1(3,r) = wrpr
					else if (xcapr >= -0.5) then

!  			clenshaw*s original recurrence.

						factor = 2.0*xcapr
						bjp1 = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							bjp2 = bjp1
							bjp1 = bj
							bj = work2(1,j) - bjp2 + factor*bjp1
						end do
						wrpr = w(r)*(0.5*work2(1,1)-bjp1+0.5*factor*bj)
						work1(3,r) = wrpr
					else

!  			gentleman*s modified recurrence.

						factor = 2.0*(1.0+xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 1, i
							j = j - 1
							dj = work2(1,j) - dj + factor*bj
							bj = dj - bj
						end do
						wrpr = w(r)*(0.5*work2(1,1)-dj+0.5*factor*bj)
						work1(3,r) = wrpr
					end if
				end if

!  			the coefficient  ci  of the  i th  orthogonal polynomial and
!  			the
!  			coefficients  alpip1  and  beta i  in the
!  			three-term recurrence relation for the orthogonal
!  			polynomials are computed.

				wrprsq = wrpr**2
				di = di + wrprsq
				ci = ci + wrpr*work1(1,r)
				alpip1 = alpip1 + wrprsq*xcapr
			end do
			ci = ci/di
			if (iplus1 /= 1) betai = di/dim1
			alpip1 = 2.0*alpip1/di

!  		the weighted residuals  work1(1, r)  (r = 1, 2, ..., m)  for
!  		degree  i  are computed, together with their sum of squares,
!  		sigmai.

			sigmai = 0.0
			do r = 1, m
				epsr = work1(1,r) - ci*work1(3,r)
				work1(1,r) = epsr
				sigmai = sigmai + epsr**2
			end do

!  		the root mean square residual  s(i + 1)  for degree  i  is
!  		theoretically undefined if  m = i + 1  (the condition for the
!  		polynomial to pass exactly through the data points).  should
!  		this
!  		case arise the r.m.s. residual is set to zero.

			if (iplus1 < m) then
				df = m - iplus1
				s(iplus1) = sqrt(sigmai/df)
			else
				s(iplus1) = 0.0_IDP
			end if

!  		the chebyshev coefficients  a(i + 1, 1), a(i + 1, 2), ...,
!  		a(i + 1, i + 1)  together with the coefficients
!  		work2(1, 1), work2(1, 2), ..., work2(1, i + 1),   in the
!  		chebyshev-series representation of the  i th  orthogonal
!  		polynomial are computed.

			do jplus1 = 1, iplus1
				jplus2 = jplus1 + 1
				pij = work2(1,jplus1)
				a(iplus1,jplus1) = a(i,jplus1) + ci*pij
				if (jplus2 > kplus1) exit
				work2(1,jplus1) = work2(1,jplus2) + work2(2,jplus1) - alpip1*pij - betai*work2(2,jplus2)
				work2(2,jplus2) = pij
			end do
			if (iplus1 < kplus1) then
				dim1 = di
				i = iplus1
			end if
		end do
		ifail = 0

	end subroutine e02ade

	subroutine e02afe(nplus1,f,a,ifail)

		use param
		implicit none

!  	nag library subroutine  e02afe

!  	e02afe  computes the coefficients of a polynomial,
!  	in its chebyshev-series form, which interpolates
!  	(passes exactly through) data at a special set of
!  	points.  least-squares polynomial approximations
!  	can also be obtained.

!  	clenshaw method with modifications due to reinsch
!  	and gentleman.

!  	uses nag library routines  p01aae  and  x01aae.
!  	uses basic external function  sin.

!  	started - 1973.
!  	completed - 1976.
!  	author - mgc and jgh.

!  	nag copyright 1975
!  	mark 5 release
!  	mark 5b revised  ier-73
!  	mark 11.5(f77) revised. (sept 1985.)


!  	.. parameters ..
		character(len=6), parameter :: srname="e02afe"
!  	.. scalar arguments ..
		integer :: ifail, nplus1
!  	.. array arguments ..
		real(IDP), dimension(nplus1) :: a, f
!  	.. local scalars ..
		real(IDP) :: bk, bkp1, bkp2, dk, f0, factor, fli, fln, halffn, piby2n
		integer :: i, ierror, iplus1, j, k, krev, n, n2, nless1
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec
!  	.. data statements ..
		real(IDP) :: pi=3.14159265358979323846_IDP

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic sin
!  	.. executable statements ..
		ierror = 0
		if (nplus1 < 2) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
		else if (nplus1 == 2) then
			a(1) = f(1) + f(2)
			a(2) = 0.5*(f(1)-f(2))
			ifail = 0
		else
			n = nplus1 - 1
			fln = n
			n2 = 2*n
			nless1 = n - 1
			piby2n = 0.5*pi/fln
			f0 = f(1)
			halffn = 0.5*f(nplus1)
			do iplus1 = 1, nplus1
				i = iplus1 - 1
				k = nplus1
				j = 3*i
				if (j > n2) then

!  			gentleman*s modified recurrence.

					fli = n - i
					factor = 4.0*(sin(piby2n*fli))**2
					dk = halffn
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						dk = f(k) - dk + factor*bk
						bk = dk - bk
					end do
					a(iplus1) = (f0-2.0*dk+factor*bk)/fln
				else if (j >= n) then

!  			clenshaw*s original recurrence.

					fli = n - 2*i
					factor = 2.0*sin(piby2n*fli)
					bkp1 = 0.0_IDP
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						bkp2 = bkp1
						bkp1 = bk
						bk = f(k) - bkp2 + factor*bkp1
					end do
					a(iplus1) = (f0-2.0*bkp1+factor*bk)/fln
				else

!  			reinsch*s modified recurrence.

					fli = i
					factor = 4.0*(sin(piby2n*fli))**2
					dk = halffn
					bk = halffn
					do krev = 1, nless1
						k = k - 1
						dk = f(k) + dk - factor*bk
						bk = bk + dk
					end do
 					a(iplus1) = (f0+2.0*dk-factor*bk)/fln
				end if
			end do
			a(nplus1) = 0.5*a(nplus1)
			ifail = 0
		end if

	end subroutine e02afe

	subroutine e02age(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 9 revised. ier-314 (sep 1981).
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 14c revised. ier-878 (nov 1990).

!  	**************************************************************

!  	npl algorithms library routine confit

!  	created  20/8/1979      updated  16/5/80     release no. 00/05

!  	author... gerald t anthony.
!  	national physical laboratory,
!  	teddington, middlesex tw11 0lw,
!  	england.

!  	**************************************************************

!  	e02age calls aewe01 to determine a polynomial mu(x) which
!  	interpolates the given constraints and a polynomial nu(x)
!  	which has value zero where a constrained value is specified.
!  	it then calls adze02 to fit y-mu(x) as a polynomial in x
!  	with factor nu(x). finally the coefficients of mu are added
!  	to these fits to give the coefficients of the constrained
!  	fits to y. all polynomials are expressed in chebyshev
!  	series form

!  	input parameters
!  		m        the number of data points to be fitted
!  		kplus1   fits with up to kplus1 coefficients are required
!  		nrows    first dimension of array a where coefficients are
!  					to be stored
!  		xmin,    end points of the range of the
!  		xmax     independent variable
!  		x, y, w  arrays of data values of the independent variable,
!  					dependent variable and weight, respectively
!  		mf       number of x values at which a constraint
!  					is specified
!  		xf       array of values of the independent
!  					variable at which constraints are
!  					specified
!  		yf       array of specified values and derivatives of the
!  					dependent variable in the order
!  						y1, y1 deriv, y1 2nd deriv,...., y2,....
!  		lyf      dimension of array yf
!  		ip       integer array of degrees of derivatives
!  					specified at each point xf

!  	output parameters
!  		a        on exit, 2 parameter array containing the
!  					coefficients of the chebyshev series
!  					representation of the fits, a(i+1, j+1)
!  					contains the coefficient of tj in the fit
!  					of degree i, i = n,n+1,...,k,  j =
!  					0,1,...,i  where n = np1 - 1
!  		s        on exit, array containing the r.m.s. residual for
!  					each degree of fit from n to k
!  		np1      on exit, contains n + 1, where n is the
!  					total number of interpolation conditions

!  		ifail    failure indicator
!  					 0 - successful termination
!  					 1 - at least one of the following conditions
!  						  has been violated
!  						  lyf    at least n
!  						  lwrk   at least 2*n + 2 + the larger of
!  									4*m + 3*kplus1 and 8*np1 +
!  									5*imax + mf - 3 where imax =
!  									1 + max(ip(i))
!  						  liwrk  at least 2*mf + 2
!  						  kplus1 at least np1
!  						  m      at least 1
!  						  nrows  at least kplus1
!  						  mf     at least 1
!  					 2 - for some i, ip(i) is less than 0
!  					 3 - xmin is not strictly less than xmax
!  						  or for some i, xf(i) is not in range
!  						  xmin to xmax or the xf(i) are not
!  						  distinct
!  					 4 - for some i, x(i) is not in range xmin to xmax
!  					 5 - the x(i) are not non-decreasing
!  					 6 - the number of distinct values of x(i) with
!  						  non-zero weight is less than kplus1 - np1
!  					 7 - aewe01 has failed to converge, ie
!  						  the constraint cannot be satisfied
!  						  with sufficient accuracy

!  	workspace parameters
!  		wrk      real workspace array
!  		lwrk     dimension of wrk.   lwrk must be at least
!  					2*n + 2 + the larger of
!  					4*m + 3*kplus1 and 8*np1 + 5*imax + mf - 3
!  					where imax = 1 + max(ip(i))
!  		iwrk     integer workspace array
!  		liwrk    dimension of iwrk.   liwrk must be at least
!  					2*mf + 2

!  	.. parameters ..
		character(len=6), parameter :: srname="e02age"
!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, kplus1, liwrk, lwrk, lyf, m, mf, np1, nrows
!  	.. array arguments ..
		real(IDP), dimension(nrows,kplus1) :: a
		real(IDP), dimension(kplus1) :: s
		real(IDP), dimension(m) :: w, x, y
		real(IDP), dimension(lwrk) :: wrk
		real(IDP), dimension(mf) :: xf
		real(IDP), dimension(lyf) :: yf
		integer, dimension(mf) :: ip
		integer, dimension(liwrk) :: iwrk
!  	.. local scalars ..
		real(IDP) :: amuj, xi, xmu
		integer :: i, ierror, im1, imax, ipi, iymux, j, lw, mdist, n, nanu, neps, nser, nwrk, nwrk1, nwrk2
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
!  	.. external subroutines ..
			subroutine aewe01(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, itmax, itmin, liwrk, lwrk, m, n, np1
				real(IDP), dimension(n) :: a, y
				real(IDP), dimension(np1) :: b
				real(IDP), dimension(lwrk) :: wrk
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
				integer, dimension(liwrk) :: iwrk
			end subroutine aewe01
			subroutine adze02(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: ifail, inup1, kall, kplus1, mfirst, mlast, mtot, ndv, nrows
				real(IDP), dimension(ndv,nrows,kplus1) :: a
				real(IDP), dimension(ndv,mlast) :: eps, y
				real(IDP), dimension(inup1) :: nu
				real(IDP), dimension(ndv,kplus1) :: s
				real(IDP), dimension(kplus1) :: serr
				real(IDP), dimension(mlast) :: w, x
				real(IDP), dimension(2,mtot) :: work1
				real(IDP), dimension(2,kplus1) :: work2
			end subroutine adze02
			subroutine e02ake(np1,xmin,xmax,a,ia1,la,x,result,ifail)
				use param
				implicit none
				real(IDP) :: result, x, xmax, xmin
				integer :: ia1, ifail, la, np1
				real(IDP), dimension(la) :: a
			end subroutine e02ake
		end interface

!  	.. data statements ..
		real(IDP) :: one=1.0_IDP, zero=0.0_IDP
!  	.. executable statements ..
		if (mf < 1) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		imax = 0
		np1 = 1
		do i = 1, mf
			ipi = ip(i) + 1
			if (ipi < 1) then
				ierror = 2
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			if (ipi > imax) imax = ipi
			np1 = np1 + ipi
		end do
		n = np1 - 1
		if (lyf < n .or. liwrk < 2*mf+2) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		i = 4*m + 3*kplus1
		lw = 8*np1 + 5*imax + mf - 3
		if (lw < i) lw = i
		lw = lw + 2*np1
		if (lw > lwrk .or. np1 > kplus1 .or. m < 1) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		nanu = np1 + 1
		nwrk = nanu + np1
		if (xmax <= xmin) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (xf(1) > xmax .or. xf(1) < xmin) then
			ierror = 3
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (mf > 1) then
			do i = 2, mf
				xi = xf(i)
				if (xi > xmax .or. xi < xmin) then
					ierror = 3
					ifail = p01abe(ifail,ierror,srname,0,p01rec)
					return
				end if
				im1 = i - 1
				do j = 1, im1
					if (xi == xf(j)) then
						ierror = 3
						ifail = p01abe(ifail,ierror,srname,0,p01rec)
						return
					end if
				end do
			end do
		end if
		xmu = xmin
		if (x(1) == xmin) xmu = xmax
		mdist = 0
		do i = 1, m
			xi = x(i)
			if (xi == xmu .or. w(i) == zero) cycle
			do j = 1, mf
				if (xi == xf(j)) exit
			end do
			if (j > mf) mdist = mdist + 1
			xmu = xi
		end do
		if (mdist < kplus1-n) then
			ierror = 6
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		iwrk(1) = 1
		ierror = 1
		call aewe01(mf,xmin,xmax,xf,yf,ip,n,np1,5,20,wrk,wrk(nanu),wrk(nwrk),lwrk-nwrk+1,iwrk,liwrk,ierror)
		if (ierror /= 0) then
			ierror = 7
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		do i = 1, m
			ierror = 1
			call e02ake(n,xmin,xmax,wrk,1,np1,x(i),xmu,ierror)
			if (ierror /= 0) then
				ierror = 4
				ifail = p01abe(ifail,ierror,srname,0,p01rec)
				return
			end if
			iymux = nwrk + i - 1

!  		store y - mu(x) at ith data point

			wrk(iymux) = y(i) - xmu
		end do
		nwrk1 = nwrk + m
		nwrk2 = nwrk1 + 2*m
		nser = nwrk2 + 2*kplus1
		neps = nser + kplus1
		ierror = 1
		call adze02(1,m,m,kplus1,nrows,1,1,x,wrk(nwrk),w,xmin,xmax,np1,wrk(nanu),wrk(nwrk1),wrk(nwrk2),a,s,wrk(nser),wrk(neps),ierror)
		if (ierror == 0) then
			do j = 1, n
				amuj = wrk(j)
				do  i = np1, kplus1
					a(i,j) = a(i,j) + amuj
				end do
			end do
		else
			ierror = ierror + 3
			if (ierror == 8) ierror = 1
		end if
		ifail = p01abe(ifail,ierror,srname,0,p01rec)

	end subroutine e02age

	subroutine e02ake(np1,xmin,xmax,a,ia1,la,x,result,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine tval1c

!  	created 19/3/79    updated 6/7/79    release no. 00/04.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1      np1 = n + 1. n is the degree of the
!  					chebyshev series
!  		xmin     minimum value of x
!  		xmax     maximum value of x
!  		a        the array where the coefficients are stored
!  		ia1      the address increment of a
!  		la       dimension of a
!  		x        unnormalized argument in the range (xmin, xmax)

!  	output parameters
!  		result   value of the summation
!  		ifail    error indicator

!  	np1 chebyshev coefficients a0, a1, ..., an, are
!  	stored in the array a in positions 1, 1+ia1, 1+2*ia1, ...,
!  	1+n*ia1, where n = np1 - 1.
!  	ia1 must not be negative.
!  	la must be at least equal to 1 + n*ia1.
!  	the value of the polynomial of degree n
!  	a0t0(xcap)/2 + a1t1(xcap) + a2t2(xcap) + + ... + antn(xcap),
!  	is calculated for the argument xcap, where xcap is
!  	the normalized value of x in the range (xmin, xmax),
!  	storing it in result.
!  	unless the routine detects an error, ifail contains
!  	zero on exit.
!  	ifail = 1 indicates at least one of the restrictions on
!  		input parameters is violated - ie
!  	np1  >  0
!  	ia1  >=  0
!  	la  >=  1 + n * ia1
!  	xmin  <  xmax
!  	ifail = 2 indicates that
!  	x does not satisfy the restriction xmin  <=  x  <=  xmax.
!  	the recurrence relation by clenshaw, modified by reinsch
!  	and gentleman, is used.

!  	.. parameters ..
		character(len=6), parameter :: srname="e02ake"
!  	.. scalar arguments ..
		real(IDP) :: result, x, xmax, xmin
		integer :: ia1, ifail, la, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
!  	.. local scalars ..
		real(IDP) :: xcap
		integer :: ierror
!  	.. local arrays ..
		character(len=1), dimension(1) :: p01rec

		interface
!  	.. external functions ..
			integer function p01abe(ifail,ierror,srname,nrec,rec)
				implicit none
				integer :: ierror, ifail, nrec
				character(*) :: srname
				character(*), dimension(:) :: rec
			end function p01abe
			subroutine akye02(xmin,xmax,x,xcap)
				use param
				implicit none
				real(IDP) :: x, xcap, xmax, xmin
			end subroutine akye02
			subroutine akze02(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine akze02
		end interface

!  	.. executable statements ..
		if (np1 < 1 .or. ia1 < 1 .or. la < 1+(np1-1)*ia1 .or. xmax <= xmin) then
			ierror = 1
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		if (x > xmax .or. x < xmin) then
			ierror = 2
			ifail = p01abe(ifail,ierror,srname,0,p01rec)
			return
		end if
		ierror = 0
		call akye02(xmin,xmax,x,xcap)
		call akze02(np1,a,ia1,la,xcap,result)
		ifail = p01abe(ifail,ierror,srname,0,p01rec)

	end subroutine e02ake

	subroutine adze02(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)

		use param
		implicit none

!  	mark 7 release. nag copyright 1978.
!  	mark 8 revised. ier-228 (apr 1980).
!  	mark 11.5(f77) revised. (sept 1985.)

!  	adze02  computes weighted least-squares polynomial
!  	approximations to an arbitrary set of data points,
!  	with, if required, several sets of values of the
!  	dependent variable.

!  	forsythe-clenshaw method with modifications due to
!  	reinsch and gentleman.

!  	started - 1973.
!  	completed - 1978.
!  	authors - mgc and gta.

!  	work1  and  work2  are workspace areas.
!  	work1(1, r)  contains the value of  x(r)  transformed
!  	to the range  -1  to  +1.
!  	work1(2, r)  contains the weighted value of the current
!  	orthogonal polynomial (of degree  i)  at the  r th
!  	data point.
!  	work2(1, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the current orthogonal polynomial
!  	(of degree  i).
!  	work2(2, j)  contains the coefficient of the chebyshev
!  	polynomial of degree  j - 1  in the chebyshev-series
!  	representation of the previous orthogonal polynomial
!  	(of degree  i - 1).

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, inup1, kall, kplus1, mfirst, mlast, mtot, ndv, nrows
!  	.. array arguments ..
		real(IDP), dimension(ndv,nrows,kplus1) :: a
		real(IDP), dimension(ndv,mlast) :: eps, y
		real(IDP), dimension(inup1) :: nu
		real(IDP), dimension(ndv,kplus1) :: s
		real(IDP), dimension(kplus1) :: serr
		real(IDP), dimension(mlast) :: w, x
		real(IDP), dimension(2,mtot) :: work1
		real(IDP), dimension(2,kplus1) :: work2
!  	.. local scalars ..
		real(IDP) :: alpip1, betai, bj, bjp1, bjp2, cil, d, df, di, dim1, dj, epslr, factor, pij, sigmai, wr, wrpr, wrprsq, x1, &
						 xcapr, xm
		integer :: i, ii, im1, inu, iplus1, iplus2, j, jplus1, jplus2, jrev, k, l, m, mdist, mr, r
		logical :: wnz
!  	.. local arrays ..
		real(IDP), dimension(10) :: ci
!  	.. intrinsic functions ..
!  	intrinsic sqrt
!  	.. executable statements ..
		k = kplus1 - 1
		inu = inup1 - 1

!  	test the validity of the data.

!  	check input parameters.

		m = mlast - mfirst + 1
		i = kplus1 - inu
		if (mfirst < 1 .or. inup1 < 1 .or. kplus1 < inup1 .or. m < i .or. ndv < 1 .or. (kall /= 1 .and. kall /= 0)) then
			ifail = 5
			return
		end if

!  	check that the values of x(r) are non-decreasing and
!  	determine the number (mdist) of distinct values of x(r)
!  	with non-zero weight

		mdist = 1
		if (w(mfirst) == 0.0_IDP) mdist = 0
		l = mfirst + 1
		if (l <= mlast) then
			wnz = w(mfirst)  /=  0.0_IDP
			do r = l, mlast
				if (x(r) < x(r-1)) then
					ifail = 2
					return
				end if
				if (x(r) > x(r-1)) wnz = .false.
				if (w(r) == 0.0_IDP .or. wnz) cycle
				mdist = mdist + 1
				wnz = .true.
			end do
		end if

!  	check that xmin < xmax and that xmin and xmax span the data
!  	x values.

		if (xmin > x(mfirst) .or. xmax < x(mlast) .or. xmin >= xmax) then
			ifail = 1
			return
		end if

!  	if the number of distinct values of  x(r)  with non-zero
!  	weight is less than the number of independent coefficients
!  	in the fit of maximum degree  k  there is no unique
!  	polynomial
!  	approximation of that degree.

		l = k - inu
		if (mdist <= l) then
			ifail = 3
			return
		end if

!  	check that  nrows  has been set sufficiently large.

		if (kall == 1 .and. nrows < kplus1) then
			ifail = 5
			return
		end if
		if (inup1 > 1) then

!  	normalize the forcing factor so that its leading coefficient
!  	is unity, checking that this coefficient was not zero.

			di = nu(inup1)
			if (di == 0.0_IDP) then
				ifail = 4
				return
			end if
			do i = 1, inup1
				work2(1,i) = nu(i)/di
				work2(2,i) = 0.0_IDP
			end do
		end if

		x1 = xmin
		xm = xmax
		d = xm - x1

!  	the initial values of eps(l,r) (l = 1,2,....ndv and r =
!  	mfirst, mfirst+1,....mlast) of the weighted residuals and
!  	the values work1(1,r)(r=1,2...m) of the normalized
!  	independent variable are computed. n.b. work1(1,r) is
!  	computed from the expression below rather than the more
!  	natural form   (2.0*x(r) - x1 - xm)/d
!  	since the former guarantees the computed value to differ from
!  	the true value by at most  4.0*machine accuracy,  whereas the
!  	latter has no such guarantee.

!  	mdist is now used to record the total number of data points
!  	with non-zero weight.

		mdist = 0
		do r = mfirst, mlast
			wr = w(r)
			if (wr /= 0.0_IDP) mdist = mdist + 1
			mr = r - mfirst + 1
			do l = 1, ndv
				eps(l,r) = wr*y(l,r)
			end do
			work1(1,mr) = ((x(r)-x1)-(xm-x(r)))/d
		end do
		im1 = inu*kall + 1
		betai = 0.0_IDP
		do jplus1 = 1, kplus1
			serr(jplus1) = 0.0_IDP
			do l = 1, ndv
				a(l,im1,jplus1) = 0.0_IDP
			end do
		end do
		do iplus1 = inup1, kplus1

!  		set starting values for degree  i.

			ii = (iplus1-1)*kall + 1
			iplus2 = iplus1 + 1
			if (iplus1 < kplus1) then
				if (kall == 1) then
					do jplus1 = iplus2, kplus1
						do l = 1, ndv
							a(l,ii,jplus1) = 0.0_IDP
						end do
					end do
				end if
				work2(1,iplus2) = 0.0_IDP
				work2(2,iplus2) = 0.0_IDP
			end if
			alpip1 = 0.0_IDP
			di = 0.0_IDP
			do l = 1, ndv
				ci(l) = 0.0_IDP
			end do
			work2(1,iplus1) = 1.0_IDP
			if (kplus1 > 1) work2(2,1) = work2(1,2)
			do r = mfirst, mlast
				if (w(r) == 0.0_IDP) cycle
				mr = r - mfirst + 1
				xcapr = work1(1,mr)

!  			the weighted value work1(2, r)  of the orthogonal polynomial
!  			of degree i at x = x(r) is computed by recurrence from its
!  			chebyshev-series representation.

				if (iplus1 == 1) then
					wrpr = w(r)*0.5*work2(1,1)
					work1(2,mr) = wrpr
				else
					j = iplus2
					if (xcapr > 0.5_IDP) then

!  					reinsch*s modified recurrence.

						factor = 2.0*(1.0-xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							dj = work2(1,j) + dj - factor*bj
							bj = bj + dj
						end do
						wrpr = w(r)*(0.5*work2(1,1)+dj-0.5*factor*bj)
						work1(2,mr) = wrpr
					else if (xcapr >= -0.5_IDP) then

!  					clenshaw*s original recurrence.

						factor = 2.0*xcapr
						bjp1 = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							bjp2 = bjp1
							bjp1 = bj
							bj = work2(1,j) - bjp2 + factor*bjp1
						end do
						wrpr = w(r)*(0.5*work2(1,1)-bjp1+0.5*factor*bj)
						work1(2,mr) = wrpr
					else

!  					gentleman*s modified recurrence.

						factor = 2.0*(1.0+xcapr)
						dj = 0.0_IDP
						bj = 0.0_IDP
						do jrev = 2, iplus1
							j = j - 1
							dj = work2(1,j) - dj + factor*bj
							bj = dj - bj
						end do
						wrpr = w(r)*(0.5*work2(1,1)-dj+0.5*factor*bj)
						work1(2,mr) = wrpr
					end if
				end if

!  			the coefficients ci(l) of the i th orthogonal polynomial
!  			l=1,2....ndv and the coefficients alpip1 and betai in the
!  			three-term recurrence relation for the orthogonal
!  			polynomials are computed.

				wrprsq = wrpr**2
				di = di + wrprsq
				do l = 1, ndv
					ci(l) = ci(l) + wrpr*eps(l,r)
				end do
				alpip1 = alpip1 + wrprsq*xcapr
			end do
			do l = 1, ndv
				ci(l) = ci(l)/di
			end do
			if (iplus1 /= inup1) betai = di/dim1
			alpip1 = 2.0*alpip1/di

!  		the weighted residuals eps(l,r)(l=1,2....ndv and r=mfirst,
!  		mfirst+1....mlast) for degree i are computed, together
!  		with their sum of squares, sigmai

			df = mdist - (iplus1-inu)
			do l = 1, ndv
				cil = ci(l)
				sigmai = 0.0_IDP
				do r = mfirst, mlast
					if (w(r) == 0.0_IDP) cycle
					mr = r - mfirst + 1
					epslr = eps(l,r) - cil*work1(2,mr)
					eps(l,r) = epslr
					sigmai = sigmai + epslr**2
				end do

!  			the root mean square residual  s(l, i + 1)  for degree  i
!  			is theoretically undefined if  m = i + 1 - inu  (the
!  			condition for the polynomial to pass exactly through the
!  			data points). should this case arise the r.m.s. residual
!  			is set to zero.

				if (df <= 0.0_IDP) s(l,iplus1) = 0.0_IDP
				if (df > 0.0_IDP) s(l,iplus1) = sqrt(sigmai/df)
			end do

!  		the chebyshev coefficients a(l, i+1, 1), a(l, i+1, 2)....
!  		a(l, i+1, i+1) in the polynomial approximation of degree i
!  		to each set of values of the independent variable
!  		(l=1,2,...,ndv) together with the coefficients
!  		work2(1, 1), work2(1, 2), ..., work2(1, i + 1),   in the
!  		chebyshev-series representation of the  (i + 1) th
!  		orthogonal polynomial are computed.

			do jplus1 = 1, iplus1
				jplus2 = jplus1 + 1
				pij = work2(1,jplus1)
				serr(jplus1) = serr(jplus1) + pij**2/di
				do l = 1, ndv
					a(l,ii,jplus1) = a(l,im1,jplus1) + ci(l)*pij
				end do
				if (jplus1 == kplus1) exit
				work2(1,jplus1) = work2(1,jplus2) + work2(2,jplus1) - alpip1*pij - betai*work2(2,jplus2)
				work2(2,jplus2) = pij
			end do
  			if (iplus1 < kplus1) then
				dim1 = di
				im1 = ii
			end if
		end do
		do iplus1 = 1, kplus1
			serr(iplus1) = 1.0/sqrt(serr(iplus1))
		end do
		ifail = 0

	end subroutine adze02

	subroutine aeue01(m,x,ip,np1,b,w)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 13 revised. use of mark 12 x02 functions (apr 1988).

!  	*******************************************************

!  	npl algorithms library routine q0poly

!  	created 02 05 80.  updated 13 05 80.  release 00/08

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england.

!  	*******************************************************

!  	aeue01.  an algorithm to determine the chebyshev series
!  	representation of a zeroizing polynomial  q0(x),
!  	i.e. a polynomial which takes on zero values (and
!  	possibly zero derivative values) at specified points

!  	input parameters
!  		m        number of distinct x-values.
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		ip       highest order of derivative at each x-value
!  		np1      n + 1,  where  n = number of zeros (including
!  						those of derivatives) to be taken on by
!  						q0(x).  n = m + ip(1) + ip(2) + ... + ip(m).

!  	output parameters
!  		b        chebyshev coefficients of  q0(x)

!  	workspace parameters
!  		w        workspace

!  	.. scalar arguments ..
		integer :: m, np1
!  	.. array arguments ..
		real(IDP), dimension(np1) :: b, w
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
!  	.. local scalars ..
		real(IDP) :: ai, anbig, ct, eps, factor, ovfl, ri, sfac, test, unfl, xi, xtrema
		integer :: i, i2, ifail, ip1, k, l, n, nu

		interface
!  	.. external functions ..
			function x02ame(x)
				use param
				implicit none
				real(IDP) :: x, x02ame
			end function x02ame
!  	.. external subroutines ..
			subroutine e02afe(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine e02afe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, log, sin
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, sxtnth=0.0625_IDP, one=1.0_IDP, two=2.0_IDP, sxteen=16.0_IDP, pi=3.14159265358979323846_IDP
!  	.. executable statements ..

!  	evaluate  q0(x)  at the extrema of the chebyshev polynomial
!  	of degree  n,  employing dynamic scaling to avoid the
!  	possibility of overflow or underflow

		ovfl = sxtnth/x02ame(sxtnth)
		unfl = sxtnth*ovfl
		eps = epsilon(sxtnth)
		n = np1 - 1
		factor = 2*n
		factor = pi/factor
		do i = 1, np1
			w(i) = one
		end do
		do k = 1, m
			ip1 = ip(k) + 1
			xi = x(k)
			do l = 1, ip1
				anbig = zero
				i2 = n + 2
				do i = 1, np1
					i2 = i2 - 2
					ri = i2
					xtrema = sin(factor*ri)
					ai = w(i)*(xtrema-xi)
					w(i) = ai
					if (abs(ai) > anbig) anbig = abs(ai)
				end do
				do
					sfac = one
					if (anbig > ovfl) sfac = sxtnth
					if (anbig < unfl) sfac = sxteen
					if (sfac == one) exit
					anbig = anbig*sfac
					do i = 1, np1
						w(i) = w(i)*sfac
					end do
				end do
			end do
		end do
		ct = ovfl
		do i = 1, np1
			test = abs(w(i))/anbig
			if (test <= eps) w(i) = zero
			if (test > eps .and. test < ct) ct = test
		end do
		ct = ct*anbig
		sfac = one
		do
		 	if (ct < one) exit
			ct = ct*sxtnth
			sfac = sfac*sxtnth
		end do
		do i = 1, np1
			w(i) = w(i)*sfac
		end do

!  	determine the chebyshev representation of  q0(x)

		call e02afe(np1,w,b,ifail)

!  	set the leading coefficient of  q0(x)  to an
!  	exact power of  2

		ai = b(np1)
		ct = log(abs(ai))/log(two)
		nu = ct
		sfac = two**nu
		b(np1) = sfac
		sfac = sfac/ai
		do i = 1, n
			b(i) = b(i)*sfac
		end do

	end subroutine aeue01

	subroutine aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 13 revised. use of mark 12 x02 functions (apr 1988).

!  	*******************************************************

!  	npl algorithms library routine presid

!  	created 18 02 80.  updated 14 05 80.  release 00/28

!  	author ... maurice g. cox.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aeve01.  forms performance indices and residuals
!  	for a polynomial approximation  p(x)  to a set of data
!  	which may contain derivative values.  also indicates
!  	whether the polynomial defined by the specified
!  	coefficients is a better approximation than the best
!  	so far.

!  	input parameters
!  		withq0   true if zeroizing polynomial, else false
!  		withpi   if true, performance indices produced under
!  						all circumstances.  otherwise, produced
!  						only if  p(x)  is an improvement
!  		m        number of x-values.  all distinct
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        x-values.  normalized to  (-1, 1)
!  		n        number of y-values
!  		y        values and derivatives of dependent variable
!  		ip       highest order of derivative at each x-value
!  		imax     one greater than largest value of  ip
!  		a        chebyshev coefficients of  p(x)
!  		la       dimension of  a  and  adif.
!  						 >=  n + 1 if withq0 is true,
!  						 >=  n     otherwise
!  		it       iteration number

!  	input/output parameters
!  		rnmbst   2-norms of residuals corresponding to the
!  						best polynomial so far and its derivatives

!  	output parameters
!  		rnm      2-norms of residuals corresponding to
!  						p(x)  and its derivatives
!  		improv   true if  p(x)  is an improvement, else false
!  		adif     chebyshev coefficients of  (imax - 1)-st
!  						derivative of  p(x)
!  		res      residuals corresponding to y-values
!  	 * pmax     largest performance index
!  	 * pindex   performance indices

!  	note.  the parameters marked  *  are provided only
!  			 if  improv  is true

!  	key local variables
!  		asumax   for current value of  l,  maximum value
!  						over  j = 0, 1, ..., l  of sum of moduli
!  						of chebyshev coefficients of derivative
!  						of order  j  of  p(x)
!  		eps      relative machine precision
!  		iy       location in  y  of current specified
!  						derivative value of order  l
!  		nl       number of specified derivative values
!  						of order  l
!  		nterms   number of terms in chebyshev series
!  						representation of current derivative
!  						of  p(x)
!  		res2nm   2-norm of the residuals corresponding to the
!  						specified derivative values of order  l - 1

!  	.. scalar arguments ..
		real(IDP) :: pmax, xmax, xmin
		integer :: imax, it, la, m, n
		logical :: improv, withpi, withq0
!  	.. array arguments ..
		real(IDP), dimension(la) :: a, adif
		real(IDP), dimension(imax) :: pindex, rnm, rnmbst
		real(IDP), dimension(n) :: res, y
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
!  	.. local scalars ..
		real(IDP) :: absres, asum, asumax, eps, p, plarge, resmax, rscale, t
		integer :: i, ia, iy, l, lm1, nl, nterms
		logical :: imp

		interface
!  	.. external subroutines ..
			subroutine ahze02(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)
				use param
				implicit none
				real(IDP) :: patm1, xmax, xmin
				integer :: ia1, iadif1, la, ladif, np1
				real(IDP), dimension(la) :: a
				real(IDP), dimension(ladif) :: adif
			end subroutine ahze02
			subroutine akze02(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine akze02
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, sqrt
!  	.. data statements ..
		real(IDP) :: half=0.5_IDP, zero=0.0_IDP, one=1.0_IDP, mltplr=8.0_IDP
!  	.. executable statements ..
		pmax = zero
		eps = epsilon(pmax)
		nterms = n
		if (withq0) nterms = nterms + 1
		asumax = zero
		do i = 1, nterms
			adif(i) = a(i)
		end do
		nterms = nterms + 1
		do l = 1, imax
			nterms = nterms - 1

!  		determine sum of moduli  asum  of chebyshev coefficients
!  		of derivative of order  l - 1  of  p(x),  and update
!  		asumax

			asum = half*abs(adif(1))
			if (nterms > 1) then
				do ia = 2, nterms
					asum = asum + abs(adif(ia))
				end do
			end if
			if (asum > asumax) asumax = asum

!  		pindex(l)  is used temporarily to hold the
!  		value of  asumax

			pindex(l) = asumax
			iy = l
			nl = 0

!  		to reduce the possibility of underflow and overflow,
!  		compute  res2nm  as  resmax*sqrt(rscale),  where
!  		resmax  and  rscale  are updated for each residual
!  		corresponding to a specified derivative value of
!  		order  l - 1.  at any stage,  resmax  holds the
!  		modulus of the residual of maximum magnitude so far,
!  		and  rscale  the sum so far of the squares of the
!  		residuals, each scaled by  resmax.

			resmax = zero
			rscale = one
			do i = 1, m

!  			skip if no derivative value of order  l - 1  is
!  			specified at  i-th  x-value

				if (ip(i)+1 >= l) then
					nl = nl + 1

!  				evaluate  p,  the  (l - 1)-st  derivative of  p(x)
!  				at  x = x(i)

					call akze02(nterms,adif,1,la,x(i),p)

!  				save residual corresponding to this value

					if (withq0) res(iy) = -p
					if ( .not. withq0) res(iy) = y(iy) - p

!  				update  resmax  and  rscale

					absres = abs(res(iy))
					if (absres /= zero .and. absres <= resmax) rscale = rscale + (absres/resmax)**2
					if (absres /= zero .and. absres > resmax) rscale = rscale*(resmax/absres)**2 + one
					if (absres /= zero .and. absres > resmax) resmax = absres
				end if
				iy = iy + ip(i) + 1
			end do
			rnm(l) = resmax*sqrt(rscale)

!  		form coefficients in chebyshev series representation
!  		of  l-th  derivative of  p(x)  from those of
!  		(l - 1)-st derivative

			if (l < imax) call ahze02(nterms,xmin,xmax,adif,1,la,t,adif,1,la)
		end do

!  	if not on zero-th iteration,
!  	detect whether there has been an improvement ...

		imp = (it == 0)
		if (it /= 0) then
			do l = 1, imax
				if (rnm(l) < rnmbst(l)) imp = .true.
			end do
		end if
		if (imp .or. withpi) then

!  	... and, if so, or if zero-th iteration, or if
!  	they are required anyway, form the performance
!  	indices corresponding to the improved
!  	approximation

			plarge = zero
			do l = 1, imax

!  			nl  is the number of derivative values of order  l - 1

				nl = 0
				lm1 = l - 1
				do i = 1, m
					if (ip(i) >= lm1) nl = nl + 1
				end do
				t = nl

!  			retrieve the value of asumax

				asumax = pindex(l)
				if (asumax /= zero) pindex(l) = rnm(l)/(mltplr*eps*asumax*sqrt(t))
				if (asumax == zero) pindex(l) = zero
				if (pindex(l) > plarge) plarge = pindex(l)
			end do
			pmax = plarge
		end if
		improv = imp

	end subroutine aeve01

	subroutine aewe01(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine pntrpa

!  	created 20 12 79.  updated 13 05 80.  release 00/47

!  	authors ... gerald t. anthony, maurice g. cox,
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aewe01. a routine, without checks, which determines and
!  	refines a polynomial interpolant  q(x)  to data which
!  	may contain derivatives, and an associated zeroizing
!  	polynomial  q0(x).

!  		m        number of distinct x-values
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        independent variable values (distinct)
!  		y        values and derivatives of
!  						dependent variable.
!  		ip       highest order of derivative at each x-value.
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		np1      value of  n + 1
!  		itmin,
!  		itmax    minimum and maximum number of iterations to be
!  						performed.
!  		iwrk(1)  see workspace (and associated
!  						dimension) parameters

!  	output parameters
!  		a        chebyshev coefficients of  q(x)
!  		b        chebyshev coefficients of  q0(x)

!  	workspace (and associated dimension) parameters
!  		wrk      real workspace array.  the first imax elements
!  						contain, on exit, performance indices for
!  						the interpolating polynomial, and the next
!  						n  elements the computed residuals
!  		lwrk     dimension of wrk. lwrk must be at least
!  						8*n + 5*imax + m + 5, where
!  						imax is one more than the largest element
!  						of the array ip.
!  		iwrk     integer workspace array.
!  						the first element of this array is used
!  						as an input parameter (which is destroyed
!  						on exit).  the zeroizing polynomial  q0(x)
!  						is constructed or not according to whether
!  						iwrk(1)  is non-zero or zero.
!  		liwrk    dimension of iwrk.  at least 2*m + 2.

!  	failure indicator parameter
!  		ifail    failure indicator
!  						0 - successful termination
!  						1 - iteration limit in deriving  q(x)
!  						2 - divergent iteration in deriving  q(x)
!  						3 - iteration limit in deriving  q0(x)
!  						4 - divergent iteration in deriving  q0(x)

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: ifail, itmax, itmin, liwrk, lwrk, m, n, np1
!  	.. array arguments ..
		real(IDP), dimension(n) :: a, y
		real(IDP), dimension(np1) :: b
		real(IDP), dimension(lwrk) :: wrk
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip
		integer, dimension(liwrk) :: iwrk
!  	.. local scalars ..
		real(IDP) :: pmax
		integer :: i, iadif, iatrl, ibdif, ibtrl, ic, id, ida, idb, ierror, iftau, ilocx, ilocy, imax, initq, initq0, ipiq, ipiq0, &
					  iptrl, ires, irnm, irtrnm, iw, ix
		logical :: withq0

		interface
!  	.. external subroutines ..
			subroutine aeye01(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w, &
									adif,da,rnm,rtrlnm,locx,locy,ifail)
				use param
				implicit none
				real(IDP) :: pmax, xmax, xmin
				integer :: ifail, imax, itmax, itmin, la, m, n, nit, np1
				logical :: withq0
				real(IDP), dimension(la) :: a, adif, atrial, da
				real(IDP), dimension(n) :: c, d, ftau, res, y
				real(IDP), dimension(imax) :: pindex, ptrial, rnm, rtrlnm
				real(IDP), dimension(np1) :: w
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx, locy
			end subroutine aeye01
			subroutine akye02(xmin,xmax,x,xcap)
				use param
				implicit none
				real(IDP) :: x, xcap, xmax, xmin
			end subroutine akye02
		end interface

!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP
!  	.. executable statements ..
		ierror = 0
		imax = 0
		do i = 1, m
			if (ip(i) > imax) imax = ip(i)
		end do
		imax = imax + 1

!  	indicate whether  q0(x)  is required

		withq0 = (iwrk(1) /= 0)

!  	treat the case n = 1 separately.

		if (n == 1) then
			a(1) = two*y(1)
			if (withq0) then
				call akye02(xmin,xmax,x(1),b(1))
				b(1) = -two*b(1)
				b(2) = one
			end if

!  	set to zero the numbers of iterations taken, the (sole)
!  	residual and the values of the performance indices
!  	in this special case, and finish

			iwrk(1) = 0
			if (withq0) iwrk(2) = 0
			wrk(1) = zero
			wrk(2) = zero
			if (withq0) wrk(3) = zero
		end if

!  	transform the  x*s  to  (-1, 1)

		ix = 2*imax + n
		do i = 1, m
			ix = ix + 1
			call akye02(xmin,xmax,x(i),wrk(ix))
		end do
		if (n == 1) then
			ifail = ierror
			return
		end if
		if (withq0) then

!  		workspace allocation for call to  aeye01  ...

			ires = imax + 1
			ipiq0 = ires + n
			ix = ipiq0 + imax
			ibtrl = ix + m
			iptrl = ibtrl + np1
			iftau = iptrl + imax
			ic = iftau + n
			id = ic + n
			iw = id + np1
			ibdif = iw + np1
			idb = ibdif + np1
			irnm = idb + np1
			irtrnm = irnm + imax
			initq0 = 2
			ilocx = initq0 + 1
			ilocy = ilocx + m

!  		... to determine  q0(x)

			call aeye01(.true.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,b,np1,wrk(ires),pmax,wrk(ipiq0),iwrk(initq0), &
							wrk(ibtrl),wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(ibdif),wrk(idb),wrk(irnm),wrk(irtrnm), &
							iwrk(ilocx),iwrk(ilocy),ierror)

!  		re-code failure indicator

			if (ierror /= 0) ierror = ierror + 2
			if (ierror /= 0) then
				ifail = ierror
				return
			end if
		end if

!  	workspace allocation for call to  aeye01  ...

		ipiq = 1
		ires = ipiq + imax
		ix = ires + n + imax
		iatrl = ix + m
		iptrl = iatrl + n
		iftau = iptrl + imax
		ic = iftau + n
		id = ic + n
		iw = id + n
		iadif = iw
		ida = iadif + np1
		irnm = ida + np1
		irtrnm = irnm + imax
		initq = 1
		ilocx = initq + 2
		ilocy = ilocx + m

!  	... to determine  q(x)

		call aeye01(.false.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,n,wrk(ires),pmax,wrk(ipiq),iwrk(initq),wrk(iatrl), &
						wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(iadif),wrk(ida),wrk(irnm),wrk(irtrnm),iwrk(ilocx),iwrk(ilocy), &
						ierror)

		ifail = ierror

	end subroutine aewe01

	subroutine aexe01(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine divdif

!  	created 17 07 79.  updated 14 05 80.  release 00/08

!  	authors ... maurice g. cox and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aexe01.  an algorithm to determine the next coefficient
!  	in the newton form of an interpolating polynomial

!  	input parameters
!  		m        number of distinct x-values.
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		ip       highest order of derivative at each x-value
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		locx     pointers to x-values in constructing
!  						newton form of polynomial
!  		!  		newton coefficients determined so far
!  		n!  	  number of newton coefficients determined so far
!  		xnew     element of  x  associated with new
!  						newton coefficient
!  		ixnext   number of x-values so far incorporated
!  						(including  xnew)
!  		ynew     scaled derivative value corresponding to
!  						xnew  and  nordp1
!  		nordp1   one plus order of derivative
!  						associated with  ynew

!  	input/output parameters
!  		d        elements in previous, and then new, upward
!  						sloping diagonal of divided difference table

!  	output parameters
!  		cnew     new newton coefficient generated

!  	.. scalar arguments ..
		real(IDP) :: cnew, xnew, ynew
		integer :: ixnext, m, n, nc, nordp1
!  	.. array arguments ..
		real(IDP), dimension(n) :: c, d
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx
!  	.. local scalars ..
		real(IDP) :: dif
		integer :: ic, is, ix, k, locxi
!  	.. executable statements ..
		ic = nc - nordp1 + 1
		d(1) = ynew
		if (ixnext /= 1) then
			is = 0
			ix = 0
			do k = 1, ic
				if (k > is) then
					ix = ix + 1
					locxi = locx(ix)
					is = is + ip(locxi) + 1
					dif = x(locxi) - xnew
				end if
				if (nordp1 == 1) d(k+1) = (c(k)-d(k))/dif
				if (nordp1 > 1) d(k+1) = (d(k+1)-d(k))/dif
			end do
		end if
		cnew = d(ic+1)

	end subroutine aexe01

	subroutine aeye01(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w,adif,da, &
							rnm,rtrlnm,locx,locy,ifail)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine refh

!  	created 17 07 79.  updated 14 05 80.  release 00/47

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england

!  	*******************************************************

!  	aeye01.  a routine to approximate and then refine an
!  	interpolating polynomial  q(x)  or a zeroizing
!  	polynomial  q0(x)  in its chebyshev representation

!  	input parameters
!  		withq0   true if zeroizing polynomial, else false
!  		m        the number of distinct data points.
!  		x        array containing the distinct x-values,
!  						normalized if necessary to (-1, 1).
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  	 * y        array containing values and derivatives of
!  						the dependent variable.
!  		ip       array specifying the highest order of
!  						derivative at each x-value.
!  		imax     one more than the largest element of the
!  						array ip.
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).
!  		np1      value of  n + 1
!  		itmin,
!  		itmax    the lower and upper limits on the iterative
!  						process.

!  	output (and associated dimension) parameters
!  		a        chebyshev coefficients of polynomial
!  		la       dimension of  a.
!  						 >=  n  if interpolating polynomial
!  						 >=  n + 1  if zeroizing polynomial
!  		res      residuals of polynomial
!  		pmax     largest performance index
!  		pindex   performance indices
!  		nit      number of iterations taken

!  	workspace parameters
!  		atrial   trial values of the chebyshev coefficients
!  		ptrial   performance indices corresponding to  atrial
!  	 * ftau     scaled values of  y.  if  y(i)  is the
!  						value of an  r-th  derivative, then
!  						((xmax - xmin)/2)**r/(factorial r)
!  						times  y(i)  is the value of  ftau(i)
!  	 * !  		coefficients in newton form of polynomial
!  	 * d        intermediate divided difference values
!  	 **w        values of correction polynomial at
!  						chebyshev extrema
!  		adif     chebyshev coefficients of a derivative of
!  						an approximation to the polynomial
!  		da       chebyshev coefficients of a
!  						correction polynomial
!  		rnm      residual norms corresponding to  a
!  		rtrlnm   residual norms corresponding to  atrial
!  	 * locx     pointers to x-values in constructing
!  						newton form of polynomial
!  	 * locy     pointers to y-values corresponding to x-values

!  	failure indicator parameter
!  		ifail    failure indicator
!  						0 - successful termination
!  						1 - iteration limit exceeded
!  						2 - iteration divergent

!  			notes.  (1) the elements of the arrays marked  *  are
!  							not accessed if  withq0  is  true.
!  					  (2) the elements of the array marked  **  is
!  							not accessed if  withq0  is  false.

!  	.. scalar arguments ..
		real(IDP) :: pmax, xmax, xmin
		integer :: ifail, imax, itmax, itmin, la, m, n, nit, np1
		logical :: withq0
!  	.. array arguments ..
		real(IDP), dimension(la) :: a, adif, atrial, da
		real(IDP), dimension(n) :: c, d, ftau, res, y
		real(IDP), dimension(imax) :: pindex, ptrial, rnm, rtrlnm
		real(IDP), dimension(np1) :: w
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx, locy
!  	.. local scalars ..
		real(IDP) :: amax, atrnrm, danrm, pmxtrl, scale
		integer :: i, ierror, it, itemp, itmxp1, itp1, l, nfref, npilt1, nterms
		logical :: improv, withpi, zeroda

		interface
!  	.. external subroutines ..
			subroutine aeue01(m,x,ip,np1,b,w)
				use param
				implicit none
				integer :: m, np1
				real(IDP), dimension(np1) :: b, w
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
			end subroutine aeue01
			subroutine aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)
				use param
				implicit none
				real(IDP) :: pmax, xmax, xmin
				integer :: imax, it, la, m, n
				logical :: improv, withpi, withq0
				real(IDP), dimension(la) :: a, adif
				real(IDP), dimension(imax) :: pindex, rnm, rnmbst
				real(IDP), dimension(n) :: res, y
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
			end subroutine aeve01
			subroutine aeze01(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: m, n
				real(IDP), dimension(n) :: a, c, d, ftau, y
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx, locy
			end subroutine aeze01
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, log
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, one=1.0_IDP, sxteen=16.0_IDP
!  	.. executable statements ..
		ierror = 0

!  	number of terms in polynomial

		nterms = n
		if (withq0) nterms = nterms + 1

!  	number of performance indices less than one

		npilt1 = 0

!  	indicate that performance indices are to be produced
!  	only if an improvement is obtained

		withpi = .false.

!  	indicate that the fine refinement stage has not yet started

		nfref = -2

!  	set residuals initially equal to specified y-values
!  	(if  q(x)  required) or zero (if  q0(x)  required)

		do i = 1, n
			if ( .not. withq0) res(i) = y(i)
			if (withq0) res(i) = zero
		end do

!  	initialize trial chebyshev coefficients

		do i = 1, nterms
			atrial(i) = zero
		end do

!  	commence iterative refinement

		itmxp1 = itmax + 1
		do itp1 = 1, itmxp1

!  		it  is the actual iteration number,  it = 0
!  		corresponding to the first estimate of the polynomial

			it = itp1 - 1

!  		determine chebyshev coefficients  da  of polynomial
!  		approximately satisfying the conditions in  res

			if (withq0 .and. it == 0) call aeue01(m,x,ip,np1,da,w)
			if ( .not. withq0 .or. (withq0 .and. it > 0)) call aeze01(m,xmin,xmax,x,res,ip,n,da,locx,locy,ftau,d,c)

!  		skip test for divergence if on zero-th iteration

			if (it > 0) then

!  		determine the norms of  da  and (the previous)  atrial

				danrm = half*abs(da(1))
				atrnrm = half*abs(atrial(1))
				if (n > 1) then
					do i = 2, n
						danrm = danrm + abs(da(i))
						atrnrm = atrnrm + abs(atrial(i))
					end do
				end if
				if (withq0) atrnrm = atrnrm + abs(atrial(np1))

!  		assume divergence if the norm of  da  is not
!  		less than that of  atrial  ...

				if (danrm >= atrnrm) ierror = 2
				if (danrm >= atrnrm) exit
			end if

!  		... otherwise determine new trial approximation

			zeroda = .true.
			do i = 1, n
				atrial(i) = atrial(i) + da(i)
				if (da(i) /= zero) zeroda = .false.
			end do
			if (withq0 .and. it == 0) then
				atrial(np1) = da(np1)
				if (da(np1) /= zero) zeroda = .false.
			end if

!  		determine residuals, performance indices and
!  		largest performance index corresponding to
!  		trial coefficients  atrial

			call aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,atrial,la,it,rnm,rtrlnm,improv,adif,res,pmxtrl,ptrial)

!  		set dummy, non-zero, value of  pmxtrl  if no
!  		improvement, otherwise it is undefined

			if ( .not. improv) pmxtrl = sxteen

!  		if on first iteration, or if the largest performance
!  		index is zero, or if all components of  da  are zero,
!  		take the trial set of coefficients and performance
!  		indices as the best (so far)

			if (it == 0 .or. pmxtrl == zero .or. zeroda) then
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
				pmax = pmxtrl
			end if

!  		finish if largest performance index is zero
!  		or if all components of  da  are zero
!  		(i.e. no further improvement is possible)

			if (pmxtrl == zero .or. zeroda) exit

!  		indicate whether the fine refinement stage has commenced
!  		(i.e. for the first time all performance indices are
!  		less than one)

			if (nfref == -2 .and. pmxtrl < one) nfref = -1

!  		branch according to whether the process is in the
!  		fine refinement stage  (nfref  >=  0)  or not
!  		(nfref  ==  -1)

			if (nfref < -1) then

!  		the process is in the course refinement phase.
!  		update the coefficients and the corresponding
!  		norms and performance indices if
!  			(i)  there has been an improvement in (at
!  				  least) one of the residual norms, and
!  			(ii) the number of performance indices
!  				  that are less than one has not
!  				  increased compared with those of
!  				  the best polynomial so far.

				if ( .not. improv) cycle
				itemp = 0
				do l = 1, imax
					if (ptrial(l) < one) itemp = itemp + 1
				end do
 				if (itemp < npilt1) cycle
				npilt1 = itemp
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
				pmax = pmxtrl
				cycle
			end if

!  		the process is in the fine refinement phase.
!  		update the coefficients and the corresponding
!  		norms and performance indices if
!  			(i)  there has been an improvement in (at
!  				  least) one of the residual norms, and
!  			(ii) the largest performance index is less
!  				  than the largest of that of the best
!  				  polynomial so far.
!  		increment the number of fine refinements (the number
!  		of refinements since the first occasion when all
!  		performance indices were less than unity), exiting
!  		if as many as  itmin  fine refinements have been
!  		performed

			if (improv .and. pmxtrl < pmax) then
				do i = 1, nterms
					a(i) = atrial(i)
				end do
				do l = 1, imax
					rnm(l) = rtrlnm(l)
					pindex(l) = ptrial(l)
				end do
 				pmax = pmxtrl
			end if
			nfref = nfref + 1
			if (nfref >= itmin) exit
		end do

!  	the process has not succeeded in reducing all the
!  	performance indices to less than unity

		if (itp1 > itmxp1) then
			ierror = 1
			it = itmax
		end if

!  	number of iterations actually taken

		nit = it
		if (withq0) then

!  	in the case of  q0(x),  scale its coefficients by
!  	an integral power of  16  such that the largest
!  	coefficient is of order unity

			amax = zero
			do i = 1, np1
				if (abs(a(i)) > amax) amax = abs(a(i))
			end do
			if (amax > zero) then
				i = log(amax)/log(sxteen)
				scale = sxteen**(-i)
				do i = 1, np1
					a(i) = scale*a(i)
				end do
			end if
		end if

!  	return residuals corresponding to selected coefficients

		withpi = .true.
		call aeve01(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnm,rtrlnm,improv,adif,res,pmax,pindex)
		ifail = ierror

	end subroutine aeye01

	subroutine aeze01(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	*******************************************************

!  	npl algorithms library routine qpoly

!  	created 02 05 80.  updated 14 05 80.  release 00/09

!  	authors ... gerald t. anthony, maurice g. cox
!  					j. geoffrey hayes and michael a. singer.
!  	national physical laboratory, teddington,
!  	middlesex tw11 olw, england.

!  	*******************************************************

!  	aeze01. an algorithm to determine the chebyshev series
!  	representation of a polynomial interpolant  q(x)  to
!  	arbitrary data points where derivative information may
!  	be given.

!  	input parameters
!  		m        number of distinct x-values.
!  		xmin,
!  		xmax     lower and upper endpoints of interval
!  		x        independent variable values,
!  						normalized to  (-1, 1)
!  		y        values and derivatives of dependent variable
!  		ip       highest order of derivative at each x-value
!  		n        number of interpolating conditions.
!  						n = m + ip(1) + ip(2) + ... + ip(m).

!  	output parameters
!  		a        chebyshev coefficients of  q(x)

!  	workspace parameters
!  		locx     pointers to x-values in constructing
!  						newton form of polynomial
!  		locy     pointers to y-values corresponding to x-values
!  		ftau     scaled values of  y.  if  y(i)  is the
!  						value of an  r-th  derivative, then
!  						((xmax - xmin)/2)**r/(factorial r)
!  						times  y(i)  is the value of  ftau(i)
!  		d        intermediate divided difference values
!  		!  		newton coefficients of  q(x)

!  	.. scalar arguments ..
		real(IDP) :: xmax, xmin
		integer :: m, n
!  	.. array arguments ..
		real(IDP), dimension(n) :: a, c, d, ftau, y
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx, locy
!  	.. local scalars ..
		real(IDP) :: cmin, cnew, factor, ri, rj, s, scale, sfac, v, xch
		integer :: i, i2, ic, icmin, ifail, iftau, ip1, isave, iy, j, jmax, k, krev, l, lmax, locxi, locxj, locxk, locyi, nc

		interface
!  	.. external subroutines ..
			subroutine aexe01(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)
				use param
				implicit none
				real(IDP) :: cnew, xnew, ynew
				integer :: ixnext, m, n, nc, nordp1
				real(IDP), dimension(n) :: c, d
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx
			end subroutine aexe01
			subroutine e02afe(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine e02afe
		end interface

!  	.. intrinsic functions ..
!  	intrinsic abs, sin
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP, pi=3.14159265358979323846_IDP
!  	.. executable statements ..
		scale = (xmax-xmin)/two

!  	initialize x- and y-pointers

		iy = 0
		do i = 1, m
			locx(i) = i
			iy = iy + 1
			locy(i) = iy
			ftau(iy) = y(iy)
			jmax = ip(i)
			if (jmax == 0) cycle

!  		form the scaled derivatives, i.e. an  r-th  derivative
!  		value is divided by  factorial r  and multiplied
!  		by the  r-th  power of  (xmax - xmin)/2

			sfac = one
			do j = 1, jmax
				iy = iy + 1
				rj = j
				sfac = sfac*scale/rj
				ftau(iy) = y(iy)*sfac
			end do
		end do

!  	form successive upward sloping diagonals of
!  	the divided difference table

		nc = 0
	outer:do j = 1, m

!  			choose each x-value in turn to make the corresponding
!  			newton coefficient as small in magnitude as possible

				do i = j, m
					locxi = locx(i)
					locyi = locy(locxi)
					call aexe01(m,x,ip,n,locx,c,nc,x(locxi),j,ftau(locyi),1,cnew,d)
					if (i > j .and. abs(cnew) >= abs(cmin)) cycle
					cmin = cnew
					icmin = i
				end do
				c(nc+1) = cmin
				isave = locx(j)
				locxj = locx(icmin)
				locx(icmin) = isave
				locx(j) = locxj

!  			calculate the resulting newton coefficient (i.e.
!  			repeat the above computation, but only in the case
!  			leading to the smallest new coefficient)

				iftau = locy(locxj) - 1
				ip1 = ip(locxj) + 1
	inner:	do i = 1, ip1
					iftau = iftau + 1
					call aexe01(m,x,ip,n,locx,c,nc,x(locxj),j,ftau(iftau),i,c(nc+1),d)
					nc = nc + 1
					if (nc == n) exit outer
				end do inner
			end do outer

!  	evaluate  q(x)  (from its newton form) at the extrema
!  	of the chebyshev polynomial of degree  n - 1  ...

		factor = 2*n - 2
		factor = pi/factor
		i2 = n + 1
		do i = 1, n
			i2 = i2 - 2
			ri = i2
			xch = sin(factor*ri)
			s = c(n)
			ic = n
			k = m + 1
			do krev = 1, m
				k = k - 1
				locxk = locx(k)
				lmax = ip(locxk) + 1
				if (k == m) lmax = lmax - 1
				if (lmax <= 0) cycle
				v = xch - x(locxk)
				do l = 1, lmax
					ic = ic - 1
					s = s*v + c(ic)
				end do
			end do
			d(i) = s
		end do

!  	... in order to determine the coefficients in its
!  	chebyshev representation

		call e02afe(n,d,a,ifail)

	end subroutine aeze01

	subroutine ahze02(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine auxcdf

!  	created 1/5/79    updated 23/1/80     release no. 00/03

!  	authors.. gerald t anthony, maurice g cox, j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1    = n+1 where n is degree of given polynomial
!  		xmin   lower limit of range of x
!  		xmax   upper limit of range of x
!  		a      coefficients a0, a1,...an of the given polynomial
!  		ia1       are stored in array a in positions 1, 1+ia1,...
!  					 1+n*ia1, respectively
!  		la     the declared dimension of array a

!  	output parameters
!  		patm1  the value of the given polynomial at xmin
!  		adif   the coefficients of the derivative polynomial
!  		iadif1    are returned in array adif in positions
!  					 1, 1+iadif1,...1+(n-1)*iadif1
!  		ladif  the declared dimension of array adif

!  	differentiate the series with coefficients a of degree n
!  	(i.e. np1 coefficients) to obtain the series with coefficients
!  	adif of degree n-1. also set next higher coefficient to zero.

!  	.. scalar arguments ..
		real(IDP) :: patm1, xmax, xmin
		integer :: ia1, iadif1, la, ladif, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
		real(IDP), dimension(ladif) :: adif
!  	.. local scalars ..
		real(IDP) :: ptemp, r, sclftr, u, v, w
		integer :: i, n, na, nadif
!  	.. data statements ..
		real(IDP) :: two=2.0_IDP
!  	.. executable statements ..
		u = 0.0_IDP
		v = u
		sclftr = two/(xmax-xmin)
		n = np1 - 1
		nadif = n*iadif1 + 1
		ptemp = u
		if (n > 0) then
			na = n*ia1 + 1
			do i = 1, n
				r = np1 - i
				w = u + two*r*a(na)
				ptemp = a(na) - ptemp

!  			store coeff formed previous time round. first time round
!  			store zero as coeff of degree n.

				adif(nadif) = sclftr*v
				u = v
				v = w
				na = na - ia1
				nadif = nadif - iadif1
			end do
		end if
		adif(nadif) = sclftr*v
		patm1 = a(1)/two - ptemp

	end subroutine ahze02

	subroutine akye02(xmin,xmax,x,xcap)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine nrmlze

!  	created 5/5/78    updated 11/12/78    release no. 00/02.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		xmin     minimum value of x
!  		xmax     maximum value of x
!  		x        unnormalized argument in range (xmin, xmax)

!  	output parameter
!  		xcap     normalized value of x

!  	a value of x is given, such that
!  	xmin  <=  x  <=  xmax.
!  	xcap is calculated so that -1  <=  x  <=  +1.

!  	this form for xcap ensures that the computed value has a
!  	very small absolute error.

!  	.. scalar arguments ..
		real(IDP) :: x, xcap, xmax, xmin
!  	.. executable statements ..
		xcap = ((x-xmin)-(xmax-x))/(xmax-xmin)

	end subroutine akye02

	subroutine akze02(np1,a,ia1,la,xcap,result)

		use param
		implicit none

!  	mark 8 release. nag copyright 1979.
!  	mark 11.5(f77) revised. (sept 1985.)

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	npl data fitting library routine tval1

!  	created 9/5/78    updated 6/4/79    release no. 00/07.

!  	authors.. gerald t anthony, maurice g cox, betty curtis
!  	and j geoffrey hayes.
!  	national physical laboratory
!  	teddington, middlesex, england.

!  	* * * * * * * * * * * * * * * * * * * * * * * * * * *

!  	input parameters
!  		np1      np1 = n + 1. n is the degree of the
!  					chebyshev series
!  		a        the array where the coefficients are stored
!  		ia1      the address increment of a
!  		la       dimension of a
!  		xcap     normalized argument of the polynomial

!  	output parameter
!  		result   value of the summation

!  	np1 chebyshev coefficients a0, a1, ..., an, are
!  	stored in the array a in positions 1, 1+ia1, 1+2*ia1, ...,
!  	1+n*ia1, where n = np1 - 1.
!  	ia1 must not be negative.
!  	la must be at least equal to 1 + n*ia1.
!  	the argument xcap is assumed to lie in the range
!  	-1  <=  xcap  <=  +1.
!  	the value of the polynomial of degree n
!  	a0t0(xcap)/2 + a1t1(xcap) + a2t2(xcap) + + ... + antn(xcap),
!  	is calculated for the argument xcap storing it in result.
!  	the recurrence relation by clenshaw, modified by reinsch
!  	and gentleman, is used.

!  	.. scalar arguments ..
		real(IDP) :: result, xcap
		integer :: ia1, la, np1
!  	.. array arguments ..
		real(IDP), dimension(la) :: a
!  	.. local scalars ..
		real(IDP) :: aj, bj, cj, factor, sum
		integer :: j, jrev, n
!  	.. data statements ..
		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, two=2.0_IDP
!  	.. executable statements ..
		if (np1 == 1) then
			sum = half*a(1)
		else
			n = np1 - 1
			aj = zero
			bj = zero
			j = 1 + np1*ia1
			if (xcap > half) then

!  		reinschs modified recurrence.

				factor = two - (xcap+xcap)

!  		bracketing necessary in order to avoid errors

				do jrev = 1, n
					j = j - ia1
					aj = a(j) + aj - bj*factor
					bj = aj + bj
				end do
				sum = half*a(1) + aj - half*factor*bj
			else if (xcap >= -half) then

!  		clenshaws original recurrence.

				factor = xcap + xcap
				do jrev = 1, n
					j = j - ia1
					cj = bj
					bj = aj
					aj = a(j) - cj + bj*factor
				end do
				sum = half*a(1) - bj + half*factor*aj
			else

!  		gentlemans modified recurrence.

				factor = two + (xcap+xcap)

!  		bracketing necessary so as to avoid errors

				do jrev = 1, n
					j = j - ia1
					aj = a(j) - aj + bj*factor
					bj = aj - bj
				end do
				sum = half*a(1) - aj + half*factor*bj
			end if
		end if
		result = sum

	end subroutine akze02

	integer function p01abe(ifail,ierror,srname,nrec,rec)

		use param
		implicit none

!  	mark 11.5(f77) release. nag copyright 1986.
!  	mark 13 revised. ier-621 (apr 1988).
!  	mark 13b revised. ier-668 (aug 1988).

!  	p01abe is the error-handling routine for the nag library.

!  	p01abe either returns the value of ierror through the routine
!  	name (soft failure), or terminates execution of the program
!  	(hard failure). diagnostic messages may be output.

!  	if ierror = 0 (successful exit from the calling routine),
!  	the value 0 is returned through the routine name, and no
!  	message is output

!  	if ierror is non-zero (abnormal exit from the calling routine),
!  	the action taken depends on the value of ifail.

!  	ifail =  1: soft failure, silent exit (i.e. no messages are
!  					output)
!  	ifail = -1: soft failure, noisy exit (i.e. messages are output)
!  	ifail =-13: soft failure, noisy exit but standard messages from
!  					p01abe are suppressed
!  	ifail =  0: hard failure, noisy exit

!  	for compatibility with certain routines included before mark 12
!  	p01abe also allows an alternative specification of ifail in which
!  	it is regarded as a decimal integer with least significant digits
!  	cba. then

!  	a = 0: hard failure  a = 1: soft failure
!  	b = 0: silent exit   b = 1: noisy exit

!  	except that hard failure now always implies a noisy exit.

!  	s.hammarling, m.p.hooper and j.j.du croz, nag central office.

!  	.. scalar arguments ..
		integer :: ierror, ifail, nrec
		character(*) :: srname
!  	.. array arguments ..
		character(*), dimension(:) :: rec
!  	.. local scalars ..
		integer :: i, nerr
		character(len=72) :: mess
!  	.. intrinsic functions ..
!  	intrinsic  				 abs, mod

		interface
!  	.. external subroutines ..
			subroutine abzp01
			end subroutine abzp01
			subroutine x04aae(i,nerr)
				implicit none
				integer :: i, nerr
			end subroutine x04aae
			subroutine x04bae(nout,rec)
				implicit none
				integer :: nout
				character(*) :: rec
			end subroutine x04bae
		end interface

!  	.. executable statements ..
		if (ierror /= 0) then
!  	   abnormal exit from calling routine
			if (ifail == -1 .or. ifail == 0 .or. ifail == -13 .or. (ifail > 0 .and. mod(ifail/10,10) /= 0)) then
!  			noisy exit
				call x04aae(0,nerr)
				do i = 1, nrec
					call x04bae(nerr,rec(i))
				end do
				if (ifail /= -13) then
					write (mess,'(" ** abnormal exit from nag library routine ",a,": ifail"," =",i6)') srname, ierror
					call x04bae(nerr,mess)
					if (abs(mod(ifail,10)) /= 1) then
!  					hard failure
						call x04bae(nerr," ** nag hard failure - execution terminated")
						call abzp01
					else
!  					soft failure
						call x04bae(nerr," ** nag soft failure - control returned")
					end if
				end if
			end if
		end if
		p01abe = ierror

	end function p01abe

	subroutine abzp01

!  	mark 11.5(f77) release. nag copyright 1986.

!  	terminates execution when a hard failure occurs.

!  	******************** implementation note ********************
!  	the following stop statement may be replaced by a call to an
!  	implementation-dependent routine to display a message and/or
!  	to abort the program.
!  	*************************************************************
!  	.. executable statements ..
		stop

	end subroutine abzp01

	subroutine x04aae(i,nerr)

		implicit none

!  	mark 7 release. nag copyright 1978
!  	mark 7c revised ier-190 (may 1979)
!  	mark 11.5(f77) revised. (sept 1985.)
!  	mark 14 revised. ier-829 (dec 1989).
!  	if i = 0, sets nerr to current error message unit number
!  	(stored in nerr1).
!  	if i = 1, changes current error message unit number to
!  	value specified by nerr.

!  	.. scalar arguments ..
		integer :: i, nerr
!  	.. local scalars ..
		integer, save :: nerr1=0
!  	.. executable statements ..
		if (i == 0) nerr = nerr1
		if (i == 1) nerr1 = nerr

	end subroutine x04aae

	subroutine x04bae(nout,rec)

		implicit none

!  	mark 11.5(f77) release. nag copyright 1986.

!  	x04bae writes the contents of rec to the unit defined by nout.

!  	trailing blanks are not output, except that if rec is entirely
!  	blank, a single blank character is output.
!  	if nout < 0, i.e. if nout is not a valid fortran unit identifier,
!  	then no output occurs.

!  	.. scalar arguments ..
		integer :: nout
		character(*) :: rec
!  	.. local scalars ..
		integer :: i
!  	.. intrinsic functions ..
!  	intrinsic  	    len
!  	.. executable statements ..
		if (nout >= 0) then
!  	   remove trailing blanks
			do i = len(rec), 2, -1
				if (rec(i:i) /= " ") exit
			end do
!  	   write record to external file
			write (nout,'(a)') rec(1:i)
		end if

	end subroutine x04bae

	function x02ame(x)

		use param
		implicit none

!  	mark 12 release. nag copyright 1986.

!  	returns the 'safe range' parameter
!  	i.e. the smallest positive model number z such that
!  	for any x which satisfies x >= z and x <= 1/z
!  	the following can be computed without overflow, underflow or other
!  	error

!  	   -x
!  	   1.0/x
!  	   sqrt(x)
!  	   log(x)
!  	   exp(log(x))
!  	   y**(log(x)/log(y)) for any y

!  	x is a dummy argument

		real(IDP) :: x, x02ame

!  	.. intrinsic functions ..
!  	intrinsic  			 sqrt
!  	.. executable statements ..
		x02ame = sqrt(2.0_IDP)*tiny(x)

	end function x02ame