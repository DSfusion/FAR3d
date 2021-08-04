	subroutine fitter(yeq,yeqr,yfar,ir,m)

		use param
		use domain
		implicit none

		integer :: ir,m,ior
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
			subroutine cwlspa(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)
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
			end subroutine cwlspa
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
		end interface


!	vmec equilibrium arrays are fitted using constrained weighted least squares polynomial approximations for the first ir grid points and
!	cubic splines for the rest, in order to assure the correct behavior at the origin.

		ior=(3*ir)/10
		irp1=ir+1
		flux2(1)=0.0_IDP
		do i=1,ir
			flux2(i+1)=rfar(i)*rfar(i)
		end do

		do i=1,ir
			temp(i+1)=yfar(i)/rfar(i)**m
		end do
		temp(1)=5.*temp(5)-4.*temp(6)
		ave = 0.0_IDP
	 	do i=1,irp1
			ave = ave + abs(temp(i))
		end do
		ave = ave/irp1
		do i=1,irp1
			w(i) = 10.0_IDP
			if (abs(temp(i)) > 10.*ave) w(i) = ave/abs(temp(i))
		end do
		do i=1,ior
			w(i)=1.e-10_IDP
		end do
		pmin   = flux2(1)
		pmax   = flux2(irp1)
		xf(1)  = flux2(irp1)
		yf(1)  = temp(irp1)
		aux    = yfar(irp1)/rfar(irp1)**m
		yf(2)  = (aux-temp(ir))/(2.*flux2(2))
		yf(3)  = (aux-2.*temp(irp1)+temp(ir))/(flux2(2)*flux2(2))
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
		call cwlspa(ir,kplus1,nrows,pmin,pmax,flux2,temp,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

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

	end subroutine fitter

	subroutine cwlspa(m,kplus1,nrows,xmin,xmax,x,y,w,mf,xf,yf,lyf,ip,a,s,np1,wrk,lwrk,iwrk,liwrk,ifail)

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

		real(IDP) :: amuj, xi, xmu, xcap
		integer :: i, ierror, im1, imax, ipi, iymux, j, lw, mdist, n, nanu, neps, nser, nwrk, nwrk1, nwrk2

		interface
			subroutine pidazp(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)
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
			end subroutine pidazp
			subroutine wlspad(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)
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
			end subroutine wlspad
			subroutine sumchs(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine sumchs
		end interface

		real(IDP) :: zero=0.0_IDP

		if (mf < 1) then
			write(*,'("*** cwlspa: no x-values at which a constraint is specified ***")')
			stop
		end if
		imax = 0
		np1 = 1
		do i = 1, mf
			ipi = ip(i) + 1
			if (ipi < 1) then
				write(*,'("*** cwlspa: ip(",i3,") < 0 ***")') i
				stop
			end if
			if (ipi > imax) imax = ipi
			np1 = np1 + ipi
		end do
		n = np1 - 1
		if (lyf < n .or. liwrk < 2*mf+2) then
			write(*,'("*** cwlspa: insufficient work storage ***")')
			stop
		end if
		i = 4*m + 3*kplus1
		lw = 8*np1 + 5*imax + mf - 3
		if (lw < i) lw = i
		lw = lw + 2*np1
		if (lw > lwrk .or. np1 > kplus1 .or. m < 1) then
			write(*,'("*** cwlspa: insufficient work storage ***")')
			stop
		end if
		nanu = np1 + 1
		nwrk = nanu + np1
		if (xmax <= xmin) then
			write(*,'("*** cwlspa: xmax <= xmin ***")')
			stop
		end if
		if (xf(1) > xmax .or. xf(1) < xmin) then
			write(*,'("*** cwlspa: xf(1) > xmax or xf(1) < xmin ***")')
			stop
		end if
		if (mf > 1) then
			do i = 2, mf
				xi = xf(i)
				if (xi > xmax .or. xi < xmin) then
					write(*,'("*** cwlspa: xf(",i3,") > xmax or xf(",i3,") < xmin ***")') i,i
					stop
				end if
				im1 = i - 1
				do j = 1, im1
					if (xi == xf(j)) then
						write(*,'("*** cwlspa: xf(",i3,") = xf(",i3,") ***")') i,j
						stop
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
			write(*,'("*** cwlspa: number of distinct values of x(i) with non-zero weight =",i3 " <")')
			write(*,'("            number of coefficients - total number of interpolation conditions ***")')
			stop
		end if
		iwrk(1) = 1
		ierror = 1
		call pidazp(mf,xmin,xmax,xf,yf,ip,n,np1,5,20,wrk,wrk(nanu),wrk(nwrk),lwrk-nwrk+1,iwrk,liwrk,ierror)
		if (ierror /= 0) then
			if (ierror == 1) then
				write(*,'("*** pidazp: iteration limit in deriving q(x) ***")')
			else if (ierror == 2) then
				write(*,'("*** pidazp: diverging iteration in deriving q(x) ***")')
			else if (ierror == 3) then
				write(*,'("*** pidazp: iteration limit in deriving q0(x) ***")')
			else if (ierror == 4) then
				write(*,'("*** pidazp: diverging iteration in deriving q0(x) ***")')
			end if
			stop
		end if
		do i = 1, m
			xcap = ((x(i)-xmin)-(xmax-x(i)))/(xmax-xmin)
			call sumchs(n,wrk,1,np1,xcap,xmu)
			iymux = nwrk + i - 1
			wrk(iymux) = y(i) - xmu
		end do
		nwrk1 = nwrk + m
		nwrk2 = nwrk1 + 2*m
		nser = nwrk2 + 2*kplus1
		neps = nser + kplus1
		ierror = 1
		call wlspad(1,m,m,kplus1,nrows,1,1,x,wrk(nwrk),w,xmin,xmax,np1,wrk(nanu),wrk(nwrk1),wrk(nwrk2),a,s,wrk(nser),wrk(neps),ierror)
		if (ierror /= 0) then
			write(*,'("*** wlspad: ifail =",i2," ***")') ierror
			stop
		end if
		do j = 1, n
			amuj = wrk(j)
			do  i = np1, kplus1
				a(i,j) = a(i,j) + amuj
			end do
		end do

	end subroutine cwlspa

	subroutine pidazp(m,xmin,xmax,x,y,ip,n,np1,itmin,itmax,a,b,wrk,lwrk,iwrk,liwrk,ifail)

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

		real(IDP) :: pmax
		integer :: i, iadif, iatrl, ibdif, ibtrl, ic, id, ida, idb, ierror, iftau, ilocx, ilocy, imax, initq, initq0, ipiq, ipiq0, &
			   iptrl, ires, irnm, irtrnm, iw, ix
		logical :: withq0

		interface
			subroutine aipozp(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w, &
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
			end subroutine aipozp
		end interface

		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP

		ierror = 0
		imax = 0
		do i = 1, m
			if (ip(i) > imax) imax = ip(i)
		end do
		imax = imax + 1

		withq0 = (iwrk(1) /= 0)

		if (n == 1) then
			a(1) = two*y(1)
			if (withq0) then
				b(1) = ((x(1)-xmin)-(xmax-x(1)))/(xmax-xmin)
				b(1) = -two*b(1)
				b(2) = one
			end if

			iwrk(1) = 0
			if (withq0) iwrk(2) = 0
			wrk(1) = zero
			wrk(2) = zero
			if (withq0) wrk(3) = zero
		end if

		ix = 2*imax + n
		do i = 1, m
			ix = ix + 1
			wrk(ix) = ((x(i)-xmin)-(xmax-x(i)))/(xmax-xmin)
		end do
		if (n == 1) then
			ifail = ierror
			return
		end if
		if (withq0) then

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

			call aipozp(.true.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,b,np1,wrk(ires),pmax,wrk(ipiq0),iwrk(initq0), &
				    wrk(ibtrl),wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(ibdif),wrk(idb),wrk(irnm),wrk(irtrnm), &
				    iwrk(ilocx),iwrk(ilocy),ierror)

			if (ierror /= 0) ierror = ierror + 2
			if (ierror /= 0) then
				ifail = ierror
				return
			end if
		end if

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

		call aipozp(.false.,m,wrk(ix),xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,n,wrk(ires),pmax,wrk(ipiq),iwrk(initq),wrk(iatrl), &
			    wrk(iptrl),wrk(iftau),wrk(ic),wrk(id),wrk(iw),wrk(iadif),wrk(ida),wrk(irnm),wrk(irtrnm),iwrk(ilocx),iwrk(ilocy), &
			    ierror)

		ifail = ierror

	end subroutine pidazp

	subroutine aipozp(withq0,m,x,xmin,xmax,y,ip,imax,n,np1,itmin,itmax,a,la,res,pmax,pindex,nit,atrial,ptrial,ftau,c,d,w,adif,da, &
							rnm,rtrlnm,locx,locy,ifail)

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

		real(IDP) :: amax, atrnrm, danrm, pmxtrl, scale
		integer :: i, ierror, it, itemp, itmxp1, itp1, l, nfref, npilt1, nterms
		logical :: improv, withpi, zeroda

		interface
			subroutine chszp(m,x,ip,np1,b,w)
				use param
				implicit none
				integer :: m, np1
				real(IDP), dimension(np1) :: b, w
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip
			end subroutine chszp
			subroutine pirpa(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)
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
			end subroutine pirpa
			subroutine chspi(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)
				use param
				implicit none
				real(IDP) :: xmax, xmin
				integer :: m, n
				real(IDP), dimension(n) :: a, c, d, ftau, y
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx, locy
			end subroutine chspi
		end interface

		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, one=1.0_IDP, sxteen=16.0_IDP

		ierror = 0

		nterms = n
		if (withq0) nterms = nterms + 1

		npilt1 = 0

		withpi = .false.

		nfref = -2

		do i = 1, n
			if ( .not. withq0) res(i) = y(i)
			if (withq0) res(i) = zero
		end do

		do i = 1, nterms
			atrial(i) = zero
		end do

		itmxp1 = itmax + 1
		do itp1 = 1, itmxp1

			it = itp1 - 1

			if (withq0 .and. it == 0) call chszp(m,x,ip,np1,da,w)
			if ( .not. withq0 .or. (withq0 .and. it > 0)) call chspi(m,xmin,xmax,x,res,ip,n,da,locx,locy,ftau,d,c)

			if (it > 0) then

				danrm = half*abs(da(1))
				atrnrm = half*abs(atrial(1))
				if (n > 1) then
					do i = 2, n
						danrm = danrm + abs(da(i))
						atrnrm = atrnrm + abs(atrial(i))
					end do
				end if
				if (withq0) atrnrm = atrnrm + abs(atrial(np1))

				if (danrm >= atrnrm) ierror = 2
				if (danrm >= atrnrm) exit
			end if

			zeroda = .true.
			do i = 1, n
				atrial(i) = atrial(i) + da(i)
				if (da(i) /= zero) zeroda = .false.
			end do
			if (withq0 .and. it == 0) then
				atrial(np1) = da(np1)
				if (da(np1) /= zero) zeroda = .false.
			end if

			call pirpa(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,atrial,la,it,rnm,rtrlnm,improv,adif,res,pmxtrl,ptrial)

			if ( .not. improv) pmxtrl = sxteen

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

			if (pmxtrl == zero .or. zeroda) exit

			if (nfref == -2 .and. pmxtrl < one) nfref = -1

			if (nfref < -1) then

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

		if (itp1 > itmxp1) then
			ierror = 1
			it = itmax
		end if

		nit = it
		if (withq0) then

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

		withpi = .true.
		call pirpa(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnm,rtrlnm,improv,adif,res,pmax,pindex)
		ifail = ierror

	end subroutine aipozp

	subroutine chszp(m,x,ip,np1,b,w)

		use param
		implicit none

		integer :: m, np1
		real(IDP), dimension(np1) :: b, w
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip

		real(IDP) :: ai, anbig, ct, eps, factor, ovfl, ri, sfac, test, unfl, xi, xtrema
		integer :: i, i2, ifail, ip1, k, l, n, nu

		interface
			subroutine cpchs(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine cpchs
		end interface

		real(IDP) :: zero=0.0_IDP, sxtnth=0.0625_IDP, one=1.0_IDP, two=2.0_IDP, sxteen=16.0_IDP, pi=3.14159265358979323846_IDP

		ovfl = sxtnth/(sqrt(2.0_IDP)*tiny(sxtnth))
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

		call cpchs(np1,w,b,ifail)

		ai = b(np1)
		ct = log(abs(ai))/log(two)
		nu = ct
		sfac = two**nu
		b(np1) = sfac
		sfac = sfac/ai
		do i = 1, n
			b(i) = b(i)*sfac
		end do

	end subroutine chszp

	subroutine pirpa(withq0,withpi,m,xmin,xmax,x,n,y,ip,imax,a,la,it,rnmbst,rnm,improv,adif,res,pmax,pindex)

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

		real(IDP) :: absres, asum, asumax, eps, p, plarge, resmax, rscale, t
		integer :: i, ia, iy, l, lm1, nl, nterms
		logical :: imp

		interface
			subroutine difchs(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)
				use param
				implicit none
				real(IDP) :: patm1, xmax, xmin
				integer :: ia1, iadif1, la, ladif, np1
				real(IDP), dimension(la) :: a
				real(IDP), dimension(ladif) :: adif
			end subroutine difchs
			subroutine sumchs(np1,a,ia1,la,xcap,result)
				use param
				implicit none
				real(IDP) :: result, xcap
				integer :: ia1, la, np1
				real(IDP), dimension(la) :: a
			end subroutine sumchs
		end interface

		real(IDP) :: half=0.5_IDP, zero=0.0_IDP, one=1.0_IDP, mltplr=8.0_IDP

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

			asum = half*abs(adif(1))
			if (nterms > 1) then
				do ia = 2, nterms
					asum = asum + abs(adif(ia))
				end do
			end if
			if (asum > asumax) asumax = asum

			pindex(l) = asumax
			iy = l
			nl = 0

			resmax = zero
			rscale = one
			do i = 1, m

				if (ip(i)+1 >= l) then
					nl = nl + 1

					call sumchs(nterms,adif,1,la,x(i),p)

					if (withq0) res(iy) = -p
					if ( .not. withq0) res(iy) = y(iy) - p

					absres = abs(res(iy))
					if (absres /= zero .and. absres <= resmax) rscale = rscale + (absres/resmax)**2
					if (absres /= zero .and. absres > resmax) rscale = rscale*(resmax/absres)**2 + one
					if (absres /= zero .and. absres > resmax) resmax = absres
				end if
				iy = iy + ip(i) + 1
			end do
			rnm(l) = resmax*sqrt(rscale)

			if (l < imax) call difchs(nterms,xmin,xmax,adif,1,la,t,adif,1,la)
		end do

		imp = (it == 0)
		if (it /= 0) then
			do l = 1, imax
				if (rnm(l) < rnmbst(l)) imp = .true.
			end do
		end if
		if (imp .or. withpi) then

			plarge = zero
			do l = 1, imax

!  			nl  is the number of derivative values of order  l - 1

				nl = 0
				lm1 = l - 1
				do i = 1, m
					if (ip(i) >= lm1) nl = nl + 1
				end do
				t = nl

				asumax = pindex(l)
				if (asumax /= zero) pindex(l) = rnm(l)/(mltplr*eps*asumax*sqrt(t))
				if (asumax == zero) pindex(l) = zero
				if (pindex(l) > plarge) plarge = pindex(l)
			end do
			pmax = plarge
		end if
		improv = imp

	end subroutine pirpa

	subroutine chspi(m,xmin,xmax,x,y,ip,n,a,locx,locy,ftau,d,c)

		use param
		implicit none

		real(IDP) :: xmax, xmin
		integer :: m, n
		real(IDP), dimension(n) :: a, c, d, ftau, y
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx, locy

		real(IDP) :: cmin, cnew, factor, ri, rj, s, scale, sfac, v, xch
		integer :: i, i2, ic, icmin, ifail, iftau, ip1, isave, iy, j, jmax, k, krev, l, lmax, locxi, locxj, locxk, locyi, nc

		interface
			subroutine ncnfip(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)
				use param
				implicit none
				real(IDP) :: cnew, xnew, ynew
				integer :: ixnext, m, n, nc, nordp1
				real(IDP), dimension(n) :: c, d
				real(IDP), dimension(m) :: x
				integer, dimension(m) :: ip, locx
			end subroutine ncnfip
			subroutine cpchs(nplus1,f,a,ifail)
				use param
				implicit none
				integer :: ifail, nplus1
				real(IDP), dimension(nplus1) :: a, f
			end subroutine cpchs
		end interface

		real(IDP) :: zero=0.0_IDP, one=1.0_IDP, two=2.0_IDP, pi=3.14159265358979323846_IDP

		scale = (xmax-xmin)/two

		iy = 0
		do i = 1, m
			locx(i) = i
			iy = iy + 1
			locy(i) = iy
			ftau(iy) = y(iy)
			jmax = ip(i)
			if (jmax == 0) cycle

			sfac = one
			do j = 1, jmax
				iy = iy + 1
				rj = j
				sfac = sfac*scale/rj
				ftau(iy) = y(iy)*sfac
			end do
		end do

		nc = 0
	outer:          do j = 1, m
				do i = j, m
					locxi = locx(i)
					locyi = locy(locxi)
					call ncnfip(m,x,ip,n,locx,c,nc,x(locxi),j,ftau(locyi),1,cnew,d)
					if (i > j .and. abs(cnew) >= abs(cmin)) cycle
					cmin = cnew
					icmin = i
				end do
				c(nc+1) = cmin
				isave = locx(j)
				locxj = locx(icmin)
				locx(icmin) = isave
				locx(j) = locxj
				iftau = locy(locxj) - 1
				ip1 = ip(locxj) + 1
	inner:	                do i = 1, ip1
					iftau = iftau + 1
					call ncnfip(m,x,ip,n,locx,c,nc,x(locxj),j,ftau(iftau),i,c(nc+1),d)
					nc = nc + 1
					if (nc == n) exit outer
				end do inner
			end do outer

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

		call cpchs(n,d,a,ifail)

	end subroutine chspi

	subroutine ncnfip(m,x,ip,n,locx,c,nc,xnew,ixnext,ynew,nordp1,cnew,d)

		use param
		implicit none

		real(IDP) :: cnew, xnew, ynew
		integer :: ixnext, m, n, nc, nordp1
		real(IDP), dimension(n) :: c, d
		real(IDP), dimension(m) :: x
		integer, dimension(m) :: ip, locx

		real(IDP) :: dif
		integer :: ic, is, ix, k, locxi

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

	end subroutine ncnfip

	subroutine cpchs(nplus1,f,a,ifail)

		use param
		implicit none

		integer :: ifail, nplus1
		real(IDP), dimension(nplus1) :: a, f

		real(IDP) :: bk, bkp1, bkp2, dk, f0, factor, fli, fln, halffn, piby2n
		integer :: i, ierror, iplus1, j, k, krev, n, n2, nless1

		real(IDP) :: pi=3.14159265358979323846_IDP

		ierror = 0
		if (nplus1 < 2) then
			write(*,'("*** cpchs: nplus1 =",i3," ***")') nplus1
			stop
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

	end subroutine cpchs

	subroutine wlspad(mfirst,mlast,mtot,kplus1,nrows,kall,ndv,x,y,w,xmin,xmax,inup1,nu,work1,work2,a,s,serr,eps,ifail)

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

		real(IDP) :: alpip1, betai, bj, bjp1, bjp2, cil, d, df, di, dim1, dj, epslr, factor, pij, sigmai, wr, wrpr, wrprsq, x1, &
			     xcapr, xm
		integer :: i, ii, im1, inu, iplus1, iplus2, j, jplus1, jplus2, jrev, k, l, m, mdist, mr, r
		logical :: wnz
		real(IDP), dimension(10) :: ci

		k = kplus1 - 1
		inu = inup1 - 1

		m = mlast - mfirst + 1
		i = kplus1 - inu
		if (mfirst < 1 .or. inup1 < 1 .or. kplus1 < inup1 .or. m < i .or. ndv < 1 .or. (kall /= 1 .and. kall /= 0)) then
			ifail = 5
			return
		end if

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

		if (xmin > x(mfirst) .or. xmax < x(mlast) .or. xmin >= xmax) then
			ifail = 1
			return
		end if

		l = k - inu
		if (mdist <= l) then
			ifail = 3
			return
		end if

		if (kall == 1 .and. nrows < kplus1) then
			ifail = 5
			return
		end if
		if (inup1 > 1) then

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

				if (iplus1 == 1) then
					wrpr = w(r)*0.5*work2(1,1)
					work1(2,mr) = wrpr
				else
					j = iplus2
					if (xcapr > 0.5_IDP) then

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

				if (df <= 0.0_IDP) s(l,iplus1) = 0.0_IDP
				if (df > 0.0_IDP) s(l,iplus1) = sqrt(sigmai/df)
			end do

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

	end subroutine wlspad

	subroutine difchs(np1,xmin,xmax,a,ia1,la,patm1,adif,iadif1,ladif)

		use param
		implicit none

		real(IDP) :: patm1, xmax, xmin
		integer :: ia1, iadif1, la, ladif, np1
		real(IDP), dimension(la) :: a
		real(IDP), dimension(ladif) :: adif

		real(IDP) :: ptemp, r, sclftr, u, v, w
		integer :: i, n, na, nadif

		real(IDP) :: two=2.0_IDP

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

				adif(nadif) = sclftr*v
				u = v
				v = w
				na = na - ia1
				nadif = nadif - iadif1
			end do
		end if
		adif(nadif) = sclftr*v
		patm1 = a(1)/two - ptemp

	end subroutine difchs

	subroutine sumchs(np1,a,ia1,la,xcap,result)

		use param
		implicit none

		real(IDP) :: result, xcap
		integer :: ia1, la, np1
		real(IDP), dimension(la) :: a

		real(IDP) :: aj, bj, cj, factor, sum
		integer :: j, jrev, n

		real(IDP) :: zero=0.0_IDP, half=0.5_IDP, two=2.0_IDP

		if (np1 == 1) then
			sum = half*a(1)
		else
			n = np1 - 1
			aj = zero
			bj = zero
			j = 1 + np1*ia1
			if (xcap > half) then

				factor = two - (xcap+xcap)

				do jrev = 1, n
					j = j - ia1
					aj = a(j) + aj - bj*factor
					bj = aj + bj
				end do
				sum = half*a(1) + aj - half*factor*bj
			else if (xcap >= -half) then

				factor = xcap + xcap
				do jrev = 1, n
					j = j - ia1
					cj = bj
					bj = aj
					aj = a(j) - cj + bj*factor
				end do
				sum = half*a(1) - bj + half*factor*aj
			else

				factor = two + (xcap+xcap)

				do jrev = 1, n
					j = j - ia1
					aj = a(j) - aj + bj*factor
					bj = aj - bj
				end do
				sum = half*a(1) - aj + half*factor*bj
			end if
		end if
		result = sum

	end subroutine sumchs
