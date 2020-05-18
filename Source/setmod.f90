	subroutine setmod

		use param
		use cotrol
		use domain
		use dynamo
		implicit none

		interface
			subroutine inputlist
			end subroutine inputlist	
		end interface	
		
		integer :: l,n,m,lh,mnxxx,mt,nt,ngcd,ndiv,lll,nxeq,i,lm,lq,mcnt,mcnt1,n1,nl,lp,lpp,lppn,lop,l1,mband,neq,np,nm,icou,ncou,nlin
		integer, dimension(:), allocatable :: iflag,nfm,nfmo,mcnh		

!		Set up the modes distribution and the mode couplings of the model  
		
		mmin=0
		mmax=0
		nmin=0
		nmax=0
		mmineq=0
		mmaxeq=0
		nmineq=0
		nmaxeq=0
		do l=1,lmax
			m=mm(l)
			n=nn(l)
			signl(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) signl(l)=-1
			if (n == 0 .and. m == 0) signl(l)=0
		end do
		mmax=maxval(mm(1:lmax))
		mmin=minval(mm(1:lmax))
		nmax=maxval(nn(1:lmax))
		nmin=minval(nn(1:lmax))
		mmaxx=max(mmax,abs(mmin))
		do l=1,leqmax
			m=mmeq(l)
			n=nneq(l)
			sgnleq(l)=1
			if (n < 0 .or. (n == 0 .and. m < 0)) sgnleq(l)=-1
			if (n == 0 .and. m == 0) sgnleq(l)=0
		end do
		mmaxeq=maxval(mmeq(1:leqmax))
		mmineq=minval(mmeq(1:leqmax))
		nmaxeq=maxval(nneq(1:leqmax))
		nmineq=minval(nneq(1:leqmax))
		mmaxxeq=max(mmaxeq,abs(mmineq))

		allocate (ll(mmin:mmax,nmin:nmax))

		ll=0
		do l=1,lmax
			ll(mm(l),nn(l))=l
		end do
		l0=0
		if (0 >= mmin .and. 0 <= mmax .and. 0 >= nmin .and. 0 <= nmax) l0=ll(0,0)
		if (l0 == 0) then
			write (6,'("  setmod: l0=0")')
			stop
		end if

		allocate (lleq(mmineq:mmaxeq,nmineq:nmaxeq))

		lleq=0
		do l=1,leqmax
			lleq(mmeq(l),nneq(l))=l
		end do
		leq0=0
		if (0 >= mmineq .and. 0 <= mmaxeq .and. 0 >= nmineq .and. 0 <= nmaxeq) leq0=lleq(0,0)
		if (leq0 == 0) then
			write (6,'("  setmod: leq0=0")')
			stop
		end if

!  find prime harmonics.

		lh=0
		mnxxx=1
		do l=1,lmax

			mt=mm(l)*signl(l)
			nt=nn(l)*signl(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) then
					mnxxx=max(mnxxx,nt)
					nt=1
				end if
				if (mt /= 0) then
					mnxxx=max(mnxxx,mt)
					mt=1
				end if
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mnxxx=max(mnxxx,ngcd)
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mh(lll) .and. nt == nh(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mh(lh)=mt
				nh(lh)=nt
			end if
		end do
		lhmax=lh

		lh=0
		do l=1,leqmax

			mt=mmeq(l)*sgnleq(l)
			nt=nneq(l)*sgnleq(l)

			if (mt == 0 .or. nt == 0) then
				if (nt /= 0) nt=1
				if (mt /= 0) mt=1
			else
				ngcd=1
				do ndiv=2,nt
					if(mt /= ndiv*(mt/ndiv)) cycle
					if(nt /= ndiv*(nt/ndiv)) cycle
					ngcd=ndiv
				end do
				mt=mt/ngcd
				nt=nt/ngcd
			end if

			do lll=1,lh
				if (mt == mheq(lll) .and. nt == nheq(lll)) exit
			end do
			if (lh == 0 .or. lll > lh) then
				lh=lh+1
				mheq(lh)=mt
				nheq(lh)=nt
			end if
		end do
		lheqmx=lh

		nxeq=nmaxeq
		if (nmaxeq > 0) then
			do l=1,leqmax
				if (nneq(l) > 0) nxeq=min(nxeq,nneq(l))
			end do
			do l=1,leqmax
				if (nneq(l) /= nxeq*(nneq(l)/nxeq)) then
					write (6,'("  setmod: nneq(",i2,") =",i4," is not multiple of nxeq =",i3)') l,nneq(l),nxeq
					stop
				end if
			end do
		else
			nxeq=0
		end if

		lmaxn=lmax-leqmax
		if (m0dy < 0) m0dy=0
		lmaxn=lmaxn+m0dy
		lmx=noeqn*lmaxn

		allocate (m1n(nmin:nmax))
		allocate (mrang(nmin:nmax))

		do n=nmin,nmax
			m1n(n)=mmax+1
			mrang(n)=mmin-1
			do l=1,lmaxn
				if (nn(l) /= n) cycle
				m1n(n)=min(m1n(n),mm(l))
				mrang(n)=max(mrang(n),mm(l))
			end do
			if (mrang(n) == mmin-1) then
				mrang(n)=0
			else
				mrang(n)=mrang(n)-m1n(n)+1
			endif
		end do
		if (nmin < 0 .and. nmax > 0) then
			if (nmin /= -nmax) then
				write (6,'("  setmod: nmin =",i4," nmax =",i3)') nmin,nmax
				stop
			end if			
			do n=1,nmax
				if (mrang(n)+mrang(-n) > 0 .and. (mrang(-n) /= mrang(n) .or. m1n(n) /= -(m1n(-n)+mrang(-n)-1))) then
					write (6,'("  setmod: n =",i4," mrang =",i3,i4," m1n =",i3,i5)') n,mrang(n),mrang(-n),m1n(n),m1n(-n)+mrang(-n)-1
					stop
				end if
			end do
			if (mrang(0) /= 0 .and. (mod(mrang(0),2) == 0 .or. m1n(0)+(mrang(0)-1)/2 /= 0)) then
				write (6,'("  setmod: n = 0  mrang =",i3," m1n =",i3,i5)') mrang(0),m1n(0),m1n(0)+mrang(0)-1
				stop
			end if
		end if

		lmax0=0
		do i=1,mrang(0)
			m=m1n(0)+i-1
			l=ll(m,0)
			if (l == 0) cycle
			lmax0=lmax0+1
		end do

		if (lmax0 > 0) then

			allocate (ll0(lmax0))

			lm=0
			do i=1,mrang(0)
				m=m1n(0)+i-1
				l=ll(m,0)
				if (l == 0) cycle
				lm=lm+1
				ll0(lm)=l
			end do

		end if

		allocate (lln(lmaxn))
		allocate (iflag(nmin:nmax))

		iflag=0
		nnum=0
		n1=0

		if (nocpl == 0) then

			if (nxeq > 0) then

!	helical couplings.

				allocate (nfm(nmax+1))
				allocate (nfmo(0:nmax))
				allocate (mcnh(nmax+1))
				nfmo=0

				do n=1,nmax

					if (mrang(n) == 0 .or. iflag(n) == 1) cycle
					ncou=0
					do i=0,n-1
						ncou=ncou+iflag(i)
					end do
					if (ncou > 0) then
						icou=0
						do neq=nxeq,nmaxeq,nxeq
							nm=abs(neq-n)
							icou=icou+iflag(nm)
							if (iflag(nm) == 1) n1=nfmo(nm) 
						end do
					end if
					if (ncou == 0 .or. icou == 0) then
						nnum=nnum+1
						nfm(nnum)=n
						nfmo(n)=nnum
						mcnh(nnum)=0
						n1=nnum
					end if
					mcnt1=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnh(n1)=mcnh(n1)+1
						mcnt1=mcnt1+1
					end do
					if (mcnh(nnum) == 0) then
						nnum=nnum-1
					else
						iflag(n)=1
					end if

					if (mcnt1 == 0) cycle
					nfmo(n)=n1

					do neq=nxeq,nmaxeq,nxeq

						nm=abs(neq-n)
						if (mrang(nm) > 0 .and. iflag(nm) == 0) then
							do i=1,mrang(nm)
								m=m1n(nm)+i-1
								l=ll(m,nm)
								lq=0
								if (nm >= nmineq .and. nm <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,nm)
								if (l == 0 .or. lq > m0dy) cycle
								mcnh(n1)=mcnh(n1)+1
							end do
							iflag(nm)=1
							nfmo(nm)=n1
						end if

						np=abs(neq+n)
						if (np > nmax) cycle
						if (mrang(np) > 0 .and. iflag(np) == 0) then
							do i=1,mrang(np)
								m=m1n(np)+i-1
								l=ll(m,np)
								lq=0
								if (np >= nmineq .and. np <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,np)
								if (l == 0 .or. lq > m0dy) cycle
								mcnh(n1)=mcnh(n1)+1
							end do
							iflag(np)=1
							nfmo(np)=n1
						end if

					end do

				end do

				allocate (mnumn(nnum))
				allocate (lnumn(0:nnum))

				nnd=0
				lnumn(0)=0
				lp=0
				do nl=1,nnum

					lpp=0
					do n=1,nmax
						if (nfmo(n) /= nl) cycle
						do i=1,mrang(n)
							m=m1n(n)+i-1
							l=ll(m,n)
							lq=0
							if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
							if (l == 0 .or. lq > m0dy) cycle
							lp=lp+1
							lln(lp)=l
							lpp=lpp+1
						end do
					end do
					if (nfmo(0) == nl) then
						if (nmin < 0) then
							do i=(mrang(0)+3)/2,mrang(0)
								m=m1n(0)+i-1
								l=ll(m,0)
								lq=0
								if (0 >= nmineq .and. 0 <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,0)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
							l=ll(0,0)
							lq=0
							if (0 >= nmineq .and. 0 <= nmaxeq .and. 0 >= mmineq .and. 0 <= mmaxeq) lq=lleq(0,0)
							if (l /= 0 .and. lq <= m0dy) then
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end if
							do i=(mrang(0)-1)/2,1,-1
								m=m1n(0)+i-1
								l=ll(m,0)
								lq=0
								if (0 >= nmineq .and. 0 <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,0)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
						else
							do i=2,mrang(0)
								m=m1n(n)+i-1
								l=ll(m,n)
								lq=0
								if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
							l=ll(0,0)
							lq=0
							if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
							if (l /= 0 .and. lq <= m0dy) then
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end if
						end if
					end if
					if (nmin < 0) then
						do n=-1,nmin,-1
							if (nfmo(-n) /= nl) cycle
							do i=mrang(n),1,-1
								m=m1n(n)+i-1
								l=ll(m,n)
								lq=0
								if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
								if (l == 0 .or. lq > m0dy) cycle
								lp=lp+1
								lln(lp)=l
								lpp=lpp+1
							end do
						end do
					end if
					mnumn(nl)=lpp
					lnumn(nl)=lp

				end do

			else

!	toroidal couplings

				do n=nmin,nmax
					if (mrang(n) == 0) cycle
					if (n < 0 .and. nmax > 0) cycle
					nnum=nnum+1
					mcnt=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						mcnt=mcnt+1
					end do
					if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
						do i=1,mrang(-n)
							m=m1n(-n)+i-1
							l=ll(m,-n)
							lq=0
							if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
							if (l == 0 .or. lq > m0dy) cycle
							mcnt=mcnt+1
						end do
					end if
					if (mcnt == 0) nnum=nnum-1
				end do

				allocate (mnumn(nnum))
				allocate (lnumn(0:nnum))

				nnd=0
				lnumn(0)=0
				nl=0
				lp=0
				do n=nmin,nmax
					if (mrang(n) == 0) cycle
					if (n < 0 .and. nmax > 0) cycle
					nl=nl+1
					lpp=0
					lppn=0
					do i=1,mrang(n)
						m=m1n(n)+i-1
						l=ll(m,n)
						lq=0
						if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
						if (l == 0 .or. lq > m0dy) cycle
						lp=lp+1
						lln(lp)=l
						lpp=lpp+1
					end do
					if (-n >= nmin .and. -n <= nmax .and. n /= 0) then
						do i=mrang(-n),1,-1
							m=m1n(-n)+i-1
							l=ll(m,-n)
							lq=0
							if (-n >= nmineq .and. -n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,-n)
							if (l == 0 .or. lq > m0dy) cycle
							lp=lp+1
							lln(lp)=l
							lppn=lppn+1
						end do
					end if
					if (lpp > 0) then
						mnumn(nl)=lpp
						lnumn(nl)=lp
						if (nmin == 0 .or. nmax == 0 .or. n == 0) then
							nnd=nnd+1
						else if (lppn /= lpp) then
							write (6,'("  setmod: n =",i4," mnumn =",i3,i4)') n,lpp,lppn
							stop
						end if
					else
						nl=nl-1
					end if
				end do
				nst=nnd+1
			
			end if

		else

!	no couplings, cylinder.

			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if(l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nnum=nnum+1
					mcnt=1
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) mcnt=mcnt+1
				end do
			end do

			allocate (mnumn(nnum))
			allocate (lnumn(0:nnum))

			nnd=0
			lnumn(0)=0
			nl=0
			lp=0
			do n=nmin,nmax
				if (mrang(n) == 0) cycle
				do i=1,mrang(n)
					m=m1n(n)+i-1
					l=ll(m,n)
					if (l == 0) cycle
					lop=0
					if (-m >= mmin .and. -m <= mmax .and. -n >= nmin .and. -n <= nmax) lop=ll(-m,-n)
					if (signl(l) < 0 .and. lop /= 0) cycle
					lq=0
					if (n >= nmineq .and. n <= nmaxeq .and. m >= mmineq .and. m <= mmaxeq) lq=lleq(m,n)
					if (lq > m0dy) cycle
					nl=nl+1
					lp=lp+1
					lpp=1
					lln(lp)=l
					if (lop /= 0 .and. (m /= 0 .or. n /= 0)) then
						lp=lp+1
						lpp=lpp+1
						lln(lp)=lop
					end if
					mnumn(nl)=lpp
					lnumn(nl)=lp
					if (lpp == 1) nnd=nnd+1
				end do
			end do
			nst=nnd+1
			
			write(6,'(/" number of n-values with one parity:   ",i3)') nnd
			write(6,'(" number of n-values with both parities:",i3/)') nnum-nnd

		end if

		allocate (lo(lmaxn))

		lo=0
		do l=1,lmaxn
			m=-mm(lln(l))
			n=-nn(lln(l))
			if (m < mmin .or. m > mmax .or. n < nmin .or. n > nmax) cycle
			do l1=1,lmaxn
				if (mm(lln(l1)) == m .and. nn(lln(l1)) == n) exit
			end do
			if (l1 > lmaxn) then
				lo(l)=0
			else
				lo(l)=l1
			end if
		end do

		mband=0
		do m=0,mmax
			if (ll(m,0) /= 0) mband=mband+1
		end do
		mxmband=2*mband-1
		do n=1,nmax
			mband=0
			do m=mmin,mmax
				if (ll(m,n) /= 0) mband=mband+1
			end do
			if (mband > mxmband) mxmband=mband
		end do
	end subroutine setmod