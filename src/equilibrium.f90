MODULE equilibrium

  use param
  use cotrol
  use domain
  use equil
  use dynamo
  
CONTAINS

  subroutine setmod

    use mpi
    use processor
    use var_para

    implicit none

    integer :: l,n,m,lh,mt,nt,ngcd,ndiv,lll,nxeq,i,lm,lq,mcnt,mcnt1,n1,nl,lp,lpp,lppn,lop,l1,mband,neq,np,nm,icou,ncou, &
         idel,idelp1,i_mn,ip,mjst,mjnd,iPE,mxnum,i1,ierror
    integer, dimension(1) :: PE_min
    integer, dimension(:), allocatable :: nsort,nmap,nnum_PE,nrsq_PE,nrl_PE,nvalue,iflag,nfm,nfmo,mcnh
    integer, dimension(:,:), allocatable :: nr_PE

    !  Set up the modes distribution and the mode couplings of the model  

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
    lasym=.FALSE.
    do l=1,leqmax
       m=mmeq(l)
       n=nneq(l)
       sgnleq(l)=1
       if (n < 0 .or. (n == 0 .and. m < 0)) sgnleq(l)=-1
       if (n == 0 .and. m == 0) sgnleq(l)=0
       if (sgnleq(l) == -1) lasym=.TRUE.
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
       if (myPE == 0) write (6,'("  setmod: l0=0")')
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
       if (myPE == 0) write (6,'("  setmod: leq0=0")')
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

    !  if (myPE == 0) then
    !     write(6,'(/"helicities"/"   lh   mh   nh")')
    !     write(6,'(3i5)') (lh,mh(lh),nh(lh),lh=1,lhmax)
    !     write(6,'(/" mnxxx = ",i3)') mnxxx
    !  end if

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
             if (myPE == 0) write (6,'("  setmod: nneq(",i2,") =",i4," is not multiple of nxeq =",i3)') l,nneq(l),nxeq
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
          if (myPE == 0) write (6,'("  setmod: nmin =",i4," nmax =",i3)') nmin,nmax
          stop
       end if
       do n=1,nmax
          if (mrang(n)+mrang(-n) > 0 .and. (mrang(-n) /= mrang(n) .or. m1n(n) /= -(m1n(-n)+mrang(-n)-1))) then
             if (myPE == 0) write (6,'("  setmod: n =",i4," mrang =",i3,i4," m1n =",i3,i5)') &
                  n,mrang(n),mrang(-n),m1n(n),m1n(-n)+mrang(-n)-1
             stop
          end if
       end do
       if (mrang(0) /= 0 .and. (mod(mrang(0),2) == 0 .or. m1n(0)+(mrang(0)-1)/2 /= 0)) then
          if (myPE == 0) write (6,'("  setmod: n = 0  mrang =",i3," m1n =",i3,i5)') mrang(0),m1n(0),m1n(0)+mrang(0)-1
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

          ! helical couplings.

          ! nnum = number of helical families. one for linear run
          ! nfm(nh) = first n-value of nh-family.
          ! nfmo(n) = nh-family corresponding to toroidal number n
          ! mcnh(nh) = number of modes for nh-family

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

          do n1=1,nnum

             if (mod(nfm(n1),nxeq) /= 0) cycle
             n=0
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
             if (mcnt1 > 0) nfmo(n)=n1

          end do

          allocate (mnumn(nnum))
          allocate (lnumn(0:nnum))

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

          allocate (n_st(0:numPEsm1))
          allocate (n_nd(0:numPEsm1))

          if (nnum <= numPEs) then

             numPElm1=nnum-1
             do iPE=1,nnum
                n_st(iPE-1)=iPE
                n_nd(iPE-1)=iPE
             end do
             do n=nnum,numPEsm1
                n_st(n)=1
                n_nd(n)=0
             end do
             
          else

             numPElm1=numPEsm1
             np=nnum/numPEs
             nm=mod(nnum,numPEs)
             i1=0
             do iPE=1,nm
                n_st(iPE-1)=i1+1
                n_nd(iPE-1)=i1+np+1
                i1=n_nd(iPE-1)
             end do
             do iPE=nm+1,numPEs
                n_st(iPE-1)=i1+1
                n_nd(iPE-1)=i1+np
                i1=n_nd(iPE-1)
             end do
             
          end if

          if (myPE == 0) write (6,'(/"lln="/("    ",25i5))') (lln(l),l=1,lmaxn)

       else

          ! toroidal couplings

          do n=nmin,nmax
             if (mrang(n) == 0) cycle
             if (n < 0 .and. mrang(-n) /= 0) cycle
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
          allocate (nsort(nnum))
          allocate (nvalue(nnum))
          allocate (nmap(nnum))
          allocate (n_st(0:numPEsm1))
          allocate (n_nd(0:numPEsm1))

          nnum=0

          do n=nmin,nmax
             if (mrang(n) == 0) cycle
             if (n < 0 .and. mrang(-n) /= 0) cycle
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
             if (mcnt == 0) then 
                nnum=nnum-1
             else
                nvalue(nnum)=n
                mnumn(nnum)=mcnt
             end if
          end do

          if (nnum <= numPEs) then

             numPElm1=nnum-1
             do n=1,nnum
                nmap(n)=nvalue(n)
                n_st(n-1)=n
                n_nd(n-1)=n
             end do
             do n=nnum,numPEsm1
                n_st(n)=1
                n_nd(n)=0
             end do

          else

             numPElm1=numPEsm1
             do n=1,nnum
                nsort(n)=n
             end do
             do n=2,nnum
                nl=n
                do while (mnumn(nsort(nl)) < mnumn(nsort(nl-1)))
                   np=nsort(nl)
                   nsort(nl)=nsort(nl-1)
                   nsort(nl-1)=np
                   nl=nl-1
                   if (nl == 1) exit
                end do
             end do

             nl=nnum/(2*numPEs)
             mxnum=2*nl+mod(nnum,2*numPEs)
             allocate (nnum_PE(numPEs))
             allocate (nrsq_PE(numPEs))
             allocate (nrl_PE(numPEs))
             allocate (nr_PE(mxnum,numPEs))
             nnum_PE=0
             nrsq_PE=0
             nrl_PE=0
             n=nnum
             if (nl > 0) then
                do i=1,nl
                   do iPE=numPEs,1,-1
                      nr_PE(2*i-1,iPE)=nsort(n)
                      nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
                      nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
                      n=n-1
                   end do
                   do iPE=1,numPEs
                      nr_PE(2*i,iPE)=nsort(n)
                      nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
                      nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
                      n=n-1
                   end do
                end do
                nnum_PE=2*nl
                do while (n > 0)
                   PE_min=minloc(nrsq_PE)
                   iPE=PE_min(1)
                   nnum_PE(iPE)=nnum_PE(iPE)+1
                   nr_PE(nnum_PE(iPE),iPE)=nsort(n)
                   nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
                   nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
                   n=n-1
                end do
             else if(nnum < numPEs) then
                iPE=numPEs
                do i=nnum,1,-1
                   nnum_PE(iPE)=1
                   nr_PE(1,iPE)=nsort(i)
                   nrsq_PE(iPE)=mnumn(nsort(i))*mnumn(nsort(i))
                   nrl_PE(iPE)=mnumn(nsort(i))
                   iPE=iPE-1
                end do
             else
                do iPE=numPEs,1,-1
                   nr_PE(1,iPE)=nsort(n)
                   nrsq_PE(iPE)=mnumn(nsort(n))*mnumn(nsort(n))
                   nrl_PE(iPE)=mnumn(nsort(n))
                   n=n-1
                end do
                nnum_PE=1
                do while (n > 0)
                   PE_min=minloc(nrsq_PE)
                   iPE=PE_min(1)
                   nnum_PE(iPE)=nnum_PE(iPE)+1
                   nr_PE(nnum_PE(iPE),iPE)=nsort(n)
                   nrsq_PE(iPE)=nrsq_PE(iPE)+mnumn(nsort(n))*mnumn(nsort(n))
                   nrl_PE(iPE)=nrl_PE(iPE)+mnumn(nsort(n))
                   n=n-1
                end do
             end if

             do iPE=1,numPEs
                do i=1,nnum_PE(iPE)
                   if (nvalue(nr_PE(i,iPE)) == 0) exit
                end do
                if (i <= nnum_PE(iPE)) exit
             end do
             if (iPE > numPEs) then
                if (myPE == 0) write (6,'("  setmod: No modes with n = 0")')
                iPE=1
             end if
             if (iPE > 1) then
                do i=1,max(nnum_PE(1),nnum_PE(iPE))
                   n=nr_PE(i,1)
                   nr_PE(i,1)=nr_PE(i,iPE)
                   nr_PE(i,iPE)=n
                end do
                n=nnum_PE(1)
                nnum_PE(1)=nnum_PE(iPE)
                nnum_PE(iPE)=n
             end if

             n=1
             n_st(0)=1
             n_nd(0)=nnum_PE(1)
             do i=nnum_PE(1),1,-1
                nmap(n)=nvalue(nr_PE(i,1))
                n=n+1
             end do
             do iPE=2,numPEs
                n_st(iPE-1)=n_nd(iPE-2)+1
                n_nd(iPE-1)=n_nd(iPE-2)+nnum_PE(iPE)
                do i=nnum_PE(iPE),1,-1
                   nmap(n)=nvalue(nr_PE(i,iPE))
                   n=n+1
                end do
             end do

          end if

          ! if (myPE == 0) then
          !    write(6,'(/,"n_st=",20i5)') (n_st(iPE),iPE=0,numPEsm1)
          !    write(6,'("n_nd=",20i5)') (n_nd(iPE),iPE=0,numPEsm1)
          !    write(6,'("nmap=",20i5)') (nmap(n),n=1,nnum)
          ! end if

          lnumn(0)=0
          nl=0
          lp=0
          do np=1,nnum
             n=nmap(np)
             if (mrang(n) == 0) cycle
             if (n < 0 .and. mrang(-n) /= 0) cycle
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
                   lpp=lpp+1
                   lppn=lppn+1
                end do
             end if
             if (lpp > 0) then
                mnumn(nl)=lpp
                lnumn(nl)=lp
                if (n /= 0 .and. 2*lppn /= lpp) then
                   if (myPE == 0) write (6,'("  setmod: n =",i4," mnumn =",i3,i4)') n,lpp,2*lppn
                   stop
                end if
             else
                nl=nl-1
             end if
          end do

          if (myPE == 0) write (6,'(/"lln="/("    ",25i5))') (lln(l),l=1,lmaxn)

       end if
       
    else

       ! no couplings, cylinder.

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
          end do
       end do

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

    mband=0
    do m=0,mmaxeq
       if (lleq(m,0) /= 0) mband=mband+1
    end do
    mbandeq=2*mband-1
    do n=1,nmaxeq
       mband=0
       do m=mmineq,mmaxeq
          if (lleq(m,n) /= 0) mband=mband+1
       end do
       if (mband > mbandeq) mbandeq=mband
    end do

    allocate (mj_br(0:numPEsm1))
    allocate (mj_inc(0:numPEsm1))
    allocate (mj_st(0:numPEsm1))
    allocate (mj_dl(0:numPEsm1))

    mjm1=mj-1
    mjm2=mj-2
    idel=mjm1/numPEs
    idelp1=idel+1
    i_mn=numPEs*idelp1-mjm1
    if (i_mn > 1) i_mn=i_mn-1
    do ip=0,numPEsm1
       if (ip < i_mn) then
          mjst=ip*idel
          mjnd=(ip+1)*idel+1
       else
          mjst=ip*idelp1-i_mn
          mjnd=(ip+1)*idelp1-i_mn+1
       end if
       if (ip == numPEsm1) mjnd=mj
       mj_br(ip)=mjst+1
       if (ip == 0) mj_br(ip)=mjst
       mj_inc(ip)=mjnd-mj_br(ip)
       if (ip == numPEsm1) mj_inc(ip)=mj_inc(ip)+1
       mj_st(ip)=mjst
       mj_dl(ip)=mjnd-mjst+1
       if (ip == myPE) mj_start=mjst
       if (ip == myPE) mj_end=mjnd
    end do

  end subroutine setmod

  subroutine seteq

    implicit none

    integer :: i,j
    real(IDP) :: rsq

    !  Equilibrium set up

    !  User defined thermal plasma density profile

    if (cnep(0) > 0.0_IDP) then
       denseq=cnep(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             denseq(j)=denseq(j)+cnep(i)*rsq**i
          end do
       end do
       denseq=denseq/cnep(0)
    end if

    !  User defined thermal electron plasma temperature profile
    if (ctep(0) > 0.0_IDP) then
       teeq=ctep(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             teeq(j)=teeq(j)+ctep(i)*rsq**i
          end do
       end do
       teeq=teeq/ctep(0)
       do j=0,mj
          if (teeq(j) < 1.e-4_IDP) teeq(j)=1.e-4_IDP
       end do
    end if

    !  User defined energetic particles density profile
    if (cnfp(0) > 0.0_IDP) then
       nfeq=cnfp(0)
       do j=0,mj
          rsq=r(j)*r(j)
          do i=1,10
             nfeq(j)=nfeq(j)+cnfp(i)*rsq**i
          end do
       end do
       nfeq=nfeq/cnfp(0)
       do j=0,mj
          if (nfeq(j) < 1.e-4_IDP) nfeq(j)=1.e-4_IDP
       end do
    end if

    !  User defined profile for NBI density, parameter EP_dens_on, Adens and Bdens in input  !!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (EP_dens_on .eq. 1) then
       do j=0,mj
          nfeq(j) = (.5*(1.+tanh(Adens*(Bdens-r(j)))) + 0.02)/(.5*(1.+tanh(Adens*Bdens)) + 0.02) 
          if (nfeq(j) < 1.e-4_IDP) nfeq(j)=1.e-4_IDP
       end do
    end if

    !  User defined energetic particles parallel velocity profile
    if (cvfp(0) > 0.0_IDP) then
       vfova=cvfp(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             vfova(j)=vfova(j)+cvfp(i)*rsq**i
          end do
       end do
       vfova=vfova/LcA3
       do j=0,mj
          if (vfova(j) < 1.e-4_IDP) vfova(j)=1.e-4_IDP
       end do
    end if

    if (alpha_on .eq. 1) then

       !  User defined 2nd species energetic particles density profile
       if (cnfpalp(0) > 0.0_IDP) then
          nalpeq=cnfpalp(0)
          do j=0,mj
             rsq=r(j)*r(j)
             do i=1,10
                nalpeq(j)=nalpeq(j)+cnfpalp(i)*rsq**i
             end do
          end do
          nalpeq=nalpeq/cnfpalp(0)
          do j=0,mj
             if (nalpeq(j) < 1.e-4_IDP) nalpeq(j)=1.e-4_IDP
          end do
       end if

       !  User defined profile for 2nd EP species density, parameter Alpha_dens_on, Adensalp and Bdensalp in input  !!!!!!!

       if (Alpha_dens_on .eq. 1) then
          do j=0,mj
             nalpeq(j) = (.5*(1.+tanh(Adensalp*(Bdensalp-r(j)))) + 0.02)/(.5*(1.+tanh(Adensalp*Bdensalp)) + 0.02)
             if (nalpeq(j) < 1.e-4_IDP) nalpeq(j)=1.e-4_IDP
          end do
       end if

       !  User defined 2nd species energetic particles parallel velocity profile	
       if (cvfpalp(0) > 0.0_IDP) then
          valphaova=cvfpalp(0)
          do j=1,mj
             rsq=r(j)*r(j)
             do i=1,10
                valphaova(j)=valphaova(j)+cvfpalp(i)*rsq**i
             end do
          end do
          valphaova=valphaova/LcA3alp
          do j=0,mj
             if (valphaova(j) < 1.e-4_IDP) valphaova(j)=1.e-4_IDP
          end do
       end if

    end if

    !  Ion flr components: normalized thermal electron velocity profile if no external profile

    if (cvep(0) > 0.0_IDP) then
       vtherm_elc=cvep(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             vtherm_elc(j)=vtherm_elc(j)+cvep(i)*rsq**i
          end do
       end do
       do j=0,mj
          if (vtherm_elc(j) < 1.e-4_IDP) vtherm_elc(j)=1.e-4_IDP
       end do
    end if

    !  Thermal ion toroidal velocity profile if no external profile

    if (eqvt(0) > 0.0_IDP) then
       vzt_eq=eqvt(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             vzt_eq(j)=vzt_eq(j)+eqvt(i)*rsq**i
          end do
       end do
       do j=0,mj
          if (vzt_eq(j) < 1.e-4_IDP) vzt_eq(j)=1.e-4_IDP
       end do
    end if

    !  Thermal ion poloidal velocity profile if no external profile

    if (eqvp(0) > 0.0_IDP) then
       vth_eq=eqvp(0)
       do j=1,mj
          rsq=r(j)*r(j)
          do i=1,10
             vth_eq(j)=vth_eq(j)+eqvp(i)*rsq**i
          end do
       end do
       do j=0,mj
          if (vth_eq(j) < 1.e-4_IDP) vth_eq(j)=1.e-4_IDP
       end do
    end if

    !  The equilibria main parameters are calculated	

    call vmec

  end subroutine seteq

  subroutine vmec

    ! reads in high beta equilibrium data from file fareq
    ! on unit 25.  then sets up equilibrium quantities
    ! required by far.

    use processor
    use multb_mod
    use dbydb
    use dbyd
    use fit_mod
    use scratch

    implicit none

    integer :: j,l,l1,lq,mjeqp,mjeqm1,mpol,ntor,lb0,m,m2,mm1,mm2,mm12,mm22,mm2p1,mm2p12,mp2,mp22,isg,idummy,mband,n,lbm2, &
         ir,jp,lp,nbst
    real(IDP) :: pror,dummy,twopi,twophip,one,gc,iota,iotap,pprime,qmin,qmax,p,xl,pprimel
    integer, dimension(:), allocatable :: llc,lls,lbst
    real(IDP), dimension(:), allocatable :: rbinv,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,rsb,Abst,Wbst
    real(IDP), dimension(:,:), allocatable :: rmnb,zmnb,pmnsb,sqgib,sqgb,grrb,grtb,gttb,bmodb,grrojb,grtojb,gttojb, &
         jbgrrb,jbgrtb,jbgttb
    real(IDP), dimension(:,:), allocatable :: sqgieq,sqgeq,grreq,grteq,gtteq,bmodeq,bsteq,jbgrreq,jbgrteq,jbgtteq, &
         grrojeq,grtojeq,gttojeq,lplreq,lplteq,djrojeq,djtojeq,djzojeq,omdreq,omdteq,omdzeq,omdrprpeq,omdtprpeq, &
         omdzprpeq,dbsjtojeq,dbsjzojeq,dbsjtbjeq,dgttreq,dgrrteq,dgrtteq,dgttteq,dgrrzeq, &
         dgrtzeq,dgttzeq,dgrtpeq,dgttpeq,jsqeq,bsgrteq,bsgtteq,bsqeq,bsqgtteq, &
         lplrreq,lplrteq,lplrzeq,lpltreq,lpltteq,lpltzeq,lplzreq,lplzteq,lplzzeq, &
         eildreq,eildteq,eildzeq,eildrreq,eildrteq,eildrzeq,eildtteq,eildtzeq,eildzzeq, &
         sb1,sb2,sb3,sb4,sb5,sb6
    ! character(len=1) :: t
    ! character(len=32) :: formatt='("r",9(a1,i4,"/",i4))'
    ! character(len=32) :: formatv='(1pe13.6,9(a1,1pe15.8))'

    allocate (sqg(0:mj,0:leqmax),sqgi(0:mj,0:leqmax),bmod(0:mj,0:leqmax),bst(0:mj,0:leqmax), &
         grr(0:mj,0:leqmax),grt(0:mj,0:leqmax),gtt(0:mj,0:leqmax), &
         grroj(0:mj,0:leqmax),grtoj(0:mj,0:leqmax),gttoj(0:mj,0:leqmax), &
         jbgrr(0:mj,0:leqmax),jbgrt(0:mj,0:leqmax),jbgtt(0:mj,0:leqmax),lplr(0:mj,0:leqmax), &
         lplt(0:mj,0:leqmax),lplz(0:mj,0:leqmax),djroj(0:mj,0:leqmax),djtoj(0:mj,0:leqmax),djzoj(0:mj,0:leqmax), &
         jsq(0:mj,0:leqmax),omdr(0:mj,0:leqmax),omdt(0:mj,0:leqmax),omdz(0:mj,0:leqmax), &
         dbsjtoj(0:mj,0:leqmax),dbsjzoj(0:mj,0:leqmax),dbsjtbj(0:mj,0:leqmax),dgttr(0:mj,0:leqmax), &
         dgrrt(0:mj,0:leqmax),dgrtt(0:mj,0:leqmax),dgttt(0:mj,0:leqmax),dgrrz(0:mj,0:leqmax), &
         dgrtz(0:mj,0:leqmax),dgttz(0:mj,0:leqmax),dgrtp(0:mj,0:leqmax),dgttp(0:mj,0:leqmax), &
         bsgrt(0:mj,0:leqmax),bsgtt(0:mj,0:leqmax),bsq(0:mj,0:leqmax),bsqgtt(0:mj,0:leqmax), &
         lplrr(0:mj,0:leqmax),lplrt(0:mj,0:leqmax),lplrz(0:mj,0:leqmax),lpltr(0:mj,0:leqmax),lpltt(0:mj,0:leqmax), &
         lpltz(0:mj,0:leqmax),lplzr(0:mj,0:leqmax),lplzt(0:mj,0:leqmax),lplzz(0:mj,0:leqmax))
    if (ieldamp_on == 1) allocate (eildr(0:mj,0:leqmax),eildt(0:mj,0:leqmax),eildz(0:mj,0:leqmax),eildrr(0:mj,0:leqmax), &
         eildrt(0:mj,0:leqmax),eildrz(0:mj,0:leqmax),eildtt(0:mj,0:leqmax),eildtz(0:mj,0:leqmax), &
         eildzz(0:mj,0:leqmax))
    if (trapped_on == 1) allocate (omdrprp(0:mj,0:leqmax),omdtprp(0:mj,0:leqmax),omdzprp(0:mj,0:leqmax))

    open(unit=25,file=eq_name,status='old',convert='big_endian',form='unformatted')

    read(25) nfp,lbmax,mjeqp,lasym

    if (nfp == 1) then
       ndevice(1)=" 2D equi"
       ndevice(2)="librium "
    else
       ndevice(1)=" 3D equi"
       ndevice(2)="librium "
    end if

    mjeq=mjeqp-1
    mjeqm1=mjeq-1
    lbm2=2*lbmax-1

    if (lasym) then
       allocate (mmb(lbm2),nnb(lbm2),llc(lbmax),lls(lbmax))
    else
       allocate (mmb(lbmax),nnb(lbmax))
    end if
    allocate (rfar(0:mjeq),rbinv(0:mjeq),qfar(0:mjeq),pfar(0:mjeq),phip(0:mjeq),curfar(0:mjeq),ffar(0:mjeq), &
         sfar1(0:mjeq),sfar2(0:mjeq),sfar3(0:mjeq),sfar4(0:mjeq))
    if (lasym) then
       allocate (rmnb(0:mjeq,0:lbm2),zmnb(0:mjeq,0:lbm2),sqgb(0:mjeq,0:lbm2),sqgib(0:mjeq,0:lbm2), &
            bmodb(0:mjeq,0:lbm2),grrb(0:mjeq,0:lbm2),grtb(0:mjeq,0:lbm2),gttb(0:mjeq,0:lbm2), &
            grrojb(0:mjeq,0:lbm2),grtojb(0:mjeq,0:lbm2),gttojb(0:mjeq,0:lbm2), &
            jbgrrb(0:mjeq,0:lbm2),jbgrtb(0:mjeq,0:lbm2),jbgttb(0:mjeq,0:lbm2), &
            pmnsb(0:mjeq,0:lbm2))
    else
       allocate (rmnb(0:mjeq,0:lbmax),zmnb(0:mjeq,0:lbm2),sqgb(0:mjeq,0:lbmax),sqgib(0:mjeq,0:lbmax), &
            bmodb(0:mjeq,0:lbmax),grrb(0:mjeq,0:lbmax),grtb(0:mjeq,0:lbmax),gttb(0:mjeq,0:lbmax), &
            grrojb(0:mjeq,0:lbmax),grtojb(0:mjeq,0:lbmax),gttojb(0:mjeq,0:lbmax), &
            jbgrrb(0:mjeq,0:lbmax),jbgrtb(0:mjeq,0:lbmax),jbgttb(0:mjeq,0:lbmax), &
            pmnsb(0:mjeq,0:lbmax))
    end if

    rmnb(:,0)=0.0_IDP
    zmnb(:,0)=0.0_IDP
    pmnsb(:,0)=0.0_IDP
    sqgib(:,0)=0.0_IDP
    sqgb(:,0)=0.0_IDP
    grrb(:,0)=0.0_IDP
    grtb(:,0)=0.0_IDP
    gttb(:,0)=0.0_IDP
    grrojb(:,0)=0.0_IDP
    grtojb(:,0)=0.0_IDP
    gttojb(:,0)=0.0_IDP

    read(25) (mmb(l),nnb(l),l=1,lbmax)
    read(25) (phip(j),qfar(j),curfar(j),ffar(j),pfar(j),j=1,mjeq)
    do j=1,mjeq
       read(25) (rmnb(j,l),zmnb(j,l),pmnsb(j,l),bmodb(j,l),sqgb(j,l),sqgib(j,l),l=1,lbmax)
    end do
    do j=1,mjeq
       read(25) (grrb(j,l),grtb(j,l),gttb(j,l),grrojb(j,l),grtojb(j,l),gttojb(j,l),jbgrrb(j,l),jbgrtb(j,l),jbgttb(j,l), &
            l=1,lbmax)
    end do

    if (myPE == 0) then
       open(unit=89,file="orb_data",status="unknown")
       do j=1,mjeq
          write(89,'(i3,4(2x,1pe15.8))') j,phip(j),qfar(j),curfar(j),ffar(j)
       end do
       close(unit=89)
    end if

    lb0=0
    do l=1,lbmax
       if (mmb(l) == 0 .and. nnb(l) == 0) exit
    end do
    lb0=l
    if (lb0 == 0 .or. lb0 > lbmax) then
       if (myPE == 0) write (6,'("  vmec: lb0=0")')
       stop
    end if

    if (lasym) then
       l1=lbmax
       do l=1,lbmax
          if (l == lb0) cycle
          l1=l1+1
          llc(l)=l1
          lls(l)=l1
       end do
       llc(lb0)=lb0
       lls(lb0)=0

       do j=1,mjeq
          read(25) (rmnb(j,lls(l)),zmnb(j,llc(l)),pmnsb(j,llc(l)),bmodb(j,lls(l)),sqgb(j,lls(l)),sqgib(j,lls(l)),l=1,lbmax)
       end do
       do j=1,mjeq
          read(25) (grrb(j,lls(l)),grtb(j,llc(l)),gttb(j,lls(l)),grrojb(j,lls(l)),grtojb(j,llc(l)),gttojb(j,lls(l)), &
               jbgrrb(j,lls(l)),jbgrtb(j,llc(l)),jbgttb(j,lls(l)),l=1,lbmax)
       end do
    end if

    close(unit=25)

    do l=1,lbmax
       if (nnb(l) == 0 .and. mmb(l) < 0) then
          mmb(l)=-mmb(l)
          zmnb(:,l)=-zmnb(:,l)
          pmnsb(:,l)=-pmnsb(:,l)
          grtb(:,l)=-grtb(:,l)
          grtojb(:,l)=-grtojb(:,l)
          jbgrtb(:,l)=-jbgrtb(:,l)
          if (lasym) then
             rmnb(:,lls(l))=-rmnb(:,lls(l))
             bmodb(:,lls(l))=-bmodb(:,lls(l))
             sqgb(:,lls(l))=-sqgb(:,lls(l))
             sqgib(:,lls(l))=-sqgib(:,lls(l))
             grrb(:,lls(l))=-grrb(:,lls(l))
             gttb(:,lls(l))=-gttb(:,lls(l))
             grrojb(:,lls(l))=-grrojb(:,lls(l))
             gttojb(:,lls(l))=-gttojb(:,lls(l))
             jbgrrb(:,lls(l))=-jbgrrb(:,lls(l))
             jbgttb(:,lls(l))=-jbgttb(:,lls(l))
          end if
       end if
       if (nnb(l) < 0) then
          mmb(l)=-mmb(l)
          nnb(l)=-nnb(l)
          zmnb(:,l)=-zmnb(:,l)
          pmnsb(:,l)=-pmnsb(:,l)
          grtb(:,l)=-grtb(:,l)
          grtojb(:,l)=-grtojb(:,l)
          jbgrtb(:,l)=-jbgrtb(:,l)
          if (lasym) then
             rmnb(:,lls(l))=-rmnb(:,lls(l))
             bmodb(:,lls(l))=-bmodb(:,lls(l))
             sqgb(:,lls(l))=-sqgb(:,lls(l))
             sqgib(:,lls(l))=-sqgib(:,lls(l))
             grrb(:,lls(l))=-grrb(:,lls(l))
             gttb(:,lls(l))=-gttb(:,lls(l))
             grrojb(:,lls(l))=-grrojb(:,lls(l))
             gttojb(:,lls(l))=-gttojb(:,lls(l))
             jbgrrb(:,lls(l))=-jbgrrb(:,lls(l))
             jbgttb(:,lls(l))=-jbgttb(:,lls(l))
          end if
       end if
    end do

    if (lasym) then
       l1=lbmax
       do l=1,lbmax
          if (l == lb0) cycle
          l1=l1+1
          mmb(l1)=-mmb(l)
          nnb(l1)=-nnb(l)
       end do
       lbmax=lbm2
    end if

    ! normalizations

    pi=4.0_IDP*atan(1.0_IDP)
    twopi=2.0_IDP*pi
    ffar=-ffar
    pfar=2.e-7*twopi*pfar

    bigrn=1.5*rmnb(1,lb0)-0.5*rmnb(2,lb0)
    bmodn=1.5*bmodb(1,lb0)-0.5*bmodb(2,lb0)

    twophip=2.0*phip(mjeq)
    curfar=sign(1.0_IDP,twophip)*curfar/(bigrn*bmodn)
    ffar=sign(1.0_IDP,twophip)*ffar/(bigrn*bmodn)
    bmodb=bmodb/bmodn
    sqgb=-bmodn*sqgb/bigrn
    sqgib=-bigrn*sqgib/bmodn
    grrb=bmodn*grrb/abs(twophip)
    grtb=bmodn*grtb/abs(twophip)
    gttb=bmodn*gttb/abs(twophip)
    grrojb=bigrn*grrojb/abs(twophip)
    grtojb=bigrn*grtojb/abs(twophip)
    gttojb=bigrn*gttojb/abs(twophip)
    jbgrrb=bmodn*bmodn*jbgrrb/(bigrn*abs(twophip))
    jbgrtb=bmodn*bmodn*jbgrtb/(bigrn*abs(twophip))
    jbgttb=bmodn*bmodn*jbgttb/(bigrn*abs(twophip))

    eps=sqrt(abs(twophip))/(bigrn*sqrt(bmodn))

    if (myPE == 0) write(6,'(/"R_0 =",1pe13.6," B0 =",1pe13.6," eps=",1pe13.6)') bigrn,bmodn,eps

    do j=0,mjeq
       rfar(j)=sqrt(1.0_IDP*j/mjeq)
    end do

    rbinv(0)=0.0_IDP
    do j=1,mjeq
       rbinv(j)=1.0_IDP/rfar(j)
    end do

    nnb=nnb/nfp

    mmaxb=maxval(mmb(1:lbmax))
    mminb=minval(mmb(1:lbmax))
    nmaxb=maxval(nnb(1:lbmax))
    nminb=minval(nnb(1:lbmax))
    mmaxxb=max(mmaxb,abs(mminb))

    allocate (llb(mminb:mmaxb,nminb:nmaxb))

    llb=0
    do l=1,lbmax
       llb(mmb(l),nnb(l))=l
    end do

    mband=0
    do m=0,mmaxb
       if (llb(m,0) /= 0) mband=mband+1
    end do
    mxmbandb=2*mband-1
    do n=1,nmaxb
       mband=0
       do m=mminb,mmaxb
          if (llb(m,n) /= 0) mband=mband+1
       end do
       if (mband > mxmbandb) mxmbandb=mband
    end do
    call mmblims

    nnb=nnb*nfp

    sqgib(0,:)=0.0_IDP
    sqgb(0,:)=0.0_IDP
    grrb(0,:)=0.0_IDP
    grtb(0,:)=0.0_IDP
    gttb(0,:)=0.0_IDP
    bmodb(0,:)=0.0_IDP
    grrojb(0,:)=0.0_IDP
    grtojb(0,:)=0.0_IDP
    gttojb(0,:)=0.0_IDP
    jbgrrb(0,:)=0.0_IDP
    jbgrtb(0,:)=0.0_IDP
    jbgttb(0,:)=0.0_IDP
    qfar(0)=1.5*qfar(1)-0.5*qfar(2)
    pfar(0)=1.5*pfar(1)-0.5*pfar(2)
    ffar(0)=1.5*ffar(1)-0.5*ffar(2)
    do j=1,mjeq
       curfar(j)=curfar(j)*mjeq/(j-0.5_IDP)
    end do
    curfar(0)=1.5*curfar(1)-0.5*curfar(2)

    do l=1,lbmax
       if (mmb(l) == 0) then
          sqgib(0,l)=2.*sqgib(1,l)-sqgib(2,l)
          sqgb(0,l)=1.5*sqgb(1,l)-0.5*sqgb(2,l)
          bmodb(0,l)=1.5*bmodb(1,l)-0.5*bmodb(2,l)
          rmnb(0,l)=1.5*rmnb(1,l)-0.5*rmnb(2,l)
          zmnb(0,l)=1.5*zmnb(1,l)-0.5*zmnb(2,l)
          pmnsb(0,l)=1.5*pmnsb(1,l)-0.5*pmnsb(2,l)
       endif
       if (mmb(l) == 0 .or. abs(mmb(l)) == 2) then
          grrb(0,l)=2.*grrb(1,l)-grrb(2,l)
          grtb(0,l)=2.*grtb(1,l)-grtb(2,l)
          gttb(0,l)=2.*gttb(1,l)-gttb(2,l)
          grrojb(0,l)=2.*grrojb(1,l)-grrojb(2,l)
          grtojb(0,l)=2.*grtojb(1,l)-grtojb(2,l)
          gttojb(0,l)=2.*gttojb(1,l)-gttojb(2,l)
          jbgrrb(0,l)=2.*jbgrrb(1,l)-jbgrrb(2,l)
          jbgrtb(0,l)=2.*jbgrtb(1,l)-jbgrtb(2,l)
          jbgttb(0,l)=2.*jbgttb(1,l)-jbgttb(2,l)
       endif
    end do

    bet0=2.*pfar(0)/bmodn**2
    if (bet0 /= 0.0_IDP) then
       pror=pfar(0)
       pfar=pfar/pror
    end if

    dummy=1.5*qfar(mjeq)-0.5*qfar(mjeqm1)
    do j=1,mjeqm1
       qfar(j)=0.5*(qfar(j+1)+qfar(j))
    end do
    qfar(mjeq)=dummy
    call eqsplns(qqinv,sd1,qfar,0,0,"spline",1.e-07_IDP)
    dummy=1.5*pfar(mjeq)-0.5*pfar(mjeqm1)
    do j=1,mjeqm1
       pfar(j)=0.5*(pfar(j+1)+pfar(j))
    end do
    pfar(mjeq)=dummy
    call eqsplns(preq,sd1,pfar,0,0,"spline",1.e-07_IDP)
    dummy=1.5*ffar(mjeq)-0.5*ffar(mjeqm1)
    do j=1,mjeqm1
       ffar(j)=0.5*(ffar(j+1)+ffar(j))
    end do
    ffar(mjeq)=dummy
    call eqsplns(feq,sd1,ffar,0,0,"spline",1.e-07_IDP)
    dummy=1.5*curfar(mjeq)-0.5*curfar(mjeqm1)
    do j=1,mjeqm1
       curfar(j)=0.5*(curfar(j+1)+curfar(j))
    end do
    curfar(mjeq)=dummy
    call eqsplns(cureq,sd1,curfar,0,0,"spline",1.e-07_IDP)
    do j=1,mj
       cureq(j)=cureq(j)*r(j)**2
    end do
    cureq(0)=0.0_IDP

    qq=1.0_IDP/qqinv

    allocate (R_eq(0:mj,0:lbmax),Z_eq(0:mj,0:lbmax),Phimns_eq(0:mj,0:lbmax),sqgeq(0:mj,0:lbmax),sqgieq(0:mj,0:lbmax), &
         bmodeq(0:mj,0:lbmax),bsteq(0:mj,0:lbmax),grreq(0:mj,0:lbmax),grteq(0:mj,0:lbmax),gtteq(0:mj,0:lbmax), &
         grrojeq(0:mj,0:lbmax),grtojeq(0:mj,0:lbmax),gttojeq(0:mj,0:lbmax), &
         jbgrreq(0:mj,0:lbmax),jbgrteq(0:mj,0:lbmax),jbgtteq(0:mj,0:lbmax), &
         lplreq(0:mj,0:lbmax),lplteq(0:mj,0:lbmax),djrojeq(0:mj,0:lbmax),djtojeq(0:mj,0:lbmax),djzojeq(0:mj,0:lbmax), &
         jsqeq(0:mj,0:lbmax),omdreq(0:mj,0:lbmax),omdteq(0:mj,0:lbmax),omdzeq(0:mj,0:lbmax), &
         dbsjtojeq(0:mj,0:lbmax),dbsjzojeq(0:mj,0:lbmax),dbsjtbjeq(0:mj,0:lbmax),dgttreq(0:mj,0:lbmax), &
         dgrrteq(0:mj,0:lbmax),dgrtteq(0:mj,0:lbmax),dgttteq(0:mj,0:lbmax),dgrrzeq(0:mj,0:lbmax), &
         dgrtzeq(0:mj,0:lbmax),dgttzeq(0:mj,0:lbmax),dgrtpeq(0:mj,0:lbmax),dgttpeq(0:mj,0:lbmax), &
         bsgrteq(0:mj,0:lbmax),bsgtteq(0:mj,0:lbmax),bsqeq(0:mj,0:lbmax),bsqgtteq(0:mj,0:lbmax), &
         lplrreq(0:mj,0:lbmax),lplrteq(0:mj,0:lbmax),lplrzeq(0:mj,0:lbmax),lpltreq(0:mj,0:lbmax),lpltteq(0:mj,0:lbmax), &
         lpltzeq(0:mj,0:lbmax),lplzreq(0:mj,0:lbmax),lplzteq(0:mj,0:lbmax),lplzzeq(0:mj,0:lbmax), &
         sb1(0:mj,0:lbmax),sb2(0:mj,0:lbmax),sb3(0:mj,0:lbmax),sb4(0:mj,0:lbmax),sb5(0:mj,0:lbmax),sb6(0:mj,0:lbmax))

    if (trapped_on == 1) allocate (omdrprpeq(0:mj,0:lbmax),omdtprpeq(0:mj,0:lbmax),omdzprpeq(0:mj,0:lbmax))

    sqgieq=0.0_IDP
    sqgeq=0.0_IDP
    bsteq=0.0_IDP
    grreq=0.0_IDP
    grteq=0.0_IDP
    gtteq=0.0_IDP
    bmodeq=0.0_IDP
    grrojeq=0.0_IDP
    grtojeq=0.0_IDP
    gttojeq=0.0_IDP
    jbgrreq=0.0_IDP
    jbgrteq=0.0_IDP
    jbgtteq=0.0_IDP

    ir=max(15,mjeq/10)
    do l=1,lbmax
       m=min(abs(mmb(l)),30)
       m2=m
       call eqsplns(sqgieq(:,l),sd1,sqgib(:,l),m,m2,"icsscu",1.e-07_IDP)
       dummy=1.5*sqgb(mjeq,l)-0.5*sqgb(mjeqm1,l)
       do j=1,mjeqm1
          sqgb(j,l)=0.5*(sqgb(j+1,l)+sqgb(j,l))
       end do
       sqgb(mjeq,l)=dummy
       call eqsplns(sqgeq(:,l),sd1,sqgb(:,l),m,m2,"icsscu",1.e-07_IDP)
       dummy=1.5*bmodb(mjeq,l)-0.5*bmodb(mjeqm1,l)
       do j=1,mjeqm1
          bmodb(j,l)=0.5*(bmodb(j+1,l)+bmodb(j,l))
       end do
       bmodb(mjeq,l)=dummy
       call eqsplns(bmodeq(:,l),sd1,bmodb(:,l),m,m2,"icsscu",1.e-07_IDP)
       dummy=1.5*rmnb(mjeq,l)-0.5*rmnb(mjeqm1,l)
       do j=1,mjeqm1
          rmnb(j,l)=0.5*(rmnb(j+1,l)+rmnb(j,l))
       end do
       rmnb(mjeq,l)=dummy
       call eqsplns(R_eq(:,l),sd1,rmnb(:,l),m,m2,"icsscu",1.e-07_IDP)
       dummy=1.5*zmnb(mjeq,l)-0.5*zmnb(mjeqm1,l)
       do j=1,mjeqm1
          zmnb(j,l)=0.5*(zmnb(j+1,l)+zmnb(j,l))
       end do
       zmnb(mjeq,l)=dummy
       call eqsplns(Z_eq(:,l),sd1,zmnb(:,l),m,m2,"icsscu",1.e-07_IDP)
       dummy=1.5*pmnsb(mjeq,l)-0.5*pmnsb(mjeqm1,l)
       do j=1,mjeqm1
          pmnsb(j,l)=0.5*(pmnsb(j+1,l)+pmnsb(j,l))
       end do
       pmnsb(mjeq,l)=dummy
       call eqsplns(Phimns_eq(:,l),sd1,pmnsb(:,l),m,m2,"icsscu",1.e-07_IDP)

       mm1=abs(m-1)
       mm2=abs(mm1-1)
       mm12=mm1
       mm22=mm2
       mm2p1=mm2+1
       mm2p12=mm2p1
       mp2=m+2
       mp22=mp2
       isg=sign(1,mmb(l))

       sfar1=grrb(:,l)+isg*grtb(:,l)
       sfar2=grrb(:,l)-isg*grtb(:,l)
       sfar3=gttb(:,l)-isg*grtb(:,l)
       sfar4=sfar3-sfar1
       call fitter(sd1,sd6,sfar1,ir,m)
       call fitter(sd2,sd6,sfar2,ir,mm2)
       call fitter(sd4,sd6,sfar4,ir,mp2)
       grreq(:,l)=(sd1+sd2)/2.
       grteq(:,l)=isg*(sd1-sd2)/2.
       gtteq(:,l)=(2.*sd4+3.*sd1-sd2)/2.

       sfar1=grrojb(:,l)+isg*grtojb(:,l)
       sfar2=grrojb(:,l)-isg*grtojb(:,l)
       sfar3=gttojb(:,l)-isg*grtojb(:,l)
       sfar4=sfar3-sfar1
       call fitter(sd1,sd6,sfar1,ir,m)
       call fitter(sd2,sd6,sfar2,ir,mm2)
       call fitter(sd4,sd6,sfar4,ir,mp2)
       grrojeq(:,l)=(sd1+sd2)/2.
       grtojeq(:,l)=isg*(sd1-sd2)/2.
       gttojeq(:,l)=(2.*sd4+3.*sd1-sd2)/2.

       sfar1=jbgrrb(:,l)+isg*jbgrtb(:,l)
       sfar2=jbgrrb(:,l)-isg*jbgrtb(:,l)
       sfar3=jbgttb(:,l)-isg*jbgrtb(:,l)
       sfar4=sfar3-sfar1
       call fitter(sd1,sd6,sfar1,ir,m)
       call fitter(sd2,sd6,sfar2,ir,mm2)
       call fitter(sd4,sd6,sfar4,ir,mp2)
       jbgrreq(:,l)=(sd1+sd2)/2.
       jbgrteq(:,l)=isg*(sd1-sd2)/2.
       jbgtteq(:,l)=(2.*sd4+3.*sd1-sd2)/2.

    end do
!          Write interpolated R_eq, Z_eq into geom_bzr_to_cyl file  9-19-2022
    !            Include Phimns_eq                                      8-11-2023
    if (myPE == 0) then
       open(unit=30,file='geom_bzr_to_cyl.dat',status='unknown')
       write(30,*) mj, lbmax, lasym
       do l=1,lbmax
          write(30,*) mmb(l),nnb(l)
       end do
       do j=0,mj
          write(30,*) real(j)/real(mj),1./qqinv(j)
       end do
       do l=1,lbmax
          do j=0,mj
             !	       if(lasym) then                  
             !		  write(30,'(e15.8,5(2x,e15.8))') rmncbf(l,j), zmnsbf(l,j), pmnsbf(l,j),  &
             !		           rmnsbf(l,j), zmncbf(l,j), pmncbf(l,j)
             !	       else
             !     	  	         write(30,'(e15.8,2(2x,e15.8))') R_eq(j,l), Z_eq(j,l), Phimns_eq(j,l)
             !        Option to also include the Bmns:
             write(30,'(e15.8,3(2x,e15.8))') R_eq(j,l), Z_eq(j,l),  &
                  Phimns_eq(j,l), bmodeq(j,l)
             !	       endif
          end do
       end do
       close(unit=30)
    end if
      
    ! Edge extrapolation

    if (Edge_on .eq. 1) then
       do j=edge_p,mj
          sqgieq(j,:)=sqgieq(j-2,:) + (r(j)-r(j-2))*(sqgieq(j-1,:)-sqgieq(j-2,:))/(r(j-1)-r(j-2))
          sqgeq(j,:)=sqgeq(j-2,:) + (r(j)-r(j-2))*(sqgeq(j-1,:)-sqgeq(j-2,:))/(r(j-1)-r(j-2))
          grreq(j,:)=grreq(j-2,:) + (r(j)-r(j-2))*(grreq(j-1,:)-grreq(j-2,:))/(r(j-1)-r(j-2))
          grteq(j,:)=grteq(j-2,:) + (r(j)-r(j-2))*(grteq(j-1,:)-grteq(j-2,:))/(r(j-1)-r(j-2))
          gtteq(j,:)=gtteq(j-2,:) + (r(j)-r(j-2))*(gtteq(j-1,:)-gtteq(j-2,:))/(r(j-1)-r(j-2))
          bmodeq(j,:)=bmodeq(j-2,:) + (r(j)-r(j-2))*(bmodeq(j-1,:)-bmodeq(j-2,:))/(r(j-1)-r(j-2))
          grrojeq(j,:)=grrojeq(j-2,:) + (r(j)-r(j-2))*(grrojeq(j-1,:)-grrojeq(j-2,:))/(r(j-1)-r(j-2))
          grtojeq(j,:)=grtojeq(j-2,:) + (r(j)-r(j-2))*(grtojeq(j-1,:)-grtojeq(j-2,:))/(r(j-1)-r(j-2))
          gttojeq(j,:)=gttojeq(j-2,:) + (r(j)-r(j-2))*(gttojeq(j-1,:)-gttojeq(j-2,:))/(r(j-1)-r(j-2))
          jbgrreq(j,:)=jbgrreq(j-2,:) + (r(j)-r(j-2))*(jbgrreq(j-1,:)-jbgrreq(j-2,:))/(r(j-1)-r(j-2))
          jbgrteq(j,:)=jbgrteq(j-2,:) + (r(j)-r(j-2))*(jbgrteq(j-1,:)-jbgrteq(j-2,:))/(r(j-1)-r(j-2))
          jbgtteq(j,:)=jbgtteq(j-2,:) + (r(j)-r(j-2))*(jbgtteq(j-1,:)-jbgtteq(j-2,:))/(r(j-1)-r(j-2))

          preq(j)=preq(j-2) + (r(j)-r(j-2))*(preq(j-1)-preq(j-2))/(r(j-1)-r(j-2))
          feq(j)=feq(j-2) + (r(j)-r(j-2))*(feq(j-1)-feq(j-2))/(r(j-1)-r(j-2))
          cureq(j)=cureq(j-2) + (r(j)-r(j-2))*(cureq(j-1)-cureq(j-2))/(r(j-1)-r(j-2))
       end do
    end if

    if (ext_prof == 1) call ae_profiles

    ! bst is modified for resonant modes

    allocate (lbst(lbmax))
    allocate (rsb(lbmax), Abst(lbmax), Wbst(lbmax))

    bsteq=0
    one = 1.0_IDP
    lbst=0
    Abst=0
    qmin=minval(qqinv(0:mj))
    qmax=maxval(qqinv(0:mj))
    call dbydr0(sd1,preq,0.0_IDP,-1.0_IDP,0)
    call dbydr0(sd3,qqinv,0.0_IDP,1.0_IDP,0)
    lp=0
    do l=1,lbmax
       if (mmb(l) == 0 .or. nnb(l) == 0) cycle
       p=one*nnb(l)/mmb(l)
       if ((p-qmin)*(p-qmax) >= 0.0_IDP) cycle
       lp=lp+1
       lbst(lp)=l
       do j=1,mj
          if ((qqinv(j-1)-p)*(qqinv(j)-p) <= 0.) exit
       end do
       jp=j
       rsb(lp)=r(jp-1)+(r(jp)-r(jp-1))*(qqinv(jp-1)-p)/(qqinv(jp-1)-qqinv(jp))
       ! write(6,'("m =",i5," n =",i5," rs =",1pe13.6," r(",i3,") =",1pe13.6," r(",i3,") =",1pe13.6)') mmb(l), &
       !    nnb(l),rsb(lp),jp-1,r(jp-1),jp,r(jp)
       if (jp == 1) then
          Abst(lp) = (preq(0)-preq(1))/r(1)
          iotap = (qqinv(1)-qqinv(0))/r(1)
       else if (jp < mj) then
          Abst(lp) = ((r(jp)-rsb(lp))*sd1(jp-1)+(rsb(lp)-r(jp-1))*sd1(jp))/(r(jp)-r(jp-1))
          iotap = ((r(jp)-rsb(lp))*sd3(jp-1)+(rsb(lp)-r(jp-1))*sd3(jp))/(r(jp)-r(jp-1))
       else
          Abst(lp) = (preq(mjm1)-preq(mj))/(r(mj)-r(mjm1))
          iotap = (qqinv(mj)-qqinv(mjm1))/(r(mj)-r(mjm1))
       end if
       Wbst(lp) = abs(p/(mmb(l)*iotap))
    end do
    nbst=lp
    if (myPE == 0) write(6,'(/"bst is modified for",i5," resonant modes")') nbst
    do l=1,nbst
       if (myPE == 0) write(6,'("m =",i5," n =",i5," rs =",1pe13.6," A =",1pe13.6," W =",1pe13.6)') &
            mmb(lbst(l)),nnb(lbst(l)),rsb(l),Abst(l),Wbst(l)
    end do

    sd2=rinv*sd1/2.0_IDP
    sd2(0)=(preq(0)-preq(1))/r(1)**2
    do l=1,nbst
       do j=0,mj
          pprime=sd2(j)
          iota=qqinv(j)
          gc=sqgeq(j,lbst(l))
          xl=r(j)-rsb(l)
          pprimel=pprime-0.5*rinv(j)*Abst(l)*(1-(xl/Wbst(l))**2)*exp(-0.5*(xl/Wbst(l))**2)
          bsteq(j,lbst(l))=bet0*pprimel*gc/(mmb(lbst(l))*iota-nnb(lbst(l)))
       end do
    end do

    do l=1,lbmax
       if(mmb(l) == 0 .and. nnb(l) == 0) cycle
       p=0.0
       if (mmb(l) /= 0) p=one*nnb(l)/mmb(l)
       if ((p-qmin)*(p-qmax) < 0.0_IDP) cycle
       bsteq(:,l)=bet0*sd2*sqgeq(:,l)/(mmb(l)*qqinv-nnb(l))
    end do

    sb1=sqgeq
    call multb(jsqeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
    call dbydrb(sb1,sqgieq,0.0_IDP,1.0_IDP,0)
    call multb(lplreq,gtteq,1,sb1,1,0.0_IDP,-1.0_IDP)
    call dbydthb(sb2,sqgieq,1,0.0_IDP,1.0_IDP,0)
    call multb(lplreq,grteq,-1,sb2,-1,1.0_IDP,1.0_IDP)
    call multb(lplteq,grteq,-1,sb1,1,0.0_IDP,1.0_IDP)
    call multb(lplteq,grreq,1,sb2,-1,1.0_IDP,-1.0_IDP)

    call dbydrb(sb1,sqgeq,0.0_IDP,1.0_IDP,0)
    call multb(djrojeq,sqgieq,1,sb1,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb1,sqgeq,1,0.0_IDP,1.0_IDP,0)
    call multb(djtojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)
    call dbydztb(sb1,sqgeq,1,0.0_IDP,1.0_IDP)
    call multb(djzojeq,sqgieq,1,sb1,-1,0.0_IDP,1.0_IDP)

    ! WARNING: sd1 is feq-qqinv*cureq

    sd1=feq-qqinv*cureq
    call dbydztb(sb1,sqgieq,1,0.0_IDP,1.0_IDP)
    call multb(omdteq,bsteq,-1,sb1,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb2(:,l)=rinv*cureq*sb1(:,l)
    end do
    call dbydthb(sb1,sqgieq,1,0.0_IDP,1.0_IDP,0)
    call multb(omdzeq,bsteq,-1,sb1,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb2(:,l)=sb2(:,l)-feq*sb1(:,l)
    end do
    call multb(omdreq,sqgeq,1,sb2,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       omdreq(:,l)=omdreq(:,l)/(2.0_IDP*sd1)
    end do
    do l=1,lbmax
       sb1(:,l)=sd1*sqgieq(:,l)
    end do
    call dbydrb(sb2,sb1,0.0_IDP,1.0_IDP,0)
    do l=1,lbmax
       sb1(:,l)=feq*sb2(:,l)-sd1*r*omdteq(:,l)
    end do
    call multb(omdteq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       omdteq(:,l)=omdteq(:,l)/(2.0_IDP*sd1*sd1)
    end do
    do l=1,lbmax
       sb1(:,l)=sd1*r*omdzeq(:,l)-rinv*cureq*sb2(:,l)
    end do
    call multb(omdzeq,sqgeq,1,sb1,1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       omdzeq(:,l)=omdzeq(:,l)/(2.0_IDP*sd1*sd1)
    end do

    if (trapped_on .eq. 1) then
       do l=1,lbmax
          if (mm(l) .eq. 0)  then
             sb1(:,l)=4.*pi*pi*r*omcyb*sqgeq(:,l)/(rbound)
          else
             sb1(:,l)=4.*pi*pi*r*omcyb*sqgeq(:,l)/(rbound*mm(l))
          end if
       end do
       call multb(omdrprpeq,omdreq,-1,sb1,1,0.0_IDP,1.0_IDP)
       call multb(omdtprpeq,omdteq,1,sb1,1,0.0_IDP,1.0_IDP)
       call multb(omdzprpeq,omdzeq,1,sb1,1,0.0_IDP,1.0_IDP)
    end if

    call multb(sb1,bsteq,-1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb2,sb1,-1,0.0_IDP,1.0_IDP,0)
    call multb(dbsjtojeq,sb2,1,sqgieq,1,0.0_IDP,1.0_IDP)
    call multb(dbsjtbjeq,sb2,1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydztb(sb2,sb1,-1,0.0_IDP,1.0_IDP)
    call multb(dbsjzojeq,sb2,1,sqgieq,1,0.0_IDP,1.0_IDP)

    call dbydrb(sb1,gtteq,0.0_IDP,1.0_IDP,2)
    call multb(dgttreq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb1,grreq,1,0.0_IDP,1.0_IDP,2)
    call multb(dgrrteq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb1,grteq,-1,0.0_IDP,1.0_IDP,2)
    call multb(dgrtteq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb1,gtteq,1,0.0_IDP,1.0_IDP,2)
    call multb(dgttteq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydztb(sb1,grreq,1,0.0_IDP,1.0_IDP)
    call multb(dgrrzeq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydztb(sb1,grteq,-1,0.0_IDP,1.0_IDP)
    call multb(dgrtzeq,sb1,1,sqgeq,1,0.0_IDP,1.0_IDP)
    call dbydztb(sb1,gtteq,1,0.0_IDP,1.0_IDP)
    call multb(dgttzeq,sb1,-1,sqgeq,1,0.0_IDP,1.0_IDP)
    call grparb(sb1,grteq,-1,0.0_IDP,1.0_IDP)
    call multb(dgrtpeq,sb1,1,sqgieq,1,0.0_IDP,1.0_IDP)
    call grparb(sb1,gtteq,1,0.0_IDP,1.0_IDP)
    call multb(dgttpeq,sb1,-1,sqgieq,1,0.0_IDP,1.0_IDP)

    do l=1,lbmax
       dbsjtojeq(:,l)=r*dbsjtojeq(:,l)/sd1
       dbsjzojeq(:,l)=dbsjzojeq(:,l)/sd1
       dbsjtbjeq(:,l)=r*dbsjtbjeq(:,l)/sd1
       dgttreq(:,l)=dgttreq(:,l)/sd1
       dgrrteq(:,l)=r*dgrrteq(:,l)/sd1
       dgrtteq(:,l)=r*dgrtteq(:,l)/sd1
       dgttteq(:,l)=r*dgttteq(:,l)/sd1
       dgrrzeq(:,l)=dgrrzeq(:,l)/sd1
       dgrtzeq(:,l)=dgrtzeq(:,l)/sd1
       dgttzeq(:,l)=dgttzeq(:,l)/sd1
       dgrtpeq(:,l)=dgrtpeq(:,l)/sd1
       dgttpeq(:,l)=dgttpeq(:,l)/sd1
    end do

    sb1=bsteq
    call multb(bsqeq,sb1,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    call multb(bsgrteq,grtojeq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    call multb(bsgtteq,gttojeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
    call multb(bsqgtteq,bsgtteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    call dbydrb(sb1,sqgieq,0.0_IDP,1.0_IDP,0)
    do l=1,lbmax
       sb3(:,l)=sd1*gtteq(:,l)-rinv*rinv*cureq*cureq*sqgeq(:,l)/(eps*eps)
    end do
    call multb(lplrreq,sb3,1,sb1,1,0.0_IDP,-1.0_IDP)
    call multb(sb4,gtteq,1,bsteq,-1,0.0_IDP,1.0_IDP)
    call multb(sb5,sqgeq,1,bsteq,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=feq*grteq(:,l)-r*r*qqinv*sb4(:,l)-cureq*sb5(:,l)/(eps*eps)
    end do
    call multb(lplrteq,sb3,-1,sb1,1,0.0_IDP,1.0_IDP)
    call dbydthb(sb2,sqgieq,1,0.0_IDP,1.0_IDP,0)
    call multb(lpltreq,sb3,-1,sb2,-1,0.0_IDP,1.0_IDP)
    call multb(sb6,sb5,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    call multb(sb5,grteq,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=feq*(feq*grreq(:,l)-2*r*r*qqinv*sb5(:,l))/sd1-r*r*sb6(:,l)/(eps*eps)
    end do
    call multb(sb6,sb4,-1,bsteq,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=sb3(:,l)+r*r*r*r*qqinv*qqinv*sb6(:,l)/sd1
    end do
    call multb(lpltteq,sb3,1,sb2,-1,0.0_IDP,-1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=rinv*cureq*grteq(:,l)-r*sb4(:,l)
    end do
    call multb(lplrzeq,sb3,-1,sb1,1,0.0_IDP,-1.0_IDP)
    call dbydztb(sb1,sqgieq,1,0.0_IDP,1.0_IDP)
    call multb(lplzreq,sb3,-1,sb1,-1,0.0_IDP,-1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=(feq*(rinv*cureq*grreq(:,l)-r*sb5(:,l))-r*qqinv*(cureq*sb5(:,l)-r*r*sb6(:,l)))/sd1
    end do
    call multb(lpltzeq,sb3,1,sb2,-1,0.0_IDP,1.0_IDP)
    call multb(lplzteq,sb3,1,sb1,-1,0.0_IDP,1.0_IDP)
    do l=1,lbmax
       sb3(:,l)=(rinv*cureq*(rinv*cureq*grreq(:,l)-r*sb5(:,l))-(cureq*sb5(:,l)-r*r*sb6(:,l)))/sd1
    end do
    call multb(lplzzeq,sb3,1,sb1,-1,0.0_IDP,-1.0_IDP)

    sqgi=0.0_IDP
    sqg=0.0_IDP
    bst=0.0_IDP
    grr=0.0_IDP
    grt=0.0_IDP
    gtt=0.0_IDP
    bmod=0.0_IDP
    grroj=0.0_IDP
    grtoj=0.0_IDP
    gttoj=0.0_IDP
    jbgrr=0.0_IDP
    jbgrt=0.0_IDP
    jbgtt=0.0_IDP
    lplr=0.0_IDP
    lplt=0.0_IDP
    djroj=0.0_IDP
    djtoj=0.0_IDP
    djzoj=0.0_IDP
    omdr=0.0_IDP
    omdt=0.0_IDP
    omdz=0.0_IDP
    dbsjtoj=0.0_IDP
    dbsjzoj=0.0_IDP
    dbsjtbj=0.0_IDP
    dgttr=0.0_IDP
    dgrrt=0.0_IDP
    dgrtt=0.0_IDP
    dgttt=0.0_IDP
    dgrrz=0.0_IDP
    dgrtz=0.0_IDP
    dgttz=0.0_IDP
    dgrtp=0.0_IDP
    dgttp=0.0_IDP
    jsq=0.0_IDP
    bsq=0.0_IDP
    bsgrt=0.0_IDP
    bsgtt=0.0_IDP
    bsqgtt=0.0_IDP
    lplrr=0.0_IDP
    lplrt=0.0_IDP
    lplrz=0.0_IDP
    lpltr=0.0_IDP
    lpltt=0.0_IDP
    lpltz=0.0_IDP
    lplzr=0.0_IDP
    lplzt=0.0_IDP
    lplzz=0.0_IDP

    if (trapped_on .eq. 1) then
       omdrprp=0.0_IDP
       omdtprp=0.0_IDP
       omdzprp=0.0_IDP
    end if

    do lq=1,leqmax
       do l=1,lbmax
          if(mmb(l).eq.mmeq(lq).and.nnb(l).eq.nneq(lq)) exit
       end do
       if (l > lbmax) cycle
       sqgi(:,lq)=sqgieq(:,l)
       sqg(:,lq)=sqgeq(:,l)
       bst(:,lq)=bsteq(:,l)
       grr(:,lq)=grreq(:,l)
       grt(:,lq)=grteq(:,l)
       gtt(:,lq)=gtteq(:,l)
       bmod(:,lq)=bmodeq(:,l)
       grroj(:,lq)=grrojeq(:,l)
       grtoj(:,lq)=grtojeq(:,l)
       gttoj(:,lq)=gttojeq(:,l)
       jbgrr(:,lq)=jbgrreq(:,l)
       jbgrt(:,lq)=jbgrteq(:,l)
       jbgtt(:,lq)=jbgtteq(:,l)
       lplr(:,lq)=lplreq(:,l)
       lplt(:,lq)=lplteq(:,l)
       djroj(:,lq)=djrojeq(:,l)
       djtoj(:,lq)=djtojeq(:,l)
       djzoj(:,lq)=djzojeq(:,l)
       omdr(:,lq)=omdreq(:,l)
       omdt(:,lq)=omdteq(:,l)
       omdz(:,lq)=omdzeq(:,l)
       dbsjtoj(:,lq)=dbsjtojeq(:,l)
       dbsjzoj(:,lq)=dbsjzojeq(:,l)
       dbsjtbj(:,lq)=dbsjtbjeq(:,l)
       dgttr(:,lq)=dgttreq(:,l)
       dgrrt(:,lq)=dgrrteq(:,l)
       dgrtt(:,lq)=dgrtteq(:,l)
       dgttt(:,lq)=dgttteq(:,l)
       dgrrz(:,lq)=dgrrzeq(:,l)
       dgrtz(:,lq)=dgrtzeq(:,l)
       dgttz(:,lq)=dgttzeq(:,l)
       dgrtp(:,lq)=dgrtpeq(:,l)
       dgttp(:,lq)=dgttpeq(:,l)
       jsq(:,lq)=jsqeq(:,l)
       bsq(:,lq)=bsqeq(:,l)
       bsgrt(:,lq)=bsgrteq(:,l)
       bsgtt(:,lq)=bsgtteq(:,l)
       bsqgtt(:,lq)=bsqgtteq(:,l)
       lplrr(:,lq)=lplrreq(:,l)
       lplrt(:,lq)=lplrteq(:,l)
       lplrz(:,lq)=lplrzeq(:,l)
       lpltr(:,lq)=lpltreq(:,l)
       lpltt(:,lq)=lpltteq(:,l)
       lpltz(:,lq)=lpltzeq(:,l)
       lplzr(:,lq)=lplzreq(:,l)
       lplzt(:,lq)=lplzteq(:,l)
       lplzz(:,lq)=lplzzeq(:,l)

       if (trapped_on .eq. 1) then
          omdrprp(:,lq)=omdrprpeq(:,l)
          omdtprp(:,lq)=omdtprpeq(:,l)
          omdzprp(:,lq)=omdzprpeq(:,l)
       end if

    end do

    ! lplrr=gttoj
    ! call dbydreq(lplr,lplrr,1.0_IDP,1.0_IDP,2)
    ! do l=1,leqmax
    !    lplr(:,l)=lplr(:,l)+rinv*lplrr(:,l)
    ! end do

    ! lplrt=-grtoj
    ! call dbydtheq(lplr,lplrt,-1,1.0_IDP,1.0_IDP,2)
    ! call dbydreq(lplt,lplrt,1.0_IDP,1.0_IDP,2)
    ! lplrt=2.0*lplrt

    ! lpltt=grroj
    ! call dbydtheq(lplt,lpltt,1,1.0_IDP,1.0_IDP,2)

    if (ieldamp_on == 1) then

       allocate (eildreq(0:mj,0:lbmax),eildteq(0:mj,0:lbmax),eildzeq(0:mj,0:lbmax),eildrreq(0:mj,0:lbmax), &
            eildrteq(0:mj,0:lbmax),eildrzeq(0:mj,0:lbmax),eildtteq(0:mj,0:lbmax),eildtzeq(0:mj,0:lbmax), &
            eildzzeq(0:mj,0:lbmax))

       call dbydrb(sb4,omdreq,0.0_IDP,1.0_IDP,3)
       call multb(eildreq,omdreq,-1,sb4,-1,0.0_IDP,1.0_IDP)
       call dbydthb(sb4,omdreq,-1,0.0_IDP,1.0_IDP,3)
       call multb(eildreq,omdteq,1,sb4,1,1.0_IDP,1.0_IDP)
       call dbydztb(sb4,omdreq,-1,0.0_IDP,1.0_IDP)
       call multb(eildreq,omdzeq,1,sb4,1,1.0_IDP,1.0_IDP)

       call dbydrb(sb4,omdteq,0.0_IDP,1.0_IDP,3)
       do l=1,lbmax
          sb4(:,l)=sb4(:,l)-rinv*omdteq(:,l)
          if (mmb(l) == 0 .or. abs(mmb(l)) == 2) sb4(0,l)=sb4(0,l)-rinv(1)*omdteq(1,l)
       end do
       call multb(eildteq,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
       call dbydthb(sb4,omdteq,1,0.0_IDP,1.0_IDP,3)
       call multb(eildteq,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
       call dbydztb(sb4,omdteq,1,0.0_IDP,1.0_IDP)
       call multb(eildteq,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

       call dbydrb(sb4,omdzeq,0.0_IDP,1.0_IDP,0)
       call multb(eildzeq,omdreq,-1,sb4,1,0.0_IDP,1.0_IDP)
       call dbydthb(sb4,omdzeq,1,0.0_IDP,1.0_IDP,0)
       call multb(eildzeq,omdteq,1,sb4,-1,1.0_IDP,1.0_IDP)
       call dbydztb(sb4,omdzeq,1,0.0_IDP,1.0_IDP)
       call multb(eildzeq,omdzeq,1,sb4,-1,1.0_IDP,1.0_IDP)

       sb1=omdreq
       call multb(eildrreq,sb1,-1,omdreq,-1,0.0_IDP,1.0_IDP)
       call multb(eildrteq,sb1,-1,omdteq,1,0.0_IDP,2.0_IDP)
       call multb(eildrzeq,sb1,-1,omdzeq,1,0.0_IDP,2.0_IDP)
       sb1=omdteq
       call multb(eildtteq,sb1,1,omdteq,1,0.0_IDP,1.0_IDP)
       call multb(eildtzeq,sb1,1,omdzeq,1,0.0_IDP,2.0_IDP)
       sb1=omdzeq
       call multb(eildzzeq,sb1,1,omdzeq,1,0.0_IDP,1.0_IDP)

       eildr=0.0_IDP
       eildt=0.0_IDP
       eildz=0.0_IDP
       eildrr=0.0_IDP
       eildrt=0.0_IDP
       eildrz=0.0_IDP
       eildtt=0.0_IDP
       eildtz=0.0_IDP
       eildzz=0.0_IDP

       do lq=1,leqmax
          do l=1,lbmax
             if(mmb(l).eq.mmeq(lq).and.nnb(l).eq.nneq(lq)) exit
          end do
          if (l > lbmax) cycle
          eildr(:,lq)=sd1*preq*eildreq(:,l)
          eildt(:,lq)=sd1*preq*eildteq(:,l)
          eildz(:,lq)=sd1*preq*eildzeq(:,l)
          eildrr(:,lq)=sd1*preq*eildrreq(:,l)
          eildrt(:,lq)=sd1*preq*eildrteq(:,l)
          eildrz(:,lq)=sd1*preq*eildrzeq(:,l)
          eildtt(:,lq)=sd1*preq*eildtteq(:,l)
          eildtz(:,lq)=sd1*preq*eildtzeq(:,l)
          eildzz(:,lq)=sd1*preq*eildzzeq(:,l)
       end do

       deallocate (eildreq,eildteq,eildzeq,eildrreq,eildrteq,eildrzeq,eildtteq,eildtzeq,eildzzeq)

    end if

    ! equilibrium arrays for delperpsq

    lplr=lplrr+lpltr+lplzr
    lplt=lplrt+lpltt+lplzt
    lplz=lplrz+lpltz+lplzz

    do l=1,leqmax
       lplrr(:,l)=sd1*gttoj(:,l)
    end do
    lplrr(:,leq0)=lplrr(:,leq0)-rinv*rinv*cureq*cureq/(eps*eps)
    call dbydreq(lplr,lplrr,1.0_IDP,1.0_IDP,2)
    do l=1,leqmax
       lplr(:,l)=lplr(:,l)+rinv*lplrr(:,l)
    end do

    do l=1,leqmax
       lplrt(:,l)=-feq*grtoj(:,l)+r*r*qqinv*bsgtt(:,l)+cureq*bst(:,l)/(eps*eps)
    end do
    call dbydtheq(lplr,lplrt,-1,1.0_IDP,1.0_IDP,2)
    call dbydreq(lplt,lplrt,1.0_IDP,1.0_IDP,2)
    lplrt=2.0*lplrt

    do l=1,leqmax
       lplrz(:,l)=rinv*cureq*grtoj(:,l)-r*bsgtt(:,l)
    end do
    call dbydzteq(lplr,lplrz,-1,1.0_IDP,1.0_IDP)
    call dbydreq(lplz,lplrz,1.0_IDP,1.0_IDP,3)
    do l=1,leqmax
       lplz(:,l)=lplz(:,l)+rinv*lplrz(:,l)
    end do
    lplrz=2.0*lplrz

    do l=1,leqmax
       lpltt(:,l)=(feq*feq*grroj(:,l)-r*r*qqinv*(2*feq*bsgrt(:,l)-r*r*qqinv*bsqgtt(:,l)))/sd1-r*r*bsq(:,l)/(eps*eps)
    end do
    call dbydtheq(lplt,lpltt,1,1.0_IDP,1.0_IDP,2)

    do l=1,leqmax
       lpltz(:,l)=-(rinv*feq*cureq*grroj(:,l)-r*((feq+qqinv*cureq)*bsgrt(:,l)-r*r*qqinv*bsqgtt(:,l)))/sd1
    end do
    call dbydzteq(lplt,lpltz,1,1.0_IDP,1.0_IDP)
    call dbydtheq(lplz,lpltz,1,1.0_IDP,1.0_IDP,3)
    lpltz=2.0*lpltz

    do l=1,leqmax
       lplzz(:,l)=(cureq*(rinv*rinv*cureq*grroj(:,l)-2*bsgrt(:,l))+r*r*bsqgtt(:,l))/sd1
    end do
    call dbydzteq(lplz,lplzz,1,1.0_IDP,1.0_IDP)

    ! END WARNING

    if (myPE == 0) write(0,'(" ====> Equilibria set up DONE !! ")')

    deallocate (llb,rbinv,qfar,pfar,phip,curfar,ffar,sfar1,sfar2,sfar3,sfar4,lbst,rsb,Abst,Wbst)
    if (lasym) deallocate (llc,lls)
    deallocate (rmnb,sqgb,sqgib,bmodb,grrb,grtb,gttb,grrojb,grtojb,gttojb,jbgrrb,jbgrtb,jbgttb)
    deallocate (sqgeq,sqgieq,bmodeq,bsteq,grreq,grteq,gtteq,grrojeq,grtojeq,gttojeq,jbgrreq,jbgrteq,jbgtteq, &
         lplreq,lplteq,djrojeq,djtojeq,djzojeq,jsqeq,omdreq,omdteq,omdzeq,dbsjtojeq,dbsjzojeq,dbsjtbjeq,dgttreq, &
         dgrrteq,dgrtteq,dgttteq,dgrrzeq,dgrtzeq,dgttzeq,dgrtpeq,dgttpeq,bsgrteq,bsgtteq,bsqeq,bsqgtteq, &
         lplrreq,lplrteq,lplrzeq,lpltreq,lpltteq,lpltzeq,lplzreq,lplzteq,lplzzeq,sb1,sb2,sb3,sb4,sb5,sb6)

  end subroutine vmec

  subroutine ae_profiles

    ! Set up profiles using an external source

    use mpi
    use processor
    use tools

    implicit none

    integer :: j,ic,i,l,ns0,nunit,je
    real(IDP) :: rsq,xpi,B0_e,R0_e,a_e,uion_e,kappa_e,delt_e,beta0_e,rmax_e,xx,etactte,qe,me,epsil,kblotz, &
         B_inboard,B_outboard,Rin,Rout,btor,bpol,bnrm,tecnt,ticnt,dnecnt,dnicnt,bt0,rmajr,rminr,mu0,va0,omgcya, &
         betalf,bet00,vf0,vtor0,dnnbinn_max,dnenn_max,dninn_max,tinn_max,tenn_max,tbnnn_max,pthermalnn_max, &
         ptotnn_max,dnalphann_max,talphann_max
    real(IDP), dimension(500) :: rho_e,rhosq_e,qprof,den_beam_e,den_ion_e,den_elec_e,den_imp_e,temp_beam_e,temp_ion_e, &
         temp_elec_e,pres_beam_e,pres_thermal_e,pres_equil_e,zeff_e,pol_rot_vel_e,tor_rot_vel_e,tor_rot_freq_e, &
         den_alpha_e,temp_alpha_e
    real(IDP), dimension(:), allocatable :: bspl,cspl,dspl
    character*1 :: cdum,cdum2

    !  Option to read parameters and profiles directly from a data file
    !  assumed to be ext_prof_name. The profiles are interpolated. We use SI units.
    !  Only one energetic particle species

    xpi = 4.0_IDP*atan(1.0_IDP)

    nunit = 17
    open(unit=nunit,file=ext_prof_name,status="old",form="formatted")

    read(nunit,'(a1)') cdum
    read(nunit,'(a1)') cdum
    read(nunit,*) B0_e
    read(nunit,'(a1)') cdum
    read(nunit,*) R0_e
    read(nunit,'(a1)') cdum
    read(nunit,*) a_e
    read(nunit,'(a1)') cdum
    read(nunit,*) kappa_e
    read(nunit,'(a1)') cdum
    read(nunit,*) delt_e
    read(nunit,'(a1)') cdum
    read(nunit,'(a1)') cdum
    read(nunit,'(a1)') cdum
    read(nunit,*) uion_e
    ! read(nunit,'(a25,f8.5,a7,f7.5)') cdum,beta0_e,cdum2,rmax_e
    read(nunit,'(a1)') cdum
    read(nunit,'(a1)') cdum
    read(nunit,'(a1)') cdum

    pol_rot_vel_e = 0.
    zeff_e = 0.
    tor_rot_freq_e = 0.

    if (alpha_on .eq. 0) then

       if (DIIID_u .eq. 0) then
          i=1
          do
             read(nunit,*,end=10) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                  den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                  temp_elec_e(i),pres_beam_e(i),pres_thermal_e(i), &
                  pres_equil_e(i),tor_rot_vel_e(i),pol_rot_vel_e(i)
             i=i+1
          end do
10        close(unit=nunit)
       end if

       if (DIIID_u .eq. 1 .or. DIIID_u .eq. 2) then
          i=1
          do
             read(nunit,*,end=20) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
                  den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
                  temp_elec_e(i),pres_beam_e(i),pres_thermal_e(i), &
                  pres_equil_e(i),zeff_e(i),tor_rot_freq_e(i),tor_rot_vel_e(i)
             i=i+1
          end do
20        close(unit=nunit)
       end if

    end if

    if (alpha_on .eq. 1) then
       i=1
       do
          read(nunit,*,end=30) rho_e(i),qprof(i),den_beam_e(i),den_ion_e(i),den_elec_e(i), &
               den_alpha_e(i),den_imp_e(i),temp_beam_e(i),temp_ion_e(i), &
               temp_elec_e(i),temp_alpha_e(i),pres_beam_e(i),pres_thermal_e(i), &
               pres_equil_e(i),tor_rot_vel_e(i),pol_rot_vel_e(i)
          i=i+1
       end do
30     close(unit=nunit)
    end if

    if (rho_e(1) == 0.0) then
       ns0=i-1
       rhosq_e=rho_e*rho_e
    else
       ns0=i
       do i=ns0,2,-1
          rho_e(i)=rho_e(i-1)
          qprof(i)=qprof(i-1)
          den_beam_e(i)=den_beam_e(i-1)
          den_ion_e(i)=den_ion_e(i-1)
          den_elec_e(i)=den_elec_e(i-1)
          den_imp_e(i)=den_imp_e(i-1)
          temp_beam_e(i)=temp_beam_e(i-1)
          temp_ion_e(i)=temp_ion_e(i-1)
          temp_elec_e(i)=temp_elec_e(i-1)
          pres_beam_e(i)=pres_beam_e(i-1)
          pres_thermal_e(i)=pres_thermal_e(i-1)
          pres_equil_e(i)=pres_equil_e(i-1)
          zeff_e(i)=zeff_e(i-1)
          tor_rot_vel_e(i)=tor_rot_vel_e(i-1)
          tor_rot_freq_e(i)=tor_rot_freq_e(i-1)
          pol_rot_vel_e(i)=pol_rot_vel_e(i-1)
       end do
       rho_e(1)=0.0
       rhosq_e=rho_e*rho_e
       qprof(1)=(qprof(2)*rhosq_e(3)-qprof(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       den_beam_e(1)=(den_beam_e(2)*rhosq_e(3)-den_beam_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       den_ion_e(1)=(den_ion_e(2)*rhosq_e(3)-den_ion_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       den_elec_e(1)=(den_elec_e(2)*rhosq_e(3)-den_elec_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       den_imp_e(1)=(den_imp_e(2)*rhosq_e(3)-den_imp_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       temp_beam_e(1)=(temp_beam_e(2)*rhosq_e(3)-temp_beam_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       temp_ion_e(1)=(temp_ion_e(2)*rhosq_e(3)-temp_ion_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       temp_elec_e(1)=(temp_elec_e(2)*rhosq_e(3)-temp_elec_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       pres_beam_e(1)=(pres_beam_e(2)*rhosq_e(3)-pres_beam_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       pres_thermal_e(1)=(pres_thermal_e(2)*rhosq_e(3)-pres_thermal_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       pres_equil_e(1)=(pres_equil_e(2)*rhosq_e(3)-pres_equil_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       zeff_e(1)=(zeff_e(2)*rhosq_e(3)-zeff_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       tor_rot_vel_e(1)=(tor_rot_vel_e(2)*rhosq_e(3)-tor_rot_vel_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       tor_rot_freq_e(1)=(tor_rot_freq_e(2)*rhosq_e(3)-tor_rot_freq_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       pol_rot_vel_e(1)=(pol_rot_vel_e(2)*rhosq_e(3)-pol_rot_vel_e(3)*rhosq_e(2))/ &
            (rhosq_e(3)-rhosq_e(2))
       if (alpha_on .eq. 1) then
          do i=ns0,2,-1
             den_alpha_e(i)=den_alpha_e(i-1)
             temp_alpha_e(i)=temp_alpha_e(i-1)
          end do
          den_alpha_e(1)=(den_alpha_e(2)*rhosq_e(3)-den_alpha_e(3)*rhosq_e(2))/ &
               (rhosq_e(3)-rhosq_e(2))
          temp_alpha_e(1)=(temp_alpha_e(2)*rhosq_e(3)-temp_alpha_e(3)*rhosq_e(2))/ &
               (rhosq_e(3)-rhosq_e(2))
       end if
    end if

    if (ns0 > 500) then
       if (myPE == 0) write(6,'("*** ae_profiles: ns0 =",i4," > 500. stop ***")') ns0 
       stop
    end if

    uion = uion_e                                                                    !! ion species
    bt0 = B0_e                                                                       !! magnetic field axis
    rmajr = R0_e                                                                     !! major radius
    rminr = a_e                                                                      !! minor radius
    mu0 = 4.e-7*xpi                                                                  !! vacuum magnetic permeability
    va0 = B0_e/sqrt(uion_e*mu0*den_ion_e(1)*1.e+20*1.672e-27)                        !! Alfven velocity axis
    if (DIIID_u .eq. 1) va0 = B0_e/sqrt(uion_e*mu0*den_ion_e(1)*1.e+19*1.672e-27)  
    omgcya = 1.602e-19*B0_e*R0_e/(1.672e-27*uion_e*va0)                              !! cyclotron FR
    betalf = 2.*mu0*(pres_beam_e(1)*1.e+3)/(B0_e**2)                                 !! beta EP axis
    bet00 = betath_factor*2.*mu0*(pres_thermal_e(1)*1.e+3)/(B0_e**2)                 !! beta thermal plasma axis
    vf0 = sqrt(2.0_IDP)*(sqrt(1000.*temp_beam_e(1)*1.6e-19/(spe1*1.672e-27)))/va0    !! normalized EP thermal velocity axis
    vthi = sqrt(1000.*temp_ion_e(1)*1.6e-19/(uion_e*1.672e-27))/va0                  !! normalized ion thermal velocity axis
    vthe = sqrt(1000.*temp_elec_e(1)*1.6e-19/(9.109e-31))/va0                        !! normalized electron thermal velocity axis
    vtor0 = tor_rot_vel_e(1)*1.e+5/va0                                               !! normalized toroidal velocity axis
    qe = 1.602e-19                                                                   !! electron charge (C)
    me = 9.109e-31                                                                   !! electron mass (kg)
    epsil = 8.85e-12                                                                 !! vacuum permittivity (F/m)
    kblotz = 1.602e-19                                                               !! Conversion factor eV to J

    allocate(bspl(ns0),cspl(ns0),dspl(ns0))

    if (nstres == 0) then

       if (myPE == 0) then
          write(*,'("External Profiles: omgcya = ",f8.4," betalf = ",f8.4, &
               " bet00 = ",f8.4,/,"vthi = ",f8.4," vthe = ",f8.4, &
               " vtor0 = ",f8.4)') omgcya, betalf, bet00, vthi, vthe, vtor0

          write(*,'("bt0 = ",1pe12.5," rmajr = ",0pf8.4, &
               " rminr = ",f8.4,/,"va0 = ",1pe15.7, &
               " pres_beam0 = ",e15.7)') B0_e, rmajr,rminr,va0, pres_beam_e(1)
          write(*,'("Conversion factor from code frequency to kHz: ",1pe15.7)') va0/(2000.*xpi*R0_e)  
       end if

       allocate(qprofile(0:mj),dnnbi(0:mj),dne(0:mj),tbn(0:mj),ti(0:mj),te(0:mj),pep(0:mj),pthermal(0:mj), &
            dnnbinn(0:mj),tbnnn(0:mj),pepnn(0:mj),pthermalnn(0:mj),vAlfven(0:mj),vtherm_elecP(0:mj), &
            vthermalep(0:mj),vzt_eqp(0:mj),vth_eqp(0:mj),ptot(0:mj),ptotnn(0:mj),vtherm_ionP(0:mj), &
            etann(0:mj))
       if (alpha_on .eq. 1) allocate(dnalpha(0:mj),talpha(0:mj),dnalphann(0:mj),talphann(0:mj))

       call spline(ns0,rhosq_e,qprof,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          qprofile(je) = seval(ns0,xx,rhosq_e,qprof,bspl,cspl,dspl)
       end do

       call spline(ns0,rhosq_e,den_beam_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          dnnbinn(je) = 1.e+20*seval(ns0,xx,rhosq_e,den_beam_e,bspl,cspl,dspl)
          dnnbinn(je) = max(dnnbinn(je),1e15)
       end do
       dnnbinn_max=maxval(dnnbinn(0:mj))
       dnnbi=dnnbinn/dnnbinn_max

    end if

    allocate(dnenn(0:mj),dni(0:mj),dninn(0:mj),tinn(0:mj),tenn(0:mj))

    call spline(ns0,rhosq_e,den_elec_e,bspl,cspl,dspl)
    do je=0,mj
       xx = r(je)*r(je)
       dnenn(je) = 1.e+20*seval(ns0,xx,rhosq_e,den_elec_e,bspl,cspl,dspl)
    end do
    dnenn_max=maxval(dnenn(0:mj))

    call spline(ns0,rhosq_e,den_ion_e,bspl,cspl,dspl)
    do je=0,mj
       xx = r(je)*r(je)
       dninn(je) = 1.e+20*seval(ns0,xx,rhosq_e,den_ion_e,bspl,cspl,dspl)
    end do
    dninn_max=maxval(dninn(0:mj))
    dni=dninn/dninn_max  

    call spline(ns0,rhosq_e,temp_ion_e,bspl,cspl,dspl)
    do je=0,mj
       xx = r(je)*r(je)
       tinn(je) = 1.e+3*seval(ns0,xx,rhosq_e,temp_ion_e,bspl,cspl,dspl)
    end do

    call spline(ns0,rhosq_e,temp_elec_e,bspl,cspl,dspl)
    do je=0,mj
       xx = r(je)*r(je)
       tenn(je) = 1.e+3*seval(ns0,xx,rhosq_e,temp_elec_e,bspl,cspl,dspl)
    end do
    tenn_max=maxval(tenn(0:mj))

    if (nstres == 0) then

       dne=dnenn/dnenn_max  
       vAlfven = (B0_e/sqrt(uion_e*mu0*dninn*1.672e-27_IDP))/va0
       vtherm_ionP = sqrt(tinn*1.6e-19_IDP/(uion_e*1.672e-27_IDP))/va0
       tinn_max=maxval(tinn(0:mj))
       ti=tinn/tinn_max  
       vtherm_elecP = sqrt(tenn*1.6e-19/(9.1094e-31))/va0
       te=tenn/tenn_max  

       call spline(ns0,rhosq_e,temp_beam_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          tbnnn(je) = 1.e+3*seval(ns0,xx,rhosq_e,temp_beam_e,bspl,cspl,dspl)
       end do
       tbnnn_max=maxval(tbnnn(0:mj))
       tbn=tbnnn/tbnnn_max

       if (alpha_on .eq. 1) then
          call spline(ns0,rhosq_e,den_alpha_e,bspl,cspl,dspl)
          do je=0,mj
             xx = r(je)*r(je)
             dnalphann(je) = 1.e+20*seval(ns0,xx,rhosq_e,den_alpha_e,bspl,cspl,dspl)
             dnalphann(je) = max(dnalphann(je),1e15)
          end do
          dnalphann_max=maxval(dnalphann(0:mj))
          dnalpha=dnalphann/dnalphann_max 

          call spline(ns0,rhosq_e,temp_alpha_e,bspl,cspl,dspl)
          do je=0,mj
             xx = r(je)*r(je)
             talphann(je) = 1.e+3*seval(ns0,xx,rhosq_e,temp_alpha_e,bspl,cspl,dspl)
          end do
          talphann_max=maxval(talphann(0:mj))
          talpha=talphann/talphann_max
       end if

       call spline(ns0,rhosq_e,tor_rot_vel_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          vzt_eqp(je) = seval(ns0,xx,rhosq_e,tor_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
       end do

       call spline(ns0,rhosq_e,pol_rot_vel_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          vth_eqp(je) = seval(ns0,xx,rhosq_e,pol_rot_vel_e,bspl,cspl,dspl)*1.e+3_IDP/va0
       end do

       call spline(ns0,rhosq_e,pres_thermal_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          pthermalnn(je) = seval(ns0,xx,rhosq_e,pres_thermal_e,bspl,cspl,dspl)
       end do
       pthermalnn_max=maxval(pthermalnn(0:mj))
       pthermal=pthermalnn/pthermalnn_max

       call spline(ns0,rhosq_e,pres_beam_e,bspl,cspl,dspl)
       do je=0,mj
          xx = r(je)*r(je)
          pepnn(je) = seval(ns0,xx,rhosq_e,pres_beam_e,bspl,cspl,dspl)
          ptotnn(je) = pepnn(je) + pthermalnn(je)
       end do
       ptotnn_max=maxval(ptotnn(0:mj))
       ptot=ptotnn/ptotnn_max

       if (DIIID_u .eq. 1 .or. DIIID_u .eq. 2) then

          do je=0,mj
             if (EP_dens_on .eq. 0) then
                dnnbinn(je) = dnnbinn(je)/10.
                dnnbi(je) = dnnbinn(je)/dnnbinn(1)
             end if
             dnenn(je) = dnenn(je)/10.
             dne(je) = dnenn(je)/dnenn(0)
             dninn(je) = dninn(je)/10.
             dni(je) = dninn(je)/dninn(0)
             vAlfven(je) = (B0_e/sqrt(uion_e*mu0*dninn(je)*1.672e-27_IDP))/va0
             vzt_eqp(je) = 100._IDP*vzt_eqp(je)
             vth_eqp(je) = 100._IDP*vth_eqp(je)
          end do

       end if

    end if

    if (DIIID_u .eq. 1 .or. DIIID_u .eq. 2) then
       dnenn_max=dnenn_max/10.
       dninn_max=dninn_max/10.
    end if

    coul_log = 24. - log(1.e-3*sqrt(dnenn_max)/tenn_max)                                              !! Coulomb logarithm (Te > 10 eV)
    etactte = 0.02116*uion_e*uion_e*coul_log*qe*qe*sqrt(me)/(epsil*epsil*sqrt(kblotz*kblotz*kblotz))  !! ctte resistivity

    if (nstres == 0) then

       denseq=dne
       teeq=te
       tieq=ti
       if (EP_dens_on .eq. 0) nfeq=dnnbi 
       if (EP_vel_on .eq. 0) vfova=sqrt(tbnnn*qe/(spe1*1.672e-27))/(va0*LcA3)
       vthermalep=sqrt(tbnnn*qe/(spe1*1.672e-27))

       if (q_prof_on .eq. 1) qq=qprofile
       if (deltaq .ne. 0.0) then
          qq = 1.0_IDP/qqinv + deltaq
          qqinv = 1.0_IDP/qq
       else if (deltaiota .ne. 0.0) then
          qqinv = qqinv + deltaiota
          qq = 1.0_IDP/qqinv
       end if
       etann=etactte*dninn/(dnenn*sqrt(tenn*tenn*tenn))

       vzt_eq=0.0; vth_eq=0.0   !defaults if Eq_vel_on = 0 and Eq_velp_on = 0
       if (Eq_vel_on .eq. 1) vzt_eq=vzt_eqp
       if (Eq_velp_on .eq. 1) vth_eq=vth_eqp
       if (Eq_Presseq_on .eq. 1) then
          preq=pthermal
          if (Eq_Presstot_on .eq. 1) preq=ptot
       end if

       if (alpha_on .eq. 1) then
          if (Alpha_dens_on .eq. 0) nalpeq=dnalpha 
          if (Alpha_vel_on .eq. 0) valphaova=sqrt(talphann*qe/(spe2*1.672e-27))/(va0*LcA3alp)
       end if

       if (Edge_on .eq. 1) then
          do je=edge_p,mj
             etann(je)=etann(je-2) + (r(je)-r(je-2))*(etann(je-1)-etann(je-2))/(r(je-1)-r(je-2))
          end do
       end if

       eta=etann/etann(0)

    end if

    if (ieldamp_on .eq. 1) then
       xnuelc0 = 2.89e-12*dninn_max*coul_log/(tenn_max*sqrt(tenn_max))         !! electron-ion collision FR axis (MKS)
       xnuelc0 = rmajr*xnuelc0/va0                                             !! normalized electron-ion collision FR axis 
    end if

  end subroutine ae_profiles

END MODULE equilibrium
