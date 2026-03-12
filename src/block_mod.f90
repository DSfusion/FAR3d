MODULE block_mod

  IMPLICIT NONE

CONTAINS

  subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum3,l,m,n,l1,l2,l3,l1t,l2t,lp,lpp,lpn,ibnd,j,imat,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*(rinv(j))**ith
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+xa(j)*tx(j,l3)
                   end do
                else if (ir == 1) then
                   do j=1,mjm2
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                   imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      ! else if (ivar == 4) then
                      !    amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm2))/ &
                      !                (r(mjm1)-r(mjm2))
                      !    cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm1))/ &
                      !                (r(mjm1)-r(mjm2))
                   else
                      amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                   end if
                else
                   do j=1,mjm2
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                   imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      ! else if (ivar == 4) then
                      !    amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm2))/ &
                      !                (r(mjm1)-r(mjm2))
                      !    cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm1))/ &
                      !                (r(mjm1)-r(mjm2))
                   else
                      amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine blockj

  subroutine blockjl(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(:,:) :: coef
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum3,l,m,n,l1,l2,l3,l1t,l2t,lp,lpp,lpn,ibnd,j,imat,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*coef(j,lp)*(rinv(j))**ith
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+xa(j)*tx(j,l3)
                   end do
                else if (ir == 1) then
                   do j=1,mjm2
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                   imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      ! else if (ivar == 4) then
                      !    amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm2))/ &
                      !                (r(mjm1)-r(mjm2))
                      !    cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm1))/ &
                      !                (r(mjm1)-r(mjm2))
                   else
                      amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                   end if
                else
                   do j=1,mjm2
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                   imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1,l3)
                      ! else if (ivar == 4) then
                      !    amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      !    amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm2))/ &
                      !                (r(mjm1)-r(mjm2))
                      !    cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*(r(mj)-r(mjm1))/ &
                      !                (r(mjm1)-r(mjm2))
                   else
                      amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                      bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1,l3)
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine blockjl

  subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,itz,ilp,ich,i,mnum3,l,m,n,l1,l2,l2t,lp,ibnd,j,imat,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          lp=l
          if (ich == -1) lp=lo(l)
          if (lp == 0) cycle
          l1=lp+(ieqn-1)*mnumn(i)-lnumn(i-1)
          x=coef
          do il=1,ith
             x=x*m
          end do
          do il=1,izt
             x=x*n
          end do
          if (x == 0.) cycle
          x=x*iph
          do j=1,mjm1
             xa(j)=x*(rinv(j))**ith
          end do
          if (ir == 0) then
             do j=1,mjm1
                imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                amat(imat)=amat(imat)+xa(j)*tx(j)
             end do
          else if(ir == 1) then
             do j=1,mjm2
                imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j)
                cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j)
                bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j)
             end do
             imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
             if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                amat(imat)=amat(imat)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)
                cmat(imat)=cmat(imat)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)
                bmat(imat)=bmat(imat)+wt1p(mjm1,2)*xa(mjm1)*tx(mjm1)
                ! else if (ivar == 4) then
                !    amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    amat(imat)=amat(imat)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
                !                (r(mjm1)-r(mjm2))
                !    cmat(imat)=cmat(imat)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
                !                (r(mjm1)-r(mjm2))
             else
                amat(imat)=amat(imat)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)
                cmat(imat)=cmat(imat)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)
                bmat(imat)=bmat(imat)+wt1p(mjm1,1)*xa(mjm1)*tx(mjm1)
             end if
          else
             do j=1,mjm2
                imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j)
                cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j)
                bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j)
             end do
             imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
             if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                amat(imat)=amat(imat)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)
                cmat(imat)=cmat(imat)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)
                bmat(imat)=bmat(imat)+wt2p(mjm1,2)*xa(mjm1)*tx(mjm1)
                ! else if (ivar == 4) then
                !    amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
                !    amat(imat)=amat(imat)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm2))/ &
                !                (r(mjm1)-r(mjm2))
                !    cmat(imat)=cmat(imat)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*(r(mj)-r(mjm1))/ &
                !                (r(mjm1)-r(mjm2))
             else
                amat(imat)=amat(imat)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)
                cmat(imat)=cmat(imat)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)
                bmat(imat)=bmat(imat)+wt2p(mjm1,1)*xa(mjm1)*tx(mjm1)
             end if
          end if
       end do
    end do

  end subroutine block0

  subroutine block0_dlsq(tx,ieqn,ivar,coef)

    use param
    use var_para
    use cotrol
    use domain
    use dynamo
    implicit none

    integer :: ieqn,ivar
    real(IDP) :: coef
    real(IDP), dimension(0:) :: tx
    integer :: i,l,m,l1,l2,l2t,j,mnum3,imat

    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          l1=l2t+(ieqn-1)*mnumn(i)
          do j=1,mjm1
             imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
             amat(imat)=amat(imat)-coef*m*m*rinv(j)*rinv(j)*tx(j)
          end do
          imat=l1+mnum3*(l2-1)+nskpn(i)
          if (m == 0) then
             amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))
             bmat(imat)=bmat(imat)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))
          else
             amat(imat)=amat(imat)-coef*tx(1)*(del2cm(1)+del2cp(1))
             bmat(imat)=bmat(imat)+coef*tx(1)*del2cp(1)
          end if
          do j=2,mjm2
             imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
             amat(imat)=amat(imat)-coef*tx(j)*(del2cm(j)+del2cp(j))
             cmat(imat)=cmat(imat)+coef*tx(j)*del2cm(j)
             bmat(imat)=bmat(imat)+coef*tx(j)*del2cp(j)
          end do
          imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
          if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
             ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
             !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
             amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
                  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
             cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
                  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
             ! else if (ivar == 4) then
             !    amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
             !                (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))
             !    cmat(imat)=cmat(imat)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
             !                (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))
          else
             amat(imat)=amat(imat)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))
             cmat(imat)=cmat(imat)+coef*tx(mjm1)*del2cm(mjm1)
          end if
       end do
    end do

  end subroutine block0_dlsq

  subroutine block0_dlsqnr(prf,ieqn,ivar,coef,dprfdr)

    use param
    use var_para
    use cotrol
    use domain
    use dynamo
    implicit none

    integer :: ieqn,ivar
    real(IDP), dimension(:) :: coef
    real(IDP), dimension(0:) :: prf,dprfdr
    integer :: i,l,m,l1,l2,l2t,j,mnum3,imat

    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          l1=l2t+(ieqn-1)*mnumn(i)
          do j=1,mjm1
             imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
             amat(imat)=amat(imat)-coef(i)*m*m*rinv(j)*rinv(j)*prf(j)
          end do
          imat=l1+mnum3*(l2-1)+nskpn(i)
          if (m == 0) then
             amat(imat)=amat(imat)-coef(i)*(prf(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))+ &
                  dprfdr(1)*(dc1m(1)+dc1p(1)-r(2)**2*dc1m(1)/(r(2)**2-r(1)**2)))
             bmat(imat)=bmat(imat)+coef(i)*(prf(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2)) + &
                  dprfdr(1)*(dc1p(1)-r(1)**2*dc1m(1)/(r(2)**2-r(1)**2)))
          else
             amat(imat)=amat(imat)-coef(i)*(prf(1)*(del2cm(1)+del2cp(1))+dprfdr(1)*(dc1m(1)+dc1p(1)))
             bmat(imat)=bmat(imat)+coef(i)*(prf(1)*del2cp(1)+dprfdr(1)*dc2p(1))
          end if
          do j=2,mjm2
             imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
             amat(imat)=amat(imat)-coef(i)*(prf(j)*(del2cm(j)+del2cp(j))+dprfdr(j)*(dc1m(j)+dc1p(j)))
             cmat(imat)=cmat(imat)+coef(i)*(prf(j)*del2cm(j)+dprfdr(j)*dc1m(j))
             bmat(imat)=bmat(imat)+coef(i)*(prf(j)*del2cp(j)+dprfdr(j)*dc1p(j))
          end do
          imat=l1+mnum3*(l2-1+mnum3*mjm2)+nskpn(i)
          if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
             ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
             !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
             amat(imat)=amat(imat)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+ &
                  dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1))- &
                  (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
                  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
             cmat(imat)=cmat(imat)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1)- &
                  (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
                  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))
             ! else if (ivar == 4) then
             !    amat(imat)=amat(imat)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+ &
             !                                   dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1))- &
             !                                   (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
             !                                   (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))
             !    cmat(imat)=cmat(imat)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1)- &
             !                                   (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
             !                                   (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))
          else
             amat(imat)=amat(imat)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1)))
             cmat(imat)=cmat(imat)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1))
          end if
       end do
    end do

  end subroutine block0_dlsqnr

  subroutine block_dlsq(ieqn,ivar,coef)

    use param
    use domain
    use equil
    implicit none

    integer :: ieqn,ivar,l
    real(IDP), dimension(0:) :: coef
    real(IDP), dimension(0:mj,0:leqmax) :: tx

    do l=1,leqmax
       tx(:,l)=coef*lplrr(:,l)
    end do
    call blockj(tx,1,ieqn,ivar,0,2,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplrt(:,l)
    end do
    call blockj(tx,-1,ieqn,ivar,1,1,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplrz(:,l)
    end do
    call blockj(tx,-1,ieqn,ivar,0,1,1,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lpltt(:,l)
    end do
    call blockj(tx,1,ieqn,ivar,2,0,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lpltz(:,l)
    end do
    call blockj(tx,1,ieqn,ivar,1,0,1,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplzz(:,l)
    end do
    call blockj(tx,1,ieqn,ivar,0,0,2,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplr(:,l)
    end do
    call blockj(tx,1,ieqn,ivar,0,1,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplt(:,l)
    end do
    call blockj(tx,-1,ieqn,ivar,1,0,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplz(:,l)
    end do
    call blockj(tx,-1,ieqn,ivar,0,0,1,1.0_IDP)

  end subroutine block_dlsq

  subroutine blockj_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum3,l,m,n,l1,l2,l3,l1t,l2t,lp,lpp,lpn,ibnd,j,imat,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*(rinv(j))**ith*abs(n-m*qqinv(j))
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+xa(j)*tx(j,l3)
                   end do
                else if (ir == 1) then
                   do j=1,mjm1
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt10(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt1m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt1p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                else
                   do j=1,mjm1
                      imat=l1+mnum3*(l2-1+mnum3*(j-1))+nskpn(i)
                      amat(imat)=amat(imat)+wt20(j,ibnd)*xa(j)*tx(j,l3)
                      cmat(imat)=cmat(imat)+wt2m(j,ibnd)*xa(j)*tx(j,l3)
                      bmat(imat)=bmat(imat)+wt2p(j,ibnd)*xa(j)*tx(j,l3)
                   end do
                end if
             end do
          end do
       end do
    end do

  end subroutine blockj_landau_grad_parallel

  subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum2,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l1p,l2p,lp,lpp,lpn,ibnd,j,imt1,imt2,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*(rinv(j))**ith
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+xa(j)*tx(j,l3)*xt(imt2)
                   end do
                else if (ir == 1) then
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm2
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                   imt1=l1+mnum3*mjm2+nskpxn(i)
                   imt2=l2+mnum3*mjm2+nskpxn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      ! else if (ivar == 4) then
                      !    yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      !    yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      !    yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)*(r(mj)-r(mjm2))/ &
                      !            (r(mjm1)-r(mjm2))
                      !    yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                      !            (r(mjm1)-r(mjm2))
                   else
                      yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                   end if
                else
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm2
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                   imt1=l1+mnum3*mjm2+nskpxn(i)
                   imt2=l2+mnum3*mjm2+nskpxn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      ! else if (ivar == 4) then
                      !    yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      !    yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      !    yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)*(r(mj)-r(mjm2))/ &
                      !            (r(mjm1)-r(mjm2))
                      !    yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                      !            (r(mjm1)-r(mjm2))
                   else
                      yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine b2lx

  subroutine b2lxjl(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(:,:) :: coef
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum2,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l1p,l2p,lp,lpp,lpn,ibnd,j,imt1,imt2,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*coef(j,lp)*(rinv(j))**ith
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+xa(j)*tx(j,l3)*xt(imt2)
                   end do
                else if (ir == 1) then
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm2
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                   imt1=l1+mnum3*mjm2+nskpxn(i)
                   imt2=l2+mnum3*mjm2+nskpxn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      ! else if (ivar == 4) then
                      !    yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      !    yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      !    yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)*(r(mj)-r(mjm2))/ &
                      !            (r(mjm1)-r(mjm2))
                      !    yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                      !            (r(mjm1)-r(mjm2))
                   else
                      yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                   end if
                else
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm2
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                   imt1=l1+mnum3*mjm2+nskpxn(i)
                   imt2=l2+mnum3*mjm2+nskpxn(i)
                   if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                      ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                      !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                      yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      ! else if (ivar == 4) then
                      !    yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      !    yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                      !    yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)*(r(mj)-r(mjm2))/ &
                      !            (r(mjm1)-r(mjm2))
                      !    yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                      !            (r(mjm1)-r(mjm2))
                   else
                      yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1,l3)*xt(imt2-mnum3)
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine b2lxjl

  subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,itz,ilp,ich,i,mnum2,mnum3,l,m,n,l1,l2,l1p,l2p,l2t,lp,ibnd,j,imt1,imt2,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    !$OMP PARALLEL DO PRIVATE(mnum3,l2t,l,l2,m,n,ibnd,lp,l1,x,il,j,xa,imt1,imt2)
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          lp=l
          if (ich == -1) lp=lo(l)
          if (lp == 0) cycle
          l1=lp+(ieqn-1)*mnumn(i)-lnumn(i-1)
          x=coef
          do il=1,ith
             x=x*m
          end do
          do il=1,izt
             x=x*n
          end do
          if (x == 0.) cycle
          x=x*iph
          do j=1,mjm1
             xa(j)=x*(rinv(j))**ith
          end do
          if (ir == 0) then
             do j=1,mjm1
                imt1=l1+mnum3*(j-1)+nskpxn(i)
                imt2=l2+mnum3*(j-1)+nskpxn(i)
                yt(imt1)=yt(imt1)+xa(j)*tx(j)*xt(imt2)
             end do
          else if (ir == 1) then
             imt1=l1+nskpxn(i)
             imt2=l2+nskpxn(i)
             yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1)*xt(imt2)
             yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
             do j=2,mjm2
                imt1=l1+mnum3*(j-1)+nskpxn(i)
                imt2=l2+mnum3*(j-1)+nskpxn(i)
                yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j)*xt(imt2)
                yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
                yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
             end do
             imt1=l1+mnum3*mjm2+nskpxn(i)
             imt2=l2+mnum3*mjm2+nskpxn(i)
             if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                yt(imt1)=yt(imt1)+wt10(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
                yt(imt1)=yt(imt1)+wt1m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
                ! else if (ivar == 4) then
                !    yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
                !    yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
                !    yt(imt1)=yt(imt1)+dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
                !            (r(mjm1)-r(mjm2))
                !    yt(imt1)=yt(imt1)-dc1p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                !            (r(mjm1)-r(mjm2))
             else
                yt(imt1)=yt(imt1)+wt10(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
                yt(imt1)=yt(imt1)+wt1m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
             end if
          else
             imt1=l1+nskpxn(i)
             imt2=l2+nskpxn(i)
             yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1)*xt(imt2)
             yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1)*xt(imt2+mnum3)
             do j=2,mjm2
                imt1=l1+mnum3*(j-1)+nskpxn(i)
                imt2=l2+mnum3*(j-1)+nskpxn(i)
                yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j)*xt(imt2)
                yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j)*xt(imt2-mnum3)
                yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j)*xt(imt2+mnum3)
             end do
             imt1=l1+mnum3*mjm2+nskpxn(i)
             imt2=l2+mnum3*mjm2+nskpxn(i)
             if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
                ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
                !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
                yt(imt1)=yt(imt1)+wt20(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2)
                yt(imt1)=yt(imt1)+wt2m(mjm1,2)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
                ! else if (ivar == 4) then
                !    yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
                !    yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
                !    yt(imt1)=yt(imt1)+dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2)*(r(mj)-r(mjm2))/ &
                !            (r(mjm1)-r(mjm2))
                !    yt(imt1)=yt(imt1)-dc2p(mjm1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)*(r(mj)-r(mjm1))/ &
                !            (r(mjm1)-r(mjm2))
             else
                yt(imt1)=yt(imt1)+wt20(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2)
                yt(imt1)=yt(imt1)+wt2m(mjm1,1)*xa(mjm1)*tx(mjm1)*xt(imt2-mnum3)
             end if
          end if
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine b2lx0

  subroutine b2lx0_dlsq(tx,ieqn,ivar,coef)

    use param
    use var_para
    use cotrol
    use domain
    use dynamo
    implicit none

    integer :: ieqn,ivar
    real(IDP) :: coef
    real(IDP), dimension(0:) :: tx
    integer :: i,l,m,l1,l2,l2t,j,mnum2,mnum3,imt1,imt2

    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          l1=l2t+(ieqn-1)*mnumn(i)
          do j=1,mjm1
             imt1=l1+mnum3*(j-1)+nskpxn(i)
             imt2=l2+mnum3*(j-1)+nskpxn(i)
             yt(imt1)=yt(imt1)-coef*m*m*rinv(j)*rinv(j)*tx(j)*xt(imt2)
          end do
          imt1=l1+nskpxn(i)
          imt2=l2+nskpxn(i)
          if (m == 0) then
             yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2)
             yt(imt1)=yt(imt1)+coef*tx(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2))*xt(imt2+mnum3)
          else
             yt(imt1)=yt(imt1)-coef*tx(1)*(del2cm(1)+del2cp(1))*xt(imt2)
             yt(imt1)=yt(imt1)+coef*tx(1)*del2cp(1)*xt(imt2+mnum3)
          end if
          do j=2,mjm2
             imt1=l1+mnum3*(j-1)+nskpxn(i)
             imt2=l2+mnum3*(j-1)+nskpxn(i)
             yt(imt1)=yt(imt1)-coef*tx(j)*(del2cm(j)+del2cp(j))*xt(imt2)
             yt(imt1)=yt(imt1)+coef*tx(j)*del2cm(j)*xt(imt2-mnum3)
             yt(imt1)=yt(imt1)+coef*tx(j)*del2cp(j)*xt(imt2+mnum3)
          end do
          imt1=l1+mnum3*mjm2+nskpxn(i)
          imt2=l2+mnum3*mjm2+nskpxn(i)
          if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
             ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
             !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
             yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
                  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)
             yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
                  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)
             ! else if (ivar == 4) then
             !    yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1)-del2cp(mjm1)* &
             !            (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)
             !    yt(imt1)=yt(imt1)+coef*tx(mjm1)*(del2cm(mjm1)-del2cp(mjm1)* &
             !            (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)
          else
             yt(imt1)=yt(imt1)-coef*tx(mjm1)*(del2cm(mjm1)+del2cp(mjm1))*xt(imt2)
             yt(imt1)=yt(imt1)+coef*tx(mjm1)*del2cm(mjm1)*xt(imt2-mnum3)
          end if
       end do
    end do

  end subroutine b2lx0_dlsq

  subroutine b2lx0_dlsqnr(prf,ieqn,ivar,coef,dprfdr)

    use param
    use var_para
    use cotrol
    use domain
    use dynamo
    implicit none

    integer :: ieqn,ivar
    real(IDP), dimension(:) :: coef
    real(IDP), dimension(0:) :: prf,dprfdr
    integer :: i,l,m,l1,l2,l2t,j,mnum2,mnum3,imt1,imt2

    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          l1=l2t+(ieqn-1)*mnumn(i)
          do j=1,mjm1
             imt1=l1+mnum3*(j-1)+nskpxn(i)
             imt2=l2+mnum3*(j-1)+nskpxn(i)
             yt(imt1)=yt(imt1)-coef(i)*m*m*rinv(j)*rinv(j)*prf(j)*xt(imt2)
          end do
          imt1=l1+nskpxn(i)
          imt2=l2+nskpxn(i)
          if (m == 0) then
             yt(imt1)=yt(imt1)-coef(i)*(prf(1)*(del2cm(1)+del2cp(1)-r(2)**2*del2cm(1)/(r(2)**2-r(1)**2))+ &
                  dprfdr(1)*(dc1m(1)+dc1p(1)-r(2)**2*dc1m(1)/(r(2)**2-r(1)**2)))*xt(imt2)
             yt(imt1)=yt(imt1)+coef(i)*(prf(1)*(del2cp(1)-r(1)**2*del2cm(1)/(r(2)**2-r(1)**2)) + &
                  dprfdr(1)*(dc1p(1)-r(1)**2*dc1m(1)/(r(2)**2-r(1)**2)))*xt(imt2+mnum3)
          else
             yt(imt1)=yt(imt1)-coef(i)*(prf(1)*(del2cm(1)+del2cp(1))+dprfdr(1)*(dc1m(1)+dc1p(1)))*xt(imt2)
             yt(imt1)=yt(imt1)+coef(i)*(prf(1)*del2cp(1)+dprfdr(1)*dc2p(1))*xt(imt2+mnum3)
          end if
          do j=2,mjm2
             imt1=l1+mnum3*(j-1)+nskpxn(i)
             imt2=l2+mnum3*(j-1)+nskpxn(i)
             yt(imt1)=yt(imt1)-coef(i)*(prf(j)*(del2cm(j)+del2cp(j))+dprfdr(j)*(dc1m(j)+dc1p(j)))*xt(imt2)
             yt(imt1)=yt(imt1)+coef(i)*(prf(j)*del2cm(j)+dprfdr(j)*dc1m(j))*xt(imt2-mnum3)
             yt(imt1)=yt(imt1)+coef(i)*(prf(j)*del2cp(j)+dprfdr(j)*dc1p(j))*xt(imt2+mnum3)
          end do
          imt1=l1+mnum3*mjm2+nskpxn(i)
          imt2=l2+mnum3*mjm2+nskpxn(i)
          if ((ivar > 4 .and. ivar < 8) .or. (alpha_on == 1 .and. (ivar == 8 .or. ivar == 9))) then
             ! if ((ivar == 5 .and. lln(l) == l0) .or. ivar == 6 .or. ivar == 7 .or. &
             !     (alpha_on == 1 .and. ((ivar == 8 .and. lln(l) == l0) .or. ivar == 9))) then
             yt(imt1)=yt(imt1)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+ &
                  dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1))- &
                  (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
                  (r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2)
             yt(imt1)=yt(imt1)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1)- &
                  (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
                  (r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2))*xt(imt2-mnum3)
             ! else if (ivar == 4) then
             !    yt(imt1)=yt(imt1)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+ &
             !                               dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1))- &
             !                               (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
             !                               (r(mj)-r(mjm2))/(r(mjm1)-r(mjm2)))*xt(imt2)
             !    yt(imt1)=yt(imt1)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1)- &
             !                               (prf(mjm1)*del2cp(mjm1)+dprfdr(mjm1)*dc1p(mjm1))* &
             !                               (r(mj)-r(mjm1))/(r(mjm1)-r(mjm2)))*xt(imt2-mnum3)
          else
             yt(imt1)=yt(imt1)-coef(i)*(prf(mjm1)*(del2cm(mjm1)+del2cp(mjm1))+dprfdr(mjm1)*(dc1m(mjm1)+dc1p(mjm1)))*xt(imt2)
             yt(imt1)=yt(imt1)+coef(i)*(prf(mjm1)*del2cm(mjm1)+dprfdr(mjm1)*dc1m(mjm1))*xt(imt2-mnum3)
          end if
       end do
    end do

  end subroutine b2lx0_dlsqnr

  subroutine b2lx_dlsq(ieqn,ivar,coef)

    use param
    use domain
    use equil
    implicit none

    integer :: ieqn,ivar,l
    real(IDP), dimension(0:) :: coef
    real(IDP), dimension(0:mj,0:leqmax) :: tx

    do l=1,leqmax
       tx(:,l)=coef*lplrr(:,l)
    end do
    call b2lx(tx,1,ieqn,ivar,0,2,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplrt(:,l)
    end do
    call b2lx(tx,-1,ieqn,ivar,1,1,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplrz(:,l)
    end do
    call b2lx(tx,-1,ieqn,ivar,0,1,1,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lpltt(:,l)
    end do
    call b2lx(tx,1,ieqn,ivar,2,0,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lpltz(:,l)
    end do
    call b2lx(tx,1,ieqn,ivar,1,0,1,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplzz(:,l)
    end do
    call b2lx(tx,1,ieqn,ivar,0,0,2,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplr(:,l)
    end do
    call b2lx(tx,1,ieqn,ivar,0,1,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplt(:,l)
    end do
    call b2lx(tx,-1,ieqn,ivar,1,0,0,1.0_IDP)
    do l=1,leqmax
       tx(:,l)=coef*lplz(:,l)
    end do
    call b2lx(tx,-1,ieqn,ivar,0,0,1,1.0_IDP)

  end subroutine b2lx_dlsq

  subroutine b2lx_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use var_para
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt
    real(IDP) :: coef,x
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(mj) :: xa
    integer :: ity,iph,ich,i,mnum2,mnum3,l,m,n,l1,l2,l3,l1t,l2t,l1p,l2p,lp,lpp,lpn,ibnd,j,imt1,imt2,il

    ity= 1
    if (ivar == 2 .or. ivar == 4 .or. ivar == 6 .or. ivar == 7 .or. ivar == ivalp .or. ivar > iq) ity=-1
    iph=ity**(ith+izt)*(-1)**((ith+izt+1)/2)
    ity=ity*(-1)**(ith+izt)
    ich=itx*ity
    if (ieqn == 2 .or. ieqn == 4 .or. ieqn == 6 .or. ieqn == 7 .or. ieqn == ivalp .or. ieqn > iq) ich=-ich
    do i=n_start,n_end
       mnum3=noeqn*mnumn(i)
       do l2t=1,mnumn(i)
          l=l2t+lnumn(i-1)
          l2=l2t+(ivar-1)*mnumn(i)
          m=mm(lln(l))
          n=nn(lln(l))
          ibnd=1
          if (m == 0) ibnd=2
          do l1t=1,mnumn(i)
             lp=l1t+lnumn(i-1)
             lpp=lp
             if (ich == -1) lpp=lo(lp)
             if (lpp == 0) cycle
             l1=lpp+(ieqn-1)*mnumn(i)-lnumn(i-1)
             lpn=nskp2n(i)+mnumn(i)*(l2t-1)+l1t
             do l3=1,leqmax
                if (ity ==  1 .and. itx ==  1) x=coef*cmapp(lpn,l3)
                if (ity ==  1 .and. itx == -1) x=coef*cmamp(lpn,l3)
                if (ity == -1 .and. itx == -1) x=coef*cmamm(lpn,l3)
                if (ity == -1 .and. itx ==  1) x=coef*cmapm(lpn,l3)
                do il=1,ith
                   x=x*m
                end do
                do il=1,izt
                   x=x*n
                end do
                if (x == 0.) cycle
                x=x*iph
                do j=1,mjm1
                   xa(j)=x*(rinv(j))**ith*abs(n-m*qqinv(j))
                end do
                if (ir == 0) then
                   do j=1,mjm1
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+xa(j)*tx(j,l3)*xt(imt2)
                   end do
                else if (ir == 1) then
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt10(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt1p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm1
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt10(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt1m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt1p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                else
                   imt1=l1+nskpxn(i)
                   imt2=l2+nskpxn(i)
                   yt(imt1)=yt(imt1)+wt20(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2)
                   yt(imt1)=yt(imt1)+wt2p(1,ibnd)*xa(1)*tx(1,l3)*xt(imt2+mnum3)
                   do j=2,mjm1
                      imt1=l1+mnum3*(j-1)+nskpxn(i)
                      imt2=l2+mnum3*(j-1)+nskpxn(i)
                      yt(imt1)=yt(imt1)+wt20(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2)
                      yt(imt1)=yt(imt1)+wt2m(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2-mnum3)
                      yt(imt1)=yt(imt1)+wt2p(j,ibnd)*xa(j)*tx(j,l3)*xt(imt2+mnum3)
                   end do
                end if
             end do
          end do
       end do
    end do

  end subroutine b2lx_landau_grad_parallel

END MODULE block_mod
