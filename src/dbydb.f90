MODULE dbydb

CONTAINS

  subroutine dbydztb(d,a,ltype,c1,c2)

    use param
    use domain
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d

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

END MODULE dbydb
