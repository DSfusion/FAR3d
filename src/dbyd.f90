MODULE dbyd

CONTAINS

  subroutine dbydzt(d,a,ltype,c1,c2)

    use param
    use domain
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d

    a(:,0)=0. 
    !$OMP PARALLEL DO
    do l=1,lmax
       d(:,l)=c1*d(:,l)-ltype*nn(l)*c2*a(:,l)
    end do
    !$OMP END PARALLEL DO
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

    a(:,0)=0. 
    !$OMP PARALLEL DO PRIVATE(dold,temp,m)
    do l=1,lmax
       dold=d(:,l) 
       temp=-mm(l)*ltype
       d(:,l)=temp*rinv*a(:,l) 
       d(0,l)=0.
       m=abs(mm(l))
       if (k == 1 .or. k == 2) m=abs(m-1)
       if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
       d(:,l)=c1*dold+c2*d(:,l) 
    end do
    !$OMP END PARALLEL DO
    d(:,0)=0. 

  end subroutine dbydth

  subroutine dbydr(d,a,c1,c2,k)

    use param
    use domain
    implicit none

    integer :: k,l,j,m
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d
    real(IDP), dimension(0:mj) :: dold

    a(:,0)=0. 
    !$OMP PARALLEL DO PRIVATE(dold,j,m)
    do l=1,lmax
       dold=d(:,l) 
       do j=1,mjm1
          d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
       end do
       d(0,l)=0.
       m=abs(mm(l))
       if (k == 1 .or. k == 2) m=abs(m-1)
       if (m == 1) d(0,l)=rinv(1)*a(1,l)
       d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
       d(:,l)=c1*dold+c2*d(:,l)
    end do
    !$OMP END PARALLEL DO
    d(:,0)=0. 

  end subroutine dbydr

  subroutine dbydth_par(d,a,ltype,c1,c2,k)

    use param
    use processor
    use var_para
    use domain
    implicit none

    integer :: k,ltype,l,m
    real(IDP) :: c1,c2,temp
    real(IDP), dimension(mj_start:,0:) :: a,d
    real(IDP), dimension(mj_start:mj_end) :: dold

    a(mj_start:mj_end,0)=0. 
    !$OMP PARALLEL DO PRIVATE(dold,temp,m)
    do l=1,lmax
       dold=d(mj_start:mj_end,l)
       temp=-mm(l)*ltype
       d(mj_start:mj_end,l)=temp*rinv(mj_start:mj_end)*a(mj_start:mj_end,l)
       if (myPE == 0) then
          d(0,l)=0.
          m=abs(mm(l))
          if (k == 1 .or. k == 2) m=abs(m-1)
          if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
       end if
       d(mj_start:mj_end,l)=c1*dold+c2*d(mj_start:mj_end,l) 
    end do
    !$OMP END PARALLEL DO
    d(mj_start:mj_end,0)=0. 

  end subroutine dbydth_par

  subroutine grdpar(d,a,ltype,c1,c2)

    use param
    use var_para
    use domain
    use equil
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: a,d

    a(mj_start:mj_end,0)=0. 
    !$OMP PARALLEL DO
    do l=1,lmax
       d(mj_start:mj_end,l)=c1*d(mj_start:mj_end,l)-ltype*(nn(l)-mm(l)*qqinv(mj_start:mj_end))*c2*a(mj_start:mj_end,l)
    end do
    !$OMP END PARALLEL DO
    d(mj_start:mj_end,0)=0. 

  end subroutine grdpar

  subroutine dbydr_par(d,a,c1,c2,k)

    use param
    use processor
    use var_para
    use domain
    implicit none

    integer :: k,l,j,m,mjst,mjnd
    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: a,d
    real(IDP), dimension(mj_start:mj_end) :: dold

    a(mj_start:mj_end,0)=0. 
    mjst=mj_start+1
    mjnd=mj_end-1
    !$OMP PARALLEL DO PRIVATE(dold,j,m)
    do l=1,lmax
       dold=d(mj_start:mj_end,l) 
       do j=mjst,mjnd
          d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
       end do
       if (myPE == 0) then
          d(0,l)=0.
          m=abs(mm(l))
          if (k == 1 .or. k == 2) m=abs(m-1)
          if (m == 1) d(0,l)=rinv(1)*a(1,l)
       end if
       if (myPE == numPEsm1) d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
       d(mj_start:mj_end,l)=c1*dold+c2*d(mj_start:mj_end,l)
    end do
    !$OMP END PARALLEL DO
    d(mj_start:mj_end,0)=0. 

  end subroutine dbydr_par

  subroutine dbydzt_par(d,a,ltype,c1,c2)

    use param
    use var_para
    use domain
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: a,d

    a(:,0)=0. 
    !$OMP PARALLEL DO
    do l=1,lmax
       d(mj_start:mj_end,l)=c1*d(mj_start:mj_end,l)-ltype*nn(l)*c2*a(mj_start:mj_end,l)
    end do
    !$OMP END PARALLEL DO
    d(:,0)=0. 

  end subroutine dbydzt_par

  subroutine dbydtheq(d,a,ltype,c1,c2,k)

    use param
    use domain
    implicit none

    integer :: k,ltype,l,m
    real(IDP) :: c1,c2,temp
    real(IDP), dimension(0:,0:) :: a,d
    real(IDP), dimension(0:mj) :: dold

    a(:,0)=0. 
    do l=1,leqmax
       dold=d(:,l) 
       temp=-mmeq(l)*ltype
       d(:,l)=temp*rinv*a(:,l) 
       d(0,l)=0.
       m=abs(mmeq(l))
       if (k == 1 .or. k == 2) m=abs(m-1)
       if (m == 1) d(0,l)=temp*rinv(1)*a(1,l)
       d(:,l)=c1*dold+c2*d(:,l) 
    end do
    d(:,0)=0. 

  end subroutine dbydtheq

  subroutine grpareq(d,a,ltype,c1,c2)

    use param
    use domain
    use equil
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d

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

    a(:,0)=0. 
    do l=1,leqmax
       dold=d(:,l) 
       do j=1,mjm1
          d(j,l)=dc1m(j)*(a(j-1,l)-a(j,l))+dc1p(j)*(a(j+1,l)-a(j,l))
       end do
       d(0,l)=0.
       m=abs(mmeq(l))
       if (k == 1 .or. k == 2) m=abs(m-1)
       if (m == 1) d(0,l)=rinv(1)*a(1,l)
       d(mj,l)=(a(mj,l)-a(mjm1,l))/(r(mj)-r(mjm1))
       d(:,l)=c1*dold+c2*d(:,l)
    end do
    d(:,0)=0. 

  end subroutine dbydreq

  subroutine dbydzteq(d,a,ltype,c1,c2)

    use param
    use domain
    implicit none

    integer :: ltype,l
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d

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

    dold=d 
    do j=1,mjm1
       d(j)=dc1m(j)*(a(j-1)-a(j))+dc1p(j)*(a(j+1)-a(j))
    end do
    d(0)=0.
    m=abs(mm(l))
    if (k == 1 .or. k == 2) m=iabs(m-1)
    if (m == 1) d(0)=rinv(1)*a(1)
    d(mj)=(a(mj)-a(mjm1))/(r(mj)-r(mjm1))
    d=c1*dold+c2*d 

  end subroutine dbydrl

  subroutine del2c(d,a,c1,c2,k)

    use param
    use domain
    implicit none

    integer :: k,l,m,j
    real(IDP) :: c1,c2
    real(IDP), dimension(0:,0:) :: a,d
    real(IDP), dimension(0:mj) :: dold

    a(:,0)=0. 
    !$OMP PARALLEL DO PRIVATE(dold,j,m)
    do l=1,lmax
       dold=d(0:mj,l)
       m=abs(mm(l))
       do j=1,mjm1
          d(j,l)=del2cm(j)*a(j-1,l)-(del2cm(j)+del2cp(j)+(m*rinv(j))**2)*a(j,l)+del2cp(j)*a(j+1,l)
       end do
       d(mj,l)=(d(mjm1,l)*(r(mj)-r(mj-2))-d(mj-2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mj-2))
       d(0,l)=0.
       if (k == 1 .or. k == 2) m=abs(m-1)
       if (m == 0) d(0,l)=(a(1,l)-a(0,l))*4.*rinv(1)**2
       d(0:mj,l)=c1*dold+c2*d(0:mj,l)
    end do
    !$OMP END PARALLEL DO
    d(:,0)=0. 

  end subroutine del2c

  subroutine del2cl(d,a,c1,c2,k,l)

    use param
    use domain
    implicit none

    integer :: k,l,m,j
    real(IDP) :: c1,c2
    real(IDP), dimension(0:) :: a,d
    real(IDP), dimension(0:mj) :: dold

    dold=d 
    m=abs(mm(l))
    if (k == 1 .or. k == 2) m=abs(m-1)
    do j=1,mjm1
       d(j)=del2cm(j)*a(j-1)-(del2cm(j)+del2cp(j)+(m*rinv(j))**2)*a(j)+del2cp(j)*a(j+1)
    end do
    d(mj)=(d(mjm1)*(r(mj)-r(mj-2))-d(mj-2)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mj-2))
    d(0)=0.
    if (k == 1 .or. k == 2) m=abs(m-1)
    if (m == 0) d(0)=(a(1)-a(0))*4.*rinv(1)**2
    d=c1*dold+c2*d

  end subroutine del2cl

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

    dold=d
    do j=1,mjm1
       d(j)=dc1m(j)*(a(j-1)-a(j))+dc1p(j)*(a(j+1)-a(j))
    end do
    d(0)=0.
    m=0
    if (k == 1 .or. k == 2) m=abs(m-1)
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

    dold=d
    do j=1,mjm1
       d(j)=dc2m(j)*(a(j-1)-a(j))+dc2p(j)*(a(j+1)-a(j))
    end do
    d(0)=0.
    m=0
    if (k == 1 .or. k == 3) m=m+1
    if (m == 0 .or. m == 2) d(0)=d(1)
    d(mj)=d(mjm1)
    d=c1*dold+c2*d

  end subroutine d2bydr20

END MODULE dbyd
