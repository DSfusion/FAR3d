MODULE om_mod

  USE block_mod

  IMPLICIT NONE

CONTAINS

  subroutine om(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt,l
    real(IDP) :: coef
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(0:mj,0:leqmax) :: wrk

    wrk=tx
    call blockj(wrk,itx,ieqn,ivar,ith,ir,izt+1,coef)
    do l=1,leqmax
       wrk(:,l)=r*wrk(:,l)*qqinv
    end do
    call blockj(wrk,itx,ieqn,ivar,ith+1,ir,izt,-coef)

  end subroutine om

  subroutine omc(tx,itx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: itx,ieqn,ivar,ith,ir,izt,l
    real(IDP) :: coef
    real(IDP), dimension(0:,0:) :: tx
    real(IDP), dimension(0:mj,0:leqmax) :: wrk

    wrk=tx
    call b2lx(wrk,itx,ieqn,ivar,ith,ir,izt+1,coef)
    do l=1,leqmax
       wrk(:,l)=r*wrk(:,l)*qqinv
    end do
    call b2lx(wrk,itx,ieqn,ivar,ith+1,ir,izt,-coef)

  end subroutine omc

  subroutine om0(tx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: ieqn,ivar,ith,ir,izt
    real(IDP) :: coef
    real(IDP), dimension(0:) :: tx
    real(IDP), dimension(0:mj) :: wrk

    wrk=tx
    call block0(wrk,ieqn,ivar,ith,ir,izt+1,coef)
    wrk=r*wrk*qqinv
    call block0(wrk,ieqn,ivar,ith+1,ir,izt,-coef)

  end subroutine om0

  subroutine omc0(tx,ieqn,ivar,ith,ir,izt,coef)

    use param
    use cotrol
    use domain
    use equil
    use dynamo
    use scratch
    implicit none

    integer :: ieqn,ivar,ith,ir,izt
    real(IDP) :: coef
    real(IDP), dimension(0:) :: tx
    real(IDP), dimension(0:mj) :: wrk

    wrk=tx
    call b2lx0(wrk,ieqn,ivar,ith,ir,izt+1,coef)
    wrk=r*wrk*qqinv
    call b2lx0(wrk,ieqn,ivar,ith+1,ir,izt,-coef)

  end subroutine omc0

END MODULE om_mod
