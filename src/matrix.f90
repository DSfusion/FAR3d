MODULE matrix

  USE param

  IMPLICIT NONE

  ! Block-tridiagonal solver for the toroidal families, on top of BLAS/LAPACK.

  PUBLIC  :: decbt, solbt, wtime

  ! accumulated wall-clock time (s) of the family loops around decbt and solbt
  real(IDP), PUBLIC :: t_decbt = 0., t_solbt = 0.

CONTAINS

  function wtime ()

    ! wall-clock time in seconds

    implicit none

    real(IDP) :: wtime
    integer(8) :: clk, rate

    call system_clock(clk, rate)
    wtime = real(clk,IDP)/real(rate,IDP)

  end function wtime

  subroutine decbt (m, n, a, b, c, ip, ier)

    implicit none

    integer :: m, n, ier
    integer, dimension(m,n) :: ip
    real(IDP), dimension(m,m,n) :: a, b, c

    ! block-tridiagonal matrix decomposition routine.
    ! written by a. c. hindmarsh.
    ! latest revision january 26, 1977  (ag)
    ! reference!  ucid-30150
    !             solution of block-tridiagonal systems of linear
    !             algebraic equations
    !             a.c. hindmarsh
    !             february 1977
    ! the input matrix contains three blocks of elements in each block-row,
    ! including blocks in the (1,3) and (n,n-2) block positions.
    ! decbt uses block gauss elimination, with LAPACK dgetrf/dgetrs for the
    ! diagonal blocks and BLAS dgemm for the block products. partial pivoting
    ! is done within block-rows only.
    ! input..
    !     m = order of each block.
    !     n = number of blocks in each direction of the matrix.
    !         n must be 4 or more.  the complete matrix has order m*n.
    !     a = m by m by n array containing diagonal blocks.
    !         a(i,j,k) contains the (i,j) element of the k-th block.
    !     b = m by m by n array containing the super-diagonal blocks
    !         (in b(*,*,k) for k = 1,...,n-1) and the block in the (n,n-2)
    !         block position (in b(*,*,n)).
    !     c = m by m by n array containing the subdiagonal blocks
    !         (in c(*,*,k) for k = 2,3,...,n) and the block in the
    !         (1,3) block position (in c(*,*,1)).
    !    ip = integer array of length m*n for working storage.
    ! output..
    ! a,b,c = m by m by n arrays containing the block lu decomposition
    !         of the input matrix (lu factors of the diagonal blocks in
    !         LAPACK format).
    !    ip = m by n array of pivot information.  ip(*,k) contains
    !         information for the k-th digonal block.
    !   ier = 0  if no trouble occurred, or
    !       = -1 if the input value of m or n was illegal, or
    !       = k  if a singular matrix was found in the k-th diagonal block.
    ! use solbt to solve the associated linear system.

    integer :: nm1, nm2, k, info
    real(IDP), parameter :: one = 1.0_IDP, mone = -1.0_IDP
    external :: dgetrf, dgetrs, dgemm

    ier = 0
    if (m > 0 .and. n > 3) then
       nm1 = n - 1
       nm2 = n - 2
       ! process the first block-row. -----------------------------------------
       call dgetrf (m, m, a(1,1,1), m, ip(1,1), info)
       if (info /= 0) then
          ier = 1
          return
       end if
       call dgetrs ('N', m, m, a(1,1,1), m, ip(1,1), b(1,1,1), m, info)
       call dgetrs ('N', m, m, a(1,1,1), m, ip(1,1), c(1,1,1), m, info)
       ! adjust b(*,*,2). -----------------------------------------------------
       call dgemm ('N', 'N', m, m, m, mone, c(1,1,2), m, c(1,1,1), m, one, b(1,1,2), m)
       ! main loop.  process block-rows 2 to n-1. -----------------------------
       do k = 2,nm1
          call dgemm ('N', 'N', m, m, m, mone, c(1,1,k), m, b(1,1,k-1), m, one, a(1,1,k), m)
          call dgetrf (m, m, a(1,1,k), m, ip(1,k), info)
          if (info /= 0) then
             ier = k
             return
          end if
          call dgetrs ('N', m, m, a(1,1,k), m, ip(1,k), b(1,1,k), m, info)
       end do
       ! process last block-row and return. -----------------------------------
       call dgemm ('N', 'N', m, m, m, mone, b(1,1,n), m, b(1,1,nm2), m, one, c(1,1,n), m)
       call dgemm ('N', 'N', m, m, m, mone, c(1,1,n), m, b(1,1,nm1), m, one, a(1,1,n), m)
       call dgetrf (m, m, a(1,1,n), m, ip(1,n), info)
       if (info /= 0) then
          ier = n
          return
       end if
    else
       ier = -1
    end if

  end subroutine decbt

  subroutine solbt (m, n, a, b, c, y, ip)

    implicit none

    integer :: m, n
    integer, dimension(m,n) :: ip
    real(IDP), dimension(m,m,n) :: a, b, c
    real(IDP), dimension(m,n) :: y

    ! solution of block-tridiagonal linear system.
    ! coefficient matrix must have been previously processed by decbt.
    ! m, n, a, b, c, and ip  must not have been changed since call to decbt.
    ! written by a. c. hindmarsh.
    ! input..
    !     m = order of each block.
    !     n = number of blocks in each direction of matrix.
    ! a,b,c = m by m by n arrays containing block lu decomposition
    !         of coefficient matrix from decbt.
    !    ip = m by n integer array of pivot information from decbt.
    !     y = array of length m*n containg the right-hand side vector
    !         (treated as an m by n array here).
    ! output..
    !     y = solution vector, of length m*n.
    ! solbt uses LAPACK dgetrs for the diagonal blocks and BLAS dgemv for
    ! the block-vector products.

    integer :: nm1, nm2, k, kb, info
    real(IDP), parameter :: one = 1.0_IDP, mone = -1.0_IDP
    external :: dgetrs, dgemv

    nm1 = n - 1
    nm2 = n - 2
    ! forward solution sweep. ----------------------------------------------
    call dgetrs ('N', m, 1, a(1,1,1), m, ip(1,1), y(1,1), m, info)
    do k = 2,nm1
       call dgemv ('N', m, m, mone, c(1,1,k), m, y(1,k-1), 1, one, y(1,k), 1)
       call dgetrs ('N', m, 1, a(1,1,k), m, ip(1,k), y(1,k), m, info)
    end do
    call dgemv ('N', m, m, mone, c(1,1,n), m, y(1,nm1), 1, one, y(1,n), 1)
    call dgemv ('N', m, m, mone, b(1,1,n), m, y(1,nm2), 1, one, y(1,n), 1)
    call dgetrs ('N', m, 1, a(1,1,n), m, ip(1,n), y(1,n), m, info)
    ! backward solution sweep. ---------------------------------------------
    do kb = 1,nm1
       k = n - kb
       call dgemv ('N', m, m, mone, b(1,1,k), m, y(1,k+1), 1, one, y(1,k), 1)
    end do
    call dgemv ('N', m, m, mone, c(1,1,1), m, y(1,3), 1, one, y(1,1), 1)

  end subroutine solbt

END MODULE matrix
