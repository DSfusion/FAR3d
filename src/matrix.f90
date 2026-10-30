MODULE matrix

  USE param
!$ USE omp_lib, only: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num

  IMPLICIT NONE

  ! Block-tridiagonal solver for the toroidal families. decbt replaces the
  ! diagonal blocks by their inverses (LAPACK, or the built-in dec/sol with
  ! FAR3D_LAPACK=None); solbt_all splits every block over the OpenMP threads
  ! by column panels, so the time step needs no BLAS.

  PUBLIC  :: decbt, solbt, solbt_all, bt_zero, wtime
  PRIVATE :: panel, psum
#ifdef FAR3D_USE_LAPACK
  PRIVATE :: lu_inverse
#else
  PRIVATE :: dec, sol, sol_inverse
#endif

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

  subroutine bt_zero (m, n, a)

    ! Zero the n blocks of order m in a. This is the first touch of the block
    ! arrays, which places their memory pages on the NUMA node of the touching
    ! thread. Every thread zeroes the panels it later streams in solbt_all
    ! (same split, so OMP_PROC_BIND should be set).

    implicit none

    integer :: m, n
    real(IDP), dimension(m,m,n) :: a
    integer :: k, nt, id, nw, i0, i1

    !$OMP PARALLEL PRIVATE(k, nt, id, nw, i0, i1)
    nt = 1
    id = 0
!$  nt = omp_get_num_threads()
!$  id = omp_get_thread_num()
    call panel (m, nt, id, i0, i1, nw)
    do k = 1,n
       if (i1 >= i0) a(:,i0:i1,k) = 0.
    end do
    !$OMP END PARALLEL

  end subroutine bt_zero

  subroutine panel (m, nt, id, i0, i1, nw)

    ! Rows/columns i0:i1 of the panel of thread id (of nt) in solbt and bt_zero.
    ! The m rows are split in units of 8 (one 64-byte cache line) as evenly as
    ! possible; nw threads get a non-empty panel (i1 < i0 for the others).

    implicit none

    integer :: m, nt, id, i0, i1, nw
    integer :: nu, q, r, u0

    nu = (m + 7)/8
    q = nu/nt
    r = mod(nu, nt)
    u0 = id*q + min(id, r)
    if (id < r) q = q + 1
    i0 = 8*u0 + 1
    i1 = min(m, 8*(u0 + q))
    nw = min(nt, nu)

  end subroutine panel

  subroutine psum (m, i0, nr, nw, part, s)

    ! Sum of the partial products of the nw threads for rows i0:i0+nr-1
    ! (column-panel matrix-vector products in solbt).

    implicit none

    integer :: m, i0, nr, nw, t
    real(IDP), dimension(m,0:nw-1) :: part
    real(IDP), dimension(nr) :: s

    s = part(i0:i0+nr-1,0)
    do t = 1,nw-1
       s = s + part(i0:i0+nr-1,t)
    end do

  end subroutine psum

#ifdef FAR3D_USE_LAPACK
  subroutine lu_inverse (m, a, ip, w)

    ! Replace the LU factors in a (from dgetrf, pivots ip) by the inverse matrix,
    ! computed by solving with the identity as right-hand side (dgetrs scales
    ! better over many threads than dgetri). w is m by m workspace.

    implicit none

    integer :: m, i, j, info
    integer, dimension(m) :: ip
    real(IDP), dimension(m,m) :: a, w
    external :: dgetrs

    !$OMP PARALLEL DO
    do j = 1,m
       w(:,j) = 0.
       w(j,j) = 1.
    end do
    !$OMP END PARALLEL DO
    call dgetrs ('N', m, m, a, m, ip, w, m, info)
    !$OMP PARALLEL DO
    do j = 1,m
       do i = 1,m
          a(i,j) = w(i,j)
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine lu_inverse
#else
  subroutine dec (n, a, ip, ier)

    implicit none

    integer :: n, ier
    integer, dimension(:) :: ip
    real(IDP), dimension(:,:) :: a

    !  matrix triangularization by gauss elimination with partial pivoting.
    !  input..
    !     n = order of matrix.
    !     ndim = declared first dimension of array  a.
    !     a = matrix to be triangularized.
    !  output..
    !     a(i,j), i.le.j = upper triangular factor, u .
    !     a(i,j), i.gt.j = multipliers = lower triangular factor, i - l.
    !     ip(k), k.lt.n = index of k-th pivot row.
    !     ier = 0 if matrix a is nonsingular, or k if found to be
    !           singular at stage k.
    !  row interchanges are finished in u, only partly in l.
    !  use  sol  to obtain solution of linear system.
    !  if ier .ne. 0, a is singular, sol will divide by zero.
    ! ----------------------------------------------------------------------
    !  reference!  a. c. hindmarsh, l. j. sloan, k. w. fong, and
    !              g. h. rodrigue,
    !              dec/sol!  solution of dense systems of linear
    !              algebraic equations,
    !              lawrence livermore laboratory report ucid-30137,
    !              june 1976.

    integer :: nm1, kp1, i, j, k, m
    real(IDP) :: t
    ier = 0
    if (n > 1) then
       nm1 = n - 1
       do k = 1,nm1
          kp1 = k + 1
          !  find the pivot in column k.  search rows k to n. --------------------
          m = k
          do i = kp1,n
             if (abs(a(i,k)) > abs(a(m,k))) m = i
          end do
          ip(k) = m
          !  interchange elements in rows k and m. -------------------------------
          t = a(m,k)
          if (m /= k) then
             a(m,k) = a(k,k)
             a(k,k) = t
          end if
          if (t == 0.) then
             ier = k
             return
          end if
          !  store multipliers in a(i,k), i = k+1,...,n. -------------------------
          t = 1./t
          do i = kp1,n
             a(i,k) = -a(i,k)*t
          end do
          !  apply multipliers to other columns of a. ----------------------------
          do j = kp1,n
             t = a(m,j)
             a(m,j) = a(k,j)
             a(k,j) = t
             if (t == 0.) cycle
             do i = kp1,n
                a(i,j) = a(i,j) + a(i,k)*t
             end do
          end do
       end do
    end if
    if (a(n,n) == 0.) ier = n

  end subroutine dec

  subroutine sol (n, a, b, ip)

    implicit none

    integer :: n
    integer, dimension(:) :: ip
    real(IDP), dimension(:,:) :: a
    real(IDP), dimension(:) :: b

    !  solution of linear system a*x = b using output of dec.
    !  input..
    !     n = order of matrix.
    !     ndim = declared first dimension of array  a.
    !     a = triangularized matrix obtained from dec.
    !     b = right hand side vector.
    !     ip = pivot information vector obtained from dec.
    !  do not use if dec has set ier .ne. 0.
    !  output..
    !     b = solution vector, x .

    !  reference!  a. c. hindmarsh, l. j. sloan, k. w. fong, and
    !              g. h. rodrigue,
    !              dec/sol!  solution of dense systems of linear
    !              algebraic equations,
    !              lawrence livermore laboratory report ucid-30137,
    !              june 1976.

    integer :: nm1, kp1, i, k, m, kb, km1
    real(IDP) :: t
    if (n > 1) then
       nm1 = n - 1
       !  apply row permutations and multipliers to b. ------------------------
       do k = 1,nm1
          kp1 = k + 1
          m = ip(k)
          t = b(m)
          b(m) = b(k)
          b(k) = t
          do i = kp1,n
             b(i) = b(i) + a(i,k)*t
          end do
       end do
       !  back solve. ---------------------------------------------------------
       do kb = 1,nm1
          km1 = n - kb
          k = km1 + 1
          b(k) = b(k)/a(k,k)
          t = -b(k)
          do i = 1,km1
             b(i) = b(i) + a(i,k)*t
          end do
       end do
    end if
    b(1) = b(1)/a(1,1)

  end subroutine sol

  subroutine sol_inverse (m, a, ip, w)

    ! Replace the LU factors in a (from dec, pivots ip) by the inverse matrix,
    ! computed by solving with the columns of the identity. w is m by m
    ! workspace.

    implicit none

    integer :: m, j
    integer, dimension(m) :: ip
    real(IDP), dimension(m,m) :: a, w

    !$OMP PARALLEL DO
    do j = 1,m
       w(:,j) = 0.
       w(j,j) = 1.
       call sol (m, a, w(:,j), ip)
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    do j = 1,m
       a(:,j) = w(:,j)
    end do
    !$OMP END PARALLEL DO

  end subroutine sol_inverse
#endif

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
    ! diagonal blocks and BLAS dgemm for the block products (with
    ! FAR3D_LAPACK=None: the built-in dec/sol and loops). partial pivoting
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
    !         of the input matrix, with the diagonal blocks replaced by
    !         their inverses ainv(k), c(k) by ainv(k)*c(k) for k = 2,...,n
    !         and b(n) by ainv(n)*b(n), as needed by solbt_all.
    !    ip = m by n array of pivot information.  ip(*,k) contains
    !         information for the k-th digonal block.
    !   ier = 0  if no trouble occurred, or
    !       = -1 if the input value of m or n was illegal, or
    !       = k  if a singular matrix was found in the k-th diagonal block.
    ! use solbt to solve the associated linear system.

#ifdef FAR3D_USE_LAPACK
    integer :: nm1, nm2, k, info
    real(IDP), parameter :: one = 1.0_IDP, mone = -1.0_IDP
    external :: dgetrf, dgetrs, dgemm
    real(IDP), dimension(:,:), allocatable :: work

    ier = 0
    if (m > 0 .and. n > 3) then
       nm1 = n - 1
       nm2 = n - 2
       allocate (work(m,m))
       ! process the first block-row. -----------------------------------------
       call dgetrf (m, m, a(1,1,1), m, ip(1,1), info)
       if (info /= 0) then
          ier = 1
          return
       end if
       call dgetrs ('N', m, m, a(1,1,1), m, ip(1,1), b(1,1,1), m, info)
       call dgetrs ('N', m, m, a(1,1,1), m, ip(1,1), c(1,1,1), m, info)
       call lu_inverse (m, a(1,1,1), ip(1,1), work)
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
          ! c(k) is only needed by solbt from here on: store ainv(k)*c(k), so that
          ! solbt needs a single matrix-vector phase per block-row
          call dgetrs ('N', m, m, a(1,1,k), m, ip(1,k), c(1,1,k), m, info)
          call lu_inverse (m, a(1,1,k), ip(1,k), work)
       end do
       ! process last block-row and return. -----------------------------------
       call dgemm ('N', 'N', m, m, m, mone, b(1,1,n), m, b(1,1,nm2), m, one, c(1,1,n), m)
       call dgemm ('N', 'N', m, m, m, mone, c(1,1,n), m, b(1,1,nm1), m, one, a(1,1,n), m)
       call dgetrf (m, m, a(1,1,n), m, ip(1,n), info)
       if (info /= 0) then
          ier = n
          return
       end if
       ! same for the last block-row: c(n) <- ainv(n)*c(n), b(n) <- ainv(n)*b(n)
       call dgetrs ('N', m, m, a(1,1,n), m, ip(1,n), c(1,1,n), m, info)
       call dgetrs ('N', m, m, a(1,1,n), m, ip(1,n), b(1,1,n), m, info)
       call lu_inverse (m, a(1,1,n), ip(1,n), work)
    else
       ier = -1
    end if
#else
    ! built-in version (FAR3D_LAPACK=None): dec/sol for the diagonal blocks.

    integer :: nm1, nm2, km1, j, k, l
    real(IDP), dimension(m) :: dpv
    real(IDP), dimension(:,:), allocatable :: work

    if (m > 0 .and. n > 3) then
       nm1 = n - 1
       nm2 = n - 2
       allocate (work(m,m))
       ! process the first block-row. -----------------------------------------
       call dec (m, a(:,:,1), ip(:,1), ier)
       if (ier /= 0) then
          ier = 1
          return
       end if
       !$OMP PARALLEL DO
       do j = 1,m
          call sol (m, a(:,:,1), b(:,j,1), ip(:,1))
          call sol (m, a(:,:,1), c(:,j,1), ip(:,1))
       end do
       !$OMP END PARALLEL DO
       call sol_inverse (m, a(:,:,1), ip(:,1), work)
       ! adjust b(*,*,2). -----------------------------------------------------
       !$OMP PARALLEL DO PRIVATE(dpv)
       do j = 1,m
          dpv = 0.
          do l = 1,m
             dpv = dpv + c(:,l,2)*c(l,j,1)
          end do
          b(:,j,2) = b(:,j,2) - dpv
       end do
       !$OMP END PARALLEL DO
       ! main loop.  process block-rows 2 to n-1. -----------------------------
       do k = 2,nm1
          km1 = k - 1
          !$OMP PARALLEL DO PRIVATE(dpv)
          do j = 1,m
             dpv = 0.
             do l = 1,m
                dpv = dpv + c(:,l,k)*b(l,j,km1)
             end do
             a(:,j,k) = a(:,j,k) - dpv
          end do
          !$OMP END PARALLEL DO
          call dec (m, a(:,:,k), ip(:,k), ier)
          if (ier /= 0) then
             ier = k
             return
          end if
          !$OMP PARALLEL DO
          do j = 1,m
             call sol (m, a(:,:,k), b(:,j,k), ip(:,k))
             ! c(k) is only needed by solbt from here on: store ainv(k)*c(k)
             call sol (m, a(:,:,k), c(:,j,k), ip(:,k))
          end do
          !$OMP END PARALLEL DO
          call sol_inverse (m, a(:,:,k), ip(:,k), work)
       end do
       ! process last block-row and return. -----------------------------------
       !$OMP PARALLEL DO PRIVATE(dpv)
       do j = 1,m
          dpv = 0.
          do l = 1,m
             dpv = dpv + b(:,l,n)*b(l,j,nm2)
          end do
          c(:,j,n) = c(:,j,n) - dpv
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(dpv)
       do j = 1,m
          dpv = 0.
          do l = 1,m
             dpv = dpv + c(:,l,n)*b(l,j,nm1)
          end do
          a(:,j,n) = a(:,j,n) - dpv
       end do
       !$OMP END PARALLEL DO
       call dec (m, a(:,:,n), ip(:,n), ier)
       if (ier /= 0) then
          ier = n
          return
       end if
       ! same for the last block-row: c(n) <- ainv(n)*c(n), b(n) <- ainv(n)*b(n)
       !$OMP PARALLEL DO
       do j = 1,m
          call sol (m, a(:,:,n), c(:,j,n), ip(:,n))
          call sol (m, a(:,:,n), b(:,j,n), ip(:,n))
       end do
       !$OMP END PARALLEL DO
       call sol_inverse (m, a(:,:,n), ip(:,n), work)
    else
       ier = -1
    end if
#endif

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
    ! the blocks are as stored by decbt (inverses of the diagonal blocks);
    ! the products are split over the OpenMP threads, see solbt_all.

    call solbt_all (1, [m], n, [1], [1], [1], a, b, c, y, ip)

  end subroutine solbt

  subroutine solbt_all (nf, ms, n, loca, loci, locx, a, b, c, y, ip)

    ! Solve the nf independent block-tridiagonal systems (toroidal families) f,
    ! of block order ms(f) and n block rows, factored by decbt. The blocks of
    ! family f start at a(loca(f)), b(loca(f)), c(loca(f)), the pivots at
    ! ip(loci(f)) and the right-hand side (overwritten by the solution) at
    ! y(locx(f)).

    implicit none

    integer :: nf, n
    integer, dimension(nf) :: ms, loca, loci, locx
    real(IDP), dimension(*) :: a, b, c, y
    integer, dimension(*) :: ip

    ! a(*,*,k) holds the inverse of the diagonal block (from decbt), so every
    ! operation is a matrix-vector product. The products are split over the
    ! OpenMP threads by column panels: thread id owns columns (and rows)
    ! i0:i1, computes the partial product
    ! M(:,i0:i1)*x(i0:i1) into part(:,id,ib), and after a barrier sums
    ! the partial products of all threads for its own rows i0:i1. Each thread
    ! therefore streams a contiguous part of every block, and only ever reads the
    ! entries x(i0:i1) of the input vectors that it computed itself. The two
    ! buffers part(:,:,0:1) alternate, so one barrier per product suffices.
    ! decbt stores c(k) = ainv(k)*c(k) (and b(n) = ainv(n)*b(n)), so each
    ! forward block-row is a single phase y(k) = ainv(k)*y(k) - c(k)*y(k-1).
    ! All families are solved at once: in every phase each thread computes its
    ! panel products for block row k of every family, so that one barrier per
    ! block row serves all families. The panels are the same as in bt_zero, so
    ! with OMP_PROC_BIND set every thread streams memory of its own NUMA node.

    integer :: f, k, kb, nt, id, ib, mmax, lx
    integer, dimension(nf) :: i0, i1, nr, nw
    real(IDP), dimension(:,:), allocatable :: tmp
    real(IDP), dimension(:,:,:,:), allocatable :: part

    if (nf < 1) return
    mmax = maxval(ms)
    nt = 1
!$  nt = omp_get_max_threads()
    allocate (part(mmax,0:nt-1,0:1,nf), tmp(mmax,nf))

    !$OMP PARALLEL PRIVATE(f, k, kb, nt, id, ib, lx, i0, i1, nr, nw)
    nt = 1
    id = 0
!$  nt = omp_get_num_threads()
!$  id = omp_get_thread_num()
    do f = 1,nf
       call panel (ms(f), nt, id, i0(f), i1(f), nw(f))
       nr(f) = i1(f) - i0(f) + 1
    end do
    ib = 0

    ! forward solution sweep. ----------------------------------------------
    ! y(:,1) = ainv(1)*y(:,1)
    do f = 1,nf
       if (nr(f) > 0) call pgemv1 (ms(f), nr(f), a(ia(f,i0(f),1)), y(iy(f,i0(f),1)), &
            & part(1,id,ib,f))
    end do
    !$OMP BARRIER
    do f = 1,nf
       if (nr(f) > 0) call psum (mmax, i0(f), nr(f), nw(f), part(1,0,ib,f), y(iy(f,i0(f),1)))
    end do
    ib = 1 - ib
    do k = 2,n-1
       ! y(:,k) = ainv(k)*y(:,k) - c(k)*y(:,k-1)
       do f = 1,nf
          if (nr(f) > 0) then
             call pgemv2 (ms(f), nr(f), a(ia(f,i0(f),k)), y(iy(f,i0(f),k)), &
                  & c(ia(f,i0(f),k)), y(iy(f,i0(f),k-1)), part(1,id,ib,f))
          end if
       end do
       !$OMP BARRIER
       do f = 1,nf
          if (nr(f) > 0) call psum (mmax, i0(f), nr(f), nw(f), part(1,0,ib,f), y(iy(f,i0(f),k)))
       end do
       ib = 1 - ib
    end do
    ! y(:,n) = ainv(n)*y(:,n) - c(n)*y(:,n-1) - b(n)*y(:,n-2)
    do f = 1,nf
       if (nr(f) > 0) then
          call pgemv3 (ms(f), nr(f), a(ia(f,i0(f),n)), y(iy(f,i0(f),n)), &
               & c(ia(f,i0(f),n)), y(iy(f,i0(f),n-1)), b(ia(f,i0(f),n)), &
               & y(iy(f,i0(f),n-2)), part(1,id,ib,f))
       end if
    end do
    !$OMP BARRIER
    do f = 1,nf
       if (nr(f) > 0) call psum (mmax, i0(f), nr(f), nw(f), part(1,0,ib,f), y(iy(f,i0(f),n)))
    end do
    ib = 1 - ib
    ! backward solution sweep. ---------------------------------------------
    do kb = 1,n-1
       k = n - kb
       ! y(:,k) = y(:,k) - b(k)*y(:,k+1)
       do f = 1,nf
          if (nr(f) > 0) call pgemv1 (ms(f), nr(f), b(ia(f,i0(f),k)), y(iy(f,i0(f),k+1)), &
               & part(1,id,ib,f))
       end do
       !$OMP BARRIER
       do f = 1,nf
          if (nr(f) > 0) then
             call psum (mmax, i0(f), nr(f), nw(f), part(1,0,ib,f), tmp(i0(f),f))
             lx = iy(f,0,k)
             y(lx+i0(f):lx+i1(f)) = y(lx+i0(f):lx+i1(f)) - tmp(i0(f):i1(f),f)
          end if
       end do
       ib = 1 - ib
    end do
    ! y(:,1) = y(:,1) - c(1)*y(:,3)
    do f = 1,nf
       if (nr(f) > 0) call pgemv1 (ms(f), nr(f), c(ia(f,i0(f),1)), y(iy(f,i0(f),3)), &
            & part(1,id,ib,f))
    end do
    !$OMP BARRIER
    do f = 1,nf
       if (nr(f) > 0) then
          call psum (mmax, i0(f), nr(f), nw(f), part(1,0,ib,f), tmp(i0(f),f))
          lx = iy(f,0,1)
          y(lx+i0(f):lx+i1(f)) = y(lx+i0(f):lx+i1(f)) - tmp(i0(f):i1(f),f)
       end if
    end do
    !$OMP END PARALLEL

    deallocate (part, tmp)

  contains

    ! p = a*x - c*z for the m by nr panels a and c (one pass over both)
    subroutine pgemv2 (m, nr, a, x, c, z, p)
      integer :: m, nr, j
      real(IDP), dimension(m,nr) :: a, c
      real(IDP), dimension(nr) :: x, z
      real(IDP), dimension(m) :: p
      p = a(:,1)*x(1) - c(:,1)*z(1)
      do j = 2,nr
         p = p + a(:,j)*x(j) - c(:,j)*z(j)
      end do
    end subroutine pgemv2

    ! p = a*x - c*z - b*w for the m by nr panels a, c and b (one pass)
    subroutine pgemv3 (m, nr, a, x, c, z, b, w, p)
      integer :: m, nr, j
      real(IDP), dimension(m,nr) :: a, c, b
      real(IDP), dimension(nr) :: x, z, w
      real(IDP), dimension(m) :: p
      p = a(:,1)*x(1) - c(:,1)*z(1) - b(:,1)*w(1)
      do j = 2,nr
         p = p + a(:,j)*x(j) - c(:,j)*z(j) - b(:,j)*w(j)
      end do
    end subroutine pgemv3

    ! p = b*x for the m by nr panel b
    subroutine pgemv1 (m, nr, b, x, p)
      integer :: m, nr, j
      real(IDP), dimension(m,nr) :: b
      real(IDP), dimension(nr) :: x
      real(IDP), dimension(m) :: p
      p = b(:,1)*x(1)
      do j = 2,nr
         p = p + b(:,j)*x(j)
      end do
    end subroutine pgemv1

    ! index of element (1,j,k) of the blocks of family f in a, b and c
    integer function ia (f, j, k)
      integer :: f, j, k
      ia = loca(f) + ms(f)*(j-1) + ms(f)*ms(f)*(k-1)
    end function ia

    ! index of element (i,k) of the right-hand side of family f in y
    integer function iy (f, i, k)
      integer :: f, i, k
      iy = locx(f) + (i-1) + ms(f)*(k-1)
    end function iy

  end subroutine solbt_all

END MODULE matrix
