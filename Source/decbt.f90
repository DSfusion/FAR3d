	subroutine decbt (m, n, a, b, c, ip, ier)

		use param
		implicit none

		integer :: m, n, ier
		integer, dimension(m,n) :: ip
		real(IDP), dimension(m,m,n) :: a, b, c

		interface
			subroutine dec (n, a, ip, ier)
				use param
				implicit none
				integer :: n, ier
				integer, dimension(:) :: ip
				real(IDP), dimension(:,:) :: a
			end subroutine dec
			subroutine sol (n, a, b, ip)
				use param
				implicit none
				integer :: n
				integer, dimension(:) :: ip
				real(IDP), dimension(:,:) :: a
				real(IDP), dimension(:) :: b
			end subroutine sol
		end interface

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
		! decbt uses block gauss elimination and subroutines dec and sol
		! for solution of blocks.  partial pivoting is done within
		! block-rows only.
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
		!         of the input matrix.
		!    ip = m by n array of pivot information.  ip(*,k) contains
		!         information for the k-th digonal block.
		!   ier = 0  if no trouble occurred, or
		!       = -1 if the input value of m or n was illegal, or
		!       = k  if a singular matrix was found in the k-th diagonal block.
		! use solbt to solve the associated linear system.
		! decbt calls subroutines  dec(m,m0,a,ip,ier)  and  sol(m,m0,a,y,ip)
		! for solution of m by m linear systems.

		integer :: nm1, nm2, km1, i, j, k, l
		real(IDP) :: dp

		if (m > 0 .and. n > 3) then
			nm1 = n - 1
			nm2 = n - 2
		! process the first block-row. -----------------------------------------
			call dec (m, a(:,:,1), ip(:,1), ier)
			if (ier /= 0) then
		 		ier = 1
		 		return
			end if
		 	do j = 1,m
		 		call sol (m, a(:,:,1), b(:,j,1), ip(:,1))
				call sol (m, a(:,:,1), c(:,j,1), ip(:,1))
			end do
		! adjust b(*,*,2). -----------------------------------------------------
			do j = 1,m
				do i = 1,m
  	 				dp = 0.
					do l = 1,m
						dp = dp + c(i,l,2)*c(l,j,1)
					end do
					b(i,j,2) = b(i,j,2) - dp
				end do
			end do
		! main loop.  process block-rows 2 to n-1. -----------------------------
			do k = 2,nm1
				km1 = k - 1
				do j = 1,m
					do i = 1,m
						dp = 0.
						do l = 1,m
							dp = dp + c(i,l,k)*b(l,j,km1)
						end do
						a(i,j,k) = a(i,j,k) - dp
					end do
				end do
				call dec (m, a(:,:,k), ip(:,k), ier)
				if (ier /= 0) then
	 				ier = k
					return
				end if
				do j = 1,m
					call sol (m, a(:,:,k), b(:,j,k), ip(:,k))
				end do
			end do
		! process last block-row and return. -----------------------------------
			do j = 1,m
	 			do i = 1,m
					dp = 0.
					do l = 1,m
						dp = dp + b(i,l,n)*b(l,j,nm2)
	 				end do
					c(i,j,n) = c(i,j,n) - dp
				end do
			end do
			do j = 1,m
				do i = 1,m
					dp = 0.
					do l = 1,m
						dp = dp + c(i,l,n)*b(l,j,nm1)
					end do
					a(i,j,n) = a(i,j,n) - dp
				end do
			end do
			call dec (m, a(:,:,n), ip(:,n), ier)
			if (ier /= 0) ier = n
		else
	 	 ier = -1
		end if

	end subroutine decbt

	subroutine dec (n, a, ip, ier)

		use param
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

		use param
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