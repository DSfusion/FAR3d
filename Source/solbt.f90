	subroutine solbt (m, n, a, b, c, y, ip)

		use param
		implicit none

		integer :: m, n
		integer, dimension(m,n) :: ip
		real(IDP), dimension(m,m,n) :: a, b, c
		real(IDP), dimension(m,n) :: y

		interface
			subroutine sol (n, a, b, ip)
				use param
				implicit none
				integer :: n
				integer, dimension(:) :: ip
				real(IDP), dimension(:,:) :: a
				real(IDP), dimension(:) :: b
			end subroutine sol
		end interface

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
		! solbt makes calls to subroutine sol(m,m0,a,y,ip)
		! for solution of m by m linear systems.

		integer :: nm1, nm2, km1, i, j, k, kp1, kb
		real(IDP) :: dp
		nm1 = n - 1
		nm2 = n - 2
		! forward solution sweep. ----------------------------------------------
		call sol (m, a(:,:,1), y(:,1), ip(:,1))
		do k = 2,nm1
			km1 = k - 1
			do i = 1,m
				dp = 0.
				do j = 1,m
					dp = dp + c(i,j,k)*y(j,km1)
				end do
				y(i,k) = y(i,k) - dp
			end do
			call sol (m, a(:,:,k), y(:,k), ip(:,k))
		end do
		do i = 1,m
			dp = 0.
			do j = 1,m
				dp = dp + c(i,j,n)*y(j,nm1) + b(i,j,n)*y(j,nm2)
			end do
			y(i,n) = y(i,n) - dp
		end do
		call sol (m, a(:,:,n), y(:,n), ip(:,n))
		! backward solution sweep. ---------------------------------------------
		do kb = 1,nm1
			k = n - kb
			kp1 = k + 1
			do i = 1,m
				dp = 0.
				do j = 1,m
					dp = dp + b(i,j,k)*y(j,kp1)
				end do
				y(i,k) = y(i,k) - dp
			end do
		end do
		do i = 1,m
			dp = 0.
			do j = 1,m
				dp = dp + c(i,j,1)*y(j,3)
			end do
			y(i,1) = y(i,1) - dp
		end do

	end subroutine solbt