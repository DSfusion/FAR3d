        module grandom_mod
        implicit none
        integer, parameter :: r8 = kind(1.0d0)
        integer, parameter :: r4 = kind(1.0)
        integer, parameter :: wp = r8


        interface grandom_number
        module procedure                                                
     &     grandom_number_z3,grandom_number_z2,grandom_number_z1,
     &     grandom_number_d3,grandom_number_d2,grandom_number_d1
        end interface grandom_number


        contains

        subroutine grandom_number_z3(A)
        implicit none
        complex(kind=wp), dimension(:,:,:) :: A

        integer :: j

        do j=lbound(A,3),ubound(A,3)
          call grandom_number_z2(A(:,:,j))
        enddo
        return
        end subroutine grandom_number_z3

        subroutine grandom_number_z2(A)
        implicit none
        complex(kind=wp), dimension(:,:) :: A

        integer :: j

        do j=lbound(A,2),ubound(A,2)
          call grandom_number_z1(A(:,j))
        enddo

        return
        end subroutine grandom_number_z2

        subroutine grandom_number_z1(A)
        implicit none
        complex(kind=wp), dimension(:) :: A

        real(kind=wp), dimension(size(A,1),2) :: x

          call random_number(x)
          A(:) = cmplx( x(:,1), x(:,2), kind=wp)
        return
        end subroutine grandom_number_z1




        subroutine grandom_number_d3(A)
        implicit none
        real(kind=wp), dimension(:,:,:) :: A

        integer :: j

        do j=lbound(A,3),ubound(A,3)
          call grandom_number_d2(A(:,:,j))
        enddo
        return
        end subroutine grandom_number_d3

        subroutine grandom_number_d2(A)
        implicit none
        real(kind=wp), dimension(:,:) :: A

        integer :: j

        do j=lbound(A,2),ubound(A,2)
          call grandom_number_d1(A(:,j))
        enddo

        return
        end subroutine grandom_number_d2

        subroutine grandom_number_d1(A)
        implicit none
        real(kind=wp), dimension(:) :: A

        real(kind=wp), dimension(size(A,1),2) :: x

          call random_number(x)
          A(:) = x(:,1)
        return
        end subroutine grandom_number_d1


        end module grandom_mod
