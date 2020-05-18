!!!!!!!!!!!!!!!!!!!!!	columms calculation program	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	version 1.0
!!!! WARNING: READ BEFORE RUNNING THE CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The program must be adapted to the number of poloidal modes in the simulation,
!!!! adding or removing egn_vectors(k,p,m) where k is the index of the target modes
!!!! and p the list of poloidal numbers. For example, if your simulation has 6 poloidal
!!!! numbers, you need to have 12 comumns in the output list (different parities).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	parameter definition

	module imput
		implicit none
		integer, parameter :: IDP = kind(1.d0)	
		integer :: nmat,mn_col,ns,ic,nmat_rd
                real(IDP) :: fscale,freq_max,dm_test
		real(IDP), dimension(:), allocatable :: rho,im_col,in_col,dm,dm_rd	
		real(IDP), dimension(:,:), allocatable :: egn_vectors0
		real(IDP), dimension(:,:,:), allocatable :: egn_vectors	
                integer, dimension(:), allocatable :: iacept,i_orig
	        integer :: m,j,i
	        character(len=1) :: t			
		character(len=32) :: formatw='(1pe13.6,100(a1,1pe13.6))'
	end module imput	
	
	program columms
	
    use imput

!	Format of the external file 	

    open(unit=33,file="egn_mode_asci.dat",status="old")
    read(33,*) nmat
    read(33,*) mn_col
    read(33,*) ns   

    allocate(rho(ns))
    allocate(im_col(mn_col))
    allocate(in_col(mn_col))	
    allocate(dm(nmat))	
    allocate (egn_vectors0(mn_col,ns))		
    allocate (iacept(nmat))	
    allocate(i_orig(nmat))		
    allocate (dm_rd(nmat))	
        t=char(9)


    rho(:) = 0.0	
    im_col(:) = 0.0	
    in_col(:) = 0.0	
    dm(:) = 0.0	
    dm_rd(:) = 0.0
    i_orig(:) = 0.0	
    iacept(:) = 0.0
    egn_vectors0(:,:) = 0.0	

    fscale = 1.
    freq_max = 1.e+6

    do i=1,mn_col
      read(33,*) im_col(i),in_col(i)
    end do		
	
!    print *,"im_col=",(im_col(i),i=1,mn_col)
!    print *,"in_col=",(in_col(i),i=1,mn_col)	

    ic = 0	
    do i=1,nmat
      read(33,*) dm(i)	 
      dm_test = (fscale**2)*dm(i)
      if(dm_test .lt. freq_max**2) then
         ic = ic + 1
	 iacept(i) = ic
	 dm_rd(ic) = (fscale**2)*dm(i)
	 i_orig(ic) = i
       else
	 iacept(i) = 0
       endif
    end do	 

!	print *,"dm=",(dm(i),i=1,nmat)

    nmat_rd = ic
    allocate (egn_vectors(nmat_rd,mn_col,ns))
    do i=1,ns
      read(33,*) rho(i)
    end do	
	
!	print *,"rho=",(rho(i),i=1,ns)	
 	
    do j=1,nmat
        do m = 1,ns
          do i=1,mn_col
             read(33,*) egn_vectors0(i,m)
           end do
        end do
        if(iacept(j) .ne. 0) then
          do m = 1,ns
            do i=1,mn_col
              egn_vectors(iacept(j),i,m) = egn_vectors0(i,m)
            end do
          end do
        endif
    end do	
    nmat = nmat_rd

!	print *,"egn_vectors=",((egn_vectors(1,i,m),i=1,mn_col),m=1,ns)		

!	output format and file output 	
		
	open(unit=2,file='egn_vectors',status='new',recl=2048)
	do m=1,ns	
	      write(2,formatw) rho(m),t,egn_vectors(1,1,m),t,egn_vectors(1,2,m),t,egn_vectors(1,3,m),t, &  
                         egn_vectors(1,4,m),t,egn_vectors(1,5,m),t,egn_vectors(1,6,m),t,egn_vectors(1,7,m),t, &
                         egn_vectors(1,8,m),t,egn_vectors(1,9,m),t,egn_vectors(1,10,m)
	end do
        close(2)	
	open(unit=3,file='egn_vectors2',status='new',recl=2048)
	do m=1,ns	
	      write(3,formatw) rho(m),t,egn_vectors(2,1,m),t,egn_vectors(2,2,m),t,egn_vectors(2,3,m),t, &  
                         egn_vectors(2,4,m),t,egn_vectors(2,5,m),t,egn_vectors(2,6,m),t,egn_vectors(2,7,m),t, &
                         egn_vectors(2,8,m),t,egn_vectors(2,9,m),t,egn_vectors(2,10,m)   
	end do
        close(3)	
	open(unit=4,file='egn_vectors3',status='new',recl=2048)
	do m=1,ns	
	      write(4,formatw) rho(m),t,egn_vectors(3,1,m),t,egn_vectors(3,2,m),t,egn_vectors(3,3,m),t, &  
                         egn_vectors(3,4,m),t,egn_vectors(3,5,m),t,egn_vectors(3,6,m),t,egn_vectors(3,7,m),t, &
                         egn_vectors(3,8,m),t,egn_vectors(3,9,m),t,egn_vectors(3,10,m)    
	end do
        close(4)	
	open(unit=5,file='egn_vectors4',status='new',recl=2048)
	do m=1,ns	
	      write(5,formatw) rho(m),t,egn_vectors(4,1,m),t,egn_vectors(4,2,m),t,egn_vectors(4,3,m),t, &  
                         egn_vectors(4,4,m),t,egn_vectors(4,5,m),t,egn_vectors(4,6,m),t,egn_vectors(4,7,m),t, &
                         egn_vectors(4,8,m),t,egn_vectors(4,9,m),t,egn_vectors(4,10,m)  
	end do
        close(5)
	open(unit=6,file='egn_vectors5',status='new',recl=2048)
	do m=1,ns	
	      write(6,formatw) rho(m),t,egn_vectors(5,1,m),t,egn_vectors(5,2,m),t,egn_vectors(5,3,m),t, &  
                         egn_vectors(5,4,m),t,egn_vectors(5,5,m),t,egn_vectors(5,6,m),t,egn_vectors(5,7,m),t, &
                         egn_vectors(5,8,m),t,egn_vectors(5,9,m),t,egn_vectors(5,10,m)  
	end do
        close(6)
	open(unit=7,file='egn_vectors6',status='new',recl=2048)
	do m=1,ns	
	      write(7,formatw) rho(m),t,egn_vectors(6,1,m),t,egn_vectors(6,2,m),t,egn_vectors(6,3,m),t, &  
                         egn_vectors(6,4,m),t,egn_vectors(6,5,m),t,egn_vectors(6,6,m),t,egn_vectors(6,7,m),t, &
                         egn_vectors(6,8,m),t,egn_vectors(6,9,m),t,egn_vectors(6,10,m) 
	end do
        close(7)
	open(unit=8,file='egn_vectors7',status='new',recl=2048)
	do m=1,ns	
	      write(8,formatw) rho(m),t,egn_vectors(7,1,m),t,egn_vectors(7,2,m),t,egn_vectors(7,3,m),t, &  
                         egn_vectors(7,4,m),t,egn_vectors(7,5,m),t,egn_vectors(7,6,m),t,egn_vectors(7,7,m),t, &
                         egn_vectors(7,8,m),t,egn_vectors(7,9,m),t,egn_vectors(7,10,m)  
	end do
        close(8)					

	end program columms		
		