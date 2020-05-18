!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Eigensolver tools !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!-----------------------------------------------------------------------
      subroutine mat_out(ic)
      USE param
      USE cotrol
      USE domain
      USE equil
      USE dynamo
      USE scratch
      implicit none
      integer ieqn, ivar, l, ln, lmx1, lmx2
      integer j, ic, lp, l2, l3
      integer row_cab, col_cab, offset,offset2, row, col, colA, colB, colC
      integer, dimension(:,:), allocatable :: rc_location
      integer :: icount, overdim
      real(IDP), dimension(:), allocatable :: cab, cab_all
      real(IDP) :: scl
      logical :: bin_out
!-----------------------------------------------------------------------
!      bin_out = .true.
!      if(.not.bin_out) then
!       if(ic .eq. 1) open(unit=77,file="a_matrix.dat",status="unknown",form="formatted")
!       if(ic .eq. 1) open(unit=76,file="jdqz.dat",status="unknown",form="formatted")
!       if(ic .eq. 2) open(unit=77,file="b_matrix.dat",status="unknown",form="formatted")
!      else if(bin_out) then
!       if(ic .eq. 1) open(unit=77,file="a_matrix.dat",status="unknown",form="unformatted")
!       if(ic .eq. 1) open(unit=76,file="jdqz.dat",status="unknown",form="formatted")
!       if(ic .eq. 2) open(unit=77,file="b_matrix.dat",status="unknown",form="unformatted")
!      end if

      if(ic .eq. 1) open(unit=77,file="a_matrix.dat",status="unknown")
      if(ic .eq. 1) open(unit=76,file="jdqz.dat",status="unknown")
      if(ic .eq. 2) open(unit=77,file="b_matrix.dat",status="unknown")

      write(0,'("mjm1,noeqn,lmaxn: ",i4,2x,i1,2x,i3)') mjm1,noeqn,lmaxn
      if(ic .eq. 1) scl = 1.0_IDP
      if(ic .eq. 2) scl = 1.0_IDP

      offset = noeqn*lmaxn
      offset2 = 2*offset
      if(ic .eq. 1) then
       write(76,'(i5,2x,i8,2x,i8)') mjm1,lmaxn,noeqn
        do l=1,lmaxn
	 write(76,'(2x,i4,3(3x,i4))') mm(lln(l)),nn(lln(l)),lln(l),signl(l)
!	 write(76,'(2x,i4,3x,i4)') mm(lln(l)),nn(lln(l))
	end do
       do j=1,mjm1
        write(76,*) r(j)
       end do
      end if
!
      allocate(cab(3*noeqn*lmaxn))  !temporarily stores C, A, B for each j
      overdim = 3*noeqn*noeqn*lmaxn*lmaxn*mjm1
      if(bin_out) then
       allocate(cab_all(overdim))  !temporarily stores C, A, B for all j
       allocate(rc_location(overdim,2))  !temporarily stores row and col location
      end if
!
      icount = 0
      do j=1,mjm1
       do ieqn=1,noeqn
        do lp=1,lmaxn
         do ivar=1,noeqn
          do l=1,lmaxn
            lmx1 = lp + (ieqn - 1)*lmaxn
            lmx2 = l + (ivar - 1)*lmaxn
	    col_cab = l + lmaxn*(ivar - 1)
	    cab(col_cab) = cmat(lmx1,lmx2,j)*scl
	    cab(col_cab+offset) = amat(lmx1,lmx2,j)*scl
	    cab(col_cab+offset2) = bmat(lmx1,lmx2,j)*scl
	    row = lp + lmaxn*(ieqn - 1) + (j-1)*offset
	    colC = l + lmaxn*(ivar - 1) + (j-2)*offset
	    colA = l + lmaxn*(ivar - 1) + (j-1)*offset
	    colB = l + lmaxn*(ivar - 1) + j*offset
          enddo   !l
         enddo    !ivar
	 if(j .eq. 1) then
	   l2 = offset + 1
	   l3 = 3*offset
	  else if(j .eq. mjm1) then
	   l2 = 1
	   l3 = 2*offset
	  else
	   l2 = 1
	   l3 = 3*offset
	  end if
	  do l=l2,l3	
	   row = lp + lmaxn*(ieqn - 1) + (j-1)*offset
	   col = l + (j-2)*offset
	   if(.not.bin_out) then
	    if(cab(l) .ne. 0.0) then
	     write(77,'(i12,2x,i12,2x,e15.7)') row,col,cab(l)
	     icount = icount + 1
	    end if
	   else if(bin_out) then
	    if(cab(l) .ne. 0.0) then
	     icount = icount + 1
	     cab_all(icount) = cab(l)
	     rc_location(icount,1) = row
	     rc_location(icount,2) = col
	    endif
	   endif
	  enddo   !l
        enddo     !lp
       enddo      !ieqn
      enddo       !j
!
      if (ic .eq. 1) close(unit=76)
      if(bin_out) then
       write(77) overdim
       write(77) icount
       write(77) rc_location
       write(77) cab_all
      end if
      close(unit=77)
      deallocate(cab)
      if(bin_out) then
       deallocate(cab_all)
       deallocate(rc_location)
      end if
      write(0,*) icount,overdim
      return
      end    
!-----------------------------------------------------------------------
      subroutine mat_out0(ic)
      use param
      use cotrol
      use domain
      use equil
      use dynamo
      use scratch
      implicit none
      real(IDP) :: scl
      integer ieqn, ivar, l, ln, lmx1, lmx2
      integer j, ic, lp, l2, l3, iblk, ibmn, ibmx
      integer offset, row, col, colA, colB, colC
!-----------------------------------------------------------------------
      if(ic .eq. 1) open(unit=77,file="a_matrix.dat",status="unknown")
      if(ic .eq. 1) open(unit=76,file="jdqz.dat",status="unknown")
      if(ic .eq. 2) open(unit=77,file="b_matrix.dat",status="unknown")
      write(0,'("mjm1,noeqn,lmaxn: ",i4,2x,i1,2x,i3)')mjm1,noeqn,lmaxn

      offset = noeqn*lmaxn
      if(ic .eq. 1) scl = 1.0_IDP/s
!      if(ic .eq. 1) scl = 1.0_IDP
      if(ic .eq. 2) scl = 1.0_IDP
      write(0,'("scale factor = ",e15.7)') scl
      if(ic .eq. 1) then
       write(76,'(i5,2x,i8,2x,i8)') mjm1,lmaxn,noeqn
        do l=1,lmaxn
	 write(76,'(2x,i4,3x,i4)') mm(lln(l)),nn(lln(l))
	end do
       do j=1,mjm1
        write(76,*) r(j)
       end do
      end if
!
      do j=1,mjm1
       do ieqn=1,noeqn                     !equation row within block
        do lp=1,lmaxn                      !mode row within block
	 if(j .eq. 1) then
	  ibmn = 2; ibmx = 3
	 else if(j .gt. 1 .and. j .lt. mjm1) then
	  ibmn = 1; ibmx = 3
	 else if(j .eq. mjm1) then
	  ibmn = 1; ibmx = 2
	 endif
	 do iblk=ibmn,ibmx                  !loop over major tri-diagonal blocks
          do ivar=1,noeqn                   !variable column within block
           do l=1,lmaxn                     !mode column within block
            lmx1 = lp + (ieqn - 1)*lmaxn    !block row
            lmx2 = l + (ivar - 1)*lmaxn     !block column
	    row =  lmx1 + (j-1)*offset
	    colC = lmx2 + (j-2)*offset
	    colA = lmx2 + (j-1)*offset
	    colB = lmx2 +     j*offset
	   if(iblk .eq. 1) then
            if(cmat(lmx1,lmx2,j) .ne. 0.0_IDP) write(77,'(i12,2x,i12,2x,e15.7)') row, colC, cmat(lmx1,lmx2,j)*scl
	   else if(iblk .eq. 2) then
            if(amat(lmx1,lmx2,j) .ne. 0.0_IDP) write(77,'(i12,2x,i12,2x,e15.7)') row, colA, amat(lmx1,lmx2,j)*scl
	   else if(iblk .eq. 3) then
            if(bmat(lmx1,lmx2,j) .ne. 0.0_IDP) write(77,'(i12,2x,i12,2x,e15.7)') row, colB, bmat(lmx1,lmx2,j)*scl
	   end if
           enddo   !l
          enddo    !ivar
         enddo    !iblk	  
        enddo     !lp
       enddo      !ieqn
      enddo       !j

      if (ic .eq. 1) close(unit=76)
      close(unit=77)
      return
      end subroutine mat_out0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
