        module myjdqz_mod

        implicit none

        private

        integer, parameter :: r8 = kind(1.0d0)
        integer :: m,n
        complex(kind=r8), dimension(:,:,:), allocatable :: Adiag,
     >        Asub,Asup
        complex(kind=r8), dimension(:,:,:), allocatable :: Bdiag,
     >        Bsub,Bsup
        complex(kind=r8), dimension(:,:,:), allocatable :: D,L,U
        integer, dimension(:,:), allocatable :: ipiv
	integer :: noeqn, mn_col, ns

        public :: Adiag,Asub,Asup
        public :: Bdiag,Bsub,Bsup
        public :: D,L,U,ipiv
        public :: myjdqz_setup
        public :: m,n
	public :: noeqn, mn_col, ns

        contains

         subroutine matmulij(nlines,ilist,jlist,dlist,neq,x,y)
         implicit none
         integer nlines
         integer ilist(nlines),jlist(nlines)
         double precision dlist(nlines)
         integer neq
         double complex x(neq),y(neq)

         integer ip,i,j
         double precision aij

         do ip=1,neq
           y(ip) = 0.0d0
         enddo
         do ip=1,nlines
           i = ilist(ip)
           j = jlist(ip)
           aij = dlist(ip)
           y(i) = y(i) + aij * x(j)
         enddo

         return
         end subroutine matmulij
        subroutine scatter_entries(nitem,ilist,jlist,dlist,             
     &          m,n, Cdiag,Csub,Csup)
        integer, intent(in) :: nitem
        integer, dimension(nitem),intent(in) :: ilist,jlist
        real(kind=r8), dimension(nitem) :: dlist
        integer,intent(inout) :: m,n

        complex(kind=r8), dimension(m,m,n),intent(inout)::Cdiag,
     >      Csub,Csup
        integer :: ip,ia,ja,  i,j,k
        real(kind=r8) :: aij
        logical :: is_diag, is_sub, is_sup

        do ip=1,nitem
          ia = ilist(ip)
          ja = jlist(ip)
          aij = dlist(ip)

          k = int(ia/m)  + 1
          i = mod(ia,m)
          j = mod(ja,m)
          if (i.eq.0) then
                i = m
                k = k - 1
          endif
          if (j.eq.0) then
                j = m
          endif

          is_diag = ((k-1)*m+i.eq.ia).and.((k-1)*m+j.eq.ja)
          is_sub = ((k-1)*m+i.eq.ia).and.( (k-1-1)*m+j.eq.ja)
          is_sup = ((k-1)*m+i.eq.ia).and.( (k-1+1)*m+j.eq.ja)
          if (is_diag) then
!            --------------
!            diagonal block
!            --------------
             Cdiag(i,j,k) = aij
          else if (is_sub) then
             Csub(i,j,k) = aij
          else  if (is_sup) then
             Csup(i,j,k) = aij
          else
            write(*,*) 'scatter_entries: ia,ja,i,j,k ',                 
     &            ia,ja,i,j,k
          endif
        enddo

        return
        end subroutine scatter_entries
        subroutine scatter_add_entries(nitem,ilist,jlist,dlist,         
     &          m,n, Cdiag,Csub,Csup)
        integer, intent(in) :: nitem
        integer, dimension(nitem),intent(in) :: ilist,jlist
        real(kind=r8), dimension(nitem) :: dlist
        integer,intent(inout) :: m,n

        complex(kind=r8), dimension(m,m,n),intent(inout)::Cdiag,
     >     Csub,Csup
        integer :: ip,ia,ja,  i,j,k
        real(kind=r8) :: aij
        logical :: is_diag, is_sub, is_sup

        do ip=1,nitem
          ia = ilist(ip)
          ja = jlist(ip)
          aij = dlist(ip)

          k = int(ia/m)  + 1
          i = mod(ia,m)
          j = mod(ja,m)
          if (i.eq.0) then
                i = m
                k = k - 1
          endif
          if (j.eq.0) then
                j = m
          endif

          is_diag = ((k-1)*m+i.eq.ia).and.((k-1)*m+j.eq.ja)
          is_sub = ((k-1)*m+i.eq.ia).and.( (k-1-1)*m+j.eq.ja)
          is_sup = ((k-1)*m+i.eq.ia).and.( (k-1+1)*m+j.eq.ja)
          if (is_diag) then
!            --------------
!            diagonal block
!            --------------
             Cdiag(i,j,k) = Cdiag(i,j,k) + aij
          else if (is_sub) then
             Csub(i,j,k) = Csub(i,j,k) + aij
          else  if (is_sup) then
             Csup(i,j,k) = Csup(i,j,k) + aij
          else
            write(*,*) 'scatter_add_entries: ia,ja,i,j,k ',             
     &            ia,ja,i,j,k
          endif
        enddo

        return
        end subroutine scatter_add_entries
        subroutine   readmat(fileC,nlines,ilist,jlist,dlist)
        implicit none
        character*(*) , intent(in) :: fileC
        integer,allocatable,dimension(:),intent(inout):: ilist,jlist
        real(kind=r8),allocatable, dimension(:),intent(inout) :: dlist
        integer, intent(inout) :: nlines

        integer, parameter :: iodev = 15
        integer :: i,j, ierr, ip
        real(kind=r8) :: aij

        open(iodev,file=trim(fileC),form='formatted',                   
     &             status='old',iostat=ierr)
        if (ierr.ne.0) then
                write(*,*) 'readmat: problem opening file ',fileC
                stop '** error ** '
        endif

        rewind(iodev,iostat=ierr)
! ----------------------------------
! first pass to see how many entries
! ----------------------------------
        nlines = 0
        do while (.true.)
          read(iodev, *,iostat=ierr) i,j,aij
          if (ierr.eq.0) then
                  nlines = nlines + 1
          else
                  exit
          endif
        enddo

        if (allocated(ilist)) then
                deallocate(ilist)
        endif
        if (allocated(jlist)) then
                deallocate(jlist)
        endif
        if (allocated(dlist)) then
                deallocate(dlist)
        endif
        allocate(ilist(nlines),jlist(nlines),dlist(nlines),stat=ierr)
        if (ierr.ne.0) then
                write(*,*) 'allocate(ilist),nlines,ierr ',nlines,ierr
                stop '** error ** '
        endif

        rewind(iodev)
        do ip=1,nlines
          read(iodev,*) i,j,aij
          ilist(ip) = i
          jlist(ip) = j
          dlist(ip) = aij
         enddo
         close(iodev)
         return
         end subroutine readmat
c
        subroutine convertij(nlines,ilist,jlist,dlist,                 
     &               m,n, Cdiag,Csub,Csup)
        implicit none
        integer, intent(in) :: nlines
        integer, dimension(nlines), intent(in) :: ilist,jlist
        real (kind=r8), dimension(nlines), intent(in) :: dlist
        integer, intent(inout) :: m,n

        complex (kind=r8),dimension(:,:,:),allocatable,                 
     &           intent(inout) ::     Cdiag,Csub,Csup

        integer :: nrow,ncol,neq,bandw,ierr
        logical :: isok 
        integer :: i, irow,jcol
        integer, dimension(:), allocatable :: jmax, jmin

        nrow = 0
        ncol = 0
        bandw = 0
        do i=1,nlines
          irow = ilist(i)
          jcol = jlist(i)
          nrow = max(nrow,irow)
          ncol = max(ncol,jcol)
          bandw = max( bandw, abs(irow-jcol) )
        enddo
        neq = max(nrow,ncol)
!        -----------------------------------------
!        another method to estimate the bandwidth
!        -----------------------------------------
        allocate( jmax(neq), jmin(neq) )
        jmax(1:neq) = 0
        jmin(1:neq) = neq + 1
        do i=1,nlines
          irow = ilist(i)
          jcol = jlist(i)
          jmax(irow) = max(jmax(irow), jcol)
          jmin(irow) = min(jmin(irow), jcol)
        enddo
        do irow=1,neq
           bandw = max( bandw, jmax(irow)-jmin(irow) )
        enddo
        deallocate( jmax, jmin )

!         bandw = 3*70; m = 70; neq = 3430              !for TAEFL test problem
          bandw = 3*noeqn*mn_col; m = noeqn*mn_col;     !for general TAEFL problem
	  neq = noeqn*mn_col*ns
!        ---------------
!        3*m >= bandw
!        ---------------
        m =  bandw/3
        if (m*3 .lt. bandw) then 
                m = m + 1
        endif
        n = int(neq/m)
        isok = (m*n .eq. neq)
        if (.not.isok) then
          write(*,*) 'problem finding blocksize,bandw,neq,m,n ',       
     &            bandw,neq,m,n
        endif

        if (allocated(Cdiag)) then
                deallocate(Cdiag)
        endif
        if (allocated(Csub)) then
                deallocate(Csub)
        endif
        if (allocated(Csup)) then
                deallocate(Csup)
        endif
        allocate( Cdiag(m,m,n),Csub(m,m,n),Csup(m,m,n),stat=ierr)
        Cdiag = 0
        Csub = 0
        Csup = 0
        if (ierr.ne.0) then
                write(*,*) 'alloc(Cdiag), m,n,ierr ',m,n,ierr
                stop '** error ** '
        endif

        call scatter_entries(nlines,ilist,jlist,dlist,                 
     &                        m,n, Cdiag,Csub,Csup )




        return
        end subroutine convertij
c
         subroutine old_convertij(nlines,ilist,jlist,dlist,                 
     &               m,n, Cdiag,Csub,Csup)
         implicit none
         integer, intent(in) :: nlines
         integer, dimension(nlines), intent(in) :: ilist,jlist
         real(kind=r8), dimension(nlines), intent(in) :: dlist
         integer, intent(inout) :: m,n

         complex(kind=r8),dimension(:,:,:),allocatable,intent(inout) ::    
     &             Cdiag,Csub,Csup

         integer :: nrow,ncol,neq,bandw,ierr
         logical :: isok 
         integer :: i

         nrow = 0
         ncol = 0
         bandw = 0
         do i=1,nlines
           nrow = max(nrow,ilist(i))
           ncol = max(ncol,jlist(i))
           bandw = max(bandw, abs(ilist(i)-jlist(i)))
         enddo
         neq = max(nrow,ncol)

!        ---------------
!        bandw = 2*m - 1
!        ---------------
         m =  int((bandw + 1)/2)
         n = int(neq/m)
         isok = (m*n .eq. neq)
         if (.not.isok) then
           write(*,*) 'problem finding blocksize,bandw,neq,m,n ',       
     &            bandw,neq,m,n
         endif

         if (allocated(Cdiag)) then
                 deallocate(Cdiag)
         endif
         if (allocated(Csub)) then
                 deallocate(Csub)
         endif
         if (allocated(Csup)) then
                 deallocate(Csup)
         endif
         allocate( Cdiag(m,m,n),Csub(m,m,n),Csup(m,m,n),stat=ierr)
         Cdiag = 0
         Csub = 0
         Csup = 0
         if (ierr.ne.0) then
                 write(*,*) 'alloc(Cdiag), m,n,ierr ',m,n,ierr
                 stop '** error ** '
         endif

         call scatter_entries(nlines,ilist,jlist,dlist,                 
     &                        m,n, Cdiag,Csub,Csup )




         return
         end subroutine old_convertij
c	 
         subroutine myjdqz_setup( ztarget, fileA,fileB, neq)
         use symtrd_mod, only : symtrd_factor

         implicit none
         complex(kind=r8), intent(in) :: ztarget
         character*(*), intent(in) :: fileA,fileB
         integer, intent(inout) :: neq

         integer :: nlines,ierr, info
         integer, allocatable, dimension(:) :: ilist,jlist
         real(kind=r8), allocatable, dimension(:) :: dlist
         real(kind=r8) :: time_1, time_2

         integer, parameter :: idebug = 1
         integer :: i
         double complex, allocatable, dimension(:) :: xvec,yvec,yvec2

         allocate( ilist(1),jlist(1),dlist(1))

         call cpu_time(time_1 )
         call readmat( fileA, nlines, ilist,jlist,dlist)
         call convertij( nlines, ilist,jlist,dlist,                     
     &                   m,n, Adiag,Asub,Asup)
         neq = m*n
         if (idebug.ge.2) then
           allocate( xvec(neq),yvec(neq),yvec2(neq))
           do i=1,neq
             xvec(i) = dcmplx(dble(i)/dble(neq), dble(-i)/dble(neq))
           enddo

           call matmulij( nlines,ilist,jlist,dlist, neq, xvec,yvec)
           call amul( neq, xvec,yvec2 )
           write(*,*) 'diff in amul ', maxval(abs(yvec-yvec2))
         endif


         call readmat( fileB, nlines, ilist,jlist,dlist)
         call convertij( nlines, ilist,jlist,dlist,                     
     &                   m,n, Bdiag,Bsub,Bsup)


         if (idebug.ge.2) then
           call matmulij( nlines,ilist,jlist,dlist, neq, xvec,yvec)
           call bmul(neq,xvec,yvec2)
           write(*,*) 'diff in bmul ', maxval(abs(yvec-yvec2))
           deallocate( xvec,yvec,yvec2)
         endif


         if (allocated(ilist)) then
                 deallocate(ilist)
         endif
         if (allocated(jlist)) then
                 deallocate(jlist)
         endif
         if (allocated(dlist)) then
                 deallocate(dlist)
         endif
         call cpu_time(time_2)
         write(*,*) 'reading in and setup matrix took ',                
     &              time_2-time_1,'sec'


!        --------------------
!        setup preconditioner
!        --------------------
         allocate( D(m,m,n), L(m,m,n), U(m,m,n), ipiv(m,n), stat=ierr)
         if (ierr.ne.0) then
                 write(*,*) 'allocte(D), m,n,ierr ',m,n,ierr
                 stop '** error ** '
         endif

         call cpu_time(time_1)
         Adiag = Adiag - ztarget * Bdiag;
         Asub  = Asub  - ztarget * Bsub;
         Asup  = Asup  - ztarget * Bsup;
         call symtrd_factor(m,n,Adiag,Asub,Asup,                        
     &             L,D,ipiv,U, info)
         call cpu_time(time_2)
         if (info.ne.0) then
                 write(*,*) 'symtrd_factor return info = ',info
                 stop '** error ** '
         endif
         write(*,*) 'symtrd_factor took ', time_2-time_1,'sec'


!       --------------
!       extra checks
!       --------------
!
         if (idebug.ge.2) then
           allocate( xvec(neq),yvec(neq),yvec2(neq))
           do i=1,neq
             xvec(i) = dcmplx(dble(i)/dble(neq), dble(-i)/dble(neq))
           enddo
           call amul(neq,xvec,yvec2)
           yvec = yvec2
           call precon(neq,yvec)
           write(*,*) 'diff in precon ',maxval(abs(xvec-yvec))

           xvec = yvec
           call amul(neq,xvec,yvec)
           write(*,*) 'diff in resid ',maxval(abs(yvec-yvec2))

           deallocate( xvec,yvec,yvec2)

         endif

!       ------------------------
!       reread Adiag, Asub, Asup
!       ------------------------

         allocate( ilist(1),jlist(1),dlist(1))

         call cpu_time(time_1 )
         call readmat( fileA, nlines, ilist,jlist,dlist)
         call convertij( nlines, ilist,jlist,dlist,                     
     &                   m,n, Adiag,Asub,Asup)
         neq = m*n

         if (idebug.ge.2) then
           allocate( xvec(neq),yvec(neq),yvec2(neq))
           do i=1,neq
             xvec(i) = dcmplx(dble(i)/dble(neq), dble(-i)/dble(neq))
           enddo

           call matmulij( nlines,ilist,jlist,dlist, neq, xvec,yvec)
           call amul( neq, xvec,yvec2 )
           write(*,*) 'diff in amul ', maxval(abs(yvec-yvec2))
         endif

         if (allocated(ilist)) then
                 deallocate(ilist)
         endif
         if (allocated(jlist)) then
                 deallocate(jlist)
         endif
         if (allocated(dlist)) then
                 deallocate(dlist)
         endif


         return
         end subroutine myjdqz_setup




        end module myjdqz_mod

