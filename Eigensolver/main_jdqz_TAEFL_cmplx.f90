c       5/17/17: This version was created to be specific to FAR3D
c                Mostly, it just includes some minor additions over
c                the TAEFL version to keep track of both the sin
c                and cos components of the eigenmode for plotting.
c                It writes out lln and signl and the full eigenfunction
c                for potential (not just the sin component).
c       8/20/14: This version of the code was modified to work
c                with TAEFL, using data from jdqz.dat to determine
c                the matrix bandsizes in mjdqz_mod_TAEFL.f90 (rather
c                than trying to figure that out for itself), and to
c                target a complex root, using command line input
c                only for the real part. Further tuning to be done.
c
        program testjdqz
        use myjdqz_mod, only : myjdqz_setup,noeqn,mn_col,ns
        implicit none

        integer, parameter :: r8 = kind(1.0d0)

        integer :: lwork
        double complex, allocatable, dimension(:,:) :: zwork
        double complex, allocatable, dimension(:) :: alpha,beta
        double complex, allocatable, dimension(:) :: v,Av,Bv,residu
        double complex, allocatable, dimension(:,:) :: eivec
        double complex :: ztarget

        real*8, allocatable, dimension(:,:,:,:) :: egn_vectors
        real*8, allocatable, dimension(:) :: rho, omega
        integer, allocatable, dimension(:) :: im_col, in_col, lln, signl
        integer :: i, k, jj, m, kmax0, ivar, mneg
        real*8 :: enrg1, enrg2, center_freq, growth_rate
        
        integer :: kmax
        integer :: neq, jmin,jmax,maxstep, testspace
        integer :: order, method, gmres_m, cgstab_l, maxnmv
        real(kind=r8) :: tol, lock
        logical :: wanted
        real(kind=8), external :: dznrm2
        integer :: j,ierr
        character*80 fileA, fileB, arg1, arg2
        real(kind=r8) :: time_1, time_2, elapse
        
        CALL getarg(1, arg1)
        read(arg1,'(e12.4)') center_freq
        CALL getarg(2, arg2)
        read(arg2,'(e12.4)') growth_rate
        write(*,'("freq. = ",e12.4,3x,"growth = ",e12.4)')
     >     center_freq, growth_rate
        open(unit=15,file="jdqz.dat",status="old")
        read(15,*) ns, mn_col, noeqn
c       write(*,*) ns,mn_col, noeqn
        allocate(im_col(mn_col), in_col(mn_col),
     >     lln(mn_col), signl(mn_col), stat=ierr)
        allocate(rho(ns), stat=ierr)
        do i=1,mn_col
         read(15,'(2x,i4,3(3x,i4))') im_col(i),in_col(i),lln(i),signl(i)
c        write(*,*) im_col(i), in_col(i)
        end do
        do i=1,ns
         read(15,*) rho(i)
c        write(*,*) rho(i)
        end do
        close(unit=15)
        
c        ztarget = dcmplx(0.05*center_freq,center_freq)               !for TAEFL test problem
        ztarget = dcmplx(growth_rate,center_freq)               !for TAEFL test problem
        write(*,'("eigenvalue target = ",2(e15.7,2x))') ztarget       !for TAEFL test problem
        fileA = 'a_matrix.dat'
        fileB = 'b_matrix.dat'
        call myjdqz_setup(ztarget, fileA,fileB,neq)

        
      tol = 1.d-9

!...  kmax is number of eigenvalues we want
      kmax = 7
      kmax0 = kmax

!...  jmin,jmax size of search space
      jmax = min(150,max(30,3*kmax))     !vary
      jmin = min(jmax,2*kmax)            !vary
!      maxstep = 500                      !vary
      jmax = 150
      maxstep = 300
!...   take it small to avoid missing eigensolutions
      lock = 1.d-9
!...     order =  0: nearest to target
!...     order = -1: smallest real part
!...     order =  1: largest real part
!...     order = -2: smallest complex part
!...     order =  2: largest complex part
      order = 0 
!...     method = 1: gmres(m)
!...     method = 2: cgstab(l)
      method = 1
!...     for gmres(m):
!      gmres_m = 50
      gmres_m = 99
!...     for cgstab(l):
      cgstab_l= 4
!...     maximum number of matvecs in cgstab or gmres
!      maxnmv = 500
      maxnmv = 300

!...  testspace = 3 if a reasonable value for target is known, else 
!...  set testspace = 2
!...  testspace = 1 seems to better match results in matlab
      testspace = 3
!...     Testspace 1: w = "Standard Petrov" * v (Section 3.1.1)
!...     Testspace 2: w = "Standard 'variable' Petrov" * v (Section 3.1.2)
!...     Testspace 3: w = "Harmonic Petrov" * v (Section 3.5.1)
!


      if ( method .eq. 1 ) then
!...   ---------------------------
!...   storage if GMRES(M) is used
!...   ---------------------------
         lwork =  4 +  gmres_m  + 5*jmax + 3*kmax
      else
!...    --------------------------------
!...    storage if BiCGSTAB(ell) is used
!...    --------------------------------
         lwork = 10 + 6*cgstab_l + 5*jmax + 3*kmax
      end if
!... converged eigenvectors if wanted = .true.
!...  eivec are converged eigenvectors if wanted = .true.,
!...  else converged Schur vectors
      wanted = .true.
!
      write(*,*) ns, mn_col, neq, kmax
      allocate( eivec(neq,kmax), alpha(jmax),beta(jmax),                
     &          v(neq),Av(neq),Bv(neq),residu(neq),stat=ierr )
      if (ierr.ne.0) then
              write(*,*) 'allocate(eivec), neq,kmax,ierr ',             
     &                    neq,kmax,ierr 
              stop '** error ** '
      endif


      allocate( zwork(neq,lwork),stat=ierr)
      if (ierr.ne.0) then
              write(*,*) 'allocate(zwork), neq,lwork,ierr',             
     &                     neq,lwork,ierr
              stop '** error ** '
      endif



      call cpu_time(time_1)
      call jdqz(alpha, beta, eivec, wanted, neq, ztarget, tol, 
     &     kmax, jmax, jmin,
     &     method, gmres_m, cgstab_l, maxnmv, maxstep,
     &     lock, order, testspace, zwork, lwork )
!
      call cpu_time(time_2)
      elapse = time_2 - time_1

      deallocate( zwork, stat=ierr)
      if (ierr.ne.0) then
              write(*,*) 'deallocate(zwork), ierr ',ierr
              stop '** error ** '
      endif

!
!...     Compute the norms of the residuals:
      do j = 1, kmax
         v(1:neq) = eivec(1:neq,j)
         call amul(neq, v, Av )
         call bmul(neq, v, Bv )
         residu(1:neq) = beta(j)*Av(1:neq) - alpha(j)*Bv(1:neq)
         print '("lambda(",i2,"): (",1p,e11.4,",",e11.4,                
     &        " )")', j,alpha(j)/beta(j)
         print '(a30,d13.6)', '||beta Ax - alpha Bx||:',
     &          dznrm2( neq, residu, 1 )
      end do
      write(*,10)  elapse, elapse
!
   10 format(1x,'END JDQZ AFTER ',f6.2,' SEC. CPU-TIME AND ', f6.2,
     &       ' SEC. ELAPSED TIME' )

      allocate(omega(kmax), stat=ierr)
      allocate(egn_vectors(kmax,mn_col,noeqn,ns), stat=ierr)
      do k=1,kmax
       if(beta(k) .ne. 0.) omega(k)=imag(alpha(k)/beta(k))
       jj = 0
       do i=1,ns
        do ivar=1,noeqn
         do j = 1, mn_col
          jj = jj + 1
          egn_vectors(k,j,ivar,i) = real(eivec(jj,k))
         end do
        end do
       end do
      end do

      write(*,'("Freq. Range: ",e12.4," to " e12.4)')
     >   minval(omega),maxval(omega)
      open(unit=33,file="egn_mode_asci.dat",status="unknown")
      open(unit=34,file="egn_values.dat",status="unknown")
      enrg1 = 1.; enrg2 = 1.   !dummy values for now
       write(33,*) kmax
       write(33,*) mn_col
       write(33,*) ns
       do i=1,mn_col
        write(33,'(2x,i4,3(3x,i4))') im_col(i),in_col(i),
     >           lln(i),signl(i)
       end do
       do i=1,kmax
        write(33,*) imag(alpha(i)/beta(i))
c        write(34,'(e22.14,2x,e22.14,2(2x,e15.7))')
        write(34,'(e22.14,2x,e22.14)')
     >      real(alpha(i)/beta(i)), imag(alpha(i)/beta(i))
       end do
       do i=1,ns
        write(33,*) rho(i)
       end do
c  The following reduction of the eigenmodes to magnitude of
c  cos and sin components is for plotting purposes. Will need to
c  be modified for other uses.
       do j=1,kmax
        do m = 1,ns
         do i = 1,mn_col
          mneg = mn_col - i + 1
          write(33,*) egn_vectors(j,i,2,m), egn_vectors(j,i,1,m)
          end do
         end do
       end do
       close(unit=33)

      
      stop
      end program



        subroutine precon(neq,q)
        use myjdqz_mod, only : m,n, L,D,ipiv,U
        use symtrd_mod, only : symtrd_solve
!       ----------------------------------
!       subroutine to compute q = inv(K)*q
!       ----------------------------------
        integer, parameter :: r8 = kind(1.0d0)

        integer neq
        double complex q(neq)


        integer, parameter :: nrhs = 1
        complex(kind=r8) :: Brhs(m,n,nrhs)
        integer :: ip,i,k, info

        logical :: isok

        isok = (neq .eq. m*n)
        if (.not.isok) then
                write(*,*) 'precon:invalid neq ',neq
        endif

!       -------------------------------------
!       convert q(1:(m*n)) to Brhs(1:m,1:n,1)
!       -------------------------------------
        ip = 1
!$omp   parallel do private(i,k,ip)
        do k=1,n
        do i=1,m
           ip =  i + (k-1)*m
!           Brhs(i,k,1) = dble(q(ip))
!           Brhs(i,k,2) = dimag(q(ip))
           Brhs(i,k,1) = q(ip)
        enddo
        enddo
        call symtrd_solve( m,n,nrhs, L,D,ipiv,U, Brhs,info)
        if (info.ne.0) then
                write(*,*) 'precon: symtrd_solve return info',info
        endif

!       ------------------------------------------
!       convert Brhs(1:m,1:n,1) back to q(1:(m*n))
!       ------------------------------------------

        ip = 1
!$omp   parallel do private(i,k,ip)
        do k=1,n
        do i=1,m
          ip = i + (k-1)*m
!          q(ip) = dcmplx(Brhs(i,k,1),Brhs(i,k,2))
          q(ip) = Brhs(i,k,1)
        enddo
        enddo
        return
        end subroutine precon
        subroutine amul( neq, q, r )
        use myjdqz_mod, only : m,n, Adiag,Asub,Asup
        use symtrd_mod, only : symtrd_matmul
!       -----------------------------
!       subroutine to compute r = A*q
!       -----------------------------
        integer, parameter :: r8 = kind(1.0d0)

        integer neq
        double complex q(neq),r(neq)

        integer, parameter :: nrhs = 1
        complex(kind=r8), dimension(m,n,nrhs) :: X,Y
        integer :: i,k,ip

!       ---------------------------------------
!       convert q(1:(m*n)) to X(1:m,1:n,1:nrhs)
!       ---------------------------------------
        ip = 1
!$omp parallel do private(i,k,ip)
        do k=1,n
        do i=1,m
           ip = i + (k-1)*m
!           X(i,k,1) = dble(q(ip))
!           X(i,k,2) = dimag(q(ip))
            X(i,k,1) = q(ip)
        enddo
        enddo

        call symtrd_matmul( m,n,nrhs, Adiag,Asub,Asup, X,Y )
         
!       ---------------------------------------
!       convert Y(1:m,1:n,1:nrhs) to r(1:(m*n)) 
!       ---------------------------------------

         ip = 1
!$omp parallel do private(i,k,ip)
         do k=1,n
         do i=1,m
           ip = i + (k-1)*m
!           r(ip) = dcmplx(Y(i,k,1),Y(i,k,2))
           r(ip) = Y(i,k,1)
         enddo
         enddo


         return
         end subroutine amul
        subroutine bmul( neq, q, r )
        use myjdqz_mod, only : m,n, Bdiag,Bsub,Bsup
        use symtrd_mod, only : symtrd_matmul
!       -----------------------------
!       subroutine to compute r = B*q
!       -----------------------------
        integer, parameter :: r8 = kind(1.0d0)
        integer neq
        double complex q(neq),r(neq)

        integer, parameter :: nrhs = 1
        complex(kind=8), dimension(m,n,nrhs) :: X,Y
        integer :: i,k,ip

!       ---------------------------------------
!       convert q(1:(m*n)) to X(1:m,1:n,1:nrhs)
!       ---------------------------------------
        ip = 1
!$omp parallel do private(i,k,ip)
        do k=1,n
        do i=1,m
           ip = i + (k-1)*m
!           X(i,k,1) = dble(q(ip))
!           X(i,k,2) = dimag(q(ip))
           X(i,k,1) = q(ip)
        enddo
        enddo

        call symtrd_matmul( m,n,nrhs, Bdiag,Bsub,Bsup, X,Y )
         
!       ---------------------------------------
!       convert Y(1:m,1:n,1:nrhs) to r(1:(m*n)) 
!       ---------------------------------------

         ip = 1
!$omp parallel do private(i,k,ip)
         do k=1,n
         do i=1,m
           ip = i + (k-1)*m
!           r(ip) = dcmplx(Y(i,k,1),Y(i,k,2))
           r(ip) = Y(i,k,1)
         enddo
         enddo


         return
         end subroutine bmul

