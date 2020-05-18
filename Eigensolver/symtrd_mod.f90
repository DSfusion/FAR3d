










        module symtrd_mod
        use grandom_mod
        implicit none

        contains
        subroutine assert(isok,msg,ival)
        implicit none
        logical, intent(in) :: isok
        character*(*), intent(in) :: msg
        integer, intent(in) :: ival

        if (.not.isok) then
                write(*,*) msg,ival
                stop 'stop in assertion '
        endif
        return
        end subroutine assert
        subroutine symtrd_factor( m,n, A, B,C,                          
     &                                 L, D, ipiv, U, info )
        implicit none
        integer, parameter :: r8 = kind(1.0d0)

        integer, intent(in) :: m,n
        complex (kind=r8), dimension(m,m,n), intent(in) :: A,B,C
        complex (kind=r8), dimension(m,m,n), intent(inout) :: L,D,U
        integer, dimension(m,n), intent(inout) :: ipiv
        integer, intent(inout) :: info
!
!
! factor symmetric block tridiagonal matrix
!
! [A1, C1         ]   [ I             ]   [ D1, U1         ]
! [B2, A2, C2     ] = [ L2, I         ] * [     D2, U2     ]
! [    B3, A3, C3 ]   [     L3, I     ]   [         D3, U3 ]
! [        B4, A4 ]   [         L4, I ]   [             D4 ]
!
!
        logical, parameter :: use_gemm = .true.
        integer, parameter :: idebug = 1

        complex (kind=r8), dimension(m,m) :: Dkm1,Dk,Lk,Ukm1
        integer, dimension(m) :: ipivkm1, ipivk
        integer :: k, linfo
        integer :: ld1,ld2,ld3

        integer :: mm,nn,kk, nrhs
        complex (kind=r8) :: alpha, beta


!  
!  Overall algorithm
!  D1 = A1
!  for k=2:n,
!      Lk = B(k) / D(k-1), 
!      or D(k-1)'\B(k)' -> Lk'
!
!      Uk = Ck
!      Dk = Ak - Lk*U(k-1)
!  end

        info = 0

!  -------
!  D1 = A1
!  -------
        k = 1
        Dk(1:m,1:m) = A(1:m,1:m,k)

        if (idebug.ge.2) then
                write(*,*) 'k,Dk ',k,Dk
        endif

        ld1 = size(Dk,1)
        call zgetrf(m,m,Dk,ld1,ipivk,linfo)
        if (linfo.ne.0) then

          if (idebug.ge.1) then
            write(*,*) 'symtrd_factor: getrf(D1) linfo=',linfo
          endif

          info = -1
          return
        endif

        D(1:m,1:m,k) = Dk(1:m,1:m)
        ipiv(1:m,k) = ipivk(1:m)

!      ------------
!      Uk = Ck
!      ------------
        U(1:m,1:m,1:n) = C(1:m,1:m,1:n)


!      ---------
!      Main loop
!      ---------

        do k=2,n
!      ------------------
!      D(k-1)'\Bk' -> Lk'
!      ------------------
         Dkm1(1:m,1:m) = D(1:m,1:m,k-1)
         ipivkm1(1:m) = ipiv(1:m,k-1)

         Lk(1:m,1:m) = transpose(B(1:m,1:m,k))
         nrhs = m
         ld1 = size(Dkm1,1)
         ld2 = size(Lk,1)
         call zgetrs( 'T', m,nrhs, Dkm1, ld1, ipivkm1, Lk,ld2,linfo)
         if (linfo.ne.0) then

          if (idebug.ge.1) then
            write(*,*) 'symtrd_factor: getrs, k,linfo ',k,linfo
          endif


           linfo = -k
           return
          endif

          L(1:m,1:m,k) = transpose(Lk(1:m,1:m))
          Lk(1:m,1:m) =  L(1:m,1:m,k)

          if (idebug.ge.2) then
            write(*,*) 'k,Lk ',k,Lk
          endif
!      ------------------- 
!      Dk = Ak - Lk*U(k-1)
!      ------------------- 
          Dk(1:m,1:m) = A(1:m,1:m,k)
          Ukm1(1:m,1:m) = U(1:m,1:m,k-1)

          if (use_gemm) then
             mm = m
             nn = m
             kk = m
             alpha = -1.0d0
             beta = 1.0d0
             ld1 = size(Lk,1)
             ld2 = size(Ukm1,1)
             ld3 = size(Dk,1)
             call zgemm( 'N', 'N', mm,nn,kk,                            
     &            alpha, Lk,ld1,  Ukm1, ld2,                            
     &            beta,  Dk,ld3 )
          else
            Dk(1:m,1:m) = Dk(1:m,1:m) -                                 
     &                    matmul(Lk(1:m,1:m),Ukm1(1:m,1:m))
          endif

          if (idebug.ge.2) then
                  write(*,*) 'k, Dk ',k,Dk
          endif

          ld1 = size(Dk,1)
          call zgetrf(m,m,Dk,ld1,ipivk,linfo)
          if (linfo.ne.0) then

            if (idebug.ge.1) then
              write(*,*) 'symtrd_factor: getrf(Dk), k,linfo ',k,linfo
            endif


            info = -k
            return
          endif
          D(1:m,1:m,k) = Dk(1:m,1:m)
          ipiv(1:m,k) = ipivk(1:m)

          enddo


          return
          end subroutine symtrd_factor


        subroutine symtrd_solve( m,n,nrhs, L, D, ipiv, U,               
     &           Brhs_inout, info )
        implicit none
        integer, parameter :: r8 = kind(1.0d0)

        integer, intent(in) :: m,n,nrhs


        complex (kind=r8), dimension(m,m,n),target,                     
     & 
     &          intent(in) :: L, D, U
        integer, dimension(m,n), target, intent(in) :: ipiv

        complex (kind=r8), dimension(m,n,nrhs),                         
     & 
     &          intent(inout) :: Brhs_inout
        integer, intent(inout) :: info

!   note Brhs_inout has dimension(m,n,nrhs)
!   but  internally, we use Brhs with dimension (m,nrhs,n)
!   for better efficiency

!
!
! factor of symmetric block tridiagonal matrix
!
! [A1, B2'        ]   [ I             ]   [ D1, U1         ]
! [B2, A2, B3'    ] = [ L2, I         ] * [     D2, U2     ]
! [    B3, A3, B4']   [     L3, I     ]   [         D3, U3 ]
! [        B4, A4 ]   [         L4, I ]   [             D4 ]
!
!
        integer, parameter :: idebug = 1
        logical, parameter :: use_gemm = .true.

        integer :: linfo, k
        integer :: mm,nn,kk
        complex (kind=r8) :: alpha, beta


        complex (kind=r8), dimension(m,nrhs,n),target :: Brhs


        complex (kind=r8), dimension(:,:,:), pointer :: Y,X
        complex (kind=r8), dimension(:,:), pointer ::                   
     & 
     &           yk,ykm1,xk,xkp1
        complex (kind=r8), dimension(:,:), pointer :: Dk,Uk,Lk
        integer, dimension(:), pointer :: ipivk

        integer :: ld1,ld2,ld3

! Overall algorithm related to LU factorization
!  L * (U * x) =  b
!  (1) solve L*y = b,   (2) solve U * x = y
!
        do k=1,n
          Brhs(1:m,1:nrhs,k) = Brhs_inout(1:m,k,1:nrhs)
        enddo


!  --------------------------
!  Solve L y = b
!  y1 = b1
!  for k=2:n
!      yk = bk - Lk * y(k-1)
!  end
!  --------------------------
        info = 0



        Y => Brhs(1:m,1:nrhs, 1:n)
        X => Brhs(1:m,1:nrhs, 1:n)


        do k=2,n

          ykm1 => Y(1:m,1:nrhs,k-1)
          yk => Y(1:m,1:nrhs,k)
          Lk => L(1:m,1:m,k)



          if (use_gemm) then
                  mm = m
                  nn = nrhs
                  kk = m
                  alpha = -1.0d0
                  beta = 1.0d0
                  ld1 = size(Lk,1)
                  ld2 = size(ykm1,1)
                  ld3 = size(yk,1)
                  call zgemm('N','N', mm,nn,kk,                         
     & 
     &               alpha, Lk, ld1, ykm1, ld2,                         
     &               beta,  yk, ld3 )

          else
             yk(1:m,1:nrhs) = yk(1:m,1:nrhs) -                          
     &            matmul(Lk(1:m,1:m),ykm1(1:m,1:nrhs))
          endif



        enddo


!    -------------------------------------
!   Solve U x = y
!
!    xn = Dn \ yn
!    for k=(n-1) downto 1 
!          xk = Dk \ (yk - Uk * x(k+1) )
!    end
!    -------------------------------------

        k = n


        xk => Y(1:m,1:nrhs,k)
        Dk => D(1:m,1:m,k)
        ipivk => ipiv(1:m,k)


        ld1 = size(Dk,1)
        ld2 = size(xk,1)
        call zgetrs('N',m,nrhs,Dk,ld1,ipivk,  xk,ld2,linfo)
        if (linfo.ne.0) then
          if (idebug.ge.1) then
            write(*,*) 'symtrd_solve: getrs, k,linfo ',k,linfo
          endif

          info = -k
          return
        endif



        

        do k=n-1,1,-1
!          -----------------------------
!          xk = Dk \ (yk - Uk * x(k+1) )
!          -----------------------------


           yk => Y(1:m,1:nrhs,k)
           xk => yk
           xkp1 => X(1:m,1:nrhs,k+1)
           Uk => U(1:m,1:m,k)


           if (use_gemm) then
             mm = m
             nn = nrhs
             kk = m
             alpha = -1.0d0
             beta  = 1.0d0
             ld1 = size(Uk,1)
             ld2 = size(xkp1,1)
             ld3 = size(xk,1)
             call zgemm('N','N', mm,nn,kk,                              
     & 
     &         alpha, Uk,ld1, xkp1, ld2,                                
     &         beta,  xk, ld3 )
           else
             xk(1:m,1:nrhs) = xk(1:m,1:nrhs) -                          
     &                      matmul(Uk(1:m,1:m), xkp1(1:m,1:nrhs))
           endif


           Dk => D(1:m,1:m,k)
           ipivk => ipiv(1:m,k)


           ld1 = size(Dk,1)
           ld2 = size(xk,1)
           call zgetrs('N',m,nrhs,Dk,ld1,ipivk,  xk,ld2,linfo)
           if (linfo.ne.0) then
             if (idebug.ge.1) then
               write(*,*) 'symtrd_solve: getrs, k,linfo ',k,linfo
             endif

             info = -k
             return
           endif



         enddo


         do k=1,n
           Brhs_inout(1:m,k,1:nrhs) = Brhs(1:m,1:nrhs,k)
         enddo

         return
         end subroutine symtrd_solve

        subroutine symtrd_matmul(m,n,nrhs, A, B,C,  X, Y )
        implicit none
        integer, parameter :: r8 = kind(1.0d0)
        integer, intent(in) :: m,n,nrhs


         complex (kind=r8), target,                                     
     & 
     &          dimension(m,m,n), intent(in) ::  A,B,C
         complex (kind=r8), target,                                     
     & 
     &          dimension(m,n,nrhs), intent(in) :: X
         complex (kind=r8), target,                                     
     & 
     &          dimension(m,n,nrhs), intent(inout) :: Y



        integer :: k

         complex (kind=r8), dimension(:,:),                             
     & 
     &                     pointer :: Ak, Bk, Ck
         complex (kind=r8), dimension(:,:),                             
     & 
     &                     pointer :: yk, xkm1,xk,xkp1


        logical, parameter :: use_gemm = .true.
        integer :: mm,nn,kk, ld1,ld2,ld3
         complex (kind=r8) :: alpha,beta
!       ---------
!       form Y
!
!       block diagonal part
!       Y(k) = A(k) * x(k)
!       ---------


        k = 2
        xk => X(1:m,k,1:nrhs)

        do k=1,n

          Ak => A(1:m,1:m,k)
          xk => X(1:m,k,1:nrhs)
          yk => Y(1:m,k,1:nrhs)


          if (use_gemm) then
            mm = m
            nn = nrhs
            kk = m
            alpha = 1.0d0
            beta = 0.0d0
            ld1 = size(Ak,1)
            ld2 = size(xk,1)
            ld3 = size(yk,1)
            call zgemm('N','N',mm,nn,kk,                                
     & 
     &          alpha, Ak,ld1, xk, ld2, beta, yk, ld3 )
          else
            yk(1:m,1:nrhs) = matmul(Ak(1:m,1:m), xk(1:m,1:nrhs))
          endif

!       -----------------------
!       Y(k) = Y(k) + B(k) * x(k-1)  ,  k=2:n
!       -----------------------
        if ((2.le.k).and.(k.le.n)) then


           Bk => B(1:m,1:m,k)
           xkm1 => X(1:m,k-1,1:nrhs)



           if (use_gemm) then
              mm = m
              nn = nrhs
              kk = m
              alpha = 1.0d0
              beta = 1.0d0
              ld1 = size(Bk,1)
              ld2 = size(xkm1,1)
              ld3 = size(yk,1)
              call zgemm('N','N', mm,nn,kk,                             
     & 
     &          alpha, Bk,ld1, xkm1, ld2, beta, yk, ld3 )
           else
              yk(1:m,1:nrhs) = yk(1:m,1:nrhs) +                         
     &               matmul( Bk(1:m,1:m),xkm1(1:m,1:nrhs))
           endif
          endif


!       -------------------------------------
!       Y(k) = Y(k) + C(k)*x(k+1), k=1:(n-1)
!       -------------------------------------
          if ((1.le.k).and.(k.le.(n-1)) ) then


            Ck => C(1:m,1:m,k)
            xkp1 => X(1:m,k+1,1:nrhs)


            if (use_gemm) then
              mm = m
              nn = nrhs
              kk = m
              alpha = 1.0d0
              beta = 1.0d0
              ld1 = size(Ck,1)
              ld2 = size(xkp1,1)
              ld3 = size(yk,1)
              call zgemm('N','N',mm,nn,kk,                              
     & 
     &            alpha, Ck,ld1, xkp1, ld2,                             
     &            beta,  yk, ld3 )
            else
              yk(1:m,1:nrhs) = yk(1:m,1:nrhs) +                         
     &                   matmul( Ck(1:m,1:m), xkp1(1:m,1:nrhs) )
            endif
          endif




          enddo


        return
        end subroutine symtrd_matmul
        subroutine test_symtrd(m,n,nrhs )
        implicit none
        integer, parameter :: r8 = kind(1.0d0)

        integer, parameter :: idebug = 1
        integer, intent(in) :: m,n,nrhs

!       ------------------------------------------------
!       simple driver to test symtrd_factor/symtrd_solve
!       ------------------------------------------------
         complex (kind=r8), dimension(m,m,n) ::                         
     & 
     &           A,B,C,L,D,U
         complex (kind=r8), dimension(m,n,nrhs) ::                      
     & 
     &           X,Brhs,Resid,Xsolve

        integer, dimension(m,n) :: ipiv

        integer :: i,j,k, info
        logical, parameter :: make_symmetric = .false.

         complex (kind=r8),dimension(m,m) :: Ak, Bkp1,Ck


         double precision  :: t1,t2,abserr,maxerr
!       ----------------
!       generate problem
!       ----------------
        write(*,*) 'test_symtrd: m,n,nrhs ',m,n,nrhs

        A = 0
        B = 0
        C = 0
        X = 0

        call cpu_time(t1)
        call grandom_number(A)
        call grandom_number(B)
        call grandom_number(C)
        call grandom_number(X)



! ---------------------
! make matrix symmetric
! ---------------------
        if (make_symmetric) then

        do k=1,n-1
           Bkp1(1:m,1:m) = B(1:m,1:m,k+1)
           Ck(1:m,1:m) = transpose(Bkp1(1:m,1:m))
           C(1:m,1:m,k) = Ck(1:m,1:m)
        enddo

        do k=1,n
           Ak(1:m,1:m) = A(1:m,1:m,k)
           A(1:m,1:m,k) = Ak(1:m,1:m) + transpose(Ak(1:m,1:m))
        enddo

        endif

        call cpu_time(t2)

        write(*,*) 'solver in complex type '

        write(*,*) 'make_symmetric = ', make_symmetric
        write(*,*) 'time to form A,B,X ', t2-t1



        call cpu_time(t1)
        call symtrd_matmul( m,n,nrhs, A,B,C,  X, Brhs )
        call cpu_time(t2)
        write(*,*) 'time to form Brhs ', t2-t1



!       ---------------------
!       perform factorization
!       ---------------------
        info = 0
        call cpu_time(t1)
        call symtrd_factor( m,n,A,B,C,                                  
     &                      L,D,ipiv, U, info)
        call cpu_time(t2)
        if (info.ne.0) then
                write(*,*) 'symtrd_factor return info ',info
        endif

        write(*,*) 'symtrd_factor took ', t2-t1,' sec '

!       -------------
!       perform solve
!       -------------
        Xsolve(1:m,1:n,1:nrhs) = Brhs(1:m,1:n,1:nrhs)
        info = 0
        call cpu_time(t1)
        call  symtrd_solve( m,n,nrhs, L,D,ipiv,U, Xsolve,info)
        call cpu_time(t2)

        if (info.ne.0) then
                write(*,*) 'symtrd_solve return info ',info
        endif

        write(*,*) 'symtrd_solve took ', t2-t1,' sec '

!       ----------------------
!       difference in solution
!       ----------------------
        maxerr = 0.0d0
        do k=1,nrhs
        do j=1,n
        do i=1,m
           abserr = ( abs( X(i,j,k) - Xsolve(i,j,k) ) )
           if (abs(abserr) .gt. abs(maxerr)) then
              maxerr =  abserr

              if (idebug.ge.2) then
                  write(*,*) 'i,j,k ',i,j,k,' X,Xsolve,abserr ',        
     &              X(i,j,k),Xsolve(i,j,k),abs(abserr)
              endif
           endif
        enddo
        enddo
        enddo

        write(*,*) 'difference in solution ',  maxerr

        call symtrd_matmul( m,n,nrhs, A,B,C, Xsolve, Resid)
        maxerr = 0.0d0
        do k=1,nrhs
        do j=1,n
        do i=1,m
          abserr = ( abs( Resid(i,j,k) - Brhs(i,j,k) ))
          if (abs(abserr) .gt. abs(maxerr)) then
             maxerr = abserr

             if (idebug.ge.2) then
               write(*,*) 'i,j,k ',i,j,k,' Resid,Brhs,abserr ',         
     &                   Resid(i,j,k),Brhs(i,j,k),abs(abserr)
             endif
          endif
        enddo
        enddo
        enddo

        write(*,*) 'difference in residual ', maxerr

        return
        end subroutine test_symtrd
        end module symtrd_mod
