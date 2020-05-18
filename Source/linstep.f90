	subroutine linstep

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: l

		interface
			subroutine b2lx0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine b2lx0
			subroutine b2lx(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine b2lx
			subroutine clgam(rgam,na,nb,nc,c1,c2)
				use param
				implicit none
				integer :: na,nb,nc
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: rgam
			end subroutine clgam
			subroutine solbt(m,n,a,b,c,y,ip)
				use param
				implicit none
				integer :: m,n
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
				real(IDP), dimension(m,n) :: y
			end subroutine solbt
			subroutine cnvt(idir)
				implicit none
				integer :: idir
			end subroutine cnvt
		end interface

!		In this subroutine the simulation is advanced in time 		

		yt=-yt

!  put in l.h.s. of equation


!  psi equation

		sd1=2.0_IDP
		call b2lx0(sd1,1,1,0,0,0,1.0_IDP)

!  u-zeta equation

		call b2lx0(sd1,2,4,0,0,0,1.0_IDP)

! p equation

		call b2lx0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density equation

		call b2lx0(sd1,5,5,0,0,0,1.0_IDP)

! fast ion v-parallel equation

		call b2lx0(sd1,6,6,0,0,0,1.0_IDP)
		
! thermal moment equation

		call b2lx0(sd1,7,7,0,0,0,1.0_IDP)		

        if(alpha_on .eq. 1) then 
		
! fast alpha density equation

		  call b2lx0(sd1,8,8,0,0,0,1.0_IDP)

! fast alpha v-parallel equation

		  call b2lx0(sd1,9,9,0,0,0,1.0_IDP)	

        end if		  

!		The subroutine solbt inverts the matrix		
		
		xt(:,1:mjm1)=yt(:,1:mjm1)
		call solbt(lmx,mjm1,amat,bmat,cmat,xt,ipc)

!  	transfer time advanced values to original arrays
		
		call cnvt(2)

	end subroutine linstep