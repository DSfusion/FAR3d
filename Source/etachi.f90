	subroutine etachi

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: mjp,j
		real(IDP) :: ezw,etanrm,etar,aaeta
		real(IDP), dimension(:), allocatable :: aa,bb,cc,dd

		interface
			subroutine spline(n,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP), dimension(:) :: x,y,b,c,d
			end subroutine spline
			function seval(n,u,x,y,b,c,d)
				use param
				implicit none
				integer :: n
				real(IDP) :: seval,u
				real(IDP), dimension(:) :: x,y,b,c,d
			end function seval
		end interface

!		Set up the magnetic diffusivity radial profile		
		
		mjp=mj+1
		select case (ietaeq)
			case(1)
				if(ext_prof .eq. 0) then
					eta=1.0_IDP/teeq**1.5_IDP
				endif
			case(2)
				eta=etascl 
			case(3)
				eta=eta0*(1.+(r/reta)**(2.*etalmb)) **(1./etalmb) 
		end select

	end subroutine etachi