	subroutine clgam(rgam,na,nb,nc,c1,c2)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: na,nb,nc
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: rgam

		interface
			subroutine gamc(rgam,agam,bgam,itb,n,c1,c2,k)
				use param
				implicit none
				integer :: itb,n,k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: rgam,agam,bgam
			end subroutine gamc
		end interface

!		Operator in linstart and lincheck		
		
		if (nc == 0) then
			select case (na)
				case (1)
					select case (nb)
						case (1)
							call gamc(rgam,grt,grr,1,na,c1,c2,2)			
						case (2)
							call gamc(rgam,gtt,grt,-1,na,c1,c2,2)
					end select
			end select
		else if (nc == 1) then
			select case (na)
				case (1)
					select case (nb)
						case (1)
							call gamc(rgam,grtoj,grroj,1,na,c1,c2,2)
						case (2)
							call gamc(rgam,gttoj,grtoj,-1,na,c1,c2,2)
					end select
			end select
		else if(nc == 2) then
			select case (na)
				case (1)
					select case (nb)
						case (1)
							call gamc(rgam,jbgrt,jbgrr,1,na,c1,c2,2)
						case (2)
							call gamc(rgam,jbgtt,jbgrt,-1,na,c1,c2,2)
					end select
			end select
		end if

	end subroutine clgam

	subroutine gamc(rgam,agam,bgam,itb,n,c1,c2,k)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: itb,n,k
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: rgam,agam,bgam

		interface
			subroutine dbydtheq(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydtheq
			subroutine dbydreq(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydreq
			subroutine dbydzteq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzteq
		end interface

!		Operator in linstart and lincheck			
		
		if (n > 2) then
			call dbydtheq(rgam,agam,itb,c1,c2,k)
		else
			call dbydreq(rgam,agam,c1,c2,k)
		end if
		if (n > 1) then
			call dbydzteq(rgam,bgam,itb,1.0_IDP,-c2)
		else
			call dbydtheq(rgam,bgam,itb,1.0_IDP,-c2,k)
		end if

	end subroutine gamc