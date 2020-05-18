	subroutine eqtodyn(adyn,aeq,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: adyn,aeq
		integer :: le,l
		
!		Equilibrium variables are transformed to dynamic variables		
		
		adyn=c1*adyn
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			adyn(:,l)=adyn(:,l)+c2*aeq(:,le) 
		end do
		adyn(:,0)=0. 
		aeq(:,0)=0. 

	end subroutine eqtodyn

	subroutine dyntoeq(aeq,adyn,c1,c2)

		use param
		use domain
		implicit none

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: aeq,adyn
		integer :: le,l

!		Dynamic variables are transformed to equilibrium variables				
		
		aeq=c1*aeq
		do le=1,leqmax
			l=ll(mmeq(le),nneq(le))
			if (l == 0) cycle
			aeq(:,le)=aeq(:,le)+c2*adyn(:,l) 
		end do
		aeq(:,0)=0. 
		adyn(:,0)=0. 

	end subroutine dyntoeq