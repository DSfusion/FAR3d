	program rdwr

		implicit none
		integer, parameter :: IDP = kind(1.d0)

		integer :: j,l,m,n,jmax,lmax
		integer, dimension(1) :: lmx,jmx
		real(IDP) :: theta,phimax
		real(IDP), dimension(1000) :: rp
		real(IDP), dimension(199) :: r
		real(IDP), dimension(199,24) :: phi
		real(IDP), dimension(199,12) :: phisq,phir,phii
		real(IDP), dimension(12) :: phismx
		character(len=1) :: tb
		character(len=32) :: confil

		open(unit=97,file='egn_vectors4_FR=0.75_n=3.txt',status="old")
		read(97,*) confil
		do j=1,199
			read(97,*) r(j),(phi(j,l),l=1,24)
		end do
		close(unit=97)
		do l=1,12
			do j=1,199
				phisq(j,l)=phi(j,l)*phi(j,l)+phi(j,l+12)*phi(j,l+12)
			end do
			phismx(l)=maxval(phisq(:,l))
		end do
		lmx=maxloc(phismx)
		lmax=lmx(1)
		jmx=maxloc(phisq(:,lmax))
		jmax=jmx(1)
!		print *," lmax=",lmax," jmax=",jmax," phisq=",phisq(jmax,lmax)


		theta=atan2(phi(jmax,lmax),phi(jmax,lmax+12))
		phimax=phi(jmax,lmax)*sin(theta)+phi(jmax,lmax+12)*cos(theta)
!		print *," theta=",180*theta/3.141592654," phimax=",phimax

		do l=1,12
			do j=1,199
				phir(j,l)=phi(j,l)*sin(theta)+phi(j,l+12)*cos(theta)
				phii(j,l)=-phi(j,l)*cos(theta)+phi(j,l+12)*sin(theta)
			end do
		end do

		n=3
		tb=char(9)
		open(unit=97,file='complex_egn_vectors.txt',status="new")
		write(97,'("r",24(a1,i3,"/",i1))') (tb,m,n,m=15,4,-1),(tb,m,n,m=15,4,-1) 
		do j=1,199
			write(97,'(f5.3,24(a1,1pe13.6))') r(j),(tb,phir(j,l),l=1,12),(tb,phii(j,l),l=1,12)
		end do
		close(unit=97)

	end program rdwr
