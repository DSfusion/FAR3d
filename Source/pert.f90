	subroutine pert

		use param
		use domain
		use equil
		use dynamo
		use cotrol
                use ffunctions
		implicit none

		integer :: l,l1,lp,iamm,iann,j,jj,jm,jp,ndivi,n,k
		integer, dimension(ldim) :: jsl
		real(IDP) :: pi,sq2,teps,xran,yran,p,sigwid,qp,psipp,beta,rbeta,awid,psitil,phitil,yprb,xx,psij,phij,r32,y,rr,dr,gam, &
						 rmrbar,rsl
		integer, dimension(8) :: values_r
		integer, dimension(:), allocatable :: i_r
		real(IDP) :: xle=0.01_IDP

!		Set up the perturbation over the equilibria

		if (pertscl /= 1.) then
			psi=psi*pertscl
			phi=phi*pertscl
			pr=pr*pertscl
		endif

		pi=4.0_IDP*atan(1.0_IDP)
		sq2=sqrt(2.0_IDP)
		teps = 0.025
		call date_and_time(values=values_r)
		call random_seed(size=k)
		allocate (i_r(k))
		call random_seed(get=i_r)

		do l1=1,lmaxn
			l=lln(l1)
			lp=l
			if (lo(l1) /= 0) lp=lln(lo(l1))
			if (signl(l) > 0 .or. lo(l1) == 0) then
				call random_number(xran)
				call random_number(yran)
				if (ipert == 0) then
					gammai(l)=gammai(l)*xran*cos(2.*yran*pi)
					if (lo(l1) /= 0) gammai(lp)=gammai(lp)*xran*sin(2.*yran*pi)
				else
					widthi(l)=widthi(l)*xran*cos(2.*yran*pi)
					if (lo(l1) /= 0) widthi(lp)=widthi(lp)*xran*sin(2.*yran*pi)
				end if
			end if

!  		compute rs

			iamm=mm(l)
			iann=nn(l)
			jsl(l)=0
			p=0.
			if (iann /= 0) p=iamm*1.0_IDP/iann
			do j=1,mj
				if (((1./qqinv(j-1))-p)*((1./qqinv(j))-p) <= 0.) exit
			end do
			if (j > mj) then
				rs(l)=.5
				jsl(l)=mj/2
			else
				rs(l)=r(j-1)+(r(j)-r(j-1))*((1./qqinv(j-1))-p)/((1./qqinv(j-1))-(1./qqinv(j)))
				jsl(l)=j
			end if

			if (widthi(l) == 0) cycle
			sigwid=sign(1.0_IDP,widthi(l))


			select case (ipert)
				case (0)

!  			old perturbation

					do j=1,mj
						if (rs(l) < r(j)) exit
					end do
					if (j > mj) j = mj
					qp=((1./qqinv(j))-(1./qqinv(j-1)))/(r(j)-r(j-1))
					psipp=rs(l)*qp*qqinv(j)**2
					if(iamm == 1 .and. iann == 1) then
						beta=-3.
						awid=abs(widthi(l))
						rbeta=beta*awid+rs(l)
						phitil=-abs(psipp)*s*sqrt(2.*pi)*widthi(l)**2/(32.*rs(l))
						psitil=-abs(psipp)*(widthi(l)**2/16.)*(exp(-beta*beta*.5)-.5*sqrt(2.*pi)*beta*(1.-erf(beta/sq2)))
						do j=1,mjm1
							if (rbeta < r(j)) exit
						end do
						yprb=rbeta*(-1.+1.*qqinv(j))
						do j=1,mjm1
							xx=(r(j)-rs(l))/awid
							phij=phitil*r(j)*(1.-erf(xx/sq2))
							if (r(j) < rbeta) then
								psij=(r(j)*(-1.+qqinv(j))/yprb)*psitil
							else
								psij=-abs(psipp)*(widthi(l)**2/16.)*(exp(-xx*xx*.5)-.5*sqrt(2.*pi)*xx*(1.-erf(xx/sq2)))
							end if
							phi(j,l)=phi(j,l)+phij*sigwid
						end do
					else
!						guess a perturbation.
						psitil=-(widthi(l)*.25)**2*psipp*sigwid/(rs(l)**iamm*(1.-rs(l)))
						do j=1,mjm1
							psij=psitil*r(j)**iamm*(1.-r(j))*2./(1.+exp(10.*(-1.+r(j)/rs(l))))
							phi(j,l)=phi(j,l)+psij*gammai(l)*(rs(l)-r(j))/(rs(l)+r(j))
						end do
					end if

				case (1)

!  				0/0 perturbation

					if (l == l0) then
						p=1.5
						do j=1,mj
							if(((1./qqinv(j-1))-p)*((1./qqinv(j))-p).le.0.) exit
						end do
						if (j > mj) then
							r32=.5
						else
							r32=r(j-1)+(r(j)-r(j-1))*((1./qqinv(j-1))-p)/((1./qqinv(j-1))-(1./qqinv(j)))
						end if
						phitil=0.
						phi(mj,l) = phi(mj,l) + phitil
						y=phitil
						ndivi = 1000
						do jj=1,mj
							jm=mj-jj
							jp=jm+1
							dr=(r(jp)-r(jm))/ndivi
							rr=r(jp)+.5*dr
							do n=1,ndivi
								rr=rr-dr
								y=y-dr*rr*tanh((rr-r32)/xle)
							end do
							phitil=widthi(l)*y/r32
							phi(jm,l) = phi(jm,l) + phitil
						end do
					else

!  				new perturbation

						rsl = rs(l)
						gam = teps
						if(iamm > 0) gam = teps*iamm**(-1./3.)
						do j = 1,mjm1
							rmrbar = ((r(j)-rsl)/gam)**2

!  						set perturbations

							phitil = widthi(l)*exp(-0.5*rmrbar)
							psitil = 0.0

!  						perturb phi and psi

							phi(j,l) = phi(j,l) + phitil
						end do
					end if

			end select

		end do

	end subroutine pert
