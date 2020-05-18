	subroutine cnvt(idir)

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: idir,j,l,lskp

		select case (idir)
			case (1)
				do j=1,mj
					do l=1,lmaxn
						xt(l,j)=psi(j,lln(l))
					end do
				end do
				lskp=lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=phi(j,lln(l))
					end do
				end do
				lskp=2*lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=pr(j,lln(l))
					end do
				end do
				lskp=3*lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=uzt(j,lln(l))
					end do
				end do
				lskp=4*lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=nf(j,lln(l))
					end do
				end do
				lskp=5*lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=vprlf(j,lln(l))
					end do
				end do
				lskp=6*lmaxn
				do j=1,mj
					do l=1,lmaxn
						xt(l+lskp,j)=vthprlf(j,lln(l))
					end do
				end do
				if (alpha_on == 1) then
					lskp=7*lmaxn
					do j=1,mj
						do l=1,lmaxn
							xt(l+lskp,j)=nalp(j,lln(l))
						end do
					end do
					lskp=8*lmaxn
					do j=1,mj
						do l=1,lmaxn
							xt(l+lskp,j)=vprlalp(j,lln(l))
						end do
					end do
				end if
			case (2)
				do j=1,mj
					do l=1,lmaxn
						psi(j,lln(l))=xt(l,j)
					end do
				end do
				do l=1,lmaxn
					if (mm(lln(l)) /= 0) cycle
					psi(0,lln(l))=(r(2)**2*psi(1,lln(l))-r(1)**2*psi(2,lln(l)))/(r(2)**2-r(1)**2)
					call mult(sceq1,test,-1,psi,1,0.0_IDP,1.0_IDP)
					psi(mj,lln(l))=-((r(mj)-r(mjm1))*sceq1(mj,lln(l))*mm(l)/r(mj)) - psi(mjm1,lln(l))
				end do
				lskp=lmaxn
				do j=1,mj
					do l=1,lmaxn
						phi(j,lln(l))=xt(l+lskp,j)
					end do
				end do
				do l=1,lmaxn
					if (mm(lln(l)) /= 0) cycle
					phi(0,lln(l))=(r(2)**2*phi(1,lln(l))-r(1)**2*phi(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				lskp=2*lmaxn
				do j=1,mj
					do l=1,lmaxn
						pr(j,lln(l))=xt(l+lskp,j)
					end do

				end do
				do l=1,lmaxn
					if (mm(lln(l)) /= 0) cycle
					pr(0,lln(l))=(r(2)**2*pr(1,lln(l))-r(1)**2*pr(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				lskp=3*lmaxn
				do j=1,mj
					do l=1,lmaxn
						uzt(j,lln(l))=xt(l+lskp,j)
					end do
				end do
				do l=1,lmaxn
					if (mm(lln(l)) /= 0) cycle
					uzt(0,lln(l))=(r(2)**2*uzt(1,lln(l))-r(1)**2*uzt(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				lskp=4*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						nf(j,lln(l))=xt(l+lskp,j)
					end do
				end do
				do l=1,lmaxn
					nf(mj,lln(l))=(nf(mjm1,lln(l))*(r(mj)-r(mjm2))**2-nf(mjm2,lln(l))*(r(mj)-r(mjm1))**2)/ &
						      ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
					if (mm(lln(l)) /= 0) cycle
					nf(0,lln(l))=(r(2)**2*nf(1,lln(l))-r(1)**2*nf(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				lskp=5*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						vprlf(j,lln(l))=xt(l+lskp,j)
					end do
				end do
				do l=1,lmaxn
					vprlf(mj,lln(l))=(vprlf(mjm1,lln(l))*(r(mj)-r(mjm2))**2-vprlf(mjm2,lln(l))*(r(mj)-r(mjm1))**2)/ &
							 ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
					if (mm(lln(l)) /= 0) cycle
					vprlf(0,lln(l))=(r(2)**2*vprlf(1,lln(l))-r(1)**2*vprlf(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				lskp=6*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						vthprlf(j,lln(l))=xt(l+lskp,j)
					end do
				end do
				do l=1,lmaxn
					vthprlf(mj,lln(l))=(vthprlf(mjm1,lln(l))*(r(mj)-r(mjm2))**2-vthprlf(mjm2,lln(l))*(r(mj)-r(mjm1))**2)/ &
							   ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
					if (mm(lln(l)) /= 0) cycle
					vthprlf(0,lln(l))=(r(2)**2*vthprlf(1,lln(l))-r(1)**2*vthprlf(2,lln(l)))/(r(2)**2-r(1)**2)
				end do
				if (alpha_on == 1) then
					lskp=7*lmaxn
					do j=1,mjm1
						do l=1,lmaxn
							nalp(j,lln(l))=xt(l+lskp,j)
						end do
					end do
					do l=1,lmaxn
						nalp(mj,lln(l))=(nalp(mjm1,lln(l))*(r(mj)-r(mjm2))**2-nalp(mjm2,lln(l))*(r(mj)-r(mjm1))**2)/ &
								((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
						if (mm(lln(l)) /= 0) cycle
						nalp(0,lln(l))=(r(2)**2*nalp(1,lln(l))-r(1)**2*nalp(2,lln(l)))/(r(2)**2-r(1)**2)
					end do
					lskp=8*lmaxn
					do j=1,mjm1
						do l=1,lmaxn
							vprlalp(j,lln(l))=xt(l+lskp,j)
						end do
					end do
					do l=1,lmaxn
						vprlalp(mj,lln(l))=(vprlalp(mjm1,lln(l))*(r(mj)-r(mjm2))**2-vprlalp(mjm2,lln(l))*(r(mj)-r(mjm1))**2)/ &
								   ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
						if (mm(lln(l)) /= 0) cycle
						vprlalp(0,lln(l))=(r(2)**2*vprlalp(1,lln(l))-r(1)**2*vprlalp(2,lln(l)))/(r(2)**2-r(1)**2)
					end do
				end if
		end select

	end subroutine cnvt