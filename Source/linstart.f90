	subroutine linstart

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		use scratch
		implicit none

		integer :: j,l,m,n,lp,l1,l2,l3,le,it1,it2,it3,m1,n1,is1,m2,n2,is2,m3,n3,is3,l3p,mp3,np3,mm3,nm3,sgn12,lskp,ibnd,ier
		real(IDP) :: epsq,oneos,betfc,betfc_f,betfc_alp,coef,omcyd,omcydalp,beteom,betiom,x,aux,xnuion0,abkprl, &
			     xsii,xsie,xsi2,xsi3,xsi4,ztr,zti,zir,zii,zer,zei,vther,vthir,tauie,xnuelc,xnuion
		complex(IDP) :: zetai,zetae,zetai2,zetai3,zetai4,zi,y0i,y1i,y2i,ddii,zetae2,zetae3,zetae4,ze,y0e,y1e,y2e,ddee,rei,sei, &
				stfe,stfi,reii,sei1,sei2,sei3,zetaiinv,zetaeinv,cmplx1
		real(IDP), dimension(:), allocatable :: xa
		character*1 :: tb
		
		interface
			subroutine cnvt(idir)
				implicit none
				integer :: idir
			end subroutine cnvt
			subroutine clgam(rgam,na,nb,nc,c1,c2)
				use param
				implicit none
				integer :: na,nb,nc
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: rgam
			end subroutine clgam
			subroutine om(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine om
			subroutine blockj(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj
			subroutine blockjl(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP), dimension(0:,0:) :: tx
				real(IDP), dimension(:,:) :: coef
			end subroutine blockjl
			subroutine blockj_landau_grad_parallel(tx,itx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: itx,ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:,0:) :: tx
			end subroutine blockj_landau_grad_parallel
			subroutine om0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine om0
			subroutine block0(tx,ieqn,ivar,ith,ir,izt,coef)
				use param
				implicit none
				integer :: ieqn,ivar,ith,ir,izt
				real(IDP) :: coef
				real(IDP), dimension(0:) :: tx
			end subroutine block0
			subroutine block_dlsq(ieqn,ivar,coef)
				use param
				implicit none
				integer :: ieqn,ivar
				real(IDP) :: coef
			end subroutine block_dlsq
			subroutine block_dlsq_r(ieqn,ivar,coef)
				use param
				implicit none
				integer :: ieqn,ivar
				real(IDP) :: coef
			end subroutine block_dlsq_r
			subroutine mat_out(ic)
				use param
				implicit none
				integer :: ic
			end subroutine mat_out
!			subroutine block0_dlsq(tx,ieqn,ivar,coef)
!				use param
!				implicit none
!				integer :: ieqn,ivar
!				real(IDP) :: coef
!				real(IDP), dimension(0:) :: tx
!			end subroutine block0_dlsq
			subroutine dbydth(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydth
			subroutine dbydzt(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydzt
			subroutine dbydtheq(d,a,ltype,c1,c2,k)
				use param
				implicit none
				integer :: k,ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine dbydtheq
			subroutine grpareq(d,a,ltype,c1,c2)
				use param
				implicit none
				integer :: ltype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: a,d
			end subroutine grpareq
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
			subroutine dbydr0(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine dbydr0
			subroutine d2bydr20(d,a,c1,c2,k)
				use param
				implicit none
				integer :: k
				real(IDP) :: c1,c2
				real(IDP), dimension(0:) :: a,d
			end subroutine d2bydr20
			subroutine decbt(m,n,a,b,c,ip,ier)
				use param
				implicit none
				integer :: m,n,ier
				integer, dimension(m,n) :: ip
				real(IDP), dimension(m,m,n) :: a, b, c
			end subroutine decbt
			subroutine dlstar(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine dlstar
			subroutine dlsq(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine dlsq
			subroutine dlsq_r(ss,ff,itypf,wk1,wk2,wk3,c1,c2)
				use param
				implicit none
				integer :: itypf
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: ss,ff,wk1,wk2,wk3
			end subroutine dlsq_r
			subroutine zzdisp(x,y,zzr,zzi)
				use param
				implicit none
				real(IDP) :: x,y,zzr,zzi
			end subroutine zzdisp
		end interface

		tb = char(9)
		epsq=eps*eps
		oneos=0.0_IDP
		if (s > 0.0_IDP) oneos=1.0_IDP/s
		betfc=bet0/(2.*epsq)
		betfc_f=LcA2*bet0_f/(2.*epsq)
		omcyd = omcy
		beteom=dpres*bet0/(2*epsq*omcyd)
		betiom=(1-dpres)*bet0/(2*epsq*omcyd)

		call dbydr0(qqinvp,qqinv,0.0_IDP,1.0_IDP,0)
		call dbydr0(denseqr,denseq,0.0_IDP,1.0_IDP,0)
		call dbydr0(dnfeqdr,nfeq,0.0_IDP,1.0_IDP,0)
		vfova2 = vfova*vfova
		if (alpha_on == 1) then
			betfc_alp=LcA2alp*bet0_alp/(2.*epsq)
			omcydalp = omcyalp
			valphaova2=valphaova*valphaova
			call dbydr0(dnalpeqdr,nalpeq,0.0_IDP,1.0_IDP,0)
		end if

!		Output files with the simulation main profiles
		
		call dbydreq(sceq1,sqg,0.0_IDP,1.0_IDP,0)
		sceq2=sqg
		sd2=qqinvp*r/qqinv

		if(ext_prof .eq. 0) then

		  open(unit=44,file="profiles.dat",status="unknown")
		  write(44,'("        r",a1,"      denseq",a1,"       teeq",a1,"       nfeq",a1,"     dnfeqdr",a1, &
		           "      vfova",a1,"      cureq",a1,"       feq",a1,"    pressure",a1,"       iota",a1, &
                           "        1/D",a1,"    curvature",a1,"      shear",a1,"     eta")') &                         
			   tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
		  do j=0,mj
		  write(44,'(e15.7,13(a1,e15.7))') r(j),tb,denseq(j),tb,teeq(j),tb,nfeq(j),tb,dnfeqdr(j), &
		          tb,vfova(j),tb,cureq(j),tb,feq(j),tb,preq(j),tb,qqinv(j),tb,sceq2(j,leq0), &
                          tb,sceq1(j,leq0),tb,sd2(j),tb,eta(j)
		  end do
		  close(unit=44)

		end if
		
		if(ext_prof .eq. 1 .and. alpha_on .eq. 0) then
	  	  open(unit=45,file="profiles.dat",status="unknown")
		  write(45,'("        r",a1,"      dnnbi",a1,"       dne",a1,"       dni",a1,"      tbn",a1, &
		             "        ti",a1,"        te",a1,"      vfova",a1,"      cureq",a1, &
                             "       feq",a1,"     pressure",a1,"       iota",a1,"        q",a1,"       1/D",a1, &
                             "    curvature",a1,"      shear",a1,"        eta")') &
		  	     tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
		  do j=0,mj
		   write(45,'(e15.7,16(a1,e15.7))') r(j),tb,dnnbi(j),tb,dne(j),tb,dni(j),tb,tbn(j), &
		          tb,ti(j),tb,te(j),tb,vfova(j),tb,cureq(j),tb,feq(j),tb,preq(j),tb,qqinv(j), &
			  tb,qq(j),tb,sceq2(j,leq0),tb,sceq1(j,leq0),tb,sd2(j),tb,eta(j)
		  end do
		  close(unit=45)

	  	  open(unit=46,file="profiles_ex.dat",status="unknown")
		  write(46,'("       r",a1,"   dnnbi(m^-3)",a1,"    dne(m^-3)",a1,"    dni(m^-3)",a1,"     tbn(keV)",a1, &
		             "     ti(keV)",a1,"     te(keV)",a1,"   vzt_eq(m/s)",a1,"   vth_eq(m/s)",a1,"  vthermalep(m/s)",a1,"       etann(MKS)")') &
		  	     tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
		  do j=0,mj
		   write(46,'(e15.7,10(a1,e15.7))') r(j),tb,dnnbinn(j),tb,dnenn(j),tb,dninn(j),tb,tbnnn(j), &
		          tb,tinn(j),tb,tenn(j),tb,vzt_eq(j),tb,vth_eq(j),tb,vthermalep(j),tb,eta(j)
		  end do
		  close(unit=46)
        end if		
		
		if(ext_prof .eq. 1 .and. alpha_on .eq. 1) then
	  	  open(unit=47,file="profiles.dat",status="unknown")
		  write(47,'("       r",a1,"      dnnbi",a1,"       dne",a1,"        dni",a1,"     dnalpha",a1, &
		             "      tbn",a1,"      ti",a1,"       te",a1,"      talpha",a1, &
                             "       vfova",a1,"      cureq",a1,"       feq",a1,"     pressure",a1, &
                             "       iota",a1,"        q",a1,"       1/D",a1, &
                             "          curvature",a1,"      shear",a1,"           eta")') &
		  	     tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
		  do j=0,mj
		   write(47,'(e15.7,18(a1,e15.7))') r(j),tb,dnnbi(j),tb,dne(j),tb,dni(j),tb,dnalpha(j),tb,tbn(j), &
		          tb,ti(j),tb,te(j),tb,talpha(j),tb,vfova(j),tb,cureq(j),tb,feq(j), &
                          tb,preq(j),tb,qqinv(j),tb,qq(j),tb,sceq2(j,leq0),tb,sceq1(j,leq0),tb,sd2(j),tb,eta(j)
		  end do
		  close(unit=47)

	  	  open(unit=48,file="profiles_ex.dat",status="unknown")
		  write(48,'("        r",a1,"        dnnbi(m^-3)",a1,"      dne(m^-3)",a1,"      dni(m^-3)",a1, &
		             "     dnalpha(m^-3)",a1,"     tbn(keV)",a1, "       ti(keV)",a1,"       te(keV)",a1, &
                             "      talpha(keV)",a1,"     vzt_eq(m/s)",a1,"     vth_eq(m/s)",a1,"     vthermalep(m/s)",a1,"       eta(MKS)")') &
		  	     tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
		  do j=0,mj
		   write(48,'(e15.7,9(a1,e15.7))') r(j),tb,dnnbinn(j),tb,dnenn(j),tb,dninn(j),tb,dnalphann(j),tb,tbnnn(j), &
		          tb,tinn(j),tb,tenn(j),tb,talphann(j),tb,vzt_eq(j),tb,vth_eq(j),tb,vthermalep(j),tb,etann(j)
		  end do
		  close(unit=48)
        end if		

!  calculate derivative weights to be used in blockj, block0, b2lx and b2lx0

		wt1m(1,1)=0.
		wt2m(1,1)=0.
		wt1m(1,2)=0.
		wt2m(1,2)=0.
		wt10(1,1)=-dc1m(1)-dc1p(1)
		wt20(1,1)=-dc2m(1)-dc2p(1)
		wt10(1,2)=r(2)**2*dc1m(1)/(r(2)**2-r(1)**2)-dc1m(1)-dc1p(1)
		wt20(1,2)=r(2)**2*dc2m(1)/(r(2)**2-r(1)**2)-dc2m(1)-dc2p(1)
		wt1p(1,1)=dc1p(1)
		wt2p(1,1)=dc2p(1)
		wt1p(1,2)=-r(1)**2*dc1m(1)/(r(2)**2-r(1)**2)+dc1p(1)
		wt2p(1,2)=-r(1)**2*dc2m(1)/(r(2)**2-r(1)**2)+dc2p(1)
		wt1m(2:mjm1,1)=dc1m(2:mjm1)
		wt2m(2:mjm1,1)=dc2m(2:mjm1)
		wt10(2:mjm1,1)=-dc1m(2:mjm1)-dc1p(2:mjm1)
		wt20(2:mjm1,1)=-dc2m(2:mjm1)-dc2p(2:mjm1)
		wt1p(2:mjm2,1)=dc1p(2:mjm2)
		wt2p(2:mjm2,1)=dc2p(2:mjm2)
		wt1p(mjm1,1)=0.
		wt2p(mjm1,1)=0.
		wt1m(2:mjm2,2)=dc1m(2:mjm2)
		wt2m(2:mjm2,2)=dc2m(2:mjm2)
		wt1m(mjm1,2)=dc1m(mjm1)-dc1p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt2m(mjm1,2)=dc2m(mjm1)-dc2p(mjm1)*(r(mj)-r(mjm1))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt10(2:mjm2,2)=-dc1m(2:mjm2)-dc1p(2:mjm2)
		wt20(2:mjm2,2)=-dc2m(2:mjm2)-dc2p(2:mjm2)
		wt10(mjm1,2)=-dc1m(mjm1)-dc1p(mjm1)+dc1p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt20(mjm1,2)=-dc2m(mjm1)-dc2p(mjm1)+dc2p(mjm1)*(r(mj)-r(mjm2))**2/((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
		wt1p(2:mjm2,2)=dc1p(2:mjm2)
		wt2p(2:mjm2,2)=dc2p(2:mjm2)
		wt1p(mjm1,2)=0.
		wt2p(mjm1,2)=0.

		allocate (cmamm(lmaxn,leqmax,lmaxn),cmamp(lmaxn,leqmax,lmaxn),cmapm(lmaxn,leqmax,lmaxn),cmapp(lmaxn,leqmax,lmaxn))

		cmamm=0.0_IDP
		cmamp=0.0_IDP
		cmapm=0.0_IDP
		cmapp=0.0_IDP

		allocate (amat(lmx,lmx,mjm1),bmat(lmx,lmx,mjm1),cmat(lmx,lmx,mjm1))

		amat=0.0_IDP
		bmat=0.0_IDP
		cmat=0.0_IDP

		allocate (ipc(lmx,mjm1))
		allocate (xt(lmx,mj))
		allocate (yt(lmx,mjm1))

		ipc=0
		xt=0.0_IDP
		yt=0.0_IDP

!  read start or continue values into proper arrays

		call dlstar(uzt,phi,-1,sc1,sc2,sc3,0.0_IDP,1.0_IDP)
		do l=1,lmax
			uzt(mj,l)=(uzt(mjm1,l)*(r(mj)-r(mjm2))-uzt(mjm2,l)*(r(mj)-r(mjm1)))/(r(mjm1)-r(mjm2))
		end do
		call cnvt(1)

!  calculate "c" coefficients to be used in blockj, b2lx and b2lx0

		allocate (llno(lmax))
		do l=1,lmax
			do l1=1,lmaxn
				if(lln(l1) == l) exit
			end do
			if (l1 > lmaxn) then
				llno(l)=0
			else
				llno(l)=l1
			end if
		end do

		do le=1,leqmax
			m1=mmeq(le)*sgnleq(le)
			n1=nneq(le)*sgnleq(le)

			do l2=1,lmaxn
				l=lln(l2)
				m2=mm(l)*signl(l)
				n2=nn(l)*signl(l)

				if (sgnleq(le) == 0) then

					cmapp(l2,le,l2)=1.0_IDP
					cmapm(l2,le,l2)=1.0_IDP
					l3=lo(l2)
					if (l3 > 0) then
						cmamp(l3,le,l2)=1.0_IDP
						cmamm(l3,le,l2)=1.0_IDP
					end if

				else if (signl(l) == 0) then

					do it2=-1,1,2
						m3=it2*mmeq(le)
						n3=it2*nneq(le)
						if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) then
							l3p=ll(m3,n3)
							if (l3p > 0) then
								l3=llno(l3p)
								if (l3 > 0) then
									if (it2 == 1) then
										cmapp(l3,le,l2)=1.0_IDP
										cmamp(l3,le,l2)=1.0_IDP
									else
										cmapm(l3,le,l2)=1.0_IDP
										cmamm(l3,le,l2)=1.0_IDP
									end if
								end if
							end if
						end if
					end do

				else

					do it1=-1,1,2
						do it2=-1,1,2
							it3=it1*it2
							is1=it1*sgnleq(le)
							is2=it2*signl(l)
							is3=is1*is2
							mp3=m1+m2
							np3=n1+n2
							m3=is3*it3*mp3
							n3=is3*it3*np3
							l3p=0
							if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) then
								l3p=ll(m3,n3)
								if (l3p > 0) then
									l3=llno(l3p)
									if (l3 > 0) then
										coef=0.5_IDP
										if (is1 == -1 .and. is2 == -1) coef=-0.5_IDP
										if (it2 == 1) then
											if (it1 == 1) then
												cmapp(l3,le,l2)=coef
											else
												cmamp(l3,le,l2)=coef
											end if
										else
											if (it1 == 1) then
												cmapm(l3,le,l2)=coef
											else
												cmamm(l3,le,l2)=coef
											end if
										end if
									end if
								end if
							end if

							if (n1 < n2) then
								sgn12=-1
								mm3=m2-m1
								nm3=n2-n1
							else if (n1 > n2) then
								sgn12=1
								mm3=m1-m2
								nm3=n1-n2
							else
								sgn12=sign(1,m1-m2)
								mm3=abs(m1-m2)
								nm3=0
							end if
							m3=is3*it3*mm3
							n3=is3*it3*nm3
							l3p=0
							if (m3 >= mmin .and. m3 <= mmax .and. n3 >= nmin .and. n3 <= nmax) then
								l3p=ll(m3,n3)
								if (l3p > 0) then
									l3=llno(l3p)
									if (l3 > 0) then
										coef=0.5_IDP
										if (is3 == -1) then
											if (is1*sgn12 == 1) coef=-0.5_IDP
											if (mm3 == 0 .and. nm3 == 0) coef=0.0_IDP
										end if
										if (it2 == 1) then
											if (it1 == 1) then
												cmapp(l3,le,l2)=coef
											else
												cmamp(l3,le,l2)=coef
											end if
										else
											if (it1 == 1) then
												cmapm(l3,le,l2)=coef
											else
												cmamm(l3,le,l2)=coef
											end if
										end if
									end if
								end if
							end if
						end do
					end do

				end if
			end do
		end do

!  compute rest of values into proper xt-arrays

!   Ion FLR effects auxiliary variable

		if (iflr_on == 1) then
!			call dlsq(sc1,psi,1,sc2,sc3,sc4,0.0_IDP,1.0_IDP)
			call dlsq_r(sc1,psi,1,sc2,sc3,sc4,0.0_IDP,1.0_IDP)
			lskp=(iq-1)*lmaxn
			do j=1,mjm1
				do l=1,lmaxn
					xt(l+lskp,j)=sc1(j,lln(l))
				end do
			end do
		end if

!   EP FLR effects auxiliary variables

		if (epflr_on == 1) then
			allocate (amatw(lmaxn,lmaxn,mjm1),bmatw(lmaxn,lmaxn,mjm1),cmatw(lmaxn,lmaxn,mjm1))
			allocate (xa(0:mj))
			allocate (xw(lmaxn,mjm1))
			allocate (ipcw(lmaxn,mjm1))

			amatw=0.0_IDP
			bmatw=0.0_IDP
			cmatw=0.0_IDP

			do l=1,lmaxn
				m=mm(lln(l))
				n=nn(lln(l))
				ibnd=1
				if (m == 0) ibnd=2
				do lp=1,lmaxn
					do leq=1,leqmax
						x=cmapm(lp,leq,l)
						if (x == 0.0_IDP) cycle
						xa=x
						do j=1,mjm2
							amatw(lp,l,j)=amatw(lp,l,j)+wt20(j,ibnd)*xa(j)*lplrr(j,leq)
							cmatw(lp,l,j)=cmatw(lp,l,j)+wt2m(j,ibnd)*xa(j)*lplrr(j,leq)
							bmatw(lp,l,j)=bmatw(lp,l,j)+wt2p(j,ibnd)*xa(j)*lplrr(j,leq)
						end do
						amatw(lp,l,mjm1)=amatw(lp,l,mjm1)+wt20(mjm1,1)*xa(mjm1)*lplrr(mjm1,leq)
						cmatw(lp,l,mjm1)=cmatw(lp,l,mjm1)+wt2m(mjm1,1)*xa(mjm1)*lplrr(mjm1,leq)
						bmatw(lp,l,mjm1)=bmatw(lp,l,mjm1)+wt2p(mjm1,1)*xa(mjm1)*lplrr(mjm1,leq)
						do j=1,mjm2
							amatw(lp,l,j)=amatw(lp,l,j)+wt10(j,ibnd)*xa(j)*lplr(j,leq)
							cmatw(lp,l,j)=cmatw(lp,l,j)+wt1m(j,ibnd)*xa(j)*lplr(j,leq)
							bmatw(lp,l,j)=bmatw(lp,l,j)+wt1p(j,ibnd)*xa(j)*lplr(j,leq)
						end do
						amatw(lp,l,mjm1)=amatw(lp,l,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplr(mjm1,leq)
						cmatw(lp,l,mjm1)=cmatw(lp,l,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplr(mjm1,leq)
						bmatw(lp,l,mjm1)=bmatw(lp,l,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplr(mjm1,leq)
						xa=-x*m*m*rinv*rinv
						do j=1,mjm1
							amatw(lp,l,j)=amatw(lp,l,j)+xa(j)*lpltt(j,leq)
						end do
						xa=-x*m*n*rinv
						do j=1,mjm1
							amatw(lp,l,j)=amatw(lp,l,j)+xa(j)*lpltz(j,leq)
						end do
						xa=-x*n*n
						do j=1,mjm1
							amatw(lp,l,j)=amatw(lp,l,j)+xa(j)*lplzz(j,leq)
						end do
					end do
					do leq=1,leqmax
						x=cmamp(lp,leq,l)
						if (x == 0.0_IDP) cycle
						xa=x*m*rinv
						do j=1,mjm2
							amatw(lp,l,j)=amatw(lp,l,j)+wt10(j,ibnd)*xa(j)*lplrt(j,leq)
							cmatw(lp,l,j)=cmatw(lp,l,j)+wt1m(j,ibnd)*xa(j)*lplrt(j,leq)
							bmatw(lp,l,j)=bmatw(lp,l,j)+wt1p(j,ibnd)*xa(j)*lplrt(j,leq)
						end do
						amatw(lp,l,mjm1)=amatw(lp,l,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrt(mjm1,leq)
						cmatw(lp,l,mjm1)=cmatw(lp,l,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrt(mjm1,leq)
						bmatw(lp,l,mjm1)=bmatw(lp,l,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrt(mjm1,leq)
						do j=1,mjm1
							amatw(lp,l,j)=amatw(lp,l,j)+xa(j)*lplt(j,leq)
						end do
						xa=x*n
						do j=1,mjm2
							amatw(lp,l,j)=amatw(lp,l,j)+wt10(j,ibnd)*xa(j)*lplrz(j,leq)
							cmatw(lp,l,j)=cmatw(lp,l,j)+wt1m(j,ibnd)*xa(j)*lplrz(j,leq)
							bmatw(lp,l,j)=bmatw(lp,l,j)+wt1p(j,ibnd)*xa(j)*lplrz(j,leq)
						end do
						amatw(lp,l,mjm1)=amatw(lp,l,mjm1)+wt10(mjm1,1)*xa(mjm1)*lplrz(mjm1,leq)
						cmatw(lp,l,mjm1)=cmatw(lp,l,mjm1)+wt1m(mjm1,1)*xa(mjm1)*lplrz(mjm1,leq)
						bmatw(lp,l,mjm1)=bmatw(lp,l,mjm1)+wt1p(mjm1,1)*xa(mjm1)*lplrz(mjm1,leq)
						do j=1,mjm1
							amatw(lp,l,j)=amatw(lp,l,j)+xa(j)*lplz(j,leq)
						end do
					end do
				end do
			end do
			call dlsq(sc1,phi,-1,sc2,sc3,sc4,0.0_IDP,1.0_IDP)
			call dbydzt(sc2,psi,1,0.0_IDP,1.0_IDP)
			call dbydth(sc4,psi,1,0.0_IDP,1.0_IDP,0)
			do l=1,lmax
				sc3(:,l)=r*sc4(:,l)
			end do

			if (alpha_on == 1 .and. r_epflralp > 0.0) then
				allocate (amatwalp(lmaxn,lmaxn,mjm1),bmatwalp(lmaxn,lmaxn,mjm1),cmatwalp(lmaxn,lmaxn,mjm1))
				allocate (ipcwalp(lmaxn,mjm1))
				coef=-r_epflralp*r_epflralp
				amatwalp=coef*amatw
				bmatwalp=coef*bmatw
				cmatwalp=coef*cmatw
				do l=1,lmaxn
					do j=1,mjm1
						amatwalp(l,l,j)=1.0+amatwalp(l,l,j)
					end do
				end do

				bmatwalp(:,:,mjm1)=0.0_IDP
				cmatwalp(:,:,1)=0.0_IDP

				call decbt(lmaxn,mjm1,amatwalp,bmatwalp,cmatwalp,ipcwalp,ier)
				if (ier /= 0) stop 18

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=coef*sc1(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
				lskp=(iwa-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=sc2(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
				lskp=(ix1a-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=sc3(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatwalp,bmatwalp,cmatwalp,xw,ipcwalp)
				lskp=(ix2a-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do
				deallocate (amatwalp,bmatwalp,cmatwalp,ipcwalp)
			end if

			if (r_epflr > 0.0) then
				coef=-r_epflr*r_epflr
				amatw=coef*amatw
				bmatw=coef*bmatw
				cmatw=coef*cmatw
				do l=1,lmaxn
					do j=1,mjm1
						amatw(l,l,j)=1.0+amatw(l,l,j)
					end do
				end do

				bmatw(:,:,mjm1)=0.0_IDP
				cmatw(:,:,1)=0.0_IDP

				call decbt(lmaxn,mjm1,amatw,bmatw,cmatw,ipcw,ier)
				if (ier /= 0) stop 17

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=coef*sc1(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatw,bmatw,cmatw,xw,ipcw)
				lskp=(iw-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=sc2(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatw,bmatw,cmatw,xw,ipcw)
				lskp=(ix1-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do

				do j=1,mjm1
					do l=1,lmaxn
						xw(l,j)=sc3(j,lln(l))
					end do
				end do
				call solbt(lmaxn,mjm1,amatw,bmatw,cmatw,xw,ipcw)
				lskp=(ix2-1)*lmaxn
				do j=1,mjm1
					do l=1,lmaxn
						xt(l+lskp,j)=xw(l,j)
					end do
				end do
			end if
			deallocate (amatw,bmatw,cmatw,xa,xw,ipcw)
		end if

!  put in r.h.s. of equations

!  psi equation

		sd1=1.0_IDP
		call om0(sd1,1,2,0,0,0,1.0_IDP)
		call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*sceq2(:,l)
		end do
		call blockj(sceq1,-1,1,1,1,0,0,oneos)
		call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*sceq2(:,l)
		end do
		call blockj(sceq1,1,1,1,0,1,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*grroj(:,l)
		end do
		call blockj(sceq1,1,1,1,2,0,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=-2.*eta*feq*grtoj(:,l)
		end do
		call blockj(sceq1,-1,1,1,1,1,0,oneos)
		do l=1,leqmax
			sceq1(:,l)=eta*feq*gttoj(:,l)
		end do
		call blockj(sceq1,1,1,1,0,2,0,oneos)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call blockj(sceq2,1,1,1,0,1,0,oneos)

!   Ion FLR effects

		if (iflr_on == 1) then

			if (ext_prof == 1) then
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*ti*bmod(:,l)/(dni*vtherm_elecP*(feq-qqinv*cureq))
				end do
			else
				do l=1,leqmax
					sceq5(:,l)=1.2533*iflr*iflr*ti*bmod(:,l)/(dni*vtherm_elecP*(feq-qqinv*cureq))
				end do
			end if
			call blockj_landau_grad_parallel(sceq5,1,1,iq,0,0,0,1.0_IDP)

!   Ion FLR effects auxiliary equation

!			call block_dlsq(iq,1,1.0_IDP)
			call block_dlsq_r(iq,1,1.0_IDP)
			sd1=1.0_IDP
			call block0(sd1,iq,iq,0,0,0,-1.0_IDP)

		end if

!   Two fluid terms

		if (twofl_on == 1) then
			sd1=beteom/denseq
			call om0(sd1,1,3,0,0,0,-1.0_IDP)
		end if

!  u-zeta equation

		call dbydreq(sceq1,sqg,0.0_IDP,-1.0_IDP,0)
		call blockj(sceq1,1,2,3,1,0,0,betfc)
		call dbydtheq(sceq2,sqg,1,0.0_IDP,1.0_IDP,0)
		call blockj(sceq2,-1,2,3,0,1,0,betfc)
		write(0,'(/"bet0 = ",e12.4,3x,"betf0_f = ",e12.4)') bet0, bet0_f
		write(0,'("eps = <a>/<R> = ",e12.4,3x,"S = tau_R/tau_A = ",e12.4)') eps, s
		write(0,'("bet0/eps**2 = ",e12.4,3x,"bet0_f/eps**2 = ",e12.4)') betfc, betfc_f

!  fast ion coupling

		call blockj(sceq1,1,2,5,1,0,0,betfc_f)
		call blockj(sceq2,-1,2,5,0,1,0,betfc_f)

		if (alpha_on == 1) then
			call blockj(sceq1,1,2,8,1,0,0,betfc_alp)
			call blockj(sceq2,-1,2,8,0,1,0,betfc_alp)
		end if

!  Sheared equilibrium toroidal flow velocity for u-zeta equation

		call block0(vzt_eq,2,4,0,0,1,-1.0_IDP)

		call clgam(sceq2,1,1,1,0.0_IDP,-1.0_IDP)
		call om(sceq2,-1,2,1,1,0,0,1.0_IDP)
		call grpareq(sceq1,sceq2,-1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,1,2,1,1,0,0,1.0_IDP)
		call clgam(sceq2,1,2,1,0.0_IDP,1.0_IDP)
		call om(sceq2,1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,sceq2,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,grtoj,-1,0.0_IDP,-2.0_IDP)
		call blockj(sceq1,1,2,1,1,1,0,1.0_IDP)
		call grpareq(sceq1,gttoj,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)
		call grpareq(sceq1,grroj,1,0.0_IDP,1.0_IDP)
		call blockj(sceq1,-1,2,1,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=-grtoj(:,l)
		end do
		call om(sceq1,-1,2,1,1,1,0,2.0_IDP*1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=grroj(:,l)
		end do
		call om(sceq1,1,2,1,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=gttoj(:,l)
		end do
		call om(sceq1,1,2,1,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=rinv*sceq1(:,l)
		end do
		call om(sceq2,1,2,1,0,1,0,1.0_IDP)
		call dbydr0(sd1,cureq,0.0_IDP,1.0_IDP,0)
		sd2=rinv*sd1
		call d2bydr20(sd2,cureq,1.0_IDP,-1.0_IDP,0)
		sd1=rinv*sd2/epsq
		call block0(sd1,2,1,1,0,0,1.0_IDP)
		call dbydtheq(sceq2,bst,-1,0.0_IDP,-1.0_IDP,0)
		do l=1,leqmax
			sceq1(:,l)=r*sceq2(:,l)/epsq
		end do
		call dbydreq(sceq2,sceq1,0.0_IDP,-1.0_IDP,0)
		call blockj(sceq2,1,2,1,1,0,0,1.0_IDP)
		call dbydtheq(sceq2,sceq1,1,0.0_IDP,1.0_IDP,0)
		call blockj(sceq2,-1,2,1,0,1,0,1.0_IDP)

!   Ion FLR effects

		if (iflr_on == 1) then
			coef=omegar*iflr*iflr
!			call block_dlsq(2,4,coef)
			call block_dlsq_r(2,4,coef)
		end if

!   Electron-ion Landau damping

		if (ieldamp_on == 1) then   
			coef=-(1-dpres)*bet0/(2.*epsq*omegar)
			cmplx1 = cmplx(0.0_IDP,1.0_IDP)                                               !! complex number i 
			tauie = tinn(0)/tenn(0)
			xnuelc0 = xnuelc0/omegar                                                      !! collision FR normalized to omegar
			xnuion0 = 0.0165*xnuelc0/(tauie*sqrt(uion*tauie))                             !! ion-ion collision FR axis

			allocate (eilnd(mjm1,lmaxn))
			eilnd=0.0_IDP
			do l=1,lmaxn     
				if (lln(l) == l0) cycle
				do j=1,mjm1
					vther = vthe*sqrt(teeq(j))
					vthir = vthi*sqrt(ti(j))
					if(vthir .le. 0.) vthir=.001
					if(vther .le. 0.) vther=.01
					if(teeq(j) .le. 0.) teeq(j)=.01
					tauie = tinn(j)/tenn(j)
					sd1(j) = tauie                                                !! sd1 is tau_i = T_i/T_e
					xnuelc = xnuelc0*dni(j)/(teeq(j)*sqrt(teeq(j)))               !! electron-ion collision FR
					xnuion = xnuion0*dni(j)/(ti(j)*sqrt(ti(j)))                   !! ion-ion collision FR
					abkprl = abs(nn(lln(l)) - mm(lln(l))*qqinv(j))                !! parallel gradient
					if(abkprl .lt. 1.e-50_IDP) then
						zetaiinv = (sqrt(2.0_IDP)*vthir*abkprl)/(omegar*(1. + cmplx1*xnuion))     
						zetaeinv = (sqrt(2.0_IDP)*vther*abkprl)/(omegar*(1. + cmplx1*xnuelc))
						y0i = -(1./(1.+cmplx1*xnuion))*(1. + 0.5*(zetaiinv**2))
						y1i = -(1./(1.+cmplx1*xnuion))*(1. + (zetaiinv**2))
						y2i = -(1.75/(1.+cmplx1*xnuion))*(1. +(23./14.)*(zetaiinv**2))
						y0e = -(1./(1.+cmplx1*xnuelc))*(1. + 0.5*(zetaeinv**2))
						y1e = -(1./(1.+cmplx1*xnuelc))*(1. + (zetaeinv**2))
						y2e = -(1.75/(1.+cmplx1*xnuelc))*(1. +(23./14.)*(zetaeinv**2))
						ddii = 1. + cmplx1*xnuion*y0i
						ddee = 1. + cmplx1*xnuelc*y0e
					else
						xsii = omegar/(sqrt(2.)*vthir*abkprl)
						xsie = xsii*vthir/vther
						xsi2 = xsii*xsii
						xsi3 = xsii*xsi2
						xsi4 = xsii*xsi3
						zetai = (1.+cmplx1*xnuion)*xsii
						zetae = (1.+cmplx1*xnuelc)*xsie
						zetai2 = zetai*zetai
						zetai3 = zetai*zetai2
						zetai4 = zetai*zetai3
						if(abs(zetai) .le. 10.0_IDP) then
							ztr = real(zetai)
							zti = aimag(zetai)
							call zzdisp(ztr,zti,zir,zii)
							zi = cmplx(zir,zii)
							y0i = xsii*zi
							y1i = xsii*(zetai+(0.5+zetai2)*zi)
							y2i = xsii*(1.5*zetai + zetai3 +(0.5 + zetai2 + zetai4)*zi)
						else
							y0i = -(xsii/zetai)*(1. + (0.5/zetai2))
							y1i = -(xsii/zetai)*(1. + (1./zetai2))
							y2i = -1.75*(xsii/zetai)*(1. + (23./(14.*zetai2)))
						endif
						ddii = 1. + cmplx1*xnuion*y0i
						zetae2 = zetae*zetae
						zetae3 = zetae*zetae2
						zetae4 = zetae*zetae3
						if(abs(zetae) .le. 10.0_IDP) then
							ztr = real(zetae)
							zti = aimag(zetae)
							call zzdisp(ztr,zti,zer,zei)
							ze = cmplx(zer,zei)
							y0e = xsie*ze
							y1e = xsie*(zetae+(0.5+zetae2)*ze)
							y2e = xsie*(1.5*zetae + zetae3 +(0.5 + zetae2 + zetae4)*ze)
						else
							y0e = -(xsie/zetae)*(1. + (0.5/zetae2))
							y1e = -(xsie/zetae)*(1. + (1./zetae2))
							y2e = -1.75*(xsie/zetae)*(1. + (23./(14.*zetae2)))
						endif
						ddee = 1. + cmplx1*xnuelc*y0e
					endif
					stfe = y0e/ddee
					stfi = (1. + y0i/ddii)/tauie
					reii = 1. + stfe + stfi
					rei = 1./reii
					sei1 = y2e + y2i*tauie
					stfe = (y1e*y1e)/ddee
					stfi = (y1i*y1i)/ddii
					sei2 = -cmplx1*(xnuion*stfi + xnuelc*stfe)
					sei3 = rei*(((y1e/ddee) - (y1i/ddii))**2)
					sei = sei1 + sei2 + sei3 
					eilnd(j,l) = coef*aimag(sei)	  
				end do
			end do
			call blockjl(eildrr,1,2,2,0,2,0,eilnd)
			call blockjl(eildrt,-1,2,2,1,1,0,eilnd)
			call blockjl(eildrz,-1,2,2,0,1,1,eilnd)
			call blockjl(eildtt,1,2,2,2,0,0,eilnd)
			call blockjl(eildtz,1,2,2,1,0,1,eilnd)
			call blockjl(eildzz,1,2,2,0,0,2,eilnd)
			call blockjl(eildr,1,2,2,0,1,0,eilnd)
			call blockjl(eildt,-1,2,2,1,0,0,eilnd)
			call blockjl(eildz,-1,2,2,0,0,1,eilnd)
		end if

!   Two fluid terms

		if (twofl_on == 1) then
			call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
			sd2=feq-qqinv*cureq
			do l=1,leqmax
				sceq1(:,l)=feq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=feq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=feq*sd1*jbgtt(:,l)/sd2
			end do
			call blockj(sceq1,1,2,2,3,0,0,-betiom)
			call blockj(sceq2,-1,2,2,2,1,0,2*betiom)
			call blockj(sceq3,1,2,2,1,2,0,-betiom)
			call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,3)
			call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,3)
			do l=1,leqmax
				sceq4(:,l)=sceq4(:,l)+rinv*sceq2(:,l)
			end do
			call blockj(sceq4,-1,2,2,2,0,0,-betiom)
			call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,3)
			call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,3)
			call blockj(sceq4,1,2,2,1,1,0,betiom)
			do l=1,leqmax
				sceq1(:,l)=rinv*cureq*sd1*jbgrr(:,l)/sd2
				sceq2(:,l)=rinv*cureq*sd1*jbgrt(:,l)/sd2
				sceq3(:,l)=rinv*cureq*sd1*jbgtt(:,l)/sd2
			end do
			call blockj(sceq1,1,2,2,2,0,1,betiom)
			call blockj(sceq2,-1,2,2,1,1,1,-2*betiom)
			call blockj(sceq3,1,2,2,0,2,1,betiom)
			call dbydtheq(sceq4,sceq1,1,0.0_IDP,1.0_IDP,4)
			call dbydreq(sceq4,sceq2,1.0_IDP,-1.0_IDP,4)
			call blockj(sceq4,-1,2,2,1,0,1,betiom)
			call dbydtheq(sceq4,sceq2,-1,0.0_IDP,1.0_IDP,4)
			call dbydreq(sceq4,sceq3,1.0_IDP,-1.0_IDP,4)
			do l=1,leqmax
				sceq4(:,l)=sceq4(:,l)-rinv*sceq3(:,l)
			end do
			call blockj(sceq4,1,2,2,0,1,1,-betiom)
			do l=1,leqmax
				sceq1(:,l)=rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=rinv*cureq*sd1*dgrtz(:,l)
				sceq3(:,l)=rinv*cureq*sd1*dgttz(:,l)
			end do
			call blockj(sceq1,-1,2,2,2,0,0,0.5*betiom)
			call blockj(sceq2,1,2,2,1,1,0,-betiom)
			call blockj(sceq3,-1,2,2,0,2,0,0.5*betiom)
			sd3=rinv*feq*sd1
			sd3(0)=2*feq(0)*(preq(1)-preq(0))/(r(1)*r(1))
			do l=1,leqmax
				sceq1(:,l)=sd3*dgrrt(:,l)-rinv*cureq*sd1*dgrrz(:,l)
				sceq2(:,l)=feq*sd1*dgttr(:,l)+2*sd3*jbgtt(:,l)/sd2-rinv*cureq*sd1*dgrtz(:,l)
			end do
			call dbydtheq(sceq3,sceq1,-1,0.0_IDP,1.0_IDP,2)
			call dbydreq(sceq3,sceq2,1.0_IDP,-1.0_IDP,2)
			call dbydr0(sd4,qqinv,0.0_IDP,1.0_IDP,0)
			call dbydr0(sd5,cureq,0.0_IDP,1.0_IDP,0)
			call dbydreq(sceq2,jsq,0.0_IDP,0.5_IDP,0)
			do l=1,leqmax
				sceq1(:,l)=cureq*sd1*(qqinv*(rinv*dgrtt(:,l)-dgttr(:,l))-(sd4+2*rinv*qqinv)*jbgtt(:,l)/sd2+ &
						      rinv*(dbsjtbj(:,l)-rinv*(cureq*sceq2(:,l)+sd5*jsq(:,l))/sd2)/(eps*eps))
			end do
			call dbydreq(sceq3,sceq1,1.0_IDP,-1.0_IDP,4)
			call blockj(sceq3,1,2,2,1,0,0,-0.5*betiom)
			call dbydtheq(sceq3,sceq1,1,0.0_IDP,1.0_IDP,4)
			do l=1,leqmax
				sceq1(:,l)=sd1*(feq*dgttt(:,l)-cureq*dgttz(:,l))
				sceq2(:,l)=2*sd3*(dgrtt(:,l)-jbgtt(:,l)/sd2)-sd1*(feq*dgttr(:,l)+rinv*cureq*dgrtz(:,l))
			end do
			call dbydreq(sceq5,sceq1,0.0_IDP,1.0_IDP,3)
			do l=1,leqmax
				sceq3(:,l)=sceq3(:,l)+rinv*sceq5(:,l)
			end do
			call dbydtheq(sceq3,sceq2,1,1.0_IDP,-1.0_IDP,2)
			call blockj(sceq3,-1,2,2,0,1,0,-0.5*betiom)
		end if

! diffusion term added
		call block_dlsq_r(2,4,stdifu)

! p equation

		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)
		sd2=feq/(feq-qqinv*cureq)
		sd3=cureq/(feq-qqinv*cureq)
		sd4=sd1*sd2
		call block0(sd4,3,2,1,0,0,1.0_IDP)
		sd4=-sd1*rinv*sd3
		call block0(sd4,3,2,0,0,1,1.0_IDP)

!		do l=1,leqmax
!			sceq1(:,l)=preq*djroj(:,l)
!		end do
!		call blockj(sceq1,1,3,2,1,0,0,gamma)
!		do l=1,leqmax
!			sceq1(:,l)=-preq*djtoj(:,l)
!		end do
!		call blockj(sceq1,-1,3,2,0,1,0,gamma)

		do l=1,leqmax
			sceq1(:,l)=-preq*(sd2*djtoj(:,l)-rinv*sd3*djzoj(:,l))
		end do
		call blockj(sceq1,-1,3,2,0,1,0,gamma)
		do l=1,leqmax
			sceq1(:,l)=-preq*(r*dbsjzoj(:,l)-sd2*djroj(:,l))
		end do
		call blockj(sceq1,1,3,2,1,0,0,gamma)
		do l=1,leqmax
			sceq1(:,l)=-preq*(rinv*sd3*djroj(:,l)-dbsjtoj(:,l))
		end do
		call blockj(sceq1,1,3,2,0,0,1,gamma)
		call dbydr0(sd4,sd2,0.0_IDP,1.0_IDP,0)
		sd2=preq*sd4
		call block0(sd2,3,2,1,0,0,gamma)
		call dbydr0(sd4,sd3,0.0_IDP,1.0_IDP,0)
		sd2=-rinv*preq*sd4
		call block0(sd2,3,2,0,0,1,gamma)

!  parallel thermal velocity term		

		do l=1,leqmax
			sceq1(:,l)=preq*bmod(:,l)/(feq-qqinv*cureq)
		end do
		call om(sceq1,1,3,7,0,0,0,-gamma)
		call grpareq(sceq2,bmod,1,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=preq*sceq2(:,l)/(feq-qqinv*cureq)
		end do
		call blockj(sceq1,-1,3,7,0,0,0,gamma)

!  Sheared equilibrium toroidal flow velocity for pressure equation

		call block0(vzt_eq,3,3,0,0,1,-1.0_IDP)

!   Two fluid terms

		if (twofl_on == 1) then
			sd2=gamma*betiom*preq/denseq
			sd3=sd2/(feq-qqinv*cureq)
 	 		call dbydr0(sd4,feq,0.0_IDP,1.0_IDP,0)
			sd4=sd3*sd4
			call block0(sd4,3,3,1,0,0,1.0_IDP)
 	 		call dbydr0(sd4,cureq,0.0_IDP,1.0_IDP,0)
			sd4=rinv*sd3*sd4
			call block0(sd4,3,3,0,0,1,-1.0_IDP)
			call dbydzteq(sceq1,bst,-1,0.0_IDP,1.0_IDP)
			call dbydtheq(sceq2,bst,-1,0.0_IDP,1.0_IDP,0)
			do l=1,leqmax
				sceq1(:,l)=r*sd3*sceq1(:,l)
				sceq2(:,l)=r*sd3*sceq2(:,l)
			end do
			call blockj(sceq1,1,3,3,1,0,0,-1.0_IDP)
			call blockj(sceq2,1,3,3,0,0,1,1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd2*omdr(:,l)
				sceq2(:,l)=sd2*omdt(:,l)
				sceq3(:,l)=sd2*omdz(:,l)
			end do
			call blockj(sceq1,-1,3,3,0,1,0,-2.0_IDP)
			call blockj(sceq2,1,3,3,1,0,0,-2.0_IDP)
			call blockj(sceq3,1,3,3,0,0,1,-2.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd1*sd3*grtoj(:,l)
				sceq2(:,l)=sd1*sd3*gttoj(:,l)
			end do
			call om(sceq1,-1,3,1,1,0,0,-epsq)
			call om(sceq2,1,3,1,0,1,0,epsq)
			sd4=rinv*sd1*sd3*cureq
			call block0(sd4,3,1,1,1,0,-1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=r*sd1*sd3*bst(:,l)
			end do
			call blockj(sceq1,-1,3,1,2,0,0,1.0_IDP)
			do l=1,leqmax
				sceq1(:,l)=sd1*sd2*dgrtp(:,l)
				sceq2(:,l)=sd1*sd2*dgttp(:,l)
			end do
			call blockj(sceq1,1,3,1,1,0,0,-epsq)
			call blockj(sceq2,-1,3,1,0,1,0,epsq)
			do l=1,leqmax
				sceq1(:,l)=rinv*sd1*sd3*cureq*djtoj(:,l)
				sceq2(:,l)=sd1*sd2*dbsjtoj(:,l)
			end do
			call blockj(sceq1,-1,3,1,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,3,1,1,0,0,1.0_IDP)
		end if

! diffusion term added
		call block_dlsq_r(3,3,stdifp)

!  u-zeta expression

		do l=1,leqmax
			sceq2(:,l)=denseq*jbgrt(:,l)
		end do
		call blockj(sceq2,-1,4,2,1,1,0,-2.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=denseq*jbgrr(:,l)
		end do
		call blockj(sceq2,1,4,2,2,0,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=denseq*jbgtt(:,l)
		end do
		call blockj(sceq2,1,4,2,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=rinv*sceq2(:,l)
		end do
		call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
		call clgam(sceq2,1,1,2,0.0_IDP,-1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=denseq*sceq2(:,l)-denseqr*jbgrt(:,l)
		end do
		call blockj(sceq1,-1,4,2,1,0,0,1.0_IDP)
		call clgam(sceq2,1,2,2,0.0_IDP,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=denseq*sceq2(:,l)+denseqr*jbgtt(:,l)
		end do
		call blockj(sceq1,1,4,2,0,1,0,1.0_IDP)
		sd1=-1.0_IDP
		call block0(sd1,4,4,0,0,0,1.0_IDP)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	NBI particle effects	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

!  Load 1st Omega-d terms in fast ion density equation:
		write(0,'(/"vfova2 = ",e12.4)') vfova2(1)

		do l=1,leqmax
			sceq1(:,l)=vfova2*omdr(:,l)/(epsq*omcyd)
			sceq2(:,l)=vfova2*omdt(:,l)/(epsq*omcyd)
			sceq3(:,l)=vfova2*omdz(:,l)/(epsq*omcyd)
		end do
		
		call blockj(sceq1,-1,5,5,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,5,5,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,5,5,0,0,1,-1.0_IDP)

!  Sheared equilibrium toroidal flow velocity for fast ion density equation

		call block0(vzt_eq,5,5,0,0,1,-1.0_IDP)

!  Load 1st Omega-d terms in fast ion parallel velocity equation:

		call blockj(sceq1,-1,6,6,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,6,6,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,6,6,0,0,1,-1.0_IDP)

!  Sheared equilibrium toroidal flow velocity for fast ion parallel velocity equation

		call block0(vzt_eq,6,6,0,0,1,-1.0_IDP)

!  Load 2nd Omega-d terms in fast ion density equation:
		write(0,'("nfeq = ",e12.4)') nfeq(1)

		do l=1,leqmax
			sceq1(:,l)=nfeq(:)*omdr(:,l)
			sceq2(:,l)=nfeq(:)*omdt(:,l)
			sceq3(:,l)=nfeq(:)*omdz(:,l)
		end do

		call blockj(sceq1,-1,5,2,0,1,0,-1.0_IDP)
		call blockj(sceq2,1,5,2,1,0,0,-1.0_IDP)
		call blockj(sceq3,1,5,2,0,0,1,-1.0_IDP)

! diffusion term added
		call block_dlsq_r(5,5,stdifnf)

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

		do l=1,leqmax
			sceq1(:,l)=nfeq*bmod(:,l)/(feq-qqinv*cureq)
		end do
		
		call om(sceq1,1,5,6,0,0,0,-1.0_IDP)

!    Omega* term

		sd1=dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
		sd2=dnfeqdr*feq/(feq-qqinv*cureq)
		
		call block0(sd1,5,2,0,0,1,1.0_IDP)
		call block0(sd2,5,2,1,0,0,-1.0_IDP)

!   EP FLR effects

		if (epflr_on == 1 .and. r_epflr > 0.0) then
			sd1=-dnfeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=-dnfeqdr*feq/(feq-qqinv*cureq)
			call block0(sd1,5,iw,0,0,1,1.0_IDP)
			call block0(sd2,5,iw,1,0,0,-1.0_IDP)
			sd1=epsq*omcyd*omegar*nfeq/vfova2
			call block0(sd1,5,iw,0,0,0,1.0_IDP)
		end if

!  Load remaining terms in fast ion parallel velocity equation

!    Landau closure term

		do l=1,leqmax
			sceq1(:,l)=1.414213*LcA1*vfova*bmod(:,l)/(feq-qqinv*cureq)
		end do
		call blockj_landau_grad_parallel(sceq1,1,6,6,0,0,0,1.0_IDP)

!    Parallel gradient term

		do l=1,leqmax
			sceq1(:,l)=2*LcA0*vfova2*bmod(:,l)/(nfeq*(feq-qqinv*cureq))
		end do
		call om(sceq1,1,6,5,0,0,0,-1.0_IDP)

!    Omega* term

		if(epflr_on .eq. 0 .or. r_epflr == 0.0) then 

			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*feq/(nfeq*(feq-qqinv*cureq))
			call block0(sd1,6,1,0,0,1,1.0_IDP)
			call block0(sd2,6,1,1,0,0,-1.0_IDP)
		 
		else 

!   EP FLR effects

			sd1=vfova2*dnfeqdr*rinv*cureq/(nfeq*(feq-qqinv*cureq))
			sd2=vfova2*dnfeqdr*rinv*feq/(nfeq*(feq-qqinv*cureq))
			call block0(sd1,6,ix1,0,0,0,1.0_IDP)
			call block0(sd2,6,ix2,0,0,0,-1.0_IDP)

!   EP FLR effects auxiliary equations

			sd1=1.0_IDP
			call block0(sd1,iw,iw,0,0,0,1.0_IDP)
			coef=-r_epflr*r_epflr
			call block_dlsq(iw,iw,coef)
			call block_dlsq(iw,2,-coef)

			call block0(sd1,ix1,ix1,0,0,0,1.0_IDP)
			call block_dlsq(ix1,ix1,coef)
			call block0(sd1,ix1,1,0,0,1,-1.0_IDP)

			call block0(sd1,ix2,ix2,0,0,0,1.0_IDP)
			call block_dlsq(ix2,ix2,coef)
			call block0(r,ix2,1,1,0,0,-1.0_IDP)
		 
		end if 

! diffusion term added
		call block_dlsq_r(6,6,stdifvf)

!  End of fast ion moment equations

!  Load terms of the thermal moment of energetic particles

!  Pressure gradient term

		do l=1,leqmax
			sceq1(:,l)=bet0*bmod(:,l)/(2.*denseq*(feq-qqinv*cureq))
		end do
		call om(sceq1,1,7,3,0,0,0,-1.0_IDP)
		
!  Magnetic field perturbation term		
		
		call dbydr0(sd1,preq,0.0_IDP,1.0_IDP,0)		
		do l=1,leqmax
			sceq2(:,l)=sceq1(:,l)*sd1
		end do		
		call blockj(sceq2,1,7,1,1,0,0,1.0_IDP)		
		
!  Sheared equilibrium toroidal flow velocity for the thermal moment of energetic particles equation

		call block0(vzt_eq,7,7,0,0,1,-1.0_IDP)

! diffusion term added
		call block_dlsq_r(7,7,stdifv)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Alpha particle effects	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

		if (alpha_on == 1) then

!  Load 1st Omega-d terms in fast ion density equation:

			do l=1,leqmax
				sceq1(:,l)=valphaova2*omdr(:,l)/(epsq*omcydalp)
				sceq2(:,l)=valphaova2*omdt(:,l)/(epsq*omcydalp)
				sceq3(:,l)=valphaova2*omdz(:,l)/(epsq*omcydalp)
			end do
		
			call blockj(sceq1,-1,8,8,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,8,8,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,8,8,0,0,1,-1.0_IDP)

!  Sheared equilibrium toroidal flow velocity for fast ion density equation

			call block0(vzt_eq,8,8,0,0,1,-1.0_IDP)

!  Load 1st Omega-d terms in fast ion parallel velocity equation:

			call blockj(sceq1,-1,9,9,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,9,9,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,9,9,0,0,1,-1.0_IDP)

!  Sheared equilibrium toroidal flow velocity for fast ion parallel velocity equation

			call block0(vzt_eq,9,9,0,0,1,-1.0_IDP)

!  Load 2nd Omega-d terms in fast ion density equation:

			do l=1,leqmax
				sceq1(:,l)=nalpeq(:)*omdr(:,l)
				sceq2(:,l)=nalpeq(:)*omdt(:,l)
				sceq3(:,l)=nalpeq(:)*omdz(:,l)
			end do

			call blockj(sceq1,-1,8,2,0,1,0,-1.0_IDP)
			call blockj(sceq2,1,8,2,1,0,0,-1.0_IDP)
			call blockj(sceq3,1,8,2,0,0,1,-1.0_IDP)

! diffusion term added
			call block_dlsq_r(8,8,stdifnalp)

!  Load remaining terms in fast ion density equation

!    Parallel gradient term

			do l=1,leqmax
				sceq1(:,l)=nalpeq*bmod(:,l)/(feq-qqinv*cureq)
			end do
		
			call om(sceq1,1,8,9,0,0,0,-1.0_IDP)

!    Omega* term

			sd1=dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
			sd2=dnalpeqdr*feq/(feq-qqinv*cureq)
		
			call block0(sd1,8,2,0,0,1,1.0_IDP)
			call block0(sd2,8,2,1,0,0,-1.0_IDP)
		 

!   EP FLR effects

			if (epflr_on == 1 .and. r_epflralp > 0.0) then
				sd1=-dnalpeqdr*rinv*cureq/(feq-qqinv*cureq)
				sd2=-dnalpeqdr*feq/(feq-qqinv*cureq)
				call block0(sd1,8,iwa,0,0,1,1.0_IDP)
				call block0(sd2,8,iwa,1,0,0,-1.0_IDP)
				sd1=epsq*omcydalp*omegar*nalpeq/valphaova2
				call block0(sd1,8,iwa,0,0,0,1.0_IDP)
			end if

!  Load remaining terms in fast ion parallel velocity equation

!    Landau closure term

			do l=1,leqmax
				sceq1(:,l)=1.414213*LcA1alp*valphaova*bmod(:,l)/(feq-qqinv*cureq)
			end do
			call blockj_landau_grad_parallel(sceq1,1,9,9,0,0,0,1.0_IDP)

!    Parallel gradient term

			do l=1,leqmax
				sceq1(:,l)=2*LcA0alp*valphaova2*bmod(:,l)/(nalpeq*(feq-qqinv*cureq))
			end do
			call om(sceq1,1,9,8,0,0,0,-1.0_IDP)

!    Omega* term

			if(epflr_on .eq. 0 .or. r_epflralp == 0.0) then 

				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*feq/(nalpeq*(feq-qqinv*cureq))
				call block0(sd1,9,1,0,0,1,1.0_IDP)
				call block0(sd2,9,1,1,0,0,-1.0_IDP)

			else

!   EP FLR effects

				sd1=valphaova2*dnalpeqdr*rinv*cureq/(nalpeq*(feq-qqinv*cureq))
				sd2=valphaova2*dnalpeqdr*rinv*feq/(nalpeq*(feq-qqinv*cureq))
				call block0(sd1,9,ix1a,0,0,0,1.0_IDP)
				call block0(sd2,9,ix2a,0,0,0,-1.0_IDP)

!   EP FLR effects auxiliary equations

				sd1=1.0_IDP
				call block0(sd1,iwa,iwa,0,0,0,1.0_IDP)
				coef=-r_epflralp*r_epflralp
				call block_dlsq(iwa,iwa,coef)
				call block_dlsq(iwa,2,-coef)

				call block0(sd1,ix1a,ix1a,0,0,0,1.0_IDP)
				call block_dlsq(ix1a,ix1a,coef)
				call block0(sd1,ix1a,1,0,0,1,-1.0_IDP)

				call block0(sd1,ix2a,ix2a,0,0,0,1.0_IDP)
				call block_dlsq(ix2a,ix2a,coef)
				call block0(r,ix2a,1,1,0,0,-1.0_IDP)
		 
			end if 

! diffusion term added
			call block_dlsq_r(9,9,stdifvalp)

!  End of alpha particles moment equations

		end if

		if (matrix_out) then
			call mat_out(1)
		else
			amat=-dtd2*amat
			bmat=-dtd2*bmat
			cmat=-dtd2*cmat
		end if

  		if (matrix_out) then
			amat=0.0
			bmat=0.0
			cmat=0.0
		end if

!  put in l.h.s. of equation

!  psi equation

		sd1=1.0_IDP
		call block0(sd1,1,1,0,0,0,1.0_IDP)

!  u-zeta equation

		call block0(sd1,2,4,0,0,0,1.0_IDP)

! p equation

		call block0(sd1,3,3,0,0,0,1.0_IDP)

! fast ion density moment equation

		call block0(sd1,5,5,0,0,0,1.0_IDP)

! fast ion parallel velocity moment equation

		call block0(sd1,6,6,0,0,0,1.0_IDP)
		
! thermal moment equation

		call block0(sd1,7,7,0,0,0,1.0_IDP)		

		if (alpha_on == 1) then

! fast ion density moment equation

			call block0(sd1,8,8,0,0,0,1.0_IDP)

! fast ion parallel velocity moment equation

			call block0(sd1,9,9,0,0,0,1.0_IDP)
		
		end if

		if (matrix_out) then
			call mat_out(2)
			write(0,'("==================================================")')
			write(0,'("========FAR3d Eigensolver module activated========")')
			write(0,'("==================================================")')
			write(0,'("Execute xEigen located at /FAR3d/Addon/Eigensolver")')
			write(0,'("==================================================")')
			write(0,'("==./../../xEigen frq grwth========================")')
			write(0,'("==frq : reference normalized mode frequency=======")')    
			write(0,'("==grwth : reference normalized growth rate========")')    
			write(0,'("==================================================")')    
			stop
		end if

		yt=0.0_IDP
		do j=1,mjm1
			do l1=1,lmx
				do l2=1,lmx
					if (j > 1) yt(l1,j)=yt(l1,j)+cmat(l1,l2,j)*xt(l2,j-1)
					yt(l1,j)=yt(l1,j)+amat(l1,l2,j)*xt(l2,j)
					yt(l1,j)=yt(l1,j)+bmat(l1,l2,j)*xt(l2,j+1)
				end do
			end do
		end do

!  u-dlpersq(phi)=0

		lskp=3*lmaxn
		do j=1,mjm1
			do l=1,lmaxn
				yt(l+lskp,j)=0.0_IDP
			end do
		end do

!   Ion FLR effects auxiliary equation

		if (iflr_on == 1) then
			lskp=(iq-1)*lmaxn
			do j=1,mjm1
				do l=1,lmaxn
					yt(l+lskp,j)=0.0_IDP
				end do
			end do
		end if

!   EP FLR effects auxiliary equations

		if (epflr_on == 1) then
			if (r_epflr > 0.0) then
				lskp=(iw-1)*lmaxn
				do j=1,mjm1
					do l=1,3*lmaxn
						yt(l+lskp,j)=0.0_IDP
					end do
				end do
			end if
			if (alpha_on == 1 .and. r_epflralp > 0.0) then
				lskp=(iwa-1)*lmaxn
				do j=1,mjm1
					do l=1,3*lmaxn
						yt(l+lskp,j)=0.0_IDP
					end do
				end do
			end if
		end if

		bmat(:,:,mjm1)=0.0_IDP
		cmat(:,:,1)=0.0_IDP

		call decbt(lmx,mjm1,amat,bmat,cmat,ipc,ier)
		if (ier /= 0) stop 19

	end subroutine linstart
