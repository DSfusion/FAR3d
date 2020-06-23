!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program metric_element_create.f, which is currently
!       under development by D. A. Spong of the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to periodic change
!       and improvement without notice.
!
!
      use read_boozer_mod, iota_bw=>iota_b, bmnc_bw=>bmnc_b, &
     	 rmnc_bw=>rmnc_b,zmns_bw=>zmns_b,gmnc_bw=>gmnc_b, &
     	 pmns_bw=>pmns_b,bmns_bw=>bmns_b, &
     	 rmns_bw=>rmns_b,zmnc_bw=>zmnc_b,gmns_bw=>gmns_b, &
     	 pmnc_bw=>pmnc_b,  idx_bw=>idx_b, &
	 pres_bw=>pres_b,phip_bw=>phip_b, &
     	 bvco_bw=>bvco_b,buco_bw=>buco_b,ixm_bw=>ixm_b, &
     	 ixn_bw=>ixn_b
      use stel_kinds
!
!
      implicit none

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(kind=rprec), parameter :: p5=0.5_dp, one=1.0_dp
      REAL(kind=rprec), PARAMETER :: twopi  = 6.28318530717958623_dp
      REAL(kind=rprec), PARAMETER :: mu_0   = 2.0e-7_dp*twopi
      real(kind=rprec), parameter :: zero = 0.0_dp, two = 2.0_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: hiota, hpres, &
        hphip, hjtor, hjpol
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: iotaf, jpolf, &
        jtorf, phipf, presf
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: iotapf, jpolpf, &
        jtorpf, phippf, prespf
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: xmb,xnb,xm,xn
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: radii, radii_flux
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmncbf, zmnsbf, &
        pmnsbf, bmncbf
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmnsbf, zmncbf, &
        pmncbf, bmnsbf
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmncbh,zmnsbh, &
        pmnsbh, bmncbh
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: rmnsbh,zmncbh, &
        pmncbh, bmnsbh
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: jprl_coef0, &
       jprl_coef1, jprl_coef2
       integer :: nsurf
       real(kind=rprec):: phipc, iotac, jtorc, jpolc, &
        presc, prespc, jtorpc, jpolpc, phippc, iotapc 
       real(kind=rprec) :: phis, phize, phith, rboo, rs, &
        rze, rth, zs, zze, zth, arg, ccosi, ssine, num1, &
        num2, ibf2, ibf3, rboo2, rjac2, rjac2i, gtssub, gstsub, &
        gzssub, gszsub, gtzsub, gztsub, gttsub, gsssub, &
        gzzsub, gaasup, gsasup, gtssup, gzssup, gztsup, &
        den1, cka, beta, t1, t2, t3, t4, t5, cks, &
        thetang, zetang, det, error, rjac2_avg, det_avg, &
        zboo, phiboo, xr, yr, zr
      INTEGER :: nfp, nsd, mnboz, i, j, k, mn, kz, kt, ks, mj
      REAL(kind=rprec) :: r0,b0,amin,beta0,r0max,r0min, &
       aspect, betaxis, ohs2, mat_test_diag, mat_test_offdiag, rr_hi, &
       r1,z1,r2,z2,r1f,z1f,r2f,z2f,tht
      REAL(kind=rprec) :: seval
      character arg1*40,warg1*45
      character*1 tb
      integer nargs, numargs, itheta, izeta, nznt
      integer iargc, unit_no, istat, ierr, ig, lf, is, l
      real viz_flux   !for plotting interior flux surfaces in AVS
      real surf_area_element, surf_area_total
      logical lasym, viz, test_jacob, test_upr_lowr, &
       make_stellgap_data,make_full_torus,surf_compute
!-----------------------------------------------
      tb = char(9)
      numargs = iargc()
      call getarg(1,arg1)
!
!
      warg1 = arg1
      call read_boozer_file(warg1,ierr)
       if (istat.ne.0) stop 22
       
       nfp = nfp_b
       nsd = ns_b
       aspect = aspect_b
       r0max = rmax_b 
       r0min = rmin_b
       betaxis = betaxis_b
       mnboz = mnboz_b
       lasym = lasym_b
       write(*,*) mj, lasym
!
!...   IOTA, PRES, PHIP (= -PHIP_VMEC), JTOR (=I = -BUCO_VMEC) and 
!         JPOL (=J = BVCO_VMEC) are ALL on HALF-MESH!!!   
    
       allocate (hiota(nsd), hpres(nsd), hjpol(nsd), hjtor(nsd),  &
        hphip(nsd), jprl_coef0(nsd), jprl_coef1(nsd), &
        jprl_coef2(nsd), stat=istat)
       if (istat .ne. 0) stop 23
       do k=1,nsd
       hiota(k) = iota_bw(k)
       hpres(k) = mu_0*pres_bw(k)  ! in VMEC versions > 6.00 pressure is given in pascals
       hphip(k) = -phip_bw(k)        ! toroidal fluxes have REVERSED sign respect to VMEC!!
       hjpol(k) = bvco_bw(k)
       hjtor(k) = -buco_bw(k)        ! toroidal fluxes have REVERSED sign respect to VMEC!! 
       end do
       r0 = (r0max+r0min)/two
       amin = r0/aspect
       beta0 = betaxis
!       if (beta0 .le. zero) stop 'Beta0 <= 0'
!       b0 = sqrt((two/beta0)*(1.5_dp*hpres(2)-.5_dp*hpres(3)))
       write(*,'(///)')
       write(*,*) nfp,nsd,aspect
       write(*,'(/)')
       write(*,*) r0max,r0min,betaxis,mnboz
       write(*,'(/)')
       write(*,*) r0,amin,beta0

       allocate (xm(mnboz), xn(mnboz), stat=istat)
       if (istat .ne. 0) stop 24
       do mn=1,mnboz
       xm(mn) = ixm_bw(mn)
       xn(mn) = ixn_bw(mn)
!       write(*,*) mn,xm(mn),xn(mn)              
       end do
!...   RMN, ZMN, PMN and BMN are ALL on HALF-MESH

       allocate (rmncbh(mnboz,nsd), zmnsbh(mnboz,nsd),  &
                pmnsbh(mnboz,nsd), bmncbh(mnboz,nsd), stat = istat)
       if (istat .ne. 0) stop 25
       if(lasym) allocate (rmnsbh(mnboz,nsd), zmncbh(mnboz,nsd),  &
                pmncbh(mnboz,nsd), bmnsbh(mnboz,nsd), stat = istat)
       if (istat .ne. 0) stop 25

       zmnsbh=zero; rmncbh=zero; pmnsbh=zero; bmncbh=zero
       do k=1, nsd
        do mn = 1,mnboz
         bmncbh(mn,k) = bmnc_bw(mn,k)
         rmncbh(mn,k) = rmnc_bw(mn,k)
         zmnsbh(mn,k) = zmns_bw(mn,k)
         pmnsbh(mn,k) = pmns_bw(mn,k)
	if(lasym) then                  
         bmnsbh(mn,k) = bmns_bw(mn,k)
         rmnsbh(mn,k) = rmns_bw(mn,k)
         zmncbh(mn,k) = zmnc_bw(mn,k)
         pmncbh(mn,k) = pmnc_bw(mn,k)
	endif                  
        enddo
       enddo
!      stop 27
      allocate (radii(nsd), radii_flux(nsd), stat = istat)
      do i=1,nsd
        radii(i)=sqrt(real(i-1,kind=rprec)/(nsd-1))     !  r/a= (radii_flux)**1/2
        radii_flux(i) = real(i-1,kind=rprec)/(nsd-1)    !  s = radii_flux
      enddo

!      stop 27
       if(rmncbh(1,nsd) .eq. 0.) then
        write(*,'("May need to add last surface to xbooz calculation")')
!	 stop 55
       end if
!=================
!    MIGRATE TO FULL MESH
!=================

!...  Store surface quantites on RADIAL full mesh

      allocate(iotaf(nsd), presf(nsd), jtorf(nsd), jpolf(nsd), &
              phipf(nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error'

      iotaf=zero ; presf =zero; jtorf =zero; jpolf=zero; phipf=zero

      iotaf(2:nsd-1) = p5*(hiota(2:nsd-1) + hiota(3:nsd))
      presf(2:nsd-1) = p5*(hpres(2:nsd-1) + hpres(3:nsd))
      jtorf(2:nsd-1) = p5*(hjtor(2:nsd-1) + hjtor(3:nsd))
      jpolf(2:nsd-1) = p5*(hjpol(2:nsd-1) + hjpol(3:nsd))
      phipf(2:nsd-1) = p5*(hphip(2:nsd-1) + hphip(3:nsd))

!...  Evaluate and store surface quantities derivatives on RADIAL full mesh
       
      allocate(iotapf(nsd), jpolpf(nsd), jtorpf(nsd), &
              phippf(nsd), prespf(nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error in get_ballooning_grate'

      iotapf=zero; prespf =zero; jtorpf =zero; jpolpf=zero; phippf=zero

!      ohs2 = 2.0_dp*dble(nsd-1)                       ! ds to differentiate on RADIAL half mesh
      ohs2 = dble(nsd-1)                                            ! 2 comes from interpolating to radial FULL mesh
      iotapf(2:nsd-1) = ohs2*(hiota(3:nsd) - hiota(2:nsd-1))
      prespf(2:nsd-1) = ohs2*(hpres(3:nsd) - hpres(2:nsd-1))
      jtorpf(2:nsd-1) = ohs2*(hjtor(3:nsd) - hjtor(2:nsd-1))
      jpolpf(2:nsd-1) = ohs2*(hjpol(3:nsd) - hjpol(2:nsd-1))
      phippf(2:nsd-1) = ohs2*(hphip(3:nsd) - hphip(2:nsd-1))
       
      deallocate(hiota,hpres,hjtor,hjpol,hphip)

!...  store (m,n)-descriptors
  
      allocate (xnb(mnboz), xmb(mnboz), stat=k)
      if (k .ne. 0) stop 'Allocation error in get_ballooning_grate'
      xnb=xn
      xmb=xm

!...  EVALUATE AND STORE Boozer Fourier coefficients and
!           their derivatives on RADIAL full mesh

      allocate (rmncbf(mnboz,nsd), zmnsbf(mnboz,nsd), pmnsbf(mnboz,nsd), &
     	       bmncbf(mnboz,nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error'
      	if(lasym) allocate (rmnsbf(mnboz,nsd), zmncbf(mnboz,nsd), pmncbf(mnboz,nsd), &
     	       bmnsbf(mnboz,nsd), stat=k)
      if (k .ne. 0) stop 'Allocation error'
 
      rmncbf=zero; zmnsbf=zero; pmnsbf=zero; bmncbf=zero
      rmnsbf=zero; zmncbf=zero; pmncbf=zero; bmnsbf=zero
      
      do mn = 1,mnboz
      do k = 1,nsd-1

!...   Boozer Fourier coefficients on RADIAL full mesh

         rmncbf(mn,k) = p5*(rmncbh(mn,k+1)+rmncbh(mn,k))
         zmnsbf(mn,k) = p5*(zmnsbh(mn,k+1)+zmnsbh(mn,k))
         pmnsbf(mn,k) = p5*(pmnsbh(mn,k+1)+pmnsbh(mn,k))
         bmncbf(mn,k) = p5*(bmncbh(mn,k+1)+bmncbh(mn,k))
	 
	if(lasym) then                  
         rmnsbf(mn,k) = p5*(rmnsbh(mn,k+1)+rmnsbh(mn,k))
         zmncbf(mn,k) = p5*(zmncbh(mn,k+1)+zmncbh(mn,k))
         pmncbf(mn,k) = p5*(pmncbh(mn,k+1)+pmncbh(mn,k))
         bmnsbf(mn,k) = p5*(bmnsbh(mn,k+1)+bmnsbh(mn,k))
        endif
	
       enddo
       enddo
       iotaf(1) = iotaf(2)
       open(unit=30,file='geom_bzr_to_cyl.dat',status='unknown')
		write(30,*) nsd-1, mnboz, lasym
		do l=1,mnboz
		  write(*,*) ixm_bw(l),ixn_bw(l)
		  write(30,*) ixm_bw(l),ixn_bw(l)
		end do
		do j=1,nsd-1
		 write(30,*) real(j)/real(nsd),1./iotaf(j)
		end do
		do l=1,mnboz
		 do j=1,nsd-1
	       if(lasym) then                  
		  write(30,'(e15.8,5(2x,e15.8))') rmncbf(l,j), zmnsbf(l,j), pmnsbf(l,j),  &
		           rmnsbf(l,j), zmncbf(l,j), pmncbf(l,j)
	       else
		  write(30,'(e15.8,2(2x,e15.8))') rmncbf(l,j), zmnsbf(l,j), pmnsbf(l,j)
	       endif
		 end do		 
		end do
		close(unit=30)

       call read_boozer_deallocate
       deallocate(xn, xm, rmncbh, zmnsbh, pmnsbh, bmncbh, stat=k)
     
       stop
       end

