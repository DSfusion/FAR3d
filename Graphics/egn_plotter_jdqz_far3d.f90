!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program egn_plotter.f, which is currently
!       under development by D. A. Spong of the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to periodic change
!       and improvement without notice.
!
!      To compile (Intel compiler):
!
!      ifort -c -O egn_plotter.f
!      ifort -O -o xplot_egn_nw egn_plotter.o `PKG_CONFIG_PATH=/usr/local/plplot5.9.7/lib/pkgconfig pkg-config --cflags --libs plplotd-f77`
!
!
!    Command line arguments (these are optional, if no arguments are given,
!     then the following default values are used):
!
!    1 - inum ( = 0 by default) - format i10
!        if inum = 0, then all eigenvectors are plotted, provided
!        their associated frequency (eigenvalue) < freq_max. A single
!        quantity [phi(rho,theta=0,zeta=0)] is plotted along with the
!        frequency and dominant m,n pair
!        if inum > 0, then only the inum-th eigenvector is plotted and the
!        the individual m,n components of phi are plotted. Different colored
!        curves are used to distinguish them, but at this time, these are
!        not labeled with m,n values
!
!
      module egn_dat
       real*8, DIMENSION(:,:,:), ALLOCATABLE :: egn_vectors
       real*8, DIMENSION(:), ALLOCATABLE :: dm, rho, y_phi, dm_rd
      end module egn_dat
      
      program egn_plot
      use egn_dat
      use bspline
      USE PLPLOT
      implicit none
      include "silo.inc"
      integer :: mn_col, ns, i, j, k, istat, nmat, icol, inum, kz, kt, itheta, &
       izeta, nfp, m, ic, nmat_rd, mpol, icycl
      real*8, DIMENSION(:,:,:), ALLOCATABLE :: egn_vec_extnd
      real*8, DIMENSION(:,:), ALLOCATABLE :: y_plot, egn_vectors0, &
       egn_vectors0_psi,egn_vectors_sngl, egn_vec_extnd_sngl, &
       phi_00plane
      real*8, DIMENSION(:), ALLOCATABLE :: rr, yy, egn_max_local, &
       ymn, ymx, ax, gm, wr, rad
      real*8, DIMENSION(:), ALLOCATABLE :: egn_max,es_energy,em_energy
      integer, DIMENSION(:), ALLOCATABLE :: im_col, in_col, &
       mmm, nnn, iacept, i_orig, color_set, lln, signl
      real*8 :: rho_eval, phi, ymin, ymax, thetang, zetang, twopi, &
       potential, fscale, freq_max, dm_test, dum1, egn_max_local0, &
       cutoff, viz_flux, xcol, moment_0, moment_r, moment_r2, zeta
      integer :: lbmax, mj, l, pi, jx
      real*8, dimension(:,:), allocatable :: rmnc, zmns, &
       phimns
      real*8, dimension(:,:), allocatable :: rmns, zmnc, &
       phimnc
      integer, dimension(:), allocatable :: mmb,nnb
      integer :: icc, lp, itht, mj1, ierr, ndims, npts, dbfile, err
      real, dimension(:), allocatable :: x, z, tht_p, phi_plot, &
       flux_surfs,yp,temp,qq,qcontour
      real, dimension(:,:), allocatable :: egn_vectors_sngl_intrp
      real :: sigma,nw_yfit,slp1,slpn, arg, gm_tar
      real :: curv2
      integer :: islpsw, kxord
      real*8 :: xval
!      real*8 :: bspline_mp_dbsval
      real*8, dimension(:), allocatable :: bcoef,xfit,yfit,xknot
      logical, dimension(:), allocatable :: egn_keep
      logical :: kp, b_and_w, lasym
     
      character*5 case_num
      character*39 title
      character*12 eigenvalue_r, eigenvalue_i
      character*4 :: mmm_ch, nnn_ch
      character*18 :: title1
      character*35 :: file1,file2
      character*20 :: arg1
      character*4 :: input_form
      logical :: visit_out, displacement, visit_vs_radius,  &
        i_exist
      displacement = .false.    !choice to plot displacement rather than potential
      cutoff = 1.e+20  !only put eigenmodes in curve file that are within cutoff factor of max.
      pi = 4.*atan(1.d0)
      input_form = "asci"
      b_and_w = .false.
      visit_out = .true.      
      visit_vs_radius = .true.
      write(*,*)input_form
      inum = 0
      allocate(color_set(9), stat=istat)
      color_set(1) = 15;color_set(2) = 1;color_set(3) = 3
      color_set(4) = 4;color_set(5) = 8;color_set(6) = 9
      color_set(7) = 10;color_set(8) = 12;color_set(9) = 13
      call getarg(1, arg1)
      write(*,*) arg1
      read(arg1,'(i10)') inum
      file1 = "egn_mode_asci.dat"
      file2 = "egn_values.dat"
      write(*,*) file1, file2
      twopi = 8.d0*atan(1.d0)
      itheta = 200; izeta = 200
      if(input_form .eq. "binr") then
       open(unit=33,file="egn_mode_bin.dat",form="unformatted",  &
          status="old")
      else if(input_form .eq. "asci") then
       open(unit=33,file="egn_mode_asci.dat",status="old")
      end if
      
      if(input_form .eq. "binr") then
       read(33) nmat
       read(33) mn_col
       read(33) ns
              
      else if(input_form .eq. "asci") then

       read(33,*) nmat
       read(33,*) mn_col
       read(33,*) ns

       write(*,'("# of eigenmodes available = ",i10,/,  &
       "# of m,n pairs included = ",i6,/,  &
       "# of flux surfaces = ",i4)') nmat, mn_col, ns
      endif
      open(unit=25,  &
      file="egn_values.dat",status="old")
      
       allocate(rho(ns), stat=istat)
       allocate(yy(ns), stat=istat)
       allocate(ax(ns), stat=istat)
       allocate(im_col(mn_col), stat=istat)
       allocate(in_col(mn_col), stat=istat)
       allocate(lln(mn_col), stat=istat)
       allocate(signl(mn_col), stat=istat)
       allocate(mmm(nmat), stat=istat)
       allocate(nnn(nmat), stat=istat)
       allocate(gm(nmat), stat=istat)
       allocate(wr(nmat), stat=istat)
       allocate(egn_max(nmat), stat=istat)
       allocate(dm(nmat), stat=istat)
       allocate(dm_rd(nmat), stat=istat)
       allocate(ymn(nmat), stat=istat)
       allocate(ymx(nmat), stat=istat)
       allocate(y_phi(ns), stat=istat)
       allocate(y_plot(nmat,ns), stat=istat)
       allocate(rr(ns), stat=istat)
       allocate(egn_vectors0(mn_col,ns), stat=istat)
       allocate(egn_vectors0_psi(mn_col,ns), stat=istat)
       allocate(egn_vectors_sngl(mn_col,ns), stat=istat)
       allocate(egn_vec_extnd_sngl(mn_col,ns), stat=istat)
       allocate(iacept(nmat), stat=istat)
       allocate(i_orig(nmat), stat=istat)
       allocate(es_energy(nmat), stat=istat)
       allocate(em_energy(nmat), stat=istat)
       allocate(egn_max_local(mn_col), stat=istat)
       allocate(egn_keep(mn_col), stat=istat)
        do i=1,nmat
         read(25,*) gm(i), wr(i)
        end do
      
       do i=1,mn_col
        read(33,'(2x,i4,3(3x,i4))') im_col(i),in_col(i),  &
            lln(i),signl(i)
       end do

       allocate (egn_vectors(nmat,mn_col,ns), stat=istat)

       do i=1,7
        read(33,*) gm_tar
       end do

       do i=1,ns
        read(33,*) rho(i)
	rr(i) = rho(i)
       end do
       
       do j=1,nmat
        do m = 1,ns
         do i=1,mn_col
          read(33,*,END=51) egn_vectors0(i,m), egn_vectors0_psi(i,m)
          end do
         end do
         if (inum .eq. 0) then
         do m = 1,ns
          do i=1,mn_col
           egn_vectors(j,i,m) = egn_vectors0(i,m)
          end do
         end do
	 else if(inum .ne. 0 .and. inum .eq. j) then
	  do m = 1,ns
          do i=1,mn_col
           egn_vectors_sngl(i,m) = egn_vectors0(i,m)
          end do
         end do
	end if   !if (inum .eq. 0 or 1
       end do
   51  continue
       do i=1,nmat
	egn_max(i) = -1.e+30
        do j=1,mn_col
	  do k=1,ns
	   if(abs(egn_vectors(i,j,k)) .gt. egn_max(i)) then
            egn_max(i) = abs(egn_vectors(i,j,k))
	    mmm(i) = im_col(j); nnn(i) = in_col(j)
	   end if
	   end do
          end do
	 end do

      if(inum .eq. 0) then
      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plsfam(1,1,200000)
      call plssub(2,2)
      call plinit
      call plscolor(1)
      call plcol0(15)
      call plfont(2)
       do i=1,nmat    !loop over eigenfunctions
	ymax = -1.e+30; ymin = 1.e+30
        do j=1,mn_col
	 mpol = im_col(j)
	 do k=1,ns
	  y_phi(k) = egn_vectors(i,lln(j),k)
	  phi = y_phi(k)
	  ymin = min(ymin,phi)
	  ymax = max(ymax,phi)
	 end do
	end do
	ymx(i) = ymax
	ymn(i) = ymin
	write(*,*) ymn(i),ymx(i)
       end do

      do i=1,nmat
       call plcol0(15)
       call plschr(0.d0,0.7d0)
       if(input_form.eq."asci")write(case_num,'(i5)') i
       if(input_form .eq. "binr")write(case_num,'(i5)') i
       write(mmm_ch,'(i4)') mmm(i)
       write(nnn_ch,'(i4)') nnn(i)
       write(eigenvalue_i,'(e12.3)') gm(i)
       write(eigenvalue_r,'(e12.3)') wr(i)
        title = "freq=" // eigenvalue_r //  &
        ", growth=" // eigenvalue_i
       title1 = "m =" // mmm_ch // ", n =" // nnn_ch
       ymin = ymn(i); ymax = ymx(i)
       call plenv(0.d0, 1.d0, ymin, ymax, 0, 0)
       call plschr(0.d0,0.7d0)
       call pllab('<r>/<a>','Phi_mn_real',title)
       call plschr(0.d0,0.7d0)
       call plmtex("t",-2.0d0,0.5d0,0.5d0,title1)
       do j = 1,mn_col
        do k=1,ns
         yy(k) = egn_vectors(i,lln(j),k)
        end do
        if(j .ne. 9) icycl = mod(j,9)
        icol = color_set(icycl)
        call plcol0(icol)
	call plwidth(3.d0)
        call plline(rr,yy)
	call plwidth(1.d0)
	call plcol0(15)
       end do
      end do
      call plend
!
!
      else if(inum .ne. 0) then
      if(visit_out) open(unit=21,file="egn_out.curve",status="unknown")

      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plssub(1,1)
      call plinit
      call plcol0(15)
      call plfont(2)
      ymin = 1.d+30; ymax = -1.d+30
        do j=1,mn_col
	 mpol = im_col(j)
	 if(visit_out)  &
         write(21,'("#",i4,",",i4)') im_col(j),in_col(j)
	 do k=1,ns
	  if(displacement) then
	    y_phi(k) = egn_vectors_sngl(lln(j),k)*dble(mpol)
	  else
	    y_phi(k) = egn_vectors_sngl(lln(j),k)
	  endif
	 end do
	 do k=1,ns
	  phi = y_phi(k)
	  if(visit_out .and. signl(j) .gt. 0.)  &
                   write(21,*) sqrt(rr(k)), y_phi(k)
	  if(visit_out .and. signl(j) .lt. 0.)  &
                  write(21,*) sqrt(rr(k)), -y_phi(k)
	  ymax = max(ymax,phi); ymin = min(ymin,phi)
	 end do
	end do
       if(ymin .lt. 0.d0) ymin = 1.1*ymin
       if(ymax .gt. 0.d0) ymax = 1.1*ymax
       call plenv(0.d0, 1.d0, ymin, ymax, 0, 0)
       if(dm(inum) .ge. 0.) write(eigenvalue_r,'(e12.5)') dm(inum)
       if(dm(inum) .lt. 0.) write(eigenvalue_r,'(e12.5)') dm(inum)
       write(case_num,'(i5)') inum
       if(dm(inum) .ge. 0.)   &
       title = case_num // ", omega =" // eigenvalue_r
       if(dm(inum) .lt. 0.)   &
       title = case_num // ", omgsq =" // eigenvalue_r
       call pllab('#gr','#gf#dmn#u(#gr)',title)
       do j = 1,mn_col
       write(*,*) j, mn_col, lln(j)
        do k=1,ns
	 yy(k) = egn_vectors_sngl(lln(j),k)
	end do
        icycl = mod(j,10)
        icol = color_set(icycl)
	icol = 15
        call plcol0(icol)
	call plwidth(3.d0)
        call plline(rr,yy)
	call plwidth(1.d0)
	call plcol0(15)
       end do
      call plend
      
             
        i_exist = .false.
	inquire(file='geom_bzr_to_cyl.dat',exist=i_exist)
	write(*,*) i_exist
	
	if(i_exist) then
	 write(*,'("File was present")')
	endif
	open(unit=30,file='geom_bzr_to_cyl.dat',status='old')
	read(30,*) mj, lbmax, lasym
	write(*,*) mj, lbmax, lasym
	 allocate (rmnc(lbmax,mj))
	 allocate (zmns(lbmax,mj))
	 allocate (phimns(lbmax,mj))
	 allocate (qq(mj))
	if(lasym) then
	 allocate (rmns(lbmax,mj))
	 allocate (zmnc(lbmax,mj))
	 allocate (phimnc(lbmax,mj))
	endif
	allocate (rad(mj))
	allocate (mmb(lbmax))
	allocate (nnb(lbmax))
	do l=1,lbmax
          read(30,*) mmb(l),nnb(l)
	end do
	do j=1,mj
	  read(30,*) rad(j), qq(j)
	end do
!        stop 86
	do l=1,lbmax
	 do j=1,mj
	if(lasym) then
          read(30,'(e15.8,5(2x,e15.8))') rmnc(l,j),zmns(l,j),phimns(l,j),  &
     	 rmns(l,j), zmnc(l,j),phimnc(l,j)
	else
          read(30,'(e15.8,2(2x,e15.8))') rmnc(l,j),  &
                      zmns(l,j),phimns(l,j)
	endif
	 end do		 
	end do
	close(unit=30)
	
	rmnc(1,1) = rmnc(1,2)
	open(unit=99,file="out_rmn.dat",status="unknown")
	do j=1,mj
	 write(99,'(e14.7,12(2x,e14.7))') rad(j),(rmnc(l,j), l=1,12)
	end do
	close(unit=99)
	
	
	open(unit=99,file="out_zmn.dat",status="unknown")
	do j=1,mj
	 write(99,'(e14.7,12(2x,e14.7))') rad(j),(zmns(l,j), l=1,12)
	end do
	close(unit=99)
	
	write(*,*) mj, rad(1), rad(mj)   !Boozer data radial resolution
	write(*,*) ns, rho(1)**2, rho(ns)**2   !FAR3D radial resolution
	allocate(yp(mj),xfit(ns),yfit(ns))
	allocate(egn_vectors_sngl_intrp(mn_col,mj))
!	allocate(temp(10*ns))
!	
!       sigma = 1.
!       islpsw = 3
       kxord = 3
       allocate(xknot(ns+kxord),bcoef(ns))
       do l=1,mn_col
	do j=1,ns
	 yfit(j) = real(egn_vectors_sngl(l,j))
	 xfit(j) = real(rho(j))**2
	end do
	call dbsnak(ns,xfit,kxord,xknot)
	call dbsint(ns,xfit,yfit,kxord,xknot,bcoef)
!	call curv1(ns,xfit,yfit,slp1,slpn,islpsw,yp,temp,sigma,ierr)
	if(ierr .ne. 0) write(*,*) ierr, l
	do j=1,mj
	 xval = rad(j)
	 if(xval .gt. xfit(1) .and. xval .lt. xfit(ns)) then
	  nw_yfit = dbsval(xval,kxord,xknot,ns,bcoef)
!	  nw_yfit = curv2(xval,ns,xfit,yfit,yp,sigma)
	 else
	  nw_yfit = 0.
	 end if
	 egn_vectors_sngl_intrp(l,j) = nw_yfit
	end do
       end do
!       deallocate(yp,xfit,yfit)
!       deallocate(temp)
	       
!       do l=1,mn_col
!	do j=1,mj
!         egn_vectors_sngl_intrp(l,j) = egn_vectors_sngl(l,j)	 
!	end do
!       end do
       
!      Generate data for Visit 2D plot	
		itht = 600; mj1 = mj; nfp = 1
		zeta = 0.
		allocate(x(itht*mj1),z(itht*mj1),phi_plot(itht*mj1),flux_surfs(itht*mj1))
		allocate(qcontour(itht*mj1))
		allocate(tht_p(itht))
		x(:) = 0.
		z(:) = 0.
		phi_plot(:) = 0.
		icc = 0
		do j = 1,mj1
		 do i = 1,itht
		  icc = icc + 1
		  tht_p(i) = twopi*real(i-1)/real(itht-1)
!		  if(j .eq. 2) write(*,*) tht_p(i)
                  x(icc) = 0.; z(icc) = 0.
		  do l = 1,lbmax
		   if(lasym) then
		   
		    arg = tht_p(i)*real(mmb(l)) - zeta*real(nnb(l))
		    x(icc) = x(icc) + rmnc(l,j)*cos(arg) + rmns(l,j)*sin(arg)
		    z(icc) = z(icc) + zmns(l,j)*sin(arg) + zmnc(l,j)*cos(arg)
     		   else
		   
		    arg = tht_p(i)*real(mmb(l)) - zeta*real(nnb(l))
!		    if(j .eq. 12) write(*,*) tht_p(i),mmb(l), arg
		    x(icc) = x(icc) + rmnc(l,j)*cos(arg)
     
		    z(icc) = z(icc) + zmns(l,j)*sin(arg)
     
                   endif
		  end do
		  
		 phi_plot(icc) = 0.
		 do l = 1,mn_col
		  lp = lln(l)
		  arg = real(im_col(lp))*tht_p(i)-zeta*real(in_col(lp))
		   if(signl(l) .eq. 1) then
		     phi_plot(icc) = phi_plot(icc) + egn_vectors_sngl_intrp(lp,j)*sin(arg)
		    else if(signl(l) .eq. -1) then
		     phi_plot(icc) = phi_plot(icc) + egn_vectors_sngl_intrp(lp,j)*cos(arg)
		    endif
		   end do
!		   if(i .eq. 1) write(*,*) x(icc), phi_plot(icc)
		   flux_surfs(icc) = real(rad(j))/real(rad(mj1))
		   qcontour(icc) = qq(j)
		 end do
		end do
!	       write(*,*) icc, itht*mj1
!	       write(*,*) minval(x),maxval(x)
!         stop 99
!	       write(*,*) minval(z),maxval(z)
!	       open(unit=9,file="out_xz.dat",status="unknown")
!	       do i=1,icc
!	         write(9,*) x(i), z(i)
!	       end do
!	       close(unit=9)
!    Write data out to Visit Silo data file
       ierr = dbcreate("2D_jdqz_plot.silo", 17, DB_CLOBBER,  &
        DB_LOCAL,"Comment about the data", 22, DB_PDB, dbfile)
	if(ierr .ne. 0) stop 27
        ndims = 2; npts = itht*mj
       err = dbputpm(dbfile, "pointmesh",9, ndims, x, z, DB_F77NULL,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 28
       err = dbputpv1(dbfile, "phi", 3, "pointmesh", 9, phi_plot,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 29
       err = dbputpv1(dbfile, "flux", 4, "pointmesh", 9, flux_surfs,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 30
       err = dbputpv1(dbfile, "q", 1, "pointmesh", 9, qcontour,  &
           npts, DB_FLOAT, DB_F77NULL, ierr)
	if(ierr .ne. 0) stop 30
      ierr = dbclose(dbfile)
      endif
      end program egn_plot

