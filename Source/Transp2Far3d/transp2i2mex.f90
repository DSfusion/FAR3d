subroutine i2mex_fromTranspEq(treq, nt1, ns, it_orient, ier)
 ! Access equilibrium profile data from TRANSP equilibrium data structure
  use transp_eq
  use i2mex_mod
  use ezspline
  use ezspline_obj
  implicit none

  type(transpeq), intent(in) :: treq
  integer, intent(in) :: nt1 ! no of poloidal rays + 1
  integer, intent(in) :: ns ! no of radial nodes
  integer, intent(in) :: it_orient  ! +1=clockwise; -1=counterclockwise
  integer, intent(out) :: ier

  integer iok, i

  real(i2mex_r8), dimension(:,:), allocatable :: x, z
  real(i2mex_r8), dimension(:), allocatable :: tnew, snew, psinew, pnew, gnew, qnew
  real(i2mex_r8), dimension(:,:), allocatable :: xnew, znew
  type(ezspline1) :: f_spl1
  type(ezspline2) :: f_spl2

  real(i2mex_r8) :: x0, x1, x2, z0, z1, z2, det
  integer ns_in, nt_in

  ier = 0

  ! Set up i2mex equilibrium mapping
  CALL eqm_select('mex2eqs (preliminary)',1)
  CALL eqi_jacheck_maxvar_Set(0.0000001d0)
  i2mex_remap = 1
  i2mex_direct = 0
  i2mex_mx = 0;  i2mex_npsi = 0;  i2mex_kb = 2  ! Boozer coordinates

  ! Get the number of input radial surfaces, theta angles
  ns_in = treq%nsm1 + 1;  nt_in = treq%ntm1 + 1
  write(*,*) ' Radial dimension ns_in = ', ns_in
  write(*,*) ' Poloidal dimension nt_in = ', nt_in

  ! Copy data
  ! Reverse theta orientation if necessary
  x0 = SUM(treq%rg(0,0:treq%ntm1))/REAL(treq%ntm1)
  z0 = SUM(treq%zg(0,0:treq%ntm1))/REAL(treq%ntm1)
  x1 = treq%rg(treq%nsm1,0)
  z1 = treq%zg(treq%nsm1,0)
  x2 = treq%rg(treq%nsm1,1)
  z2 = treq%zg(treq%nsm1,1)
  det = (x1-x0)*(z2-z1)-(x2-x1)*(z1-z0)
 
  allocate(x(nt_in, ns_in), z(nt_in, ns_in))
  if( it_orient*det < 0.0_i2mex_r8 ) then
     print *,'keep theta orientation'
     x = TRANSPOSE(treq%rg(0:ns_in-1,0:nt_in-1))
     z = TRANSPOSE(treq%zg(0:ns_in-1,0:nt_in-1))
  else
     print *,'reverse theta orientation'
     do i = 1, nt_in   !nt1
        x(i,:) = treq%rg(0:ns_in-1, nt_in-i)
        z(i,:) = treq%zg(0:ns_in-1, nt_in-i)
     enddo
  endif

  ! interpolate on finer grid and apply precision/units conversions
  PRINT *,'Spline interpolating to ',ns,' radial, ',nt1-1,' poloidal.'
  ALLOCATE(snew(ns), tnew(nt1))
  ALLOCATE(psinew(ns), pnew(ns), qnew(ns), gnew(ns))
  ALLOCATE(xnew(nt1, ns), znew(nt1, ns))
  snew = (/ (real(i-1,i2mex_r8)/real(ns-1,i2mex_r8), i=1, ns) /)
  tnew = i2mex_twopi_r8* (/ (real(i-1)/real(nt1-1), i=1, nt1) /)  !the

  ! psi
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, treq%psis(0:treq%nsm1), iok)
  call EZspline_interp(f_spl1, ns, snew, psinew, iok)
  call EZspline_free(f_spl1, iok)

  ! pressure
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, 4.0e-7_i2mex_r8*i2mex_pi_r8 * treq%ps(0:treq%nsm1), iok)
  call EZspline_interp(f_spl1, ns, snew, pnew, iok)
  call EZspline_free(f_spl1, iok)

  ! q
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, treq%qs(0:treq%nsm1), iok)
  call EZspline_interp(f_spl1, ns, snew, qnew, iok)
  call EZspline_free(f_spl1, iok)

  ! g
  call EZspline_init(f_spl1, ns_in, (/0,0/), iok)
  call EZspline_setup(f_spl1, treq%fs(0:treq%nsm1), iok)
  call EZspline_interp(f_spl1, ns, snew, gnew, iok)
  call EZspline_free(f_spl1, iok)
 
  ! x
  call EZspline_init(f_spl2, nt_in, ns_in, (/-1,-1/), (/0,0/), iok)
  call EZspline_setup(f_spl2, x, iok)
  call EZspline_interp(f_spl2, nt1, ns, tnew, snew, xnew, iok)
  call EZspline_free(f_spl2, iok)
 
  ! z
  call EZspline_init(f_spl2, nt_in, ns_in, (/-1,-1/), (/0,0/), iok)
  call EZspline_setup(f_spl2, z, iok)
  call EZspline_interp(f_spl2, nt1, ns, tnew, snew, znew, iok)
  call EZspline_free(f_spl2, iok)
  PRINT *,'Splines initialized; remap=',i2mex_remap

  ! To remap poloidally, set i2mex_remap=1
  ! Jacobian ~ x^mx /(|grad psi|^npsi B^kb)
  !PRINT *,'remap = ',i2mex_remap
  if(i2mex_remap /= 0) then
     call i2mex_poloidalRemap(i2mex_mx, i2mex_npsi, i2mex_kb, &
          & nt1, ns, tnew, psinew, gnew, xnew, znew, ier)
  endif

  PRINT *,'Initializing i2mex'
  call i2mex_init('transpeq', &
       nt1, ns, tnew, &
       psinew, pnew, gnew, qnew, &
       xnew, znew, iok)
  PRINT *,'i2mex_init returned ',iok
  call i2mex_error(iok)

  !  dmc-- I think with TRANSP runs it is better to update the q profile
  !  from Psi(s) with a numerical integration-- done here:

  !PRINT *,'Calling gen_q...'
  !call eqm_gen_q('q_i2mex',i2mex_o%id_q,iok)
  !if(iok.ne.0) iok=6
  !call i2mex_error(iok)

  DEALLOCATE(snew, tnew)
  DEALLOCATE(psinew, pnew, qnew, gnew)
  DEALLOCATE(xnew, znew)
  DEALLOCATE(x, z)

end subroutine i2mex_fromTranspEq
