! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   VERSION 2.2
!
!   f90 VERSION
!
!   This library contains routines for B-spline interpolation in
!   one, two, and three dimensions. Part of the routines are based
!   on the book by Carl de Boor: A practical guide to Splines (Springer,
!   New-York 1978) and have the same calling sequence and names as
!   the corresponding routines from the IMSL library. For documen-
!   tation see the additional files. NOTE: The results in the demo
!   routines may vary slightly on different architectures.
!
!   by W. Schadow 12/04/99
!   last changed by W. Schadow 07/28/2000
!
!
!   Wolfgang Schadow
!   TRIUMF
!   4004 Wesbrook Mall
!   Vancouver, B.C. V6T 2A3
!   Canada
!
!   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
!
!   www  : http://www.triumf.ca/people/schadow
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   Copyright (C) 2000 Wolfgang Schadow
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the
!   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!   Boston, MA  02111-1307, USA.
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module numeric

  integer, parameter :: sgl = kind(1.0)
  integer, parameter :: dbl = kind(1.0d0)

end module numeric


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module bspline
implicit none

  integer, parameter :: idebug = 0

!
!  ------------------------------------------------------------------
!
!
!   The following routines are included:
!
!            dbsnak
!
!            dbsint
!            dbsval
!            dbsder
!            dbs1gd
!
!            dbs2in
!            dbs2dr
!            dbs2vl
!            dbs2gd
!
!            dbs3in
!            dbs3vl
!            dbs3dr
!            dbs3gd
!
!  ------------------------------------------------------------------
!

  private

  public dbsnak
  public dbsint, dbsval, dbsder, dbs1gd
  public dbs2in, dbs2dr, dbs2vl, dbs2gd
  public dbs3in, dbs3vl, dbs3dr, dbs3gd

  public dbs3vd

contains


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsnak(nx,xvec,kxord,xknot)

!
!  Compute the `not-a-knot' spline knot sequence.
!  (see de Boor p. 167)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length ndata containing the location of the
!            data points.  (input)
!   kxord  - order of the spline.  (input)
!   xknot  - array of length ndata+korder containing the knot
!            sequence.  (output)
!

    use numeric

    implicit none

    integer, intent(in) :: nx, kxord

    real(kind=dbl), dimension(nx), intent(in)        :: xvec
    real(kind=dbl), dimension(nx+kxord), intent(out) :: xknot

    real(kind=dbl) :: eps
    integer        :: ix
    logical        :: first = .true.

    save first,eps


    if (first) then
       first=.false.
       eps = epsilon(1.0_dbl)
       write(6,*) "subroutine dbsnak: "
       write(6,*) "eps = ",eps
    endif

    if((kxord .lt. 0) .or. (kxord .gt. nx)) then
       write(6,*) "subroutine dbsnak: error"
       write(6,*) "0 <= kxord <= nx is required."
       write(6,*) "kxord = ", kxord, " and nx = ", nx,  " is given."
       stop
    endif

    do ix = 1, kxord
       xknot(ix) = xvec(1)
    end do

    if(mod(kxord,2) .eq. 0) then
       do ix = kxord+1, nx
          xknot(ix) = xvec(ix-kxord/2)
       end do
    else
       do ix = kxord+1, nx
          xknot(ix) = 0.5_dbl * (xvec(ix-kxord/2) + xvec(ix-kxord/2-1))
       end do
    endif

    do ix = nx+1, nx+kxord
       xknot(ix) = xvec(nx) * (1.0_dbl + eps)
    end do

  end subroutine dbsnak


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef)

!
!  Computes the spline interpolant, returning the B-spline coefficients.
!  (see de Boor p. 204)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length nx containing the data point
!            abscissas.  (input)
!   xdata  - array of length ndata containing the data point
!            ordinates.  (input)
!   kx     - order of the spline.  (input)
!            korder must be less than or equal to ndata.
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   bscoef - array of length ndata containing the B-spline
!            coefficients.  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: nx, kx
    real(kind=dbl), dimension(nx), intent(in)    :: xdata, xvec
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(out)   :: bcoef

    integer                                :: nxp1, kxm1, kpkm2, leftx, lenq
    integer                                :: ix, ik,ilp1mx, jj, iflag
    real(kind=dbl)                         :: xveci
    real(kind=dbl), dimension((2*kx-1)*nx) :: work


    nxp1  = nx + 1
    kxm1  = kx - 1
    kpkm2 = 2 * kxm1
    leftx = kx
    lenq  = nx * (kx + kxm1)

    do ix = 1, lenq
       work(ix) = 0.0_dbl
    end do

    do  ix = 1, nx
       xveci  = xvec(ix)
       ilp1mx = min0(ix+kx,nxp1)
       leftx   = max0(leftx,ix)
       if (xveci .lt. xknot(leftx)) goto 998
30     if (xveci .lt. xknot(leftx+1)) go to 40
       leftx = leftx + 1
       if (leftx .lt. ilp1mx) go to 30
       leftx = leftx - 1
       if (xveci .gt. xknot(leftx+1)) goto 998
40     call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
       jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
       do ik = 1, kx
          jj       = jj + kpkm2
          work(jj) = bcoef(ik)
       end do
    end do

    call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbsint: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do ix = 1, nx
       bcoef(ix) = xdata(ix)
    end do

    call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)

    return

998 write(6,*) "subroutine dbsint:"
    write(6,*) "xknot(ix) <= xknot(ix+1) required."
    write(6,*) ix,xknot(ix),xknot(ix+1)

    stop

  end subroutine dbsint


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsval(x,kx,xknot,nx,bcoef)

!
!  Evaluates a spline, given its B-spline representation.
!
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsval - value of the spline at x.  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: nx, kx
    real(kind=dbl)                               :: dbsval
    real(kind=dbl)                               :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef

    integer                       :: il, ik, ix, leftx
    real(kind=dbl)                :: save1, save2
    real(kind=dbl), dimension(kx) :: work, dl, dr

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0

    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbsval:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix,xknot(ix),xknot(ix+1)
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if(leftx .eq. 0) then
       write(6,*) "subroutine dbsval:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       stop
    endif

    do ik = 1, kx-1
       work(ik) = bcoef(leftx+ik-kx)
       dl(ik)   = x - xknot(leftx+ik-kx)
       dr(ik)   = xknot(leftx+ik) - x
    end do

    work(kx)  = bcoef(leftx)
    dl(kx)    = x - xknot(leftx)

    do ik = 1, kx-1
       save2 = work(ik)
       do il = ik+1, kx
          save1 = work(il)
          work(il) = (dl(il) * work(il) + dr(il-ik) * save2)                  &
               &           / (dl(il) + dr(il - ik))
          save2 = save1
       end do
    end do

    dbsval = work(kx)

  end function dbsval


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsder(iderx,x,kx,xknot,nx,bcoef)

!
!  Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsder - value of the iderx-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, kx, nx
    real(kind=dbl)                               :: dbsder
    real(kind=dbl), intent(in)                   :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef

    integer                       :: ix, ik, il, leftx
    real(kind=dbl)                :: save, save1, save2, y, sum, dik
    real(kind=dbl), dimension(kx) :: work, dl, dr,bsp

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0
    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbsder:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          stop
       endif
       if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if (leftx .eq. 0) then
       write(6,*) "subroutine dbsder:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "xknot(1)     = ", xknot(1)
       write(6,*) "xknot(nx+kx) = ", xknot(nx+kx)
       write(6,*) "         x   = ", x
       stop
    endif

    if (iderx .eq. 0) then

       do ik = 1,kx-1
          work(ik) = bcoef(leftx+ik-kx)
          dl(ik)   = x - xknot(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
       end do

       work(kx)  = bcoef(leftx)
       dl(kx)    = x - xknot(leftx)

       do ik = 1,kx-1
          save2 = work(ik)
          do il = ik+1,kx
             save1 = work(il)
             work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                  &              / (dl(il) + dr(il - ik))
             save2 = save1
          end do
       end do

       dbsder = work(kx)

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

       bsp(1) = 1.0_dbl
       do ik = 1,kx-iderx-1
          dr(ik) = xknot(leftx+ik) - x
          dl(ik) = x - xknot(leftx+1-ik)
          save   = bsp(1)
          bsp(1) = 0.0_dbl
          do il = 1, ik
             y         = save / (dr(il) + dl(ik+1-il))
             bsp(il)   = bsp(il) + dr(il) * y
             save      = bsp(il+1)
             bsp(il+1) = dl(ik+1-il) * y
          end do
       end do

       do ik = 1, kx
          work(ik) = bcoef(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
          dl(ik)   = x - xknot(leftx+ik-kx)
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          save2 = work(ik)
          do il = ik+1, kx
             save1    = work(il)
             work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
             save2    = save1
          end do
       end do

       sum = 0.0_dbl

       do ix = 1, kx-iderx
          sum = sum + bsp(ix) * work(iderx+ix)
       end do

       dbsder = sum

    else
       dbsder = 0.0_dbl
    endif

  end function dbsder


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef,val)

!
!  Evaluates the derivative of a spline on a grid, given its B-spline
!  representation.
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   nxvec  - length of vector xvec.  (input)
!   xvec   - array of length nxvec containing the points at which the
!            spline is to be evaluated.  (input)
!            xvec should be strictly increasing.
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   val    - array of length nxvec containing the values of the
!            iderx-th derivative of the spline at the points in
!            xvec.  (output)
!

    use numeric

    implicit none

    integer, intent(in)                           :: iderx, nxvec, kx, nx
    real(kind=dbl), dimension(nxvec), intent(in)  :: xvec
    real(kind=dbl), dimension(nx), intent(in)     :: bcoef
    real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dbl), dimension(nxvec), intent(out) :: val

    integer                             :: i, il, ik, ix
    integer, dimension(nxvec)           :: leftx
    real(kind=dbl)                      :: dik
    real(kind=dbl), dimension(nxvec,kx) :: dl, dr, biatx, work
    real(kind=dbl), dimension(nxvec)    :: save1, save2, term

    logical :: same, next


    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next)                                                     &
               &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs1gd:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do ix = 1, nxvec
       if ((xvec(ix).lt.xknot(1)).or.(xvec(ix).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs1gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(ix)
          stop
       endif
    end do

    if (iderx .eq. 0) then

       do ix = 1,nxvec
          biatx(ix,1) = 1._dbl
          val(ix)     = 0._dbl
       end do

       do ik = 1, kx-1
          do ix = 1, nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dbl
          end do

          do il = 1, ik
             do ix = 1,nxvec
                term(ix)     = biatx(ix,il)                                   &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do ik = 1, kx
          do ix = 1, nxvec
             val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik)
          end do
       end do

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

       do ix = 1, nxvec
          biatx(ix,1) = 1._dbl
          val(ix)     = 0._dbl
       end do

       do ik = 1, kx-iderx-1
          do ix = 1, nxvec
             dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix)    = biatx(ix,1)
             biatx(ix,1) = 0.0_dbl
             do il = 1, ik
                term(ix)       = save1(ix)                                    &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix)
                save1(ix)      = biatx(ix,il+1)
                biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix)
             end do
          end do
       end do

       do ik = 1, kx
          do ix = 1, nxvec
             work(ix,ik) = bcoef(leftx(ix)+ik-kx)
             dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx)
          end do
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          do ix = 1, nxvec
             save2(ix) = work(ix,ik)
             do il = ik+1, kx
                save1(ix)   = work(ix,il)
                work(ix,il) = dik * (work(ix,il) - save2(ix))                 &
                     &                 /(dl(ix,il) + dr(ix,il-ik))
                save2(ix)   = save1(ix)
             end do
          end do
       end do

       do i = 1, kx-iderx
          do ix = 1, nxvec
             val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i)
          end do
       end do

    else

       do ix = 1, nxvec
          val(ix) = 0.0_dbl
       end do

    endif

  end subroutine dbs1gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)
!
! This routine is equivalent to the routine dbsder, but it does not
! check the parameters!!!
!
! Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   leftx  - number of the intervall of xknot that includes x
!   dbsdca - value of the ideriv-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, kx, nx
    real(kind=dbl)                               :: dbsdca
    real(kind=dbl), intent(in)                   :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef

    integer                       :: i, ik, il, leftx
    real(kind=dbl)                :: save, save1, save2, y, sum, dik
    real(kind=dbl), dimension(kx) :: work, dl, dr,bsp


    if (iderx .eq. 0) then

       do ik = 1, kx-1
          work(ik) = bcoef(leftx+ik-kx)
          dl(ik)   = x - xknot(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
       end do

       work(kx)  = bcoef(leftx)
       dl(kx)    = x - xknot(leftx)

       do ik = 1, kx-1
          save2 = work(ik)
          do il = ik+1, kx
             save1 = work(il)
             work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                  &              / (dl(il) + dr(il - ik))
             save2 = save1
          end do
       end do

       dbsdca = work(kx)

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
       bsp(1) = 1.0_dbl
       do ik = 1,kx-iderx-1
          dr(ik) = xknot(leftx+ik) - x
          dl(ik) = x - xknot(leftx+1-ik)
          save   = bsp(1)
          bsp(1) = 0.0_dbl
          do il = 1, ik
             y         = save / (dr(il) + dl(ik+1-il))
             bsp(il)   = bsp(il) + dr(il) * y
             save      = bsp(il+1)
             bsp(il+1) = dl(ik+1-il) * y
          end do
       end do

       do ik = 1, kx
          work(ik) = bcoef(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
          dl(ik)   = x - xknot(leftx+ik-kx)
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          save2 = work(ik)
          do il = ik+1, kx
             save1    = work(il)
             work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
             save2    = save1
          end do
       end do

       sum = 0.0_dbl

       do i = 1, kx-iderx
          sum = sum + bsp(i) * work(iderx+i)
       end do

       dbsdca = sum

    else
       dbsdca = 0.0_dbl
    endif

  end function dbsdca


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs2in(nx,xvec,ny,yvec,xydata,ldf,kx,ky,xknot,yknot,bcoef)

!
!  Computes a two-dimensional tensor-product spline interpolant,
!  returning the tensor-product B-spline coefficients.
!
!    nx     - number of data points in the x-direction.  (input)
!    xvec   - array of length nx containing the data points in
!             the x-direction.  (input)
!             xdata must be strictly increasing.
!    ny     - number of data points in the y-direction.  (input)
!    yvec   - array of length ny containing the data points in
!             the y-direction.  (input)
!             ydata must be strictly increasing.
!    xydata - array of size nx by nydata containing the values to
!             be interpolated.  (input)
!             fdata(i,j) is the value at (xdata(i),ydata(j)).
!    ldf    - the leading dimension of fdata exactly as specified in
!             the dimension statement of the calling program.
!             (input)
!    kx     - order of the spline in the x-direction.  (input)
!             kxord must be less than or equal to nxdata.
!    ky     - order of the spline in the y-direction.  (input)
!             kyord must be less than or equal to nydata.
!    xknot  - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!    yknot  - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!    bcoef  - array of length nx*ny containing the
!             tensor-product B-spline coefficients.  (output)
!             bscoef is treated internally as a matrix of size nxdata
!             by nydata.
!

    use numeric

    implicit none

    integer, intent(in)                           :: nx, ny, kx, ky, ldf

    real(kind=dbl), dimension(nx), intent(in)     :: xvec
    real(kind=dbl), dimension(ny), intent(in)     :: yvec
    real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)  :: yknot
    real(kind=dbl), dimension(ldf,*), intent(in)  :: xydata
    real(kind=dbl), dimension(nx,ny), intent(out) :: bcoef

    real(kind=dbl), dimension(max(nx,ny),max(nx,ny))        :: work1
    real(kind=dbl), dimension(max(nx,ny))                   :: work2
    real(kind=dbl), dimension(max((2*kx-1)*nx,(2*ky-1)*ny)) :: work3


    call spli2d(xvec,ldf,xydata,xknot,nx,kx,ny,work2,work3,work1)
    call spli2d(yvec,ny, work1, yknot,ny,ky,nx,work2,work3,bcoef)

  end subroutine dbs2in


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine spli2d(xyvec,ld,xydata,xyknot,n,k,m,work2,work3,bcoef)

    use numeric

    implicit none


    integer, intent(in)                         :: ld, n, k, m
    real(kind=dbl), dimension(n), intent(in)    :: xyvec
    real(kind=dbl), dimension(n+k), intent(in)  :: xyknot
    real(kind=dbl), dimension(ld,m), intent(in) :: xydata
    real(kind=dbl), dimension(m,n), intent(out) :: bcoef

    real(kind=dbl), dimension(n), intent(out)         :: work2
    real(kind=dbl), dimension((2*k-1)*n), intent(out) :: work3


    integer        :: np1, km1, kpkm2, left, lenq, i, iflag, ilp1mx, j, jj
    real(kind=dbl) :: xyveci

    np1   = n + 1
    km1   = k - 1
    kpkm2 = 2 * km1
    left  = k
    lenq  = n * (k + km1)

    do i = 1,lenq
       work3(i) = 0.0_dbl
    end do

    do i = 1, n
       xyveci  = xyvec(i)
       ilp1mx = min0(i+k,np1)
       left   = max0(left,i)
       if (xyveci .lt. xyknot(left)) go to 998
30     if (xyveci .lt. xyknot(left+1)) go to 40
       left = left + 1
       if (left .lt. ilp1mx) go to 30
       left = left - 1
       if (xyveci .gt. xyknot(left+1)) go to 998
40     call bsplvb(xyknot,n+k,k,1,xyveci,left,work2)
       jj = i - left + 1 + (left - k) * (k + km1)
       do j = 1, k
          jj        = jj + kpkm2
          work3(jj) = work2(j)
       end do
    end do

    call banfac(work3,k+km1,n,km1,km1,iflag )

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbs2in: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do j = 1, m
       do i = 1, n
          work2(i) = xydata(i,j)
       end do

       call banslv(work3,k+km1,n,km1,km1,work2)

       do i = 1, n
          bcoef(j,i) = work2(i)
       end do
    end do

    return

998 write(6,*) "subroutine db2in:"
    write(6,*) "i with knot(i) <= x/y < knot(i+1) required."
    write(6,*) "knot(1)   = ", xyknot(1)
    write(6,*) "knot(n+k) = ", xyknot(n+k)
    write(6,*) "      x/y = ", xyveci

    stop

  end subroutine spli2d


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

!
!  evaluates a two-dimensional tensor-product spline, given its
!  tensor-product B-spline representation.    use numeric
!
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   bcoef  - array of length nx*ny containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny.
!   dbs2vl - value of the spline at (x,y).  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: nx, ny, kx, ky
    real(kind=dbl), intent(in)                   :: x, y
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in) :: yknot
    real(kind=dbl), dimension(nx,ny), intent(in) :: bcoef
    real(kind=dbl)                               :: dbs2vl

    integer                       :: ix, iy, iky, leftx, lefty
    real(kind=dbl), dimension(ky) :: work

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    leftx = 0

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs2vl:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if(leftx .eq. 0) then
       write(6,*) "subroutine dbs2vl:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       write(6,*)
       write(6,*) xknot
       stop
    endif

    lefty = 0

    do iy = 1, ny+ky-1
       if (yknot(iy) .gt. yknot(iy+1)) then
          write(6,*) "subroutine dbs2vl:"
          write(6,*) "yknot(iy) <= yknot(iy+1) required."
          write(6,*) iy, yknot(iy), yknot(iy+1)
          stop
       endif
       if((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) lefty = iy
    end do

    if(lefty .eq. 0) then
       write(6,*) "subroutine dbs2vl:"
       write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
       write(6,*) "yknot(iy)   = ", yknot(iy)
       write(6,*) "  y         = ", y
       write(6,*) "yknot(iy+1) = ", yknot(iy+1)
       stop
    endif

    do iky = 1, ky
       work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1,lefty-ky+iky),leftx)
    end do

    dbs2vl = dbsval(y,ky,yknot(lefty-ky+1),ky,work)

  end function dbs2vl


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

!
!  Evaluates the derivative of a two-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the derivative in the x-direction.  (input)
!   idery  - order of the derivative in the y-direction.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   bcoef  - array of length nx*ny containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny.
!   dbs2dr  - value of the (iderx,idery) derivative of the spline at
!            (x,y).  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, idery
    integer, intent(in)                          :: kx, nx, ky, ny
    real(kind=dbl)                               :: dbs2dr
    real(kind=dbl), intent(in)                   :: x, y
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in) :: yknot
    real(kind=dbl), dimension(nx,ny), intent(in) :: bcoef

    integer                       :: ix, iy, iky, nintx, ninty
    real(kind=dbl), dimension(ky) :: work

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintx = 0

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs2dr:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) nintx = ix
    end do

    if(nintx .eq. 0) then
       write(6,*) "subroutine dbs2dr:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       stop
    endif

    ninty = 0

    do iy = 1, ny+ky-1
       if (yknot(iy) .gt. yknot(iy+1)) then
          write(6,*) "subroutine dbs2dr:"
          write(6,*) "yknot(iy) <= yknot(iy+1) required."
          write(6,*) iy, yknot(iy), yknot(iy+1)
          stop
       endif
       if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) ninty = iy
    end do

    if(ninty .eq. 0) then
       write(6,*) "subroutine dbs2dr:"
       write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
       write(6,*) "y = ", y
       stop
    endif

    do iky = 1, ky
       work(iky) =  dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iky),nintx)
    end do

    dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1),ky,work)

  end function dbs2dr


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,ky,xknot,yknot,      &
       & nx,ny,bcoef,val,ldf)

!
!  Evaluates the derivative of a two-dimensional tensor-product spline,
!  given its tensor-product B-spline representation on a grid.
!
!   iderx   - order of the derivative in the x-direction.  (input)
!   idery   - order of the derivative in the y-direction.  (input)
!   nxvec   - number of grid points in the x-direction.  (input)
!   xvec    - array of length nx containing the x-coordinates at
!             which the spline is to be evaluated.  (input)
!             the points in xvec should be strictly increasing.
!   nyvec   - number of grid points in the y-direction.  (input)
!   yvec    - array of length ny containing the y-coordinates at
!             which the spline is to be evaluated.  (input)
!             the points in yvec should be strictly increasing.
!   kx      - order of the spline in the x-direction.  (input)
!   ky      - order of the spline in the y-direction.  (input)
!   xknot   - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!   yknot   - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!   nx      - number of B-spline coefficients in the x-direction.
!             (input)
!   ny      - number of B-spline coefficients in the y-direction.
!             (input)
!   bcoef   - array of length nx*ny containing the
!             tensor-product B-spline coefficients.  (input)
!             bscoef is treated internally as a matrix of size nx
!             by ny.
!   val     - value of the (iderx,idery) derivative of the spline on
!             the nx by ny grid.  (output)
!             value(i,j) contains the derivative of the spline at the
!             point (xvec(i),yvec(j)).
!   ldf     - leading dimension of value exactly as specified in the
!             dimension statement of the calling program.  (input)
!

    use numeric

    implicit none

    integer, intent(in)                           :: iderx, idery
    integer, intent(in)                           :: nxvec, nyvec
    integer, intent(in)                           :: kx, nx, ky, ny
    integer, intent(in)                           :: ldf

    real(kind=dbl), dimension(nxvec), intent(in)  :: xvec
    real(kind=dbl), dimension(nyvec), intent(in)  :: yvec
    real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)  :: yknot
    real(kind=dbl), dimension(nx,ny), intent(in)  :: bcoef
    real(kind=dbl), dimension(ldf,*), intent(out) :: val

    integer                                     :: i, ik, il, ix, iy, ikx, iky
    integer, dimension(nxvec)                   :: leftx
    integer, dimension(nyvec)                   :: lefty
    real(kind=dbl), dimension(nxvec,kx)         :: dl, dr
    real(kind=dbl), dimension(max(nxvec,nyvec)) :: save1
    real(kind=dbl), dimension(nxvec,kx)         :: biatx
    real(kind=dbl), dimension(nyvec,ky)         :: biaty
    real(kind=dbl), dimension(max(nxvec,nyvec)) :: term
    real(kind=dbl), dimension(ky)               :: work

    logical :: same,next


    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next)                                                     &
               &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do i = 1, nx+kx-1
       if (xknot(i) .gt. xknot(i+1)) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "xknot(i) <= xknot(i+1) required."
          write(6,*) i, xknot(i), xknot(i+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do i = 1, nxvec
       if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(i)
          stop
       endif
    end do

    lefty(1) = 0

    call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

    do iy = 2, nyvec
       lefty(iy) = lefty(iy-1)
       same = (yknot(lefty(iy)) .le. yvec(iy))                                &
            &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
       if(.not. same ) then
          lefty(iy) = lefty(iy) + 1
          next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
               &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
          if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
       endif
    end do

    do i = 1, ny+ky-1
       if (yknot(i) .gt. yknot(i+1)) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "yknot(i) <= yknot(i+1) required."
          write(6,*) i, yknot(i), yknot(i+1)
          write(6,*)
          write(6,*) yknot
          stop
       endif
    end do

    do i = 1, nyvec
       if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
          write(6,*) "y = ", yvec(i)
          stop
       endif
    end do

    if ((iderx .eq. 0) .and. (idery .eq. 0)) then

       do ix = 1,nxvec
          biatx(ix,1) = 1._dbl
       end do

       do ik = 1, kx-1
          do ix = 1,nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dbl
          end do

          do il = 1,ik
             do ix = 1,nxvec
                term(ix)     = biatx(ix,il)                                   &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do iy = 1, nyvec
          biaty(iy,1) = 1._dbl
       end do

       do ik = 1, ky-1
          do iy = 1, nyvec
             dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
             dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
             save1(iy) = 0._dbl
          end do

          do il = 1, ik
             do iy = 1,nyvec
                term(iy)     = biaty(iy,il)                                   &
                     &                 / (dr(iy,il) + dl(iy,ik+1-il))
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                save1(iy)    = dl(iy,ik+1-il) * term(iy)
             end do
          end do

          do iy = 1, nyvec
             biaty(iy,ik+1) = save1(iy)
          end do
       end do

       do iy = 1, nyvec
          do ix = 1, nxvec
             val(ix,iy) = 0.0_dbl
          end do
       end do

       do iky = 1, ky
          do ikx = 1, kx
             do iy = 1, nyvec
                do ix = 1, nxvec
                   val(ix,iy) = val(ix,iy)                                    &
                        & + biatx(ix,ikx) * biaty(iy,iky)                     &
                        & * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky)
                end do
             end do
          end do
       end do

    elseif (((iderx .ge. 1) .or. (idery .ge. 1))                              &
         &  .and. ( (iderx .lt. kx) .and. (idery .lt. ky))) then

       do iy = 1, nyvec
          do ix = 1, nxvec
             do iky = 1, ky
                work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,                &
                     &             bcoef(1,lefty(iy)-ky+iky),leftx(ix))
             end do
             val(ix,iy) = dbsder(idery,yvec(iy),ky,                           &
                  &              yknot(lefty(iy)-ky+1),ky,work)
          end do
       end do

    else

       do iy = 1, nyvec
          do ix = 1, nxvec
             val(ix,iy) = 0.0_dbl
          end do
       end do

    endif

  end subroutine dbs2gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kx,ky,kz,         &
       & xknot,yknot,zknot,bcoef)

!
!  Computes a three-dimensional tensor-product spline interpolant,
!  returning the tensor-product B-spline coefficients.
!
!   nx      - number of data points in the x-direction.  (input)
!   xvec    - array of length nxdata containing the data points in
!             the x-direction.  (input)
!             xdata must be increasing.
!   ny      - number of data points in the y-direction.  (input)
!   yvec    - array of length nydata containing the data points in
!             the y-direction.  (input)
!             ydata must be increasing.
!   nz      - number of data points in the z-direction.  (input)
!   zvec    - array of length nzdata containing the data points in
!             the z-direction.  (input)
!             zdata must be increasing.
!   xyzdata - array of size nx by ny by nz containing the
!             values to be interpolated.  (input)
!             xyzdata(i,j,k) contains the value at
!             (xvec(i),yvec(j),zvec(k)).
!   ldf     - leading dimension of fdata exactly as specified in the
!             dimension statement of the calling program.  (input)
!   mdf     - middle dimension of fdata exactly as specified in the
!             dimension statement of the calling program.  (input)
!   kx      - order of the spline in the x-direction.  (input)
!             kxord must be less than or equal to nxdata.
!   ky      - order of the spline in the y-direction.  (input)
!             kyord must be less than or equal to nydata.
!   kz      - order of the spline in the z-direction.  (input)
!             kzord must be less than or equal to nzdata.
!   xknot   - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!   yknot   - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!   zknot   - array of length nz+kz containing the knot
!             sequence in the z-direction.  (input)
!             zknot must be nondecreasing.
!   bcoef   - array of length nx*ny*nz containing the
!             tensor-product B-spline coefficients.  (output)
!             bscoef is treated internally as a matrix of size nx
!             by ny by nz.
!

    use numeric

    implicit none

    integer, intent(in) :: nx, ny, nz, kx, ky, kz
    integer, intent(in) :: ldf, mdf

    real(kind=dbl), dimension(nx), intent(in)         :: xvec
    real(kind=dbl), dimension(ny), intent(in)         :: yvec
    real(kind=dbl), dimension(nz), intent(in)         :: zvec
    real(kind=dbl), dimension(nx+kx), intent(in)      :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)      :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)      :: zknot
    real(kind=dbl), dimension(ldf,mdf,nz), intent(in) :: xyzdata
    real(kind=dbl), dimension(nx,ny,nz), intent(out)  :: bcoef

    integer                                :: iz
    real(kind=dbl), dimension(nx,ny,nz)    :: work1
    real(kind=dbl), dimension(nz)          :: work2
    real(kind=dbl), dimension((2*kz-1)*nz) :: work3


    call spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,work2,work3,work1,     &
         &     nx,ny,nz)

    do iz = 1, nz
       call dbs2in(nx,xvec,ny,yvec,work1(1,1,iz),nx,kx,ky,xknot,yknot,        &
            &        bcoef(1,1,iz))
    end do

  end subroutine dbs3in


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,work2,work3,       &
       & bcoef,nx,ny,nz)

    use numeric

    implicit none

    integer, intent(in)                               :: ldf, mdf, n, k, m, l
    integer, intent(in)                               :: nx, ny, nz
    real(kind=dbl), dimension(n), intent(in)          :: xyzvec
    real(kind=dbl), dimension(n+k), intent(in)        :: xyzknot
    real(kind=dbl), dimension(ldf,mdf,*), intent(in)  :: xyzdata
    real(kind=dbl), dimension(nx,ny,nz), intent(out)  :: bcoef
    real(kind=dbl), dimension(n), intent(out)         :: work2
    real(kind=dbl), dimension((2*k-1)*n), intent(out) :: work3

    integer        :: np1, km1, kpkm2, left, lenq, i, ilp1mx, j, jj, iflag, in
    real(kind=dbl) :: xyzveci


    np1   = n + 1
    km1   = k - 1
    kpkm2 = 2 * km1
    left  = k
    lenq  = n * (k + km1)

    do i = 1, lenq
       work3(i) = 0._dbl
    end do

    do i = 1, n
       xyzveci = xyzvec(i)
       ilp1mx  = min0(i+k,np1)
       left    = max0(left,i)
       if (xyzveci .lt. xyzknot(left)) go to 998
30     if (xyzveci .lt. xyzknot(left+1)) go to 40
       left = left + 1
       if (left .lt. ilp1mx) go to 30
       left = left - 1
       if (xyzveci .gt. xyzknot(left+1)) go to 998
40     call bsplvb(xyzknot,n+k,k,1,xyzveci,left,work2)
       jj = i - left + 1 + (left - k) * (k + km1)
       do j = 1, k
          jj    = jj + kpkm2
          work3(jj) = work2(j)
       end do
    end do

    call banfac(work3,k+km1,n,km1,km1,iflag)

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbs3in: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do j = 1, l
       do i = 1, m
          do in = 1, n
             work2(in) = xyzdata(i,j,in)
          end do

          call banslv(work3,k+km1,n,km1,km1,work2)

          do in = 1, n
             bcoef(i,j,in) = work2(in)
          end do

       end do
    end do

    return

998 write(6,*) "subroutine db3in:"
    write(6,*) "i with knot(i) <= x/y/z < knot(i+1) required."
    write(6,*) "knot(1)   = ", xyzknot(1)
    write(6,*) "knot(n+k) = ", xyzknot(n+k)
    write(6,*) "    x/y/z = ", xyzveci

    stop

  end subroutine spli3d


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)

!
!  Evaluates a three-dimensional tensor-product spline, given its
!  tensor-product B-spline representation.
!
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   z      - z-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   dbs3vl - value of the spline at (x,y,z).  (output)
!

    use numeric

    implicit none

    integer, intent(in)                             :: nx, ny, nz, kx, ky, kz
    real(kind=dbl), intent(in)                      :: x, y, z
    real(kind=dbl), dimension(nx+kx), intent(in)    :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)    :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)    :: zknot
    real(kind=dbl), dimension(nx,ny,nz), intent(in) :: bcoef
    real(kind=dbl)                                  :: dbs3vl

    integer                       :: iz, nintz
    real(kind=dbl), dimension(kz) :: work

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintz = 0

    do iz = 1, nz+kz-1
       if (zknot(iz) .gt. zknot(iz + 1)) then
          write(6,*) "subroutine dbs3vl:"
          write(6,*) "zknot(iz) <= zknot(iz+1) required."
          write(6,*) iz, zknot(iz), zknot(iz+1)
          stop
       endif
       if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
    end do

    if(nintz .eq. 0) then
       write(6,*) "subroutine dbs3vl:"
       write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
       write(6,*) "zknot(iz)   = ", zknot(iz)
       write(6,*) "  z         = ", z
       write(6,*) "zknot(iz+1) = ", zknot(iz+1)
       stop
    endif

    do iz = 1, kz
       work(iz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef(1,1,nintz-kz+iz))
    end do

    dbs3vl = dbsval(z,kz,zknot(nintz-kz+1),kz,work)

  end function dbs3vl


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,xknot,yknot,zknot,         &
       & nx,ny,nz,bcoef)

!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   z      - z-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   dbs3dr - value of the (iderx,idery,iderz) derivative of the
!            spline at (x,y,z).  (output)
!

    use numeric

    implicit none

    integer, intent(in)                              :: iderx, idery, iderz
    integer, intent(in)                              :: nx, ny, nz, kx, ky, kz
    real(kind=dbl), intent(in)                       :: x, y, z
    real(kind=dbl), dimension(nx+kx), intent(in)     :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)     :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)     :: zknot
    real(kind=dbl), dimension(nx,ny,nz), intent(in)  :: bcoef
    real(kind=dbl)                                   :: dbs3dr

    integer                       :: iz, nintz
    real(kind=dbl), dimension(kz) :: work

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintz = 0

    do iz = 1, nz+kz-1
       if (zknot(iz) .gt. zknot(iz + 1)) then
          write(6,*) "subroutine dbs3vl:"
          write(6,*) "zknot(iz) <= zknot(iz+1) required."
          write(6,*) iz, zknot(iz), zknot(iz+1)
          stop
       endif
       if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
    end do

    if(nintz .eq. 0) then
       write(6,*) "subroutine dbs3dr:"
       write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
       write(6,*) "zknot(iz)   = ", zknot(iz)
       write(6,*) "  z         = ", z
       write(6,*) "zknot(iz+1) = ", zknot(iz+1)
       stop
    endif

    do iz = 1, kz
       work(iz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,             &
            &        bcoef(1,1,nintz-kz+iz))
    end do

    dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1),kz,work)

  end function dbs3dr


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,zvec,       &
       & kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,val,ldf,mdf)

!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation on a grid.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   nx     - number of grid points in the x-direction.  (input)
!   xvec   - array of length nx containing the x-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in xvec should be strictly increasing.
!   ny     - number of grid points in the y-direction.  (input)
!   yvec   - array of length ny containing the y-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   nz     - number of grid points in the z-direction.  (input)
!   zvec   - array of length nz containing the z-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   val    - array of size nx by ny by nz containing the values of
!            the (iderx,idery,iderz) derivative of the spline on the
!            nx by ny by nz grid.  (output)
!            value(i,j,k) contains the derivative of the spline at
!            the point (xvec(i), yvec(j), zvec(k)).
!   ldf    - leading dimension of value exactly as specified in the
!            dimension statement of the calling program.  (input)
!   mdf    - middle dimension of value exactly as specified in the
!            dimension statement of the calling program.  (input)
!

    use numeric

    implicit none

    integer, intent(in)                               :: iderx, idery, iderz
    integer, intent(in)                               :: nxvec, nyvec, nzvec
    integer, intent(in)                               :: kx, nx, ky, ny, kz, nz
    integer, intent(in)                               :: ldf,mdf

    real(kind=dbl), dimension(nxvec), intent(in)      :: xvec
    real(kind=dbl), dimension(nyvec), intent(in)      :: yvec
    real(kind=dbl), dimension(nzvec), intent(in)      :: zvec
    real(kind=dbl), dimension(nx+kx), intent(in)      :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)      :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)      :: zknot
    real(kind=dbl), dimension(nx,ny,nz), intent(in)   :: bcoef
    real(kind=dbl), dimension(ldf,mdf,*), intent(out) :: val

    integer                                           :: i, ik, il, ix, iy, iz
    integer                                           :: ikx, iky, ikz
    integer, dimension(nxvec)                         :: leftx
    integer, dimension(nyvec)                         :: lefty
    integer, dimension(nzvec)                         :: leftz
    real(kind=dbl), dimension(nxvec,kx)               :: biatx
    real(kind=dbl), dimension(nyvec,ky)               :: biaty
    real(kind=dbl), dimension(nzvec,kz)               :: biatz
    real(kind=dbl), dimension(max(nxvec,nyvec,nzvec)) :: term, save1

    real(kind=dbl), dimension(max(nxvec,nyvec,nzvec), max(kx,ky,kz)) :: dl, dr

    logical :: same,next


    do i = 1, nx+kx-1
       if (xknot(i) .gt. xknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "xknot(i) <= xknot(i+1) required."
          write(6,*) i, xknot(i), xknot(i+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do i = 1, nxvec
       if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(i)
!debug-begin
          write(6,*) 'i,xknot(1),xknot(nx+kx) ',i,xknot(1),xknot(nx+kx)
!debug-end
     
          stop
       endif
    end do

    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next) call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do i = 1, ny+ky-1
       if (yknot(i) .gt. yknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "yknot(i) <= yknot(i+1) required."
          write(6,*) i, yknot(i), yknot(i+1)
          write(6,*)
          write(6,*) yknot
          stop
       endif
    end do

    do i = 1, nyvec
       if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
          write(6,*) "y = ", yvec(i)
          stop
       endif
    end do

    lefty(1) = 0

    call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

    do iy = 2, nyvec
       lefty(iy) = lefty(iy-1)
       same = (yknot(lefty(iy)) .le. yvec(iy))                                &
            &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
       if(.not. same ) then
          lefty(iy) = lefty(iy) + 1
          next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
               &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
          if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
       endif
    end do

    do i = 1,nz+kz-1
       if (zknot(i) .gt. zknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "zknot(i) <= zknot(i+1) required."
          write(6,*) i, zknot(i), zknot(i+1)
          write(6,*)
          write(6,*) zknot
          stop
       endif
    end do

    do i = 1, nzvec
       if ((zvec(i).lt.zknot(1)).or.(zvec(i).gt.zknot(nz+kz))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
          write(6,*) "z = ", zvec(i)
          stop
       endif
    end do

    leftz(1) = 0

    call huntn(zknot,nz+kz,kz,zvec(1),leftz(1))

    do iz = 2, nzvec
       leftz(iz) = leftz(iz-1)
       same = (zknot(leftz(iz)) .le. zvec(iz))                                &
            &        .and. (zvec(iz) .le. zknot(leftz(iz)+1))
       if(.not. same ) then
          leftz(iz) = leftz(iz) + 1
          next      = (zknot(leftz(iz)) .le. zvec(iz))                        &
               &           .and. (zvec(iz) .le. zknot(leftz(iz)+1))
          if (.not. next) call huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz))
       endif
    end do

    if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

       do ix = 1, nxvec
          biatx(ix,1) = 1.0_dbl
       end do

       do ik = 1, kx-1
          do ix = 1, nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dbl
          end do

          do il = 1, ik
             do ix = 1, nxvec
                term(ix)     = biatx(ix,il) / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do iy = 1, nyvec
          biaty(iy,1) = 1.0_dbl
       end do

       do ik = 1, ky-1
          do iy = 1, nyvec
             dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
             dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
             save1(iy) = 0._dbl
          end do

          do il = 1,ik
             do iy = 1,nyvec
                term(iy)     = biaty(iy,il) / (dr(iy,il) + dl(iy,ik+1-il))
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                save1(iy)    = dl(iy,ik+1-il) * term(iy)
             end do
          end do

          do iy = 1,nyvec
             biaty(iy,ik+1) = save1(iy)
          end do
       end do

       do iz = 1,nzvec
          biatz(iz,1) = 1.0_dbl
       end do

       do ik = 1, kz-1
          do iz = 1, nzvec
             dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz)
             dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik)
             save1(iz) = 0._dbl
          end do

          do il = 1, ik
             do iz = 1, nzvec
                term(iz)     = biatz(iz,il) / (dr(iz,il) + dl(iz,ik+1-il))
                biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
                save1(iz)    = dl(iz,ik+1-il) * term(iz)
             end do
          end do

          do iz = 1, nzvec
             biatz(iz,ik+1) = save1(iz)
          end do
       end do

       do iz = 1,nzvec
          do iy = 1,nyvec
             do ix = 1,nxvec
                val(ix,iy,iz) = 0.0_dbl
             end do
          end do
       end do

       do ikz = 1, kz
          do iky = 1, ky
             do ikx = 1, kx
                do iz = 1, nzvec
                   do iy = 1, nyvec
                      do ix = 1, nxvec
                         val(ix,iy,iz) = val(ix,iy,iz)                        &
                              &  + biatx(ix,ikx) * biaty(iy,iky)              &
                              &  * biatz(iz,ikz)                              &
                              &  * bcoef(leftx(ix)-kx+ikx,                    &
                              &          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
                      end do
                   end do
                end do
             end do
          end do
       end do

    else

       do iz = 1, nzvec
          do iy = 1, nyvec
             do ix = 1, nxvec
                val(ix,iy,iz) = dbs3dr0(iderx,idery,iderz,xvec(ix),           &
                     &  yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,               &
                     &  zknot,nx,ny,nz,bcoef,                                 &
                     &  leftx(ix),lefty(iy),leftz(iz))
             end do
          end do
       end do

    endif

  end subroutine dbs3gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

    use numeric

    implicit none

    integer, intent(in) :: n, jhigh, index, left

    real(kind=dbl), intent(in)                    :: x
    real(kind=dbl), dimension(n), intent(in)      :: t
    real(kind=dbl), dimension(jhigh), intent(out) :: biatx

    integer                          :: j = 1
    integer                          :: i, jp1
    real(kind=dbl)                   :: saved, term
    real(kind=dbl), dimension(jhigh) :: dl, dr


    if (index .eq. 1) then
       j = 1
       biatx(1) = 1.0_dbl
       if (j .ge. jhigh) return
    end if

20  jp1 = j + 1

    dr(j) = t(left+j) - x
    dl(j) = x - t(left+1-j)
    saved = 0._dbl

    do i = 1, j
       term     = biatx(i) / (dr(i) + dl(jp1-i))
       biatx(i) = saved + dr(i) * term
       saved    = dl(jp1-i) * term
    end do

    biatx(jp1) = saved
    j          = jp1

    if (j .lt. jhigh) go to 20

  end subroutine bsplvb


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

    use numeric

    implicit none

    integer, intent(in)                                  :: nroww,nrow
    integer, intent(in)                                  :: nbandl,nbandu
    integer, intent(out)                                 :: iflag
    real(kind=dbl), dimension(nroww,nrow), intent(inout) :: w

    real(kind=dbl) :: pivot, factor
    integer        :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k


    iflag  = 1
    middle = nbandu + 1
    nrowm1 = nrow - 1

    if (nrowm1 .lt. 0) goto 999
    if (nrowm1 .eq. 0) goto 900
    if (nrowm1 .gt. 0) goto 10

10  if (nbandl .gt. 0) go to 30

    do i = 1, nrowm1
       if (w(middle,i) .eq. 0._dbl) go to 999
    end do

    go to 900

30  if (nbandu .gt. 0) go to 60

    do i = 1, nrowm1
       pivot = w(middle,i)
       if(pivot .eq. 0._dbl) go to 999
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do
    end do

    return

60  do i = 1, nrowm1
       pivot = w(middle,i)
       if (pivot .eq. 0._dbl) go to 999
       jmax = min0(nbandl,nrow - i)
       do j = 1,jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do

       kmax = min0(nbandu,nrow - i)

       do k = 1, kmax
          ipk    = i + k
          midmk  = middle - k
          factor = w(midmk,ipk)
          do j = 1, jmax
             w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)                  &
                  &              * factor
          end do
       end do
    end do

900 if (w(middle,nrow) .ne. 0._dbl) return
999 iflag = 2

  end subroutine banfac


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

    use numeric

    implicit none

    integer, intent(in)                               :: nroww,nrow
    integer, intent(in)                               :: nbandl,nbandu
    real(kind=dbl), dimension(nroww,nrow), intent(in) :: w
    real(kind=dbl), dimension(nrow), intent(inout)    :: b

    integer :: middle, nrowm1, jmax, i, j

    middle = nbandu + 1
    if (nrow .eq. 1) goto 99
    nrowm1 = nrow - 1
    if (nbandl .eq. 0) goto 30

    do i = 1, nrowm1
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          b(i+j) = b(i+j) - b(i) * w(middle+j,i)
       end do
    end do

30  if (nbandu .gt. 0)  goto 50

    do i = 1, nrow
       b(i) = b(i) / w(1,i)
    end do

    return

50  do i = nrow, 2, -1
       b(i) = b(i)/w(middle,i)
       jmax = min0(nbandu,i-1)
       do j = 1, jmax
          b(i-j) = b(i-j) - b(i) * w(middle-j,i)
       end do
    end do

99  b(1) = b(1) / w(middle,1)

  end subroutine banslv


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine huntn(xx,n,kord,x,jlo)

    use numeric

    implicit none

    integer, intent(in)                      :: n, kord
    real(kind=dbl), intent(in)               :: x
    real(kind=dbl), dimension(n), intent(in) :: xx

    integer, intent(inout)                     :: jlo

    integer :: max, null, jhi, jm, inc

!
!     works only for B-Splines (order n)
!

    max  = n - kord
    null = kord

    if (jlo.le.null.or.jlo.gt.max) then
       jlo = null
       jhi = max+1
       goto 30
    endif

    inc = 1

    if (x .ge. xx(jlo)) then
10     jhi = jlo + inc
       if (jhi .gt. max) then
          jhi = max + 1
       else if (x .ge. xx(jhi)) then
          jlo = jhi
          inc = inc + inc
          goto 10
       endif
    else
       jhi = jlo
20     jlo = jhi - inc
       if (jlo .le. null) then
          jlo = null
       else if (x .lt. xx(jlo)) then
          jhi = jlo
          inc = inc + inc
          goto 20
       endif
    endif

30  if (jhi-jlo.eq.1) return

    jm = (jhi + jlo) / 2
    if (x .gt. xx(jm)) then
       jlo = jm
    else
       jhi = jm
    endif

    goto 30

  end subroutine huntn


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function dbs2dr0(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef,nintx,ninty)

!
!  Evaluates the derivative of a two-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the derivative in the x-direction.  (input)
!   idery  - order of the derivative in the y-direction.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   bcoef  - array of length nx*ny containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny.
!   dbs2dr0  - value of the (iderx,idery) derivative of the spline at
!            (x,y).  (output)
!



!      --------------------------------------------------------------------
!EFD   same as dbs2dr except nintx,ninty is passed in and no error checking
!      --------------------------------------------------------------------

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, idery
    integer, intent(in)                          :: kx, nx, ky, ny
    real(kind=dbl)                               :: dbs2dr0
    real(kind=dbl), intent(in)                   :: x, y
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in) :: yknot
    real(kind=dbl), dimension(nx,ny), intent(in) :: bcoef

    integer, intent(in)                          :: nintx,ninty

    integer                       :: ix, iy, iky 
    real(kind=dbl), dimension(ky) :: work

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

!EFD    nintx = 0
!EFD
!EFD    do ix = 1, nx+kx-1
!EFD       if (xknot(ix) .gt. xknot(ix+1)) then
!EFD          write(6,*) "subroutine dbs2dr0:"
!EFD          write(6,*) "xknot(ix) <= xknot(ix+1) required."
!EFD          write(6,*) ix, xknot(ix), xknot(ix+1)
!EFD          stop
!EFD       endif
!EFD       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) nintx = ix
!EFD    end do
!EFD
!EFD    if(nintx .eq. 0) then
!EFD       write(6,*) "subroutine dbs2dr0:"
!EFD       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
!EFD       write(6,*) "x = ", x
!EFD       stop
!EFD    endif
!EFD
!EFD    ninty = 0
!EFD
!EFD    do iy = 1, ny+ky-1
!EFD       if (yknot(iy) .gt. yknot(iy+1)) then
!EFD          write(6,*) "subroutine dbs2dr0:"
!EFD          write(6,*) "yknot(iy) <= yknot(iy+1) required."
!EFD          write(6,*) iy, yknot(iy), yknot(iy+1)
!EFD          stop
!EFD       endif
!EFD       if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) ninty = iy
!EFD    end do
!EFD
!EFD    if(ninty .eq. 0) then
!EFD       write(6,*) "subroutine dbs2dr0:"
!EFD       write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
!EFD       write(6,*) "y = ", y
!EFD       stop
!EFD    endif

    do iky = 1, ky
       work(iky) =  dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iky),nintx)
    end do

         dbs2dr0 = dbsder0(idery,y,ky,yknot(ninty-ky+1),ky,work)

  end function dbs2dr0

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs3dr0(iderx,idery,iderz,x,y,z,kx,ky,kz,xknot,yknot,zknot,         &
       & nx,ny,nz,bcoef,   nintx,ninty,nintz)
!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   z      - z-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   dbs3dr0 - value of the (iderx,idery,iderz) derivative of the
!            spline at (x,y,z).  (output)
!

    use numeric

    implicit none

    integer, intent(in)                              :: iderx, idery, iderz
    integer, intent(in)                              :: nx, ny, nz, kx, ky, kz
    real(kind=dbl), intent(in)                       :: x, y, z
    real(kind=dbl), dimension(nx+kx), intent(in)     :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)     :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)     :: zknot
    real(kind=dbl), dimension(nx,ny,nz), intent(in)  :: bcoef
    real(kind=dbl)                                   :: dbs3dr0


!EFD dbs3dr but without checking and nintx,ninty,nintz passed in

    integer, intent(in)                              :: nintx,ninty,nintz

    integer                       :: iz 
    real(kind=dbl), dimension(kz) :: work
    logical                       :: isok
    logical, parameter            :: need_check = .false.

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

!EFD    nintz = 0
!EFD
!EFD    do iz = 1, nz+kz-1
!EFD       if (zknot(iz) .gt. zknot(iz + 1)) then
!EFD          write(6,*) "subroutine dbs3vl:"
!EFD          write(6,*) "zknot(iz) <= zknot(iz+1) required."
!EFD          write(6,*) iz, zknot(iz), zknot(iz+1)
!EFD          stop
!EFD       endif
!EFD       if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
!EFD    end do
!EFD
!EFD    if(nintz .eq. 0) then
!EFD       write(6,*) "subroutine dbs3dr0:"
!EFD       write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
!EFD       write(6,*) "zknot(iz)   = ", zknot(iz)
!EFD       write(6,*) "  z         = ", z
!EFD       write(6,*) "zknot(iz+1) = ", zknot(iz+1)
!EFD       stop
!EFD    endif

!EFD    if (need_check) then
!EFD
!EFD    isok = (1.le.nintx).and.(nintx.lt.nx+kx).and.                        &
!EFD       &   (xknot(nintx).le.x).and.(x.le.xknot(nintx+1))
!EFD    if (.not.isok) then
!EFD        write(6,*) 'nintx, x ',nintx,x
!EFD        write(6,*) 'xknot(nintx),xknot(nintx+1)',                        &
!EFD       &            xknot(nintx),xknot(nintx+1)
!EFD    endif
!EFD
!EFD    isok = (1.le.ninty).and.(ninty.lt.ny+ky).and.                        &
!EFD       &   (yknot(ninty).le.y).and.(y.le.yknot(ninty+1))
!EFD    if (.not.isok) then
!EFD        write(6,*) 'ninty, y ',ninty,y
!EFD        write(6,*) 'yknot(ninty),yknot(ninty+1)',                        &
!EFD       &            yknot(ninty),yknot(ninty+1)
!EFD    endif
!EFD
!EFD    isok = (1.le.nintz).and.(nintz.lt.nz+kz).and.                        &
!EFD       &   (zknot(nintz).le.z).and.(z.le.zknot(nintz+1))
!EFD    if (.not.isok) then
!EFD        write(6,*) 'nintz, z ',nintz,z
!EFD        write(6,*) 'zknot(nintz),zknot(nintz+1)',                        &
!EFD       &            zknot(nintz),zknot(nintz+1)
!EFD    endif
!EFD
!EFD
!EFD    endif

   
    do iz = 1, kz
       work(iz) = dbs2dr0(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,             &
             &        bcoef(1,1,nintz-kz+iz),    nintx,ninty)
    end do

     dbs3dr0 = dbsder0(iderz,z,kz,zknot(nintz-kz+1),kz,work)

  end function dbs3dr0


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsder0(iderx,x,kx,xknot,nx,bcoef )

!
!  Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsder0 - value of the iderx-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, kx, nx
    real(kind=dbl)                               :: dbsder0
    real(kind=dbl), intent(in)                   :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef


!     ------------------------------------
!EFD  no error checking
!     ------------------------------------


    integer                       :: ix, ik, il, leftx
    real(kind=dbl)                :: save, save1, save2, y, sum, dik
    real(kind=dbl), dimension(kx) :: work, dl, dr,bsp
    logical                       :: isok
    integer                       :: nxpkxm1

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

!EFD    leftx = 0
!EFD    do ix = 1,nx+kx-1
!EFD       if (xknot(ix) .gt. xknot(ix+1)) then
!EFD          write(6,*) "subroutine dbsder0:"
!EFD          write(6,*) "xknot(ix) <= xknot(ix+1) required."
!EFD          stop
!EFD       endif
!EFD       if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
!EFD    end do
!EFD
!EFD    if (leftx .eq. 0) then
!EFD       write(6,*) "subroutine dbsder0:"
!EFD       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
!EFD       write(6,*) "xknot(1)     = ", xknot(1)
!EFD       write(6,*) "xknot(nx+kx) = ", xknot(nx+kx)
!EFD       write(6,*) "         x   = ", x
!EFD       stop
!EFD    endif



    leftx = 0
    nxpkxm1 = nx+kx-1
   select case(nxpkxm1)
   case (1)
    do ix=1,1
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (2)
    do ix=1,2
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (3)
    do ix=1,3
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (4)
    do ix=1,4
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (5)
    do ix=1,5
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (6)
    do ix=1,6
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (7)
    do ix=1,7
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (8)
    do ix=1,8
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (9)
    do ix=1,9
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case (10)
    do ix=1,10
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   case default
    do ix=1,nxpkxm1
       if ((xknot(ix).le.x).and.(x.le.xknot(ix+1)) ) then
          leftx = ix
          exit
       endif
    enddo
   end select


  dbsder0 =  dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)

  return
  end function dbsder0



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs3vd(iderx,idery,iderz,nxyz,xvec,yvec,zvec,                   &
       & kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,val)

!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation on a grid.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   xvec   - array of length nx containing the x-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in xvec should be strictly increasing.
!   yvec   - array of length ny containing the y-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   zvec   - array of length nz containing the z-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   val    - array of size nxyz containing the values of
!            the (iderx,idery,iderz) derivative of the spline (output)
!            value(ijk) contains the derivative of the spline at
!            the point (xvec(ijk), yvec(ijk), zvec(ijk)).
!

    use numeric

    implicit none

    integer, intent(in)                               :: iderx, idery, iderz
    integer, intent(in)                               :: nxyz
    integer, intent(in)                               :: kx, nx, ky, ny, kz, nz

    real(kind=dbl), dimension(nxyz), intent(in)      :: xvec
    real(kind=dbl), dimension(nxyz), intent(in)      :: yvec
    real(kind=dbl), dimension(nxyz), intent(in)      :: zvec
    real(kind=dbl), dimension(nx+kx), intent(in)      :: xknot
    real(kind=dbl), dimension(ny+ky), intent(in)      :: yknot
    real(kind=dbl), dimension(nz+kz), intent(in)      :: zknot
    real(kind=dbl), dimension(nx,ny,nz), intent(in)   :: bcoef
    real(kind=dbl), dimension(nxyz), intent(out)      :: val

    integer                                           :: i, ik, il, ix, iy, iz
    integer                                           :: ikx, iky, ikz
    integer, dimension(nxyz)                         :: leftx
    integer, dimension(nxyz)                         :: lefty
    integer, dimension(nxyz)                         :: leftz
    real(kind=dbl), dimension(nxyz,kx)               :: biatx
    real(kind=dbl), dimension(nxyz,ky)               :: biaty
    real(kind=dbl), dimension(nxyz,kz)               :: biatz
    real(kind=dbl), dimension(nxyz) :: term, save1

    real(kind=dbl), dimension(nxyz, max(kx,ky,kz)) :: dl, dr


    logical :: use_huntn = .false.
    logical :: use_linear = .false.
    logical :: use_vsearch = .false.
    logical :: use_hash = .true.
    logical :: isok, is_uniform
    integer :: ileft, ijk, istart,iend
    integer :: t1,t2,count_rate
    real(kind=dbl) :: delta
    real(kind=dbl), parameter  :: tol=1.0d-7

    integer :: nintx,ninty,nintz
    real(kind=dbl),dimension(kz) :: workz
    real(kind=dbl),dimension(ky) :: worky

    logical, parameter :: use_Axyz = .false.
    real(kind=dbl),dimension(kx,ky,kz) :: Axyz
    real(kind=dbl),dimension(ky,kz) :: Ayz
    real(kind=dbl),dimension(kz) :: Az

!   -----------------------------
!   double check knot positions
!   -----------------------------
   call system_clock(t1,count_rate)
   isok = all(xknot(2:(nx+kx)) .ge. xknot(1:(nx+kx-1)))
   if (.not.isok) then
    do i = 1, nx+kx-1
       if (xknot(i) .gt. xknot(i+1)) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "xknot(i) <= xknot(i+1) required."
          write(6,*) i, xknot(i), xknot(i+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
     end do
    endif


   isok = all(yknot(2:(ny+ky)) .ge. yknot(1:(ny+ky-1)))
   if (.not.isok) then
    do i = 1, ny+ky-1
       if (yknot(i) .gt. yknot(i+1)) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "yknot(i) <= yknot(i+1) required."
          write(6,*) i, yknot(i), yknot(i+1)
          write(6,*)
          write(6,*) yknot
          stop
       endif
     end do
    endif


   isok = all(zknot(2:(nz+kz)) .ge. zknot(1:(nz+kz-1)))
   if (.not.isok) then
    do i = 1, nz+kz-1
       if (zknot(i) .gt. zknot(i+1)) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "zknot(i) <= zknot(i+1) required."
          write(6,*) i, zknot(i), zknot(i+1)
          write(6,*)
          write(6,*) zknot
          stop
       endif
     end do
    endif

   call system_clock(t2,count_rate)
   if (idebug.ge.2) then
   write(*,*) 'time for 1st check ',real(t2-t1)/real(count_rate)
   endif

!   ----------------------------
!   find the location of indices
!   ----------------------------

   call system_clock(t1,count_rate)

   if (use_huntn) then
   do i=1,nxyz
     ileft = 0
     call huntn(xknot,nx+kx,kx,xvec(i),ileft)
     leftx(i) = ileft
   enddo
   do i=1,nxyz    
     ileft = 0
     call huntn(yknot,ny+ky,ky,yvec(i),ileft)
     lefty(i) = ileft
   enddo
   do i=1,nxyz    
     ileft = 0
     call huntn(zknot,nz+kz,kz,zvec(i),ileft)
     leftz(i) = ileft
   enddo
   else if (use_linear) then

    leftx(1:nxyz) = 0
    lefty(1:nxyz) = 0
    leftz(1:nxyz) = 0
    do ijk=1,nxyz
       do ix=1,nx+kx-1
         if ((xknot(ix).le.xvec(ijk)).and.(xvec(ijk).le.xknot(ix+1)))then
            leftx(ijk) = ix
            exit
         endif
       enddo

       do iy=1,ny+ky-1
         if ((yknot(iy).le.yvec(ijk)).and.(yvec(ijk).le.yknot(iy+1)))then
            lefty(ijk) = iy
            exit
         endif
       enddo

       do iz=1,nz+kz-1
         if ((zknot(iz).le.zvec(ijk)).and.(zvec(ijk).le.zknot(iz+1)))then
            leftz(ijk) = iz
            exit
         endif
       enddo
     enddo
    else if (use_vsearch) then
     if (idebug.ge.2) then
       write(*,*) 'using vsearch '
     endif

       call vsearch(nx+kx,xknot, nxyz, xvec, leftx )
       call vsearch(ny+ky,yknot, nxyz, yvec, lefty )
       call vsearch(nz+kz,zknot, nxyz, zvec, leftz )
    else if (use_hash) then

! ---------------------------------------------------
! geometric hashing works only if the grid is uniform
! ---------------------------------------------------

       istart = kx+1
       iend = nx
       delta = xknot(istart+1)-xknot(istart)
       is_uniform = all((abs(xknot((istart+1):iend)-                       &
     &                  (delta+xknot(istart:(iend-1))))).lt.tol)
       if (is_uniform) then
          call vhash(xknot,nx,kx,nxyz,xvec,leftx)
       else
          call vsearch(nx+kx,xknot, nxyz,xvec,leftx)
       endif



       istart = ky+1
       iend = ny
       delta = yknot(istart+1)-yknot(istart)
       is_uniform = all((abs(yknot((istart+1):iend)-                       &
     &                  (delta+yknot(istart:(iend-1))))).lt.tol)
       if (is_uniform) then
          call vhash(yknot,ny,ky,nxyz,yvec,lefty)
       else
          call vsearch(ny+ky,yknot, nxyz,yvec,lefty)
       endif


       istart = kz+1
       iend = nz
       delta = zknot(istart+1)-zknot(istart)
       is_uniform = all((abs(zknot((istart+1):iend)-                       &
     &                  (delta+zknot(istart:(iend-1))))).lt.tol)
       if (is_uniform) then
          call vhash(zknot,nz,kz,nxyz,zvec,leftz)
       else
          call vsearch(nz+kz,zknot, nxyz,zvec,leftz)
       endif

    endif



   call system_clock(t2,count_rate)
   if (idebug.ge.2) then
   write(*,*) 'time for index search ',real(t2-t1)/real(count_rate)
   endif



!  ------------
!  double check
!  ------------
   call system_clock(t1,count_rate)

             
   isok = .false.
   if (all((1.le.leftx(1:nxyz)).and.(leftx(1:nxyz).lt.nx+kx))) then
      isok =  all((xknot( leftx(1:nxyz) ) .le. xvec(1:nxyz)).and.             &
          &   (xvec(1:nxyz).le.xknot( 1+leftx(1:nxyz))) )
   endif
   if (.not.isok) then
       do i=1,nxyz
         ileft = leftx(i)
         isok = (1.le.ileft).and.(ileft.lt.nx+kx)
         if (isok) then
           isok= (xknot( ileft ).le. xvec(i)).and.                             &
     &           (xvec(i).le.xknot(ileft+1))
         endif

         if (.not.isok) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(i)
          write(6,*) "i,ileft ",i,ileft
          write(6,*) 'nx,kx ',nx,kx
          do ix=1,nx+kx
            write(*,*) 'ix,xknot(ix) ',ix,xknot(ix)
          enddo
          stop '** error in dbs3vd ** '
         endif
       enddo
    endif

   isok = .false.
   if (all((1.le.lefty(1:nxyz)).and.(lefty(1:nxyz).lt.ny+ky))) then
      isok =  all((yknot( lefty(1:nxyz) ) .le. yvec(1:nxyz)).and.             &
     &            (yvec(1:nxyz).le.yknot( 1+lefty(1:nxyz))) )
   endif
   if (.not.isok) then
       do i=1,nxyz
         ileft = lefty(i)
         isok = (1.le.ileft).and.(ileft.lt.ny+ky)
         if (isok) then
           isok=(yknot( ileft ).le. yvec(i)).and.                             &
     &          (yvec(i).le.yknot(ileft+1))
         endif

         if (.not.isok) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
          write(6,*) "y = ", yvec(i)
          write(6,*) "i,ileft ",i,ileft
          write(6,*) 'ny,ky ',ny,ky
          do iy=1,ny+ky
            write(*,*) 'iy,yknot(iy) ',iy,yknot(iy)
          enddo
          stop '** error in dbs3vd ** '
         endif
       enddo
    endif

   isok = .false.
   if (all((1.le.leftz(1:nxyz)).and.(leftz(1:nxyz).lt.nz+kz))) then
     isok =  all((zknot( leftz(1:nxyz) ) .le. zvec(1:nxyz)).and.             &
          &   (zvec(1:nxyz).le.zknot( 1+leftz(1:nxyz))) )
   endif
   if (.not.isok) then
       do i=1,nxyz
         ileft = leftz(i)
         isok = (1.le.ileft).and.(ileft.lt.nz+kz)
         if (isok) then
           isok= (zknot( ileft ).le. zvec(i)).and.                             &
     &           (zvec(i).le.zknot(ileft+1))
         endif

         if (.not.isok) then
          write(6,*) "subroutine dbs3vd:"
          write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
          write(6,*) "z = ", zvec(i)
          write(6,*) "i,ileft ",i,ileft
          write(6,*) 'nz,kz ',nz,kz
          do iz=1,nz+kz
            write(*,*) 'iz,zknot(iz) ',iz,zknot(iz)
          enddo
          stop '** error in dbs3vd ** '
         endif
       enddo
    endif

   call system_clock(t2,count_rate)
   if (idebug.ge.2) then
   write(*,*) 'time for 2nd check ',real(t2-t1)/real(count_rate)
   endif


    if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

       do ix = 1, nxyz
          biatx(ix,1) = 1.0_dbl
       end do
       do ix=1,nxyz
         nintx = leftx(ix)
         do ik=1,kx-1
             dr(ix,ik) = xknot(nintx+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(nintx+1-ik)
         enddo
       enddo
       

       do ik = 1, kx-1
          do ix = 1, nxyz
             save1(ix) = 0._dbl
          end do

          do il = 1, ik
             do ix = 1, nxyz
                term(ix)     = biatx(ix,il) / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxyz
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do iy = 1, nxyz
          biaty(iy,1) = 1.0_dbl
       end do
       do iy = 1, nxyz
         ninty = lefty(iy)
         do ik = 1, ky-1
             dr(iy,ik) = yknot(ninty+ik) - yvec(iy)
             dl(iy,ik) = yvec(iy) - yknot(ninty+1-ik)
         enddo
       end do
   

       do ik = 1, ky-1
          do iy = 1, nxyz
             save1(iy) = 0._dbl
          end do

          do il = 1,ik
             do iy = 1,nxyz
                term(iy)     = biaty(iy,il) / (dr(iy,il) + dl(iy,ik+1-il))
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                save1(iy)    = dl(iy,ik+1-il) * term(iy)
             end do
          end do

          do iy = 1,nxyz
             biaty(iy,ik+1) = save1(iy)
          end do
       end do

       do iz = 1,nxyz
          biatz(iz,1) = 1.0_dbl
       end do
       do iz = 1, nxyz
         nintz = leftz(iz)
         do ik = 1, kz-1
             dr(iz,ik) = zknot(nintz+ik) - zvec(iz)
             dl(iz,ik) = zvec(iz) - zknot(nintz+1-ik)
         enddo
       enddo

       do ik = 1, kz-1
          do iz = 1, nxyz
             save1(iz) = 0._dbl
          end do

          do il = 1, ik
             do iz = 1, nxyz
                term(iz)     = biatz(iz,il) / (dr(iz,il) + dl(iz,ik+1-il))
                biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
                save1(iz)    = dl(iz,ik+1-il) * term(iz)
             end do
          end do

          do iz = 1, nxyz
             biatz(iz,ik+1) = save1(iz)
          end do
       end do


       val(1:nxyz) = 0.0_dbl

!EFD    do ikz = 1, kz
!EFD       do iky = 1, ky
!EFD          do ikx = 1, kx
!EFD             do ijk=1,nxyz
!EFD                      ix = ijk
!EFD                      iy = ijk
!EFD                      iz = ijk

!EFD                      val(ijk) = val(ijk)                                  &
!EFD                           &  + biatx(ix,ikx) * biaty(iy,iky)              &
!EFD                           &  * biatz(iz,ikz)                              &
!EFD                           &  * bcoef(leftx(ix)-kx+ikx,                    &
!EFD                           &          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
!EFD             end do
!EFD          end do
!EFD       end do
!EFD    end do


       call system_clock(t1,count_rate)

       if (use_Axyz) then


       do ijk=1,nxyz
          ix = ijk
          iy = ijk
          iz = ijk
          nintx = leftx(ix)
          ninty = lefty(iy)
          nintz = leftz(iz)

          do ikz=1,kz
          do iky=1,ky
          do ikx=1,kx
             Axyz(ikx,iky,ikz) = bcoef(nintx-kx+ikx,ninty-ky+iky,nintz-kz+ikz)
          enddo
          enddo
          enddo

          do ikz=1,kz
          do iky=1,ky
              Ayz(iky,ikz) = 0.0d0
              do ikx=1,kx
                 Ayz(iky,ikz) = Ayz(iky,ikz) + Axyz(ikx,iky,ikz)*biatx(ix,ikx)
              enddo
          enddo
          enddo

          do ikz=1,kz
              Az(ikz) = 0.0d0
              do iky=1,ky
                Az(ikz) = Az(ikz) + Ayz(iky,ikz) * biaty(iy,iky)
              enddo
          enddo

          do ikz=1,kz
             val(ijk) = val(ijk) + Az(ikz)*biatz(iz,ikz)
          enddo

       enddo

       else


       do ijk=1,nxyz
	 ix = ijk
	 iy = ijk
	 iz = ijk
         nintx = leftx(ix)
         ninty = lefty(iy)
         nintz = leftz(iz)
         do ikz = 1, kz
         do iky = 1, ky
         do ikx = 1, kx

	 val(ijk) = val(ijk) +                                            &
     &    bcoef(nintx-kx+ikx,ninty-ky+iky,nintz-kz+ikz)*biatx(ix,ikx)*    &
     &         biaty(iy,iky)*biatz(iz,ikz)


         end do
         end do
         end do
       end do

       endif
    

       call system_clock(t2,count_rate)
       if (idebug.ge.2) then
       write(*,*) 'time in val loop ',real(t2-t1)/real(count_rate)
       endif

    else


       do ijk=1,nxyz
          ix = ijk
          iy = ijk
          iz = ijk

          val(ijk) = dbs3dr0(iderx,idery,iderz,xvec(ix),               &
                     &  yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,        &
                     &  zknot,nx,ny,nz,bcoef,                          &
                     &  leftx(ix),lefty(iy),leftz(iz))
       enddo


    endif
    return
  end subroutine dbs3vd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	subroutine vsearch( n, x, m, vxc, vipos )
	use numeric
	implicit none
	integer, intent(in) :: n
	real(kind=dbl), intent(in), dimension(n) :: x
        integer, intent(in) :: m
	real(kind=dbl), intent(in), dimension(m) :: vxc
	integer, intent(inout), dimension(m) :: vipos

!
!	find  ipos such that x(ipos) <= xc <= x(ipos+1)
!
	integer :: ipos,i,j, ilo,ihi,imid, nit
	integer, dimension(1) :: ipos1
	intrinsic :: log
	real(kind=dbl) :: xc
	integer, parameter :: k64=64*1024
	integer, parameter :: k128=128*1024
	integer, parameter :: k256=256*1024
	integer, parameter :: k512=512*1024
	integer, parameter :: m1=1024*1024
	integer, parameter :: m2=2*m1
	integer, parameter :: m4=4*m1
	integer, parameter :: m8=8*m1
	integer, parameter :: m16=16*m1
	integer, parameter :: m32=32*m1
	integer, parameter :: m64=64*m1

	integer itable(25)
        save
	data itable /2,4,8,16,32,                                            &
     &    64,128,256,512,1024,                                               &
     &    2048, 4096, 8192,16382,  32768,                                    &
     &    k64,k128,k512,m1,m2,                                               &
     &    m4,m8,m16,m32,m64/
       
	ipos1 = minloc(  itable - n, itable .ge. n )
	nit = ipos1(1) + 1

        do j=1,m
           xc = vxc(j)
	   ilo = 1
	   ihi = n
           if ((xc .le. x(1)).or.(xc.ge.x(n))) then
		ipos = 0
	   else
!		bineary search
		ihi = n
		ilo = 1
		ipos = 0
		do i=1,nit
		  imid = (ihi + ilo)/2
		  if ((x(imid).le.xc).and.(xc.le.x(imid+1))) then
                     ipos = imid
		     exit
		  else if (xc.le.x(imid)) then
		      ihi = imid
		  else if (xc.ge.x(imid)) then
			ilo = imid
		  endif
		enddo
	     endif
	   vipos(j) = ipos
	  enddo

	return
	end subroutine vsearch


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine vhash(zknot,nz,kz,nxyz,zvec,leftz)
	use numeric
	implicit none
        integer, intent(in) :: nz,kz,nxyz
        real(kind=dbl),dimension(nz+kz),intent(in) :: zknot
        real(kind=dbl),dimension(nxyz), intent(in) :: zvec
        integer, dimension(nxyz),intent(inout) :: leftz

        integer :: ijk,istart,iend,i,im1,ip1,ileft
        real(kind=dbl) :: z,dz

        istart = kz+1
        iend = nz
        dz = zknot(istart+1)-zknot(istart)

        do ijk=1,nxyz
           z = zvec(ijk)
!          -----------------------------------
!          uniform grid in zknot(istart:iend)
!
!          zknot(i) <= z <= zknot(i+1)
!          zknot(istart) + (i-istart)*dz <= z <= zknot(istart) + (i+1-istart)*dz
!          (i-istart) <= (z-zknot(istart))/dz <= i+1-istart
!
!          -----------------------------------
           ileft = int((z-zknot(istart))/dz) + istart
           leftz(ijk)  = max(istart,min(iend-1,ileft) )
        enddo
        do ijk=1,nxyz
           i = leftz(ijk)
           z = zvec(ijk)

           if (z.lt.zknot(i)) then
              if ((zknot(istart-1).le.z).and.(z.le.zknot(istart))) then
                  leftz(ijk) = istart-1
              else
                 im1 = max(istart,min(iend-1,i-1))
                 if ((zknot(im1).le.z).and.(z.le.zknot(im1+1))) then
                     leftz(ijk) = im1
                 endif
              endif
           else if  (z.gt.zknot(i+1)) then
              if ((zknot(iend).le.z).and.(z.le.zknot(iend+1))) then
                  leftz(ijk) = iend
              else
                   ip1 = max(istart,min(iend-1,i+1))
                   if ((zknot(ip1).le.z).and.(z.le.zknot(ip1+1))) then
                      leftz(ijk) = ip1
                   endif
              endif
           endif
         enddo
         return
         end subroutine vhash
           



end module bspline
