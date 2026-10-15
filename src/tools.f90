MODULE tools

  USE param

  IMPLICIT NONE

  PRIVATE :: uertst
  PUBLIC  :: findf, zeroin, spline, seval, icsscu, quadq, erf, zzdisp, wzdisp, numinc, elapsed_time

CONTAINS

  function findf(f)

    use findrf

    implicit none

    real(IDP) :: f,findf,xn,xndx,fnt,sum

    xn=n0-nt+1
    xndx=xn*dx
    fnt=f**nt
    if (f /= 1.) sum=f*(1.-fnt)/(1.-f)
    if (f == 1.) sum=nt
    findf=dx*sum+xndx*fnt-d

  end function findf

  function zeroin(ax,bx,tol)

    implicit none

    real(IDP) :: ax,bx,tol,zeroin

    !  a zero of the function  f(x)  is computed in the interval ax,bx .

    !  input..

    !  ax     left endpoint of initial interval
    !  bx     right endpoint of initial interval
    !  f      function subprogram which evaluates f(x) for any x in
    !     the interval  ax,bx
    !  tol    desired length of the interval of uncertainty of the
    !     final result ( .ge. 0.0)


    !  output..

    !  zeroin abcissa approximating a zero of  f  in the interval ax,bx


    !  it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
    !  without  a  check.  zeroin  returns a zero  x  in the given interval
    !  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
    !  is the relative machine precision.
    !  this function subprogram is a slightly  modified  translation  of
    !  the algol 60 procedure  zero  given in  richard brent, algorithms for
    !  minimization without derivatives, prentice - hall, inc. (1973).


    real(IDP) :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s

    !  compute eps, the relative machine precision

    eps = 1.0
    tol1 = 2.0
    do while (tol1 > 1.0)
       eps = eps/2.0
       tol1 = 1.0 + eps
    end do

    !  initialization

    a = ax
    b = bx
    fa = findf(a)
    fb = findf(b)

    !  begin step

    c = a
    fc = fa
    d = b - a
    e = d
    do while (.true.)
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if

       !  convergence test

       tol1 = 2.0*eps*abs(b) + 0.5*tol
       xm = .5*(c - b)
       if (abs(xm) <= tol1) exit
       if (fb == 0.0) exit

       !  is bisection necessary

       if (abs(e) < tol1 .or. abs(fa) <= abs(fb)) then

          !  bisection

          d = xm
          e = d

       else

          !  is quadratic interpolation possible

          if (a == c) then

             !  linear interpolation

             s = fb/fa
             p = 2.0*xm*s
             q = 1.0 - s

          else

             !  inverse quadratic interpolation

             q = fa/fc
             r = fb/fc
             s = fb/fa
             p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
             q = (q - 1.0)*(r - 1.0)*(s - 1.0)

          end if

          !  adjust signs

          if (p > 0.0) q = -q
          p = abs(p)

          !  is interpolation acceptable

          if ((2.0*p) >= (3.0*xm*q - abs(tol1*q)) .or. p >= abs(0.5*e*q)) then
             d = xm
             e = d
          else
             e = d
             d = p/q
          end if

       end if

       !  complete step

       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1, xm)
       end if
       fb = findf(b)
       if ((fb*(fc/abs(fc))) > 0.0) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
    end do

    !  done

    zeroin = b

  end function zeroin

  subroutine spline (n, x, y, b, c, d)

    implicit none

    integer :: n
    real(IDP), dimension(:) :: x, y, b, c, d

    !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
    !  for a cubic interpolating spline

    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

    !  for  x(i) .le. x .le. x(i+1)

    !  input..

    !  n = the number of data points or knots (n.ge.2)
    !  x = the abscissas of the knots in strictly increasing order
    !  y = the ordinates of the knots

    !  output..

    !  b, c, d  = arrays of spline coefficients as defined above.

    !  using  p  to denote differentiation,

    !  y(i) = s(x(i))
    !  b(i) = sp(x(i))
    !  c(i) = spp(x(i))/2
    !  d(i) = sppp(x(i))/6  (derivative from the right)

    !  the accompanying function subprogram  seval  can be used
    !  to evaluate the spline.


    integer :: nm1, ib, i
    real(IDP) :: t

    nm1 = n-1

    if ( n > 2 ) then

       !  set up tridiagonal system

       !  b = diagonal, d = offdiagonal, c = right hand side.

       d(1) = x(2) - x(1)
       c(2) = (y(2) - y(1))/d(1)
       do i = 2, nm1
          d(i) = x(i+1) - x(i)
          b(i) = 2.*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i)
          c(i) = c(i+1) - c(i)
       end do

       !  end conditions.  third derivatives at  x(1)  and  x(n)
       !  obtained from divided differences

       b(1) = -d(1)
       b(n) = -d(n-1)
       c(1) = 0.
       c(n) = 0.
       if ( n > 3 ) then
          c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
          c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
          c(1) = c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
       end if

       !  forward elimination

       do i = 2, n
          t = d(i-1)/b(i-1)
          b(i) = b(i) - t*d(i-1)
          c(i) = c(i) - t*c(i-1)
       end do

       !  back substitution

       c(n) = c(n)/b(n)
       do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       end do

       !  c(i) is now the sigma(i) of the text

       !  compute polynomial coefficients

       b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
       do i = 1, nm1
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
          d(i) = (c(i+1) - c(i))/d(i)
          c(i) = 3.*c(i)
       end do
       c(n) = 3.*c(n)
       d(n) = d(n-1)

    else if (n == 2) then

       b(1) = (y(2)-y(1))/(x(2)-x(1))
       c(1) = 0.
       d(1) = 0.
       b(2) = b(1)
       c(2) = 0.
       d(2) = 0.

    end if

  end subroutine spline

  function seval(n, u, x, y, b, c, d)

    implicit none

    integer :: n
    real(IDP) :: seval,u
    real(IDP), dimension(:) :: x, y, b, c, d

    !  this subroutine evaluates the cubic spline function

    !  seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3

    !  where  x(i) .lt. u .lt. x(i+1), using horner's rule

    !  if  u .lt. x(1) then  i = 1  is used.
    !  if  u .ge. x(n) then  i = n  is used.

    !  input..

    !  n = the number of data points
    !  u = the abscissa at which the spline is to be evaluated
    !  x,y = the arrays of data abscissas and ordinates
    !  b,c,d = arrays of spline coefficients computed by spline

    !  if  u  is not in the same interval as the previous call, then a
    !  binary search is performed to determine the proper interval.

    integer :: j, k
    real(IDP) :: dx
    integer, save :: i = 1

    if ( i >= n ) i = 1
    if ( u < x(i) .or. u > x(i+1) ) then

       !  binary search

       i = 1
       j = n+1
       do while ( j > i+1 )
          k = (i+j)/2
          if ( u < x(k) ) then
             j = k
          else
             i = k
          end if
       end do

    end if

    !  evaluate spline

    dx = u - x(i)
    seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function seval

  subroutine icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)

    !  icsscu-------s-------library 2---------------------------------------

    !  function            - cubic spline data smoothing
    !  usage               - call icsscu(x,f,df,nx,sm,y,c,ic,wk,ier)
    !  parameters   x      - vector of length nx containing the abscissae
    !                    of the nx data points (x(i),f(i)) i=1,...,
    !                    nx (input). x must be ordered so that
    !                    x(i) .lt. x(i+1).
    !           f      - vector of length nx containing the ordinates
    !                    (or function values) of the nx data points
    !                    (input).
    !           df     - vector of length nx (input).
    !                    df(i) is the relative weight of data
    !                    point i (see parameter sm below).
    !           nx     - number of elements in x, f, df, and y (input).
    !                    nx must be .ge. 2.
    !           sm     - a non-negative number which controls the
    !                    extent of smoothing (input). the spline
    !                    function s is determined such that the
    !                    sum from 1 to nx of
    !                    ((s(x(i))-f(i))/df(i))**2 .le. sm,
    !                    where equality holds unless s describes
    !                    a straight line.
    !           y,c    - spline coefficients (output). y is a vector
    !                    of length nx. c is an nx-1 by 3 matrix.
    !                    the value of the spline approximation
    !                    at t is
    !                    s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
    !                    where x(i) .le. t .lt. x(i+1) and
    !                    d = t-x(i).
    !           ic     - row dimension of matrix c in the calling
    !                    program (input). ic must be .ge. nx-1.
    !           wk     - work area vector of length greater than or
    !                    equal to 7*nx+14.
    !           ier    - error parameter
    !                    terminal error
    !                    ier = 129. ic is less than nx-1.
    !                    ier = 130, nx is less than 2.
    !                    ier = 131, input abscissae are not ordered
    !                              so that x(1) .lt. x(2) ... .lt. x(nx).
    !  precision           - single
    !  req'd imsl routines - uertst
    !  language            - fortran

    !  latest revision     - february 16, 1976
    !                    dec

    implicit none

    integer :: nx,ic,ier,m2,np1,i,np3,j
    real(IDP) :: sm,p,h,f2,ff,g,onedh,e,hmg
    real(IDP), dimension(:) :: x,f,df,y
    real(IDP), dimension(:,:) :: c
    real(IDP), dimension(:,:) :: wk

    !  check error conditions
    ier = 0
    if (ic < nx-1) then
       ier = 129
       call uertst(ier,"icsscu")
    else if (nx < 2) then
       ier = 130
       call uertst(ier,"icsscu")
    else
       !  set up working areas
       m2 = nx+2
       np1 = nx+1
       wk(1,1) = 0.0
       wk(1,2) = 0.0
       wk(2,np1) = 0.0
       wk(3,m2) = 0.0
       wk(3,np1) = 0.0
       wk(6,1) = 0.0
       wk(6,2) = 0.0
       wk(6,m2) = 0.0
       wk(6,np1) = 0.0
       p = 0.0
       h = x(2)-x(1)
       if (h <= 0.0) then
          ier = 131
          call uertst(ier,"icsscu")
       else
          f2 = -sm
          ff = (f(2)-f(1))/h
          if (nx > 2) then
             do i=3,nx
                g = h
                h = x(i)-x(i-1)
                if (h <= 0.0) then
                   ier = 131
                   call uertst(ier,"icsscu")
                   return
                end if
                onedh = 1.0/h
                e = ff
                ff = (f(i)-f(i-1))*onedh
                y(i) = ff-e
                wk(4,i) = 2.*(g+h)/3.
                wk(5,i) = h/3.0
                wk(3,i) = df(i-2)/g
                wk(1,i) = df(i)*onedh
                wk(2,i) = -df(i-1)/g-df(i-1)*onedh
             end do
             do i=3,nx
                c(i-1,1) = wk(1,i)*wk(1,i)+wk(2,i)*wk(2,i)+wk(3,i)*wk(3,i)
                c(i-1,2) = wk(1,i)*wk(2,i+1)+wk(2,i)*wk(3,i+1)
                c(i-1,3) = wk(1,i)*wk(3,i+2)
             end do
          end if
          !  next iteration
          do while (.true.)
             if (nx > 2) then
                do i=3,nx
                   wk(2,i-1) = ff*wk(1,i-1)
                   wk(3,i-2) = g*wk(1,i-2)
                   wk(1,i) = 1.0/(p*c(i-1,1)+wk(4,i)-ff*wk(2,i-1)-g*wk(3,i-2))
                   wk(6,i) = y(i)-wk(2,i-1)*wk(6,i-1)-wk(3,i-2)*wk(6,i-2)
                   ff = p*c(i-1,2)+wk(5,i)-h*wk(2,i-1)
                   g = h
                   h = c(i-1,3)*p
                end do
                np3 = nx+3
                do i=3,nx
                   j = np3-i
                   wk(6,j) = wk(1,j)*wk(6,j)-wk(2,j)*wk(6,j+1)-wk(3,j)*wk(6,j+2)
                end do
             end if
             e = 0.0
             h = 0.0
             !  compute u and accumulate e
             do i=2,nx
                g = h
                h = (wk(6,i+1)-wk(6,i))/(x(i)-x(i-1))
                hmg = h-g
                wk(7,i) = hmg*df(i-1)*df(i-1)
                e = e+wk(7,i)*hmg
             end do
             g = -h*df(nx)*df(nx)
             wk(7,np1) = g
             e = e-g*h
             g = f2
             f2 = e*p*p
             if (f2 >= sm .or. f2 <= g) exit
             ff = 0.0
             h = (wk(7,3)-wk(7,2))/(x(2)-x(1))
             if (nx > 2) then
                do i=3,nx
                   g = h
                   h = (wk(7,i+1)-wk(7,i))/(x(i)-x(i-1))
                   g = h-g-wk(2,i-1)*wk(1,i-1)-wk(3,i-2)*wk(1,i-2)
                   ff = ff+g*wk(1,i)*g
                   wk(1,i) = g
                end do
             end if
             h = e-p*ff
             if (h <= 0.0) exit
             !  update the lagrange multiplier p for the next iteration
             p = p+(sm-f2)/((sqrt(sm/e)+p)*h)
          end do
          !  if e less than or equal to s, compute the coefficients and return.
          np1 = nx-1
          do i=1,np1
             y(i) = f(i)-p*wk(7,i+1)
             c(i,2) = wk(6,i+1)
             wk(1,i) = y(i)
          end do
          wk(1,nx) = f(nx)-p*wk(7,nx+1)
          y(nx) = wk(1,nx)
          do i=2,nx
             h = x(i)-x(i-1)
             c(i-1,3) = (wk(6,i+1)-c(i-1,2))/(h+h+h)
             c(i-1,1) = (wk(1,i)-y(i-1))/h-(h*c(i-1,3)+c(i-1,2))*h
          end do
       end if
    end if

  end subroutine icsscu

  subroutine uertst(ier,name)

    !  uertst---------------library 2---------------------------------------

    !  function            - error message generation
    !  usage               - call uertst(ier,name)
    !  parameters   ier    - error parameter. type + n  where
    !                       type= 128 implies terminal error
    !                              64 implies warning with fix
    !                              32 implies warning
    !                       n   = error code relevant to calling routine
    !           name   - input scalar (double precision on dec)
    !                    containing the name of the calling routine
    !                    as a 6-character literal string.
    !  language            - fortran

    !  latest revision     - october 1,1975
    !                    dec

    implicit none

    integer :: ier,ier1,ier2,i
    character(len=6) :: name
    integer, dimension(4) :: ibit=(/32,64,128,0/)
    character(len=5), dimension(4,4) :: ityp
    !  character(len=5), dimension(4,4) :: ityp= (/"warni","ng   ","     ","     ","warni","ng(wi","th fi","x)   ", &
    !                                             "termi","nal  ","     ","     ","non-d","efine","d    ","     "/)

    ityp(:,1)=(/"warni","ng   ","     ","     "/)
    ityp(:,2)=(/"warni","ng(wi","th fi","x)   "/)
    ityp(:,3)=(/"termi","nal  ","     ","     "/)
    ityp(:,4)=(/"non-d","efine","d    ","     "/)

    ier2=ier
    if (ier2 < ibit(1)) then
       !  non-defined
       ier1=4
    else if (ier2 < ibit(2)) then
       !  warning
       ier1=1
    else if (ier2 < ibit(3)) then
       !  warning(with fix)
       ier1=2
    else
       !  terminal
       ier1=3
    end if
    !  extract 'n'
    ier2=ier2-ibit(ier1)
    !  print error message
    write (6,'(" *** i m s l(uertst) ***  ",4a5,2x,a6,2x,i2," (ier = ",i3,")")') (ityp(i,ier1),i=1,4),name,ier2,ier

  end subroutine uertst

  subroutine quadq(n,x,y,x0,b,c,d,result)

    implicit none

    integer :: n
    real(IDP) :: x0,result
    real(IDP), dimension(:) :: x,y,b,c,d

    !  sub. quadq

    !  a fortran sub. for computing an approximation
    !  to the integral of a fcn. y(x), given values of
    !  y(x) at the n distinct points  x(1).lt.x(2)...lt.x(n).

    !  n.b. this routine has been modified (from quad2) so that
    !  y(x) satisfies the constraint d**2 y/d x**2 at x=x(1)=0.

    !  the sub. uses a cubic spline fcn. to
    !  interpolate y(x). the resultant spline is then
    !  integrated, from x=x(1) to x=x0, to obtain
    !  the required approximation.
    !  the coefficients of the interpolating spline
    !  are obtained using an in-line modified version
    !  of the sub. "spline" which is published in
    !  the book 'computer methods for mathematical
    !  computations' by forsythe,malcolm, and moler.

    !  input
    !  n.........the number of data points
    !  x.........the array of data points
    !  y.........the array of fcn. values
    !  x0........the upper limit of integration
    !  b,c,d.....work arrays


    !  output
    !  result....the approximation to the integral of y(x) from
    !           x = x(1)  to  x = x0.


    !  p.w.gaffney   24th april 1978.


    real(IDP) :: sum,t,dj,alpha,a1,a3,a4,a20,b1,di
    integer :: nm1,ib,i,iup,j

    !  check that the upper limit of integration is sensible

    if (x0 < x(1)) return

    !  compute the coefficients of the interpolating spline

    if (n > 2) then

       nm1 = n-1

       !  set up tridiagonal system

       !  b = diagonal, d = offdiagonal, c = right hand side.

       d(1) = x(2) - x(1)
       c(2) = (y(2) - y(1))/d(1)
       do i = 2, nm1
          d(i) = x(i+1) - x(i)
          b(i) = 2.*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i)
          c(i) = c(i+1) - c(i)
       end do

       !  end conditions.  third derivatives at  x(1)  and  x(n)
       !  obtained from divided differences

       b(1) = 1.
       d(1) = 0.
       b(n) = -d(nm1)
       c(1) = 0.
       c(n) = 0.
       if (n > 3) then
          c(n) = c(nm1)/(x(n)-x(n-2)) - c(n-2)/(x(nm1)-x(n-3))
          c(n) = -c(n)*d(nm1)**2/(x(n)-x(n-3))
       end if

       !  forward elimination

       do i = 2, n
          t = d(i-1)/b(i-1)
          b(i) = b(i) - t*d(i-1)
          c(i) = c(i) - t*c(i-1)
       end do

       !  back substitution

       c(n) = c(n)/b(n)
       do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       end do

       !  c(i) is now the sigma(i) of the text

       !  compute second derivative array c

       do i = 1,n
          c(i) = 3.*c(i)
       end do

       !  compute the integral of y(x) from x=x(1) to x=x0

       result = 0.0
       iup = nm1
       if (x0 < x(n)) then

          !  compute j such that x(j).le.x.lt.x(j+1)

          do j=1,nm1
             if (x0 < x(j+1)) exit
          end do

          !  first compute the integral from x=x(j) to x=x0

          iup = j-1
          dj = x(j+1)-x(j)
          alpha = (x0 - x(j))/dj
          a1 = 0.5*alpha*alpha
          a20=alpha-a1
          a3=a1*(a1-1.0)
          a4=alpha*alpha*(alpha-1.-0.5*a1)
          result=dj*(a1*y(j+1)+a20*y(j))+dj**3*(a3*c(j+1)+a4*c(j))/3.
       end if

       !  then compute the integral from x=x(1) to x=x(j)

       sum = 0.0
       do i=1,iup
          di = x(i+1)-x(i)
          sum = sum + 0.5*di*(y(i)+y(i+1)) - (c(i)+c(i+1))*di**3/12.
       end do
       result = sum + result

    else if (n == 2) then

       if (x0 < x(n)) then
          b1 = (y(2)-y(1))/(x(2)-x(1))
          result = (x0-x(1))*(y(1)+0.5*b1*(x0-x(1)))
       else
          result = 0.5*(x(2)-x(1))*(y(2)+y(1))
       end if

    end if

  end subroutine quadq

  ! subroutine quanc8(a,b,abserr,relerr,result,errest,nofun,flag)

  !   implicit none

  !   real(IDP) :: a, b, abserr, relerr, result, errest, flag
  !   integer :: nofun

  ! !  estimate the integral of fun(x) from a to b
  ! !  to a user provided tolerance.
  ! !  an automatic adaptive routine based on
  ! !  the 8-panel newton-cotes rule.

  ! !  input ..

  ! !  fun     the name of the integrand function subprogram fun(x).
  ! !  a       the lower limit of integration.
  ! !  b       the upper limit of integration. (b may be less than a)
  ! !  relerr  a relative error tolerance. (should be non-negative)
  ! !  abserr  an absolute error tolerance. (should be non-negative)

  ! !  output ..

  ! !  result  an approximation to the integral hopefully satisfying the
  ! !          least stringent of the two error tolerances.
  ! !  errest  an estimate of the magnitude of the actual error.
  ! !  nofun   the number of function values used in calculation of result.
  ! !  flag    a reliability indicator.  if flag is zero, then result
  ! !          probably satisfies the error tolerance.  if flag is
  ! !          xxx.yyy , then  xxx = the number of intervals which have
  ! !          not converged and  0.yyy = the fraction of the interval
  ! !          left to do when the limit on  nofun  was approached.

  !   integer :: levmin,levmax,levout,nomax,nofin,lev,nim,i,j
  !   real(IDP) :: w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp,qprev,qnow,qdiff,qleft,esterr,tolerr
  !   real(IDP), dimension(13) :: qright
  !   real(IDP), dimension(16) :: f,x
  !   real(IDP), dimension(8,30) :: fsave,xsave


  ! !  ***   stage 1 ***   general initialization
  ! !  set constants.

  !   levmin = 1
  !   levmax = 30
  !   levout = 6
  !   nomax = 5000
  !   nofin = nomax - 8*(levmax-levout+2**(levout+1))

  ! !  trouble when nofun reaches nofin

  !   w0 =   3956.0_IDP / 14175.0_IDP
  !   w1 =  23552.0_IDP / 14175.0_IDP
  !   w2 =  -3712.0_IDP / 14175.0_IDP
  !   w3 =  41984.0_IDP / 14175.0_IDP
  !   w4 = -18160.0_IDP / 14175.0_IDP

  ! !  initialize running sums to zero.

  !   flag = 0.0
  !   result = 0.0
  !   cor11  = 0.0
  !   errest = 0.0
  !   area   = 0.0
  !   nofun = 0
  !   if (a == b) return

  ! !  ***   stage 2 ***   initialization for first interval

  !   lev = 0
  !   nim = 1
  !   x0 = a
  !   x(16) = b
  !   qprev  = 0.0
  !   f0 = fun(x0)
  !   stone = (b - a) / 16.0
  !   x(8)  =  (x0  + x(16)) / 2.0
  !   x(4)  =  (x0  + x(8))  / 2.0
  !   x(12) =  (x(8)  + x(16)) / 2.0
  !   x(2)  =  (x0  + x(4))  / 2.0
  !   x(6)  =  (x(4)  + x(8))  / 2.0
  !   x(10) =  (x(8)  + x(12)) / 2.0
  !   x(14) =  (x(12) + x(16)) / 2.0
  !   do j = 2, 16, 2
  !      f(j) = fun(x(j))
  !   end do
  !   nofun = 9

  ! !  ***   stage 3 ***   central calculation
  ! !  requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16.
  ! !  calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area.

  !   do while (lev > 0 .or. nofun == 9)

  !      x(1) = (x0 + x(2)) / 2.0
  !      f(1) = fun(x(1))
  !      do j = 3, 15, 2
  !         x(j) = (x(j-1) + x(j+1)) / 2.0
  !         f(j) = fun(x(j))
  !      end do
  !      nofun = nofun + 8
  !      step = (x(16) - x0) / 16.0
  !      qleft  =  (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6)) + w3*(f(3)+f(5))  +  w4*f(4)) * step
  !      qright(lev+1)=(w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14)) + w3*(f(11)+f(13)) + w4*f(12)) * step
  !      qnow = qleft + qright(lev+1)
  !      qdiff = qnow - qprev
  !      area = area + qdiff

  ! !    ***   stage 4 *** interval convergence test

  !      esterr = abs(qdiff) / 1023.0
  !      tolerr = max(abserr,relerr*abs(area)) * (step/stone)

  !      if (lev < levmin .or. (lev < levmax .and. nofun <= nofin .and. esterr > tolerr)) then

  ! !       ***   stage 5   ***   no convergence
  ! !       locate next interval.

  !         nim = 2*nim
  !         lev = lev+1

  ! !       store right hand elements for future use.

  !         do i = 1, 8
  !            fsave(i,lev) = f(i+8)
  !            xsave(i,lev) = x(i+8)
  !         end do

  ! !       assemble left hand elements for immediate use.

  !         qprev = qleft
  !         do i = 1, 8
  !            j = -i
  !            f(2*j+18) = f(j+9)
  !            x(2*j+18) = x(j+9)
  !         end do

  !      else

  !         if (lev >= levmax) then

  ! !          current level is levmax.

  !            flag = flag + 1.0

  !         else if (nofun > nofin) then

  ! !          ***   stage 6   ***   trouble section
  ! !          number of function values is about to exceed limit.


  !            nofin = 2*nofin
  !            levmax = levout
  !            flag = flag + (b - x0) / (b - a)

  !         end if

  ! !       ***   stage 7   ***   interval converged
  ! !       add contributions into running sums.

  !         result = result + qnow
  !         errest = errest + esterr
  !         cor11  = cor11  + qdiff / 1023.0

  ! !       locate next interval.

  !         do while (nim /= 2*(nim/2))
  !            nim = nim/2
  !            lev = lev-1
  !         end do
  !         nim = nim + 1

  !         if (lev > 0) then

  ! !       assemble elements required for the next interval.

  !            qprev = qright(lev)
  !            x0 = x(16)
  !            f0 = f(16)
  !            do i = 1, 8
  !               f(2*i) = fsave(i,lev)
  !               x(2*i) = xsave(i,lev)
  !            end do

  !         end if

  !      end if

  !   end do

  ! !  ***   stage 8   ***   finalize and return

  !   result = result + cor11

  ! !  make sure errest not less than roundoff level.

  !   if (errest /= 0.0) then
  !      do
  !         temp = abs(result) + errest
  !         if (temp /= abs(result)) exit
  !         errest = 2.0*errest
  !      end do
  !   end if

  ! contains

  !   function fun(rr)
  !      use param
  !      use equil
  !      implicit none
  !      real(IDP) :: fun,rr
  !      fun=2.*rr*exp(-100.*(rr-rsbar)**2)/(q0*(1.+(rr/ro)**2)**2)
  !   end function fun

  ! end subroutine quanc8

  function erf(x)

    implicit none

    real(IDP) :: erf,x,sign_x,xh,one,four,pi,dn,du,ai,u,fac,arg,aih
    integer :: n,np,i

    !  this routine calculates the standard error fun..

    if (x == 0.0) then
       erf=0.0
    else
       sign_x=sign(1.0_IDP,x)
       xh=x
       if (abs(x) > 5.0_IDP) then
          erf=sign_x
       else
          if (x < 0.0) x=-x
          one=1.0
          four=4.0
          pi=atan(one)*four
          n=1000
          dn=n
          du=x/dn
          if (x > 3.0_IDP) du=3.0/dn
          np=n+1
          ai=0.0
          u=0.0

          do i=1,np
             fac=1.0
             if (i == 1 .or. i == np) fac=0.5
             arg=exp(-u*u)
             ai=ai+arg*fac
             u=u+du
          end do
          aih=ai
          erf=ai*du*2.0/sqrt(pi)*sign_x
          if (x > 3.0_IDP) then

             du=(x-3.0)/dn
             u=3.0
             ai=0.0
             do i=1,np
                fac=1.0
                if (i == 1 .or. i == np) fac=0.5
                arg=exp(-u*u)
                ai=ai+arg*fac
                u=u+du
             end do
             erf=ai*du*2.0/sqrt(pi)*sign_x+erf
          end if
          x=xh

       end if
    end if

  end function erf

  subroutine zzdisp(x,y,zzr,zzi)

    implicit none

    real(IDP) :: x,y,zzr,zzi,x1,y1,wzr1,wzi1,a,b,abr,abi

    x1=abs(x)
    y1=abs(y)
    call wzdisp(x1,y1,wzr1,wzi1)
    if(y.lt.0.0_IDP) then
       a=2.0_IDP*x1*y1
       b=-(x1*x1-y1*y1)
       abr=2.0_IDP*exp(b)*cos(a)
       abi=-2.0_IDP*exp(b)*sin(a)
       wzr1=abr-wzr1
       wzi1=abi-wzi1
    end if
    if(x*y.lt.0.0_IDP) wzi1=-wzi1
    zzr=-1.7724538509055_IDP*wzi1
    zzi= 1.7724538509055_IDP*wzr1

  end subroutine zzdisp

  subroutine wzdisp(x,y,re,im)

    implicit none

    real(IDP) :: im,lambda,epsh,epsl,epsy,x,y,re,h,h2,ss,rr,ri,sr,si,tr,ti,cc,c,s
    integer :: capn, nu, nup, i, n, np1
    logical :: b

    epsh=1.e-12_IDP
    epsl = epsh; epsy = epsh
    if(y.lt.4.29_IDP .and. x.lt.5.33_IDP) then

       !  (x,y) belongs to r

       s=(1.0_IDP-y/4.29_IDP)*sqrt(1.0_IDP-x*x/28.41_IDP)
       h=1.6_IDP*s
       h2=2.0_IDP*h
       capn = int(6.0_IDP + 23.0_IDP*s + 0.5_IDP)
       nu = int(9.0_IDP + 21.0_IDP*s + 0.5_IDP)
       lambda=h2**capn
    else

       !  (x,y) belongs to q1-r

       h=0.0_IDP
       capn=0
       nu=8
    end if
    b=(h.eq.0.0_IDP .or. lambda.lt.epsl)

    !  statement (lambda.lt.epsl) covers the underflow case
    !  when h(.gt.0) is very small.

    rr=0.0_IDP
    ri=0.0_IDP
    sr=0.0_IDP
    si=0.0_IDP
    nup=nu+1
    do i=1,nup
       n=nup-i
       np1=n+1
       tr=y+h+np1*rr
       ti=x-np1*ri
       c=0.5_IDP/(tr*tr+ti*ti)
       rr=c*tr
       ri=c*ti
       if(.not.(h .gt. 0.0_IDP .and. n .le. capn)) cycle
       tr=lambda+sr
       sr=rr*tr-ri*si
       si=ri*tr+rr*si
       lambda=lambda/h2
    end do
    cc=1.12837916709551_IDP
    if(b) then
       if(y.lt.epsy) then
          re=exp(-x*x)
       else
          re=rr*cc
       end if
       im=ri*cc
    else
       if(y.lt.epsy) then
          re=exp(-x*x)
       else
          re=sr*cc
       end if
       im=si*cc
    end if

  end subroutine wzdisp

  subroutine numinc

    !  numruno=numrun
    !  increment numrun

    use cotrol

    implicit none

    integer :: i
    character(len=2), dimension(36) :: charlist=(/"0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i", &
         "j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)

    do i=1,3
       numruno(i)=numrun(i)
    end do
    do i=1,36
       if (numrun(3) == charlist(i)) exit
    end do
    if (i < 36) then
       numrun(3)=charlist(i+1)
    else
       write(6,'(" error stop in numinc. numrun(3) got to ""z"".")')
       stop
    end if

  end subroutine numinc

  subroutine elapsed_time(values_s,values_e)

    implicit none

    integer, dimension(8) :: values_s,values_e
    integer :: days, hours, minutes, seconds
    integer, dimension(12) :: days_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    if (mod(values_e(1),4) == 0) days_month(2) = 29
    days = values_e(3) - values_s(3)
    if (values_e(1) > values_s(1)) then
       days = days + 31
    else if (values_e(2) > values_s(2)) then
       days = days + days_month(values_s(2))
    end if
    hours = 24*days + values_e(5) - values_s(5)
    minutes = values_e(6) - values_s(6)
    if (minutes < 0) then
       hours = hours - 1
       minutes = minutes + 60
    end if
    seconds = values_e(7) - values_s(7)
    if (seconds < 0) then
       minutes = minutes  - 1
       seconds = seconds + 60
    end if

    write(6,'(" elapsed time:",i4," hours",i3," minutes",i3," seconds")') hours,minutes,seconds

  end subroutine elapsed_time

END MODULE tools 
