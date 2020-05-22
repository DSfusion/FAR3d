      SUBROUTINE root(t,ft,b,c,relerr,abserr,iflag)

      USE stel_kinds
      IMPLICIT NONE

!     root computes a root of the nonlinear equation f(x)=0
!     where f(x) is a continuous real function of a single real
!     variable x.  the method used is a combination of bisection
!     and the secant rule.

!     normal input consists of a continuous function f and an
!     interval (b,c) such that f(b)*f(c).le.0.0.  each iteration
!     finds new values of b and c such that the interval (b,c) is
!     shrunk and f(b)*f(c).le.0.0.  the stopping criterion is

!      abs(b-c).le.2.0*(relerr*abs(b)+abserr)

!     where relerr=relative error and abserr=absolute error are
!     input quantities.  set the flag, iflag, positive to initialize
!     the computation.  as b,c and iflag are used for both input and
!     output, they must be variables in the calling program.

!     if 0 is a possible root, one should not choose abserr=0.0.

!     the output value of b is the better approximation to a root
!     as b and c are always redefined so that abs(f(b)).le.abs(f(c)).

!     to solve the equation, root must evaluate f(x) repeatedly. this
!     is done in the calling program.  when an evaluation of f is
!     needed at t, root returns with iflag negative.  evaluate ft=f(t)
!     and call root again.  do not alter iflag.

!     when the computation is complete, root returns to the calling
!     program with iflag positive:

!      iflag=1  if f(b)*f(c).lt.0 and the stopping criterion is met.

!           =2  if a value b is found such that the computed value
!               f(b) is exactly zero.  the interval (b,c) may not
!               satisfy the stopping criterion.

!           =3  if abs(f(b)) exceeds the input values abs(f(b)),
!               abs(f(c)).   in this case it is likely that b is close
!               to a pole of f.

!           =4  if no odd order root was found in the interval.  a
!               local minimum may have been obtained.

!           =5  if too many function evaluations were made.
!               (as programmed, 500 are allowed.)

!     this code is a modification of the code  zeroin  which is completely
!     explained and documented in the text,  numerical computing:  an
!     introduction  by l. f. shampine and r. c. allen.

      integer :: iflag
      integer, save :: ic,kount
      real(rprec) :: t,ft,b,c,relerr,abserr,fb,cmb,acmb,tol,p,q
      real(rprec), save :: re,ae,a,fa,fc,fx,u,acbs
      integer, save :: mentry=1


      if (mentry == 1) then

!     if first entry then compute u

       u=epsilon(t)
       mentry = 2
      end if

      if (iflag < -1) then
       fb=ft
       if (iflag == -2) then
        fc=fa
        kount=2
        fx=max(abs(fb),abs(fc))
       else if (iflag == -3 .and. fb /= 0.0_dp) then
        kount=kount+1
        if (sign(1.0_dp,fb) == sign(1.0_dp,fc)) then
         c=a
         fc=fa
        end if
       else
        iflag=2
        return
       end if
       if (abs(fc) < abs(fb)) then

!     interchange b and c so that abs(f(b)).le.abs(f(c)).

        a=b
        fa=fb
        b=c
        fb=fc
        c=a
        fc=fa
       end if
       cmb=0.5*(c-b)
       acmb=abs(cmb)
       tol=re*abs(b)+ae

!     test stopping criterion and function count.

       if (acmb <= tol) then
        if (sign(1.0_dp,fb) == sign(1.0_dp,fc)) then
         iflag=4
        else if (abs(fb) > fx) then
         iflag=3
        else
         iflag=1
        end if
        return
       end if
       if (kount >= 500) then
        iflag=5
        return
       end if

!     calculate new iterate implicitly as b+p/q
!     where we arrange p.ge.0.  the implicit
!     form is used to prevent overflow.

       p=(b-a)*fb
       q=fa-fb
       if (p < 0.0_dp) then
        p=-p
        q=-q
       end if

!     update a, check if reduction in the size of bracketing
!     interval is satisfactory.  if not, bisect until it is.

       a=b
       fa=fb
       ic=ic+1
       if (ic > 3 .and. 8.0*acmb < acbs) then
        ic=0
        acbs=acmb
       end if
       if (ic < 4) then

!     test for too small a change.

        if (p > abs(q)*tol) then

!     root ought to be between b and (c+b)/2.

         if (p < cmb*q) then

!     use secant rule.

          b=b+p/q
         else

!     use bisection.

          b=0.5*(c+b)
         end if

        else

!     increment by tolerance.

         b=b+sign(tol,cmb)
        end if
       else

!     use bisection.

        b=0.5*(c+b)
       end if

!     have completed computation for new iterate b.

       t=b
       iflag=-3
      else if (iflag == -1) then
       fa=ft
       t=b
       iflag=-2
      else
       re=max(relerr,u)
       ae=max(abserr,0.0)
       ic=0
       acbs=abs(b-c)
       a=c
       t=a
       iflag=-1
      end if

      END SUBROUTINE root
