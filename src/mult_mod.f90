MODULE mult_mod

  PRIVATE :: cpx2cs_par, cs2cpx_par, cs2cpxeq_par
  PUBLIC  :: mult, multed, chtype

CONTAINS

  subroutine mult (f, g, gtype, h, htype, c1, c2)

    !   form the convolution of the two functions g and h, and store it in f.
    !   The logic is based on the complex exponential form of the functions,
    !   for which the convolution is

    !        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

    !   where m', n' range over all values for which both (m-m',n-n') and (m',n')
    !   are in the range of the data.

    !   The input data, however is in a cos-sin rather than complex exponential
    !   format, so it is converted to the other format as needed. The basic
    !   algorithm is: for each row (constant n) of exponential form output,
    !   find each row pair (n-n' and n') that contributes to that output row,
    !   then build that row-pair from the input functions and perform all
    !   calculations in the complex domain. When all row pairs have been
    !   considered the output row is done, and the results are mapped back into
    !   the cos-sin form in the result.

    !   when the cos-sin function is type 1, cos(mx+ny) terms are stored at
    !   ll(m,n) in the function array, sin(mx+ny) are stored at ll(-m,-n).
    !   For type -1 functions this is reversed.

    !   Since the functions represented are real, the complex exponential
    !   has symmetry through the origin: [fr,fi](-m,-n) = [fr,-fi](m,n).
    !   So only positive values of n need be computed. (negative n rows of
    !   the input functions still contribute to the convolution, however.)

    !   This is a banded convolution, since the input functions are banded.
    !   The arrays mmstart and mmend give, for each row n in the complex
    !   representation, the lower and upper bounds of m values for which
    !   data is defined. Outside these bounds the data is assumed to be zero.


    use param
    use var_para
    use domain
    implicit none

    integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
    integer :: h1idx, g1idx
    integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: g,h,f
    real(IDP), dimension(mj_end-mj_start+1,-mmaxx:mmaxx) :: fr,fi
    real(IDP), dimension(mj_end-mj_start+1,0:mxmband,0:nmax) :: gr,gi,hr,hi

    !   arrange zero column of g and h so that references ll(m,n) that
    !   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
    !   the output array

    g(mj_start:mj_end,0) = 0.0
    h(mj_start:mj_end,0) = 0.0

    !   Calculate the convolutions for one row at a time. For given row frow,
    !   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
    !   hrow takes on all values for which both it and frow - hrow have defined
    !   data. Note the sum of the row pairs hrow + grow = frow.

    !   convert operands to complex exponential format

    gr = 0.0
    gi = 0.0
    hr = 0.0
    hi = 0.0
    call cs2cpx_par (g, gtype, gr, gi)
    call cs2cpx_par (h, htype, hr, hi)

    ftype = gtype*htype

    !$OMP PARALLEL DO PRIVATE(fr,fi,fstart,fend,grow,hrow,gstart,gend,hstart,hend,midx,mstart,mend,gidx,hidx,g1idx,h1idx)
    do frow = 0, nmax
       fstart = mmstart (frow)
       fend = mmend (frow)
       if (fend < fstart) cycle
       fr = 0.0
       fi = 0.0

       do grow = frow-nmax, nmax

          !   Loop through all row pairs whose sum is frow

          hrow = frow - grow

          !   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
          !   are [fstart-gend, fend-gstart] and [hstart,hend]

          gstart = mmstart (grow)
          gend = mmend (grow)
          hstart = mmstart (hrow)
          hend = mmend (hrow)

          if (gend < gstart .or. hend < hstart) cycle
          if (fend-gstart < hstart .or. fstart-gend > hend) cycle

          !   Compute all contributions to the current frow from the current grow-hrow
          !   pair

          do midx = fstart,fend

             !   calculate the loop limits for this particular midx. The ranges are
             !   [midx-gend,midx-gstart] and [hstart,hend]

             mstart = hstart
             if (midx-gend > hstart) mstart = midx-gend
             mend = hend
             if (midx-gstart < hend) mend = midx-gstart

             if (mend < mstart) cycle

             !   Do some convolving! Since only the positive bands of the complex
             !   arrays are stored, negative band values are extracted from corresponding
             !   positive bands. Also, the bands are not stored directly at their indices,
             !   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

             !   It is possible that the following loops may be speeded up somewhat by
             !   moving the address calculations into the loop limits: i.e. let hidx
             !   be the do loop index directly.

             if (grow >= 0 .and. hrow >= 0) then

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = h1idx - mmstart(hrow)
                   gidx = g1idx - mmstart(grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
                   fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
                end do

             else if (grow < 0 .and. hrow >= 0) then

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = h1idx - mmstart(hrow)
                   gidx = -g1idx - mmstart(-grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
                   fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
                end do

             else
                !                (grow >= 0 .and. hrow < 0)

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = -h1idx - mmstart(-hrow)
                   gidx = g1idx - mmstart(grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
                   fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
                end do

             end if

          end do

       end do

       !   The current row of f has now been completed. Store it back in its
       !   packed form and loop to the next result row

       call cpx2cs_par (f, ftype, frow, fr, fi, c1, c2)

    end do
    !$OMP END PARALLEL DO

  end subroutine mult

  subroutine cpx2cs_par (f, ftype, frow, fr, fi, c1, c2)

    use param
    use var_para
    use domain
    implicit none

    integer :: ftype,frow,fstart,fend,midx

    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: f
    real(IDP), dimension(:,-mmaxx:) :: fr,fi

    !   complex exponential to cos-sin conversion

    if (frow /= 0) then

       fstart = mmstart(frow)
       fend = mmend(frow)

       if (ftype == 1) then

          do midx = fstart,fend
             f(mj_start:mj_end,ll(midx,frow)) =   c1 * f(mj_start:mj_end,ll(midx,frow))   + c2 * 2.0 * fr(:,midx)
             f(mj_start:mj_end,ll(-midx,-frow)) = c1 * f(mj_start:mj_end,ll(-midx,-frow)) - c2 * 2.0 * fi(:,midx)
          end do

       else

          do midx = fstart,fend
             f(mj_start:mj_end,ll(midx,frow))   =  c1 * f(mj_start:mj_end,ll(midx,frow))   - c2 * 2.0 * fi(:,midx)
             f(mj_start:mj_end,ll(-midx,-frow)) =  c1 * f(mj_start:mj_end,ll(-midx,-frow)) + c2 * 2.0 * fr(:,midx)
          end do

       end if

    else

       f(mj_start:mj_end,ll(0,0)) = c1 * f(mj_start:mj_end,ll(0,0)) + c2 * fr(:,0)

       do midx = 1,mmend(0)
          if (ftype == 1) then
             f(mj_start:mj_end,ll(midx,0)) =  c1 * f(mj_start:mj_end,ll(midx,0))  + c2 * 2.0 * fr(:,midx)
             f(mj_start:mj_end,ll(-midx,0)) = c1 * f(mj_start:mj_end,ll(-midx,0)) - c2 * 2.0 * fi(:,midx)
          else
             f(mj_start:mj_end,ll(midx,0)) =  c1 * f(mj_start:mj_end,ll(midx,0))  - c2 * 2.0 * fi(:,midx)
             f(mj_start:mj_end,ll(-midx,0)) = c1 * f(mj_start:mj_end,ll(-midx,0)) + c2 * 2.0 * fr(:,midx)
          end if
       end do

    end if

  end subroutine cpx2cs_par

  subroutine cs2cpx_par (f, ftype, fr, fi)

    !   This routine converts the cos-sin representation of the function
    !   f into the equivalent complex exponential form

    !   if itype = 1 (assume n >= 0)

    !      [fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
    !      [fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

    !   if itype = -1, sines and cosines ares switched in f, so

    !        [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
    !        [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

    !   The (-m,-n) values are not stored, since they are just the complex
    !   conjugate of corresponding (m,n) values.

    !   To minimize storage, each band (constant n) is stored shifted so
    !   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
    !   the next element at fr(j,1,n) etc, and similarly for fi.

    use param
    use var_para
    use domain
    implicit none

    integer :: ftype,frow,mptr,midx,lr,li

    real(IDP), dimension(mj_start:,0:) :: f
    real(IDP), dimension(:,0:,0:) :: fr,fi

    !$OMP PARALLEL DO PRIVATE(mptr,midx,lr,li)
    do frow = 1,nmax

       mptr = 0
       do midx = mmstart(frow),mmend(frow)
          if (ftype == 1) then
             lr = ll(midx,frow)
             li = ll(-midx,-frow)
          else
             lr = ll(-midx,-frow)
             li = ll(midx,frow)
          end if

          fr(:,mptr,frow) =  0.5*f(mj_start:mj_end,lr)
          fi(:,mptr,frow) = -0.5*f(mj_start:mj_end,li)

          mptr = mptr + 1
       end do

    end do
    !$OMP END PARALLEL DO

    !   Do frow = 0 as a special case, since it contains its own
    !   conjugate image

    mptr = 0
    do midx = mmstart(0), mmend(0)
       if ((ftype == 1 .and. midx >= 0) .or. (ftype == -1 .and. midx < 0)) then
          lr = ll(midx,0)
          li = ll(-midx,0)
       else
          lr = ll(-midx,0)
          li = ll(midx,0)
       end if

       if (midx > 0) then
          fr(:,mptr,0) =  0.5*f(mj_start:mj_end,lr)
          fi(:,mptr,0) = -0.5*f(mj_start:mj_end,li)
       else if (midx < 0) then
          fr(:,mptr,0) = 0.5*f(mj_start:mj_end,lr)
          fi(:,mptr,0) = 0.5*f(mj_start:mj_end,li)
       else
          fr(:,mptr,0) = f(mj_start:mj_end,lr)
          fi(:,mptr,0) = 0.0
       end if
       mptr = mptr + 1
    end do

  end subroutine cs2cpx_par

  subroutine multed (f, g, gtype, h, htype, c1, c2)

    !   form the convolution of the two functions g and h, and store it in f.
    !   The logic is based on the complex exponential form of the functions,
    !   for which the convolution is

    !        [fr,fi](m,n) = SUM {[gr,gi](m-m',n-n')*[hr,hi](m',n')}

    !   where m', n' range over all values for which both (m-m',n-n') and (m',n')
    !   are in the range of the data.

    !   This is a banded convolution, since the input functions are banded.
    !   The arrays mmstart and mmend give, for each row n in the complex
    !   representation, the lower and upper bounds of m values for which
    !   data is defined. Outside these bounds the data is assumed to be zero.


    use param
    use var_para
    use domain
    implicit none

    integer :: gtype, htype, ftype, gstart, gend, hstart, hend, fstart, fend
    integer :: h1idx, g1idx
    integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx

    real(IDP) :: c1,c2
    real(IDP), dimension(mj_start:,0:) :: h,f
    real(IDP), dimension(0:,0:) :: g
    real(IDP), dimension(mj_end-mj_start+1,-mmaxx:mmaxx) :: fr,fi
    real(IDP), dimension(mj_end-mj_start+1,0:mxmband,0:nmax) :: gr,gi,hr,hi

    !   arrange zero column of g and h so that references ll(m,n) that
    !   return zero will produce g(i,ll(m,n)) = 0.0, etc., and initialize
    !   the output array

    g(:,0) = 0.0
    h(:,0) = 0.0

    !   Calculate the convolutions for one row at a time. For given row frow,
    !   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
    !   hrow takes on all values for which both it and frow - hrow have defined
    !   data. Note the sum of the row pairs hrow + grow = frow.

    !   convert operands to complex exponential format

    gr = 0.0
    gi = 0.0
    hr = 0.0
    hi = 0.0
    call cs2cpxeq_par (g, gtype, gr, gi)
    call cs2cpx_par (h, htype, hr, hi)

    ftype = gtype*htype

    !$OMP PARALLEL DO PRIVATE(fr,fi,fstart,fend,grow,hrow,gstart,gend,hstart,hend,midx,mstart,mend,gidx,hidx,g1idx,h1idx)
    do frow = 0, nmax
       fstart = mmstart (frow)
       fend = mmend (frow)
       if (fend < fstart) cycle
       fr = 0.0
       fi = 0.0

       do grow = -nmaxeq, nmaxeq, nfp

          !   Loop through all row pairs whose sum is frow

          hrow = frow - grow
          if (hrow < -nmax .or. hrow > nmax) cycle

          !   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
          !   are [fstart-gend, fend-gstart] and [hstart,hend]

          gstart = mmstreq (grow)
          gend = mmendeq (grow)
          hstart = mmstart (hrow)
          hend = mmend (hrow)

          if (gend < gstart .or. hend < hstart) cycle
          if (fend-gstart < hstart .or. fstart-gend > hend) cycle

          !   Compute all contributions to the current frow from the current grow-hrow
          !   pair

          do midx = fstart,fend

             !   calculate the loop limits for this particular midx. The ranges are
             !   [midx-gend,midx-gstart] and [hstart,hend]

             mstart = hstart
             if (midx-gend > hstart) mstart = midx-gend
             mend = hend
             if (midx-gstart < hend) mend = midx-gstart

             if (mend < mstart) cycle

             !   Do some convolving! Since only the positive bands of the complex
             !   arrays are stored, negative band values are extracted from corresponding
             !   positive bands. Also, the bands are not stored directly at their indices,
             !   but are packed so that mmstart(n) is in fr(j,0,frow), etc.

             !   It is possible that the following loops may be speeded up somewhat by
             !   moving the address calculations into the loop limits: i.e. let hidx
             !   be the do loop index directly.

             if (grow >= 0 .and. hrow >= 0) then

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = h1idx - mmstart(hrow)
                   gidx = g1idx - mmstreq(grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
                   fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
                end do

             else if (grow < 0 .and. hrow >= 0) then

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = h1idx - mmstart(hrow)
                   gidx = -g1idx - mmstreq(-grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
                   fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
                end do

             else
                !                (grow >= 0 .and. hrow < 0)

                do h1idx = mstart,mend
                   g1idx = midx - h1idx
                   hidx = -h1idx - mmstart(-hrow)
                   gidx = g1idx - mmstreq(grow)
                   fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
                   fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
                end do

             end if

          end do

       end do

       !   The current row of f has now been completed. Store it back in its
       !   packed form and loop to the next result row

       call cpx2cs_par (f, ftype, frow, fr, fi, c1, c2)

    end do
    !$OMP END PARALLEL DO

  end subroutine multed

  subroutine cs2cpxeq_par (f, ftype, fr, fi)

    !   This routine converts the cos-sin representation of the function
    !   f into the equivalent complex exponential form

    !   The (-m,-n) values are not stored, since they are just the complex
    !   conjugate of corresponding (m,n) values.

    !   To minimize storage, each band (constant n) is stored shifted so
    !   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
    !   the next element at fr(j,1,n) etc, and similarly for fi.

    use param
    use var_para
    use domain
    implicit none

    integer :: ftype,frow,mptr,midx,lr,li,l

    real(IDP), dimension(0:,0:) :: f
    real(IDP), dimension(:,0:,0:) :: fr,fi

    if (lasym) then

       do frow = nfp, nmaxeq, nfp

          mptr = 0
          do midx = mmstreq(frow),mmendeq(frow)
             if (ftype == 1) then
                lr = lleq(midx,frow)
                li = lleq(-midx,-frow)
             else
                lr = lleq(-midx,-frow)
                li = lleq(midx,frow)
             end if
             fr(:,mptr,frow) =  0.5*f(mj_start:mj_end,lr)
             fi(:,mptr,frow) = -0.5*f(mj_start:mj_end,li)
             mptr = mptr + 1
          end do

       end do

       !   Do frow = 0 as a special case, since it contains its own
       !   conjugate image

       mptr = 0
       do midx = mmstreq(0),mmendeq(0)
          if ((ftype == 1 .and. midx >= 0) .or. (ftype == -1 .and. midx < 0)) then
             lr = lleq(midx,0)
             li = lleq(-midx,0)
          else
             lr = lleq(-midx,0)
             li = lleq(midx,0)
          end if

          if (midx > 0) then
             fr(:,mptr,0) =  0.5*f(mj_start:mj_end,lr)
             fi(:,mptr,0) = -0.5*f(mj_start:mj_end,li)
          else if (midx < 0) then
             fr(:,mptr,0) = 0.5*f(mj_start:mj_end,lr)
             fi(:,mptr,0) = 0.5*f(mj_start:mj_end,li)
          else
             fr(:,mptr,0) = f(mj_start:mj_end,lr)
             fi(:,mptr,0) = 0.0
          end if
          mptr = mptr + 1
       end do

    else

       do frow = nfp, nmaxeq, nfp

          mptr = 0
          do midx = mmstreq(frow),mmendeq(frow)
             l = lleq(midx,frow)
             if (ftype == 1) then
                fr(:,mptr,frow) =  0.5*f(mj_start:mj_end,l)
                fi(:,mptr,frow) =  0.0
             else
                fr(:,mptr,frow) =  0.0
                fi(:,mptr,frow) = -0.5*f(mj_start:mj_end,l)
             end if
             mptr = mptr + 1
          end do

       end do

       !   Do frow = 0 as a special case, since it contains its own
       !   conjugate image

       mptr = 0
       do midx = mmstreq(0), mmendeq(0)
          if (midx >= 0) then
             l = lleq(midx,0)
          else
             l = lleq(-midx,0)
          endif

          if (ftype == 1) then
             if (midx > 0) then
                fr(:,mptr,0) =  0.5*f(mj_start:mj_end,l)
                fi(:,mptr,0) =  0.0
             else if (midx < 0) then
                fr(:,mptr,0) = 0.5*f(mj_start:mj_end,l)
                fi(:,mptr,0) = 0.0
             else
                fr(:,mptr,0) = f(mj_start:mj_end,l)
                fi(:,mptr,0) = 0.0
             end if
          else
             if (midx > 0) then
                fr(:,mptr,0) =  0.0
                fi(:,mptr,0) = -0.5*f(mj_start:mj_end,l)
             else if (midx < 0) then
                fr(:,mptr,0) = 0.0
                fi(:,mptr,0) = 0.5*f(mj_start:mj_end,l)
             else
                fr(:,mptr,0) = f(mj_start:mj_end,l)
                fi(:,mptr,0) = 0.0
             end if
          end if
          mptr = mptr + 1
       end do

    end if

  end subroutine cs2cpxeq_par

  subroutine chtype(a,b)

    use param
    use var_para
    use domain
    implicit none

    integer :: l
    real(IDP), dimension(mj_start:,0:) :: a,b

    a=0.0_IDP
    !$OMP PARALLEL DO
    do l=1,lmaxn
       a(:,lln(lo(l)))=b(:,lln(l))
    end do
    !$OMP END PARALLEL DO
    a(:,0)=0.0_IDP

  end subroutine chtype

END MODULE mult_mod
