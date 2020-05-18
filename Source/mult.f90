	subroutine mult (f, g, gtype, h, htype, c1, c2)

		use param
		use domain
		implicit none

		integer :: gtype, htype,index

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f

		interface
			subroutine multc11 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc11
			subroutine multc12 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc12
			subroutine multc22 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multc22
			subroutine multr11 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr11
			subroutine multr12 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr12
			subroutine multr22 (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multr22
		end interface

		index=abs(gtype*htype)
		if (nmin < 0) then
			select case (index)
				case (1)
					call multc11 (f, g, h, gtype, htype, c1, c2)
				case (2)
					if (abs(gtype) == 2) then
						call multc12 (f, g, h, gtype, htype, c1, c2)
					else
						call multc12 (f, h, g, htype, gtype, c1, c2)
					end if
				case (4)
					call multc22 (f, g, h, gtype, htype, c1, c2)
			end select
		else
			select case (index)
				case (1)
					call multr11 (f, g, h, gtype, htype, c1, c2)
				case (2)
					if (abs(gtype) == 2) then
						call multr12 (f, g, h, gtype, htype, c1, c2)
					else
						call multr12 (f, h, g, htype, gtype, c1, c2)
					end if
				case (4)
					call multr22 (f, g, h, gtype, htype, c1, c2)
			end select
		end if

	end subroutine mult

	subroutine multeq(feq,geq,itypeg,heq,itypeh,c1,c2,f,g,h)

		use param
		use domain
		implicit none

		interface
			subroutine eqtodyn(adyn,aeq,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: adyn,aeq
			end subroutine eqtodyn
			subroutine dyntoeq(aeq,adyn,c1,c2)
				use param
				implicit none
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: aeq,adyn
			end subroutine dyntoeq
			subroutine mult(f,g,itypeg,h,itypeh,c1,c2)
				use param
				implicit none
				integer :: itypeg,itypeh
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f,g,h
			end subroutine mult
		end interface

		integer :: itypeg,itypeh
		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: feq,geq,heq
		real(IDP), dimension(0:,0:) :: f,g,h

		call eqtodyn(g,geq,0.0_IDP,1.0_IDP)
		call eqtodyn(h,heq,0.0_IDP,1.0_IDP)
		call mult(f,g,itypeg,h,itypeh,0.0_IDP,1.0_IDP)
		call dyntoeq(feq,f,c1,c2)

	end subroutine multeq

	subroutine multb (f, g, gtype, h, htype, c1, c2)

		use param
		use domain
		implicit none

		integer :: gtype, htype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f

		interface
			subroutine multcb (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multcb
			subroutine multrb (f, g, h, gtype, htype, c1, c2)
				use param
				implicit none
				integer :: gtype, htype
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: g,h,f
			end subroutine multrb
		end interface

		if (nminb < 0) then
			call multcb (f, g, h, gtype, htype, c1, c2)
		else
			call multrb (f, g, h, gtype, htype, c1, c2)
		end if

	end subroutine multb

	subroutine multcb (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxb:mmaxb) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmbandb,0:nmaxb) :: gr,gi,hr,hi

		interface
			subroutine cpx2csb (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxb:) :: fr,fi
			end subroutine cpx2csb
			subroutine cs2cpxb (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpxb
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpxb (g, gtype, gr, gi)
		call cs2cpxb (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmaxb
			fstart = mmstartb (frow)
			fend = mmendb (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = frow-nmaxb, nmaxb

				!   Loop through all row pairs whose sum is frow

				hrow = frow - grow

				!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
				!   are [fstart-gend, fend-gstart] and [hstart,hend]

				gstart = mmstartb (grow)
				gend = mmendb (grow)
				hstart = mmstartb (hrow)
				hend = mmendb (hrow)

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
							hidx = h1idx - mmstartb(hrow)
							gidx = g1idx - mmstartb(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow) + gi(:,gidx,grow)*hr(:,hidx,hrow)
						end do

					else if (grow < 0 .and. hrow >= 0) then

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = h1idx - mmstartb(hrow)
							gidx = -g1idx - mmstartb(-grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
						end do

					else
!  						(grow >= 0 .and. hrow < 0)

						do h1idx = mstart,mend
							g1idx = midx - h1idx
							hidx = -h1idx - mmstartb(-hrow)
							gidx = g1idx - mmstartb(grow)
							fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
						end do

					end if

				end do

			end do

			!   The current row of f has now been completed. Store it back in its
			!   packed form and loop to the next result row

			call cpx2csb (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multcb

	subroutine multrb (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxxb:mmaxxb) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmbandb,0:nmaxb) :: gr,gi,hr,hi

		interface
			subroutine rli2csb (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxxb:) :: fr,fi
			end subroutine rli2csb
			subroutine cs2rlib (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rlib
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rlib (g, gtype, gr, gi)

		call cs2rlib (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmaxb
				fstart = mmstartb (frow)
				fend = mmendb (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxb, nmaxb

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstartb (grow)
					gend = mmendb (grow)
					hstart = mmstartb (hrow)
					hend = mmendb (hrow)

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
								hidx = h1idx - mmstartb(hrow)
								gidx = g1idx - mmstartb(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstartb(hrow)
								gidx = -g1idx - mmstartb(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstartb(-hrow)
								gidx = g1idx - mmstartb(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2csb (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmaxb
				fstart = mmstartb (frow)
				fend = mmendb (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxb, nmaxb

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstartb (grow)
					gend = mmendb (grow)
					hstart = mmstartb (hrow)
					hend = mmendb (hrow)

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
								hidx = h1idx - mmstartb(hrow)
								gidx = g1idx - mmstartb(grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstartb(hrow)
								gidx = -g1idx - mmstartb(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstartb(-hrow)
								gidx = g1idx - mmstartb(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstartb(0)
					gidx = -mmstartb(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2csb (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmaxb
				fstart = mmstartb (frow)
				fend = mmendb (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxb, nmaxb

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstartb (grow)
					gend = mmendb (grow)
					hstart = mmstartb (hrow)
					hend = mmendb (hrow)

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
								hidx = h1idx - mmstartb(hrow)
								gidx = g1idx - mmstartb(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstartb(hrow)
								gidx = -g1idx - mmstartb(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstartb(-hrow)
								gidx = g1idx - mmstartb(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstartb(0)
					gidx = -mmstartb(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2csb (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmaxb
				fstart = mmstartb (frow)
				fend = mmendb (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxb, nmaxb

	!   Loop through all row pairs whose sum is frow

					hrow = frow - grow

	!   See if m - [grow] intersects [hrow] for any m in [frow]. The intervals
	!   are [fstart-gend, fend-gstart] and [hstart,hend]

					gstart = mmstartb (grow)
					gend = mmendb (grow)
					hstart = mmstartb (hrow)
					hend = mmendb (hrow)

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
								hidx = h1idx - mmstartb(hrow)
								gidx = g1idx - mmstartb(grow)
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstartb(hrow)
								gidx = -g1idx - mmstartb(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstartb(-hrow)
								gidx = g1idx - mmstartb(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstartb(0)
					gidx = -mmstartb(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2csb (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multrb

	subroutine cpx2csb (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmaxb:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstartb(frow)
			fend = mmendb(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,llb(midx,frow)) =   c1 * f(:,llb(midx,frow))   + c2 * 2.0 * fr(:,midx)
					f(:,llb(-midx,-frow)) = c1 * f(:,llb(-midx,-frow)) - c2 * 2.0 * fi(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,llb(midx,frow))   =  c1 * f(:,llb(midx,frow))   - c2 * 2.0 * fi(:,midx)
					f(:,llb(-midx,-frow)) =  c1 * f(:,llb(-midx,-frow)) + c2 * 2.0 * fr(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,llb(0,0)) = c1 * f(:,llb(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmendb(0)
				if (ftype > 0) then
					f(:,llb(midx,0)) =  c1 * f(:,llb(midx,0))  + c2 * 2.0 * fr(:,midx)
					f(:,llb(-midx,0)) = c1 * f(:,llb(-midx,0)) - c2 * 2.0 * fi(:,midx)
				else
					f(:,llb(midx,0)) =  c1 * f(:,llb(midx,0))  - c2 * 2.0 * fi(:,midx)
					f(:,llb(-midx,0)) = c1 * f(:,llb(-midx,0)) + c2 * 2.0 * fr(:,midx)
				end if
			end do

		end if

	end subroutine cpx2csb

	subroutine cs2cpxb (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!      	[fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!      	[fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!        [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!        [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,lr,li

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

		!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		do frow = 1,nmaxb

			mptr = 0
			do midx = mmstartb(frow),mmendb(frow)
				if (ftype == 1) then
					lr = llb(midx,frow)
					li = llb(-midx,-frow)
				else
					lr = llb(-midx,-frow)
					li = llb(midx,frow)
				end if

				fr(:,mptr,frow) =  0.5*f(:,lr)
				fi(:,mptr,frow) = -0.5*f(:,li)

				mptr = mptr + 1
			end do

		end do

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstartb(0), mmendb(0)
			if ((ftype > 0 .and. midx >= 0) .or. (ftype < 0 .and. midx < 0)) then
				lr = llb(midx,0)
				li = llb(-midx,0)
			else
				lr = llb(-midx,0)
				li = llb(midx,0)
			end if

			if (midx > 0) then
				fr(:,mptr,0) =  0.5*f(:,lr)
				fi(:,mptr,0) = -0.5*f(:,li)
			else if (midx < 0) then
				fr(:,mptr,0) = 0.5*f(:,lr)
				fi(:,mptr,0) = 0.5*f(:,li)
			else
				fr(:,mptr,0) = facfl0 * f(:,lr)
				fi(:,mptr,0) = 0.0
		  	end if
			mptr = mptr + 1
		end do

	end subroutine cs2cpxb

	subroutine rli2csb (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmaxxb:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstartb(frow)
			fend = mmendb(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,llb(midx,frow)) =   c1 * f(:,llb(midx,frow)) + c2 * 2.0 * fr(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,llb(midx,frow))  =  c1 * f(:,llb(midx,frow)) - c2 * 2.0 * fi(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,llb(0,0)) = c1 * f(:,llb(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmendb(0)
				if (ftype > 0) then
					f(:,llb(midx,0)) =  c1 * f(:,llb(midx,0)) + c2 * 2.0 * fr(:,midx)
				else
					f(:,llb(midx,0)) =  c1 * f(:,llb(midx,0)) - c2 * 2.0 * fi(:,midx)
				end if
			end do

		end if

	end subroutine rli2csb

	subroutine cs2rlib (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!         [fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!         [fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!         [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!         [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,l

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

	!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		do frow = 1,nmaxb

			mptr = 0
			do midx = mmstartb(frow),mmendb(frow)
				l = llb(midx,frow)

				if (ftype == 1) then
					fr(:,mptr,frow) =  0.5*f(:,l)
					fi(:,mptr,frow) =  0.0
				else
					fr(:,mptr,frow) =  0.0
					fi(:,mptr,frow) = -0.5*f(:,l)
				end if

				mptr = mptr + 1
			end do

		end do

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstartb(0), mmendb(0)
			if (midx >= 0) then
				l = llb(midx,0)
			else
				l = llb(-midx,0)
			endif

			if (ftype > 0) then
				if (midx > 0) then
					fr(:,mptr,0) =  0.5*f(:,l)
					fi(:,mptr,0) =  0.0
				else if (midx < 0) then
					fr(:,mptr,0) = 0.5*f(:,l)
					fi(:,mptr,0) = 0.0
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			else
				if (midx > 0) then
					fr(:,mptr,0) =  0.0
					fi(:,mptr,0) = -0.5*f(:,l)
				else if (midx < 0) then
					fr(:,mptr,0) = 0.0
					fi(:,mptr,0) = 0.5*f(:,l)
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			end if
			mptr = mptr + 1
		end do

	end subroutine cs2rlib

	subroutine multc11 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

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
!  						(grow >= 0 .and. hrow < 0)

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

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc11

	subroutine multc12 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmax
			fstart = mmstart (frow)
			fend = mmend (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

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
!  						(grow >= 0 .and. hrow < 0)

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

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc12

	subroutine multc22 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmax:mmax) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmax:) :: fr,fi
			end subroutine cpx2cs
			subroutine cs2cpx (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2cpx
		end interface

		!   arrange zero column of g and h so that references ll(m,n) that
		!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
		!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

		!   Calculate the convolutions for one row at a time. For given row frow,
		!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
		!   hrow takes on all values for which both it and frow - hrow have defined
		!   data. Note the sum of the row pairs hrow + grow = frow.

		!   convert operands to complex exponential format

		call cs2cpx (g, gtype, gr, gi)
		call cs2cpx (h, htype, hr, hi)

		ftype = gtype*htype

		do frow = 0, nmaxeq, nfp
			fstart = mmstart (frow)
			fend = mmend (frow)
			if (fend < fstart) cycle
			fr = 0.0
			fi = 0.0

			do grow = frow-nmaxeq, nmaxeq, nfp

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
!  						(grow >= 0 .and. hrow < 0)

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

			call cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		end do

	end subroutine multc22

	subroutine cpx2cs (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmax:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstart(frow)
			fend = mmend(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,ll(midx,frow)) =   c1 * f(:,ll(midx,frow))   + c2 * 2.0 * fr(:,midx)
					f(:,ll(-midx,-frow)) = c1 * f(:,ll(-midx,-frow)) - c2 * 2.0 * fi(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,ll(midx,frow))   =  c1 * f(:,ll(midx,frow))   - c2 * 2.0 * fi(:,midx)
					f(:,ll(-midx,-frow)) =  c1 * f(:,ll(-midx,-frow)) + c2 * 2.0 * fr(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,ll(0,0)) = c1 * f(:,ll(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmend(0)
				if (ftype > 0) then
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0))  + c2 * 2.0 * fr(:,midx)
					f(:,ll(-midx,0)) = c1 * f(:,ll(-midx,0)) - c2 * 2.0 * fi(:,midx)
				else
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0))  - c2 * 2.0 * fi(:,midx)
					f(:,ll(-midx,0)) = c1 * f(:,ll(-midx,0)) + c2 * 2.0 * fr(:,midx)
				end if
			end do

		end if

	end subroutine cpx2cs

	subroutine cs2cpx (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!      	[fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!      	[fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!        [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!        [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,lr,li

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

		!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		if (abs(ftype) == 1) then
		
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

					fr(:,mptr,frow) =  0.5*f(:,lr)
					fi(:,mptr,frow) = -0.5*f(:,li)

					mptr = mptr + 1
				end do

			end do
			
		end if

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstart(0), mmend(0)
			if ((ftype > 0 .and. midx >= 0) .or. (ftype < 0 .and. midx < 0)) then
				lr = ll(midx,0)
				li = ll(-midx,0)
			else
				lr = ll(-midx,0)
				li = ll(midx,0)
			end if

			if (midx > 0) then
				fr(:,mptr,0) =  0.5*f(:,lr)
				fi(:,mptr,0) = -0.5*f(:,li)
			else if (midx < 0) then
				fr(:,mptr,0) = 0.5*f(:,lr)
				fi(:,mptr,0) = 0.5*f(:,li)
			else
				fr(:,mptr,0) = facfl0 * f(:,lr)
				fi(:,mptr,0) = 0.0
		  	end if
			mptr = mptr + 1
		end do

	end subroutine cs2cpx

	subroutine multr11 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart, fend
		integer :: h1idx, g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

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
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

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
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

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
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

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
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr11

	subroutine multr12 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart,fend
		integer :: h1idx,g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

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
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

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
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

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
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmax
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = ((frow-nmax)/nfp)*nfp, nmax, nfp

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
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr12

	subroutine multr22 (f, g, h, gtype, htype, c1, c2)

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
		use domain
		implicit none

		integer :: gtype, htype, gstart, gend, hstart, hend, fstart,fend
		integer :: h1idx,g1idx
		integer :: grow, hrow, frow, midx, mstart, mend, hidx, gidx, ftype

		real(IDP) :: c1,c2
		real(IDP), dimension(0:,0:) :: g,h,f
		real(IDP), dimension(0:mj,-mmaxx:mmaxx) :: fr,fi
		real(IDP), dimension(0:mj,0:mxmband,0:nmax) :: gr,gi,hr,hi

		interface
			subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)
				use param
				use domain
				implicit none
				integer :: ftype,frow
				real(IDP) :: c1,c2
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,-mmaxx:) :: fr,fi
			end subroutine rli2cs
			subroutine cs2rli (f, ftype, fr, fi)
				use param
				implicit none
				integer :: ftype
				real(IDP), dimension(0:,0:) :: f
				real(IDP), dimension(0:,0:,0:) :: fr,fi
			end subroutine cs2rli
		end interface

	!   arrange zero column of g and h so that references ll(m,n) that
	!   return zero will produce g(i,ll(m.n)) = 0.0, etc., and initialize
	!   the output array

		g(:,0) = 0.0
		h(:,0) = 0.0

	!   Calculate the convolutions for one row at a time. For given row frow,
	!   row pairs frow - hrow ( = grow) and hrow contribute to the convolution.
	!   hrow takes on all values for which both it and frow - hrow have defined
	!   data. Note the sum of the row pairs hrow + grow = frow.

	!   convert operands to complex exponential format

		call cs2rli (g, gtype, gr, gi)

		call cs2rli (h, htype, hr, hi)

		ftype = gtype*htype

		if (gtype > 0 .and. htype > 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

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
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  						(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gr(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype > 0 .and. htype < 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

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
								fi(:,midx) = fi(:,midx) + gr(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) + gr(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) - gr(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else if (gtype < 0 .and. htype > 0) then

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

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
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fi(:,midx) = fi(:,midx) - gi(:,gidx,-grow)*hr(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fi(:,midx) = fi(:,midx) + gi(:,gidx,grow)*hr(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		else
!  			(gtype < 0 .and. htype < 0)

			do frow = 0, nmaxeq, nfp
				fstart = mmstart (frow)
				fend = mmend (frow)
				if (fend < fstart) cycle
				fr = 0.0
				fi = 0.0

				do grow = frow-nmaxeq, nmaxeq, nfp

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
								fr(:,midx) = fr(:,midx) - gi(:,gidx,grow)*hi(:,hidx,hrow)
							end do

						else if (grow < 0 .and. hrow >= 0) then

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = h1idx - mmstart(hrow)
								gidx = -g1idx - mmstart(-grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,-grow)*hi(:,hidx,hrow)
							end do

						else
!  							(grow >= 0 .and. hrow < 0)

							do h1idx = mstart,mend
								g1idx = midx - h1idx
								hidx = -h1idx - mmstart(-hrow)
								gidx = g1idx - mmstart(grow)
								fr(:,midx) = fr(:,midx) + gi(:,gidx,grow)*hi(:,hidx,-hrow)
							end do

						end if

					end do

				end do

				if (frow == 0) then
					hidx = -mmstart(0)
					gidx = -mmstart(0)
					fr(:,0) = fr(:,0) + gr(:,gidx,0)*hr(:,hidx,0)
				end if

	!   The current row of f has now been completed. Store it back in its
	!   packed form and loop to the next result row

				call rli2cs (f, ftype, frow, fr, fi, c1, c2)

			end do

		end if

	end subroutine multr22

	subroutine rli2cs (f, ftype, frow, fr, fi, c1, c2)

		use param
		use domain
		implicit none

		integer :: ftype,frow,fstart,fend,midx

		real(IDP) :: c1,c2,facrl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,-mmaxx:) :: fr,fi

	!   complex exponential to cos-sin conversion

	!   The following factor should be set to 1.0 for dtem, 2.0 for kite

		facrl0 = 1.0

		if (frow /= 0) then

			fstart = mmstart(frow)
			fend = mmend(frow)

			if (ftype > 0) then

				do midx = fstart,fend
					f(:,ll(midx,frow)) =   c1 * f(:,ll(midx,frow)) + c2 * 2.0 * fr(:,midx)
				end do

			else

				do midx = fstart,fend
					f(:,ll(midx,frow))  =  c1 * f(:,ll(midx,frow)) - c2 * 2.0 * fi(:,midx)
				end do

			end if

		else

	!   special-case n=0

			f(:,ll(0,0)) = c1 * f(:,ll(0,0)) + c2 * facrl0 * fr(:,0)

			do midx = 1,mmend(0)
				if (ftype > 0) then
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0)) + c2 * 2.0 * fr(:,midx)
				else
					f(:,ll(midx,0)) =  c1 * f(:,ll(midx,0)) - c2 * 2.0 * fi(:,midx)
				end if
			end do

		end if

	end subroutine rli2cs

	subroutine cs2rli (f, ftype, fr, fi)

	!   This routine converts the cos-sin representation of the function
	!   f into the equivalent complex exponential form

	!   if itype = 1 (assume n >= 0)

	!         [fr,fi](m,n)   = 0.5*[f(m,n),-f(-m,-n)]
	!         [fr,fi](-m,-n) = 0.5*[f(m,n), f(-m,-n)]

	!   if itype = -1, sines and cosines ares switched in f, so

	!         [fr,fi](m,n)   = 0.5*[f(-m,-n),-f(m,n)]
	!         [fr,fi](-m,-n) = 0.5*[f(-m,-n), f(m,n)]

	!   The (-m,-n) values are not stored, since they are just the complex
	!   conjugate of corresponding (m,n) values.

	!   To minimize storage, each band (constant n) is stored shifted so
	!   that its first nonzero element (mmstart(n)) is stored at fr(j,0,n),
	!   the next element at fr(j,1,n) etc, and similarly for fi.

		use param
		use domain
		implicit none

		integer :: ftype,frow,mptr,midx,l

		real(IDP) :: facfl0
		real(IDP), dimension(0:,0:) :: f
		real(IDP), dimension(0:,0:,0:) :: fr,fi

	!   The following factor should be set to 1.0 for dtem, 0.5 for kite

		facfl0 = 1.0

		if (abs(ftype) == 1) then

			do frow = 1,nmax

				mptr = 0
				do midx = mmstart(frow),mmend(frow)
					l = ll(midx,frow)

					if (ftype == 1) then
						fr(:,mptr,frow) =  0.5*f(:,l)
						fi(:,mptr,frow) =  0.0
					else
						fr(:,mptr,frow) =  0.0
						fi(:,mptr,frow) = -0.5*f(:,l)
					end if

					mptr = mptr + 1
				end do

			end do
			
		end if

	!   Do frow = 0 as a special case, since it contains its own
	!   conjugate image

		mptr = 0
		do midx = mmstart(0), mmend(0)
			if (midx >= 0) then
				l = ll(midx,0)
			else
				l = ll(-midx,0)
			endif

			if (ftype > 0) then
				if (midx > 0) then
					fr(:,mptr,0) =  0.5*f(:,l)
					fi(:,mptr,0) =  0.0
				else if (midx < 0) then
					fr(:,mptr,0) = 0.5*f(:,l)
					fi(:,mptr,0) = 0.0
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			else
				if (midx > 0) then
					fr(:,mptr,0) =  0.0
					fi(:,mptr,0) = -0.5*f(:,l)
				else if (midx < 0) then
					fr(:,mptr,0) = 0.0
					fi(:,mptr,0) = 0.5*f(:,l)
				else
					fr(:,mptr,0) = facfl0 * f(:,l)
					fi(:,mptr,0) = 0.0
				end if
			end if
			mptr = mptr + 1
		end do

	end subroutine cs2rli
