      SUBROUTINE allocate_boozer (iread)
      USE booz_params
      USE booz_persistent
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, jrad, istat1=0, istat2=0, istat3=0, iread, index
      CHARACTER(LEN=1000) :: temp
      CHARACTER(LEN=10)   :: scanset='0123456789'
c-----------------------------------------------
      IF (.not.ALLOCATED(jlist)) ALLOCATE (jlist(ns), lsurf_boz(ns), 
     1    stat=istat1)
      IF (istat1 .ne. 0) STOP 'Unable to allocate jlist/lsurf_boz'

!
!     Read in (and parse) list of radial surfaces to compute
!     Surfaces SHOULD ALL be on a single line, but IFORT doesn't write it out that way
!     so the loop below will correctly read multiple lines from the in_booz file
!
      jlist = 0
      i = 1
      READ (iread, '(a)', iostat=istat1) temp
      IF (far) istat1=1  ! All surfaces are computed for metric elements calculation
      IF (istat1 .eq. 0) THEN
         DO WHILE (istat1 .eq. 0)
            DO jrad = i, ns
               index = SCAN(temp,scanset)
               IF (index < 1) EXIT
               temp = temp(index:)
               READ(temp, *) jlist(jrad)
               IF (jlist(jrad) < 10) THEN
                  temp = temp(3:)
               ELSE IF (jlist(jrad) < 100) THEN
                  temp = temp(4:)
               ELSE
                  temp = temp(5:)
               END IF
            END DO
            READ (iread, '(a)', iostat=istat1) temp
            i = jrad
         END DO

         lsurf_boz = .FALSE.
         DO jrad = 1, ns
            i = jlist(jrad)
            IF (i.gt.1 .and. i.le.ns) lsurf_boz(i) = .TRUE.
         END DO

      ELSE IF (istat1 .ne. 0) THEN 
         WRITE(6, '(a,/,a,a,i4)')
     1    ' No jlist data was found in Boozer input file.',
     1    ' Will assume that all surfaces are needed.',
     1    ' Iostat: ', istat1
         lsurf_boz = .TRUE. 
         lsurf_boz(1) = .FALSE.

      END IF

      jsize = COUNT(lsurf_boz(1:ns))
!
!     Recompute jlist, just in case user used unordered (or repeated) indices
!
      DEALLOCATE (jlist)
      ALLOCATE (jlist(jsize), stat=istat1)
      IF (istat1 .ne. 0) STOP 'Unable to allocate jlist'

      i = 1
      DO jrad = 2, ns
         IF (lsurf_boz(jrad)) THEN
            jlist(i) = jrad
            i = i+1
         END IF
      END DO

      IF (.not.ALLOCATED(bsubumnc)) ALLOCATE(
     1    bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns),
     1    bmodmnc(mnmax_nyq,ns),
     2    rmnc(mnmax,ns), zmns(mnmax,ns), lmns(mnmax,ns),
     3    xm(mnmax), xn(mnmax),
     3    xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq),
     4    hiota(ns), phip(ns), gpsi(ns), ipsi(ns), pmns(mnmax_nyq),
     5    beta_vol(ns), pres(ns), phi(ns), buco(ns), bvco(ns),
     5    rmncb(mnboz,jsize), zmnsb(mnboz,jsize), pmnsb(mnboz,jsize), 
     6    gmncb(mnboz,jsize), bmncb(mnboz,jsize), 
     7    bmod_b(nv_boz,nu_boz), chip(ns), chi(ns), stat=istat1 )
! CRCook allocated chi (poloidal flux)

      IF (.not.ALLOCATED(sfull)) ALLOCATE(
     1    sfull(ns), scl(mnboz), xmb(mnboz), xnb(mnboz), stat=istat2)

      IF (.not.ALLOCATED(bsubumns)) ALLOCATE(
!      IF (lasym_b .AND. .not.ALLOCATED(bsubumns)) ALLOCATE(
     1    bsubumns(mnmax_nyq,ns), bsubvmns(mnmax_nyq,ns),
     1    bmodmns(mnmax_nyq,ns), 
     1    rmns(mnmax,ns), zmnc(mnmax,ns), lmnc(mnmax,ns),
     4    pmnc(mnmax_nyq),
     5    rmnsb(mnboz,jsize), zmncb(mnboz,jsize), pmncb(mnboz,jsize), 
     6    gmnsb(mnboz,jsize), bmnsb(mnboz,jsize), 
     1    stat=istat3)

      IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) THEN
          PRINT *,' problem allocating boozer memory'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2,
     1            ' istat3 = ',istat3
          STOP
      ENDIF

      END SUBROUTINE allocate_boozer
