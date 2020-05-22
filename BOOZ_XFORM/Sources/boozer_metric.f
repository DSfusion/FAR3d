      SUBROUTINE boozer_metric
      USE booz_params
      USE booz_persistent, ONLY: nu2_b, thgrd, ztgrd, scl, xm, xn, 
     1   xmb, xnb, sfull
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, js, mn, mnb, m, n, istat1=0, istat2=0,
     1   istat3=0, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lboz
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gib, grrb, grtb, gttb,
     1   grrojb, grtojb, gttojb, jbgrrb, jbgrtb, jbgttb
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   cosm_boz, sinm_boz, cosn_boz, sinn_boz,
     2   gimncb, grrmncb, grtmnsb, gttmncb, grrojmncb, grtojmnsb,
     3   gttojmncb, jbgrrmncb, jbgrtmnsb, jbgttmncb, gimnsb, grrmnsb,
     4   grtmncb, gttmnsb, grrojmnsb, grtojmncb, gttojmnsb, jbgrrmnsb,
     5   jbgrtmncb, jbgttmnsb
      CHARACTER(LEN=1) :: tb
C-----------------------------------------------

      CALL write_polcut
!
!      theta-boz = thgrd:  now uniform angle mesh in BOOZER space
!      zeta-boz  = ztgrd
!
      ALLOCATE (cosm_boz(nunv,0:mboz), sinm_boz(nunv,0:mboz),
     1   cosn_boz(nunv,0:nboz), sinn_boz(nunv,0:nboz), stat=istat1)
      IF (istat1 .ne. 0) STOP 'Deallocation error in boozer_coords'

      ALLOCATE (gib(nunv), grrb(nunv), grtb(nunv), gttb(nunv),
     1  grrojb(nunv),grtojb(nunv),gttojb(nunv),
     2  jbgrrb(nunv),jbgrtb(nunv),jbgttb(nunv), stat=istat2)
      IF (istat2 .ne. 0) STOP 'Deallocation error in boozer_coords'

      ALLOCATE (gimncb(mnboz,jsize), grrmncb(mnboz,jsize),
     1   grtmnsb(mnboz,jsize), gttmncb(mnboz,jsize),
     2   grrojmncb(mnboz,jsize), grtojmnsb(mnboz,jsize),
     3   gttojmncb(mnboz,jsize), jbgrrmncb(mnboz,jsize),
     4   jbgrtmnsb(mnboz,jsize), jbgttmncb(mnboz,jsize), stat=istat3)
      IF (istat3 .ne. 0) STOP 'Deallocation error in boozer_coords'

      IF (lasym_b) THEN
        ALLOCATE (gimnsb(mnboz,jsize), grrmnsb(mnboz,jsize),
     6   grtmncb(mnboz,jsize), gttmnsb(mnboz,jsize),
     7   grrojmnsb(mnboz,jsize), grtojmncb(mnboz,jsize),
     8   gttojmnsb(mnboz,jsize), jbgrrmnsb(mnboz,jsize),
     9   jbgrtmncb(mnboz,jsize), jbgttmnsb(mnboz,jsize), stat=istat3)
       IF (istat3 .ne. 0) STOP 'Deallocation error in boozer_coords'
      END IF

      CALL trigfunc (thgrd, ztgrd, cosm_boz, sinm_boz, cosn_boz, 
     1               sinn_boz, mboz, nboz, nunv)
!
!      COMPUTE METRIC ELEMENTS, SURFACE BY SURFACE
!
       DO js = 1, jsize
         j = jlist(js)
         CALL vcoords_gijb (rmncb, zmnsb, pmnsb, gmncb, rmnsb, zmncb,
     1     pmncb, gmnsb, xmb, xnb, cosm_boz, sinm_boz, cosn_boz,
     2     sinn_boz, mboz, nboz, mnboz, js, j, jsize, gib,
     3     grrb, grtb, gttb, grrojb, grtojb, gttojb, jbgrrb, jbgrtb,
     4     jbgttb, nunv, nfp, lasym_b)

        IF (lasym_b) THEN
         CALL boozer_gij (thgrd, ztgrd, gib, grrb, grtb, gttb, grrojb,
     1     grtojb, gttojb, jbgrrb, jbgrtb, jbgttb, xmb, xnb,
     2     gimncb(1,js), grrmncb(1,js), grtmnsb(1,js), gttmncb(1,js),
     3     grrojmncb(1,js), grtojmnsb(1,js), gttojmncb(1,js),
     4     jbgrrmncb(1,js), jbgrtmnsb(1,js), jbgttmncb(1,js),
     5     gimnsb(1,js), grrmnsb(1,js), grtmncb(1,js), gttmnsb(1,js),
     6     grrojmnsb(1,js), grtojmncb(1,js), gttojmnsb(1,js),
     7     jbgrrmnsb(1,js), jbgrtmncb(1,js), jbgttmnsb(1,js),
     8     scl, cosm_boz, sinm_boz, cosn_boz, sinn_boz, mboz, nboz,
     9     mnboz, nfp, nunv, nu2_b, nv_boz, lasym_b)
        ELSE
         CALL boozer_gij (thgrd, ztgrd, gib, grrb, grtb, gttb, grrojb,
     1     grtojb, gttojb, jbgrrb, jbgrtb, jbgttb, xmb, xnb,
     2     gimncb(1,js), grrmncb(1,js), grtmnsb(1,js), gttmncb(1,js),
     3     grrojmncb(1,js), grtojmnsb(1,js), gttojmncb(1,js),
     4     jbgrrmncb(1,js), jbgrtmnsb(1,js), jbgttmncb(1,js),
     5     gimncb(1,js), grrmncb(1,js), grtmnsb(1,js), gttmncb(1,js),
     6     grrojmncb(1,js), grtojmnsb(1,js), gttojmncb(1,js),
     7     jbgrrmncb(1,js), jbgrtmnsb(1,js), jbgttmncb(1,js),
     8     scl, cosm_boz, sinm_boz, cosn_boz, sinn_boz, mboz, nboz,
     9     mnboz, nfp, nunv, nu2_b, nv_boz, lasym_b)
        END IF
       END DO

      tb=char(9)
      open(unit=92,file="prfeq_vmec.txt")
      write(92,'("r",a1,"phip",a1,"iota",a1,"I",a1,"g",a1,"pr")')
     1  (tb,j=1,5)
      do j=2,ns
        write(92,'(f7.5,5(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),tb,
     1    phip(j),tb,hiota(j),tb,buco(j),tb,bvco(j),tb,pres(j)
      end do
      close(unit=92)

      ALLOCATE (lboz(mnmax))
      lboz = 0
      DO mn=1,mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))
         DO mnb=1,mnboz
            IF (NINT(xmb(mnb)) == m .AND. NINT(xnb(mnb)) == n) EXIT
         END DO
         IF (mnb > mnboz) STOP 'VMEC mode not included in Boozer set'
         lboz(mn) = mnb
      END DO

      open(unit=92,file="rmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,rmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="zmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,zmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="pmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,pmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="bmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,bmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,gmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gimnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,gimncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grrmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grrmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grtmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grtmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gttmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,gttmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grrojmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grrojmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grtojmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grtojmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gttojmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,gttojmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgrrmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgrrmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgrtmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgrtmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgttmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgttmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)

      IF (lasym_b) THEN
      open(unit=92,file="rmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,rmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="zmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,zmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="pmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,pmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="bmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,bmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,gmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gimns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sqrt((j-1.5_dp)/(ns-1.0)),
     1    (tb,gimnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grrmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grrmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grtmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grtmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gttmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,gttmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grrojmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grrojmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="grtojmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,grtojmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="gttojmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,gttojmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgrrmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgrrmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgrtmnc_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgrtmncb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      open(unit=92,file="jbgttmns_booz.txt",recl=8192)
      write(92,'("r",500(a1,i4,"/",i4))') (tb,nint(xm(mn)),nint(xn(mn)),
     1  mn=1,mnmax)
      do js=1,jsize
        j = jlist(js)
        write(92,'(f7.5,500(a1,1pe13.6))') sfull(j),
     1    (tb,jbgttmnsb(lboz(mn),js),mn=1,mnmax)
      end do
      close(unit=92)
      END IF

      open(unit=9,file='woutb',status='new',convert='big_endian', 
     1  form='unformatted')
      write(9) nfp,mnboz,ns,lasym_b
      write(9) (NINT(xmb(mnb)), NINT(xnb(mnb)), mnb=1,mnboz)
      write(9) (phip(j),hiota(j),buco(j),bvco(j),pres(j),j=2,ns)
      do j=1,jsize
        write(9) (rmncb(mnb,j),zmnsb(mnb,j),pmnsb(mnb,j),bmncb(mnb,j),
     1   gmncb(mnb,j),gimncb(mnb,j),mnb=1,mnboz)
      end do
      do j=1,jsize
        write(9) (grrmncb(mnb,j),grtmnsb(mnb,j),gttmncb(mnb,j),
     1   grrojmncb(mnb,j),grtojmnsb(mnb,j),gttojmncb(mnb,j),
     2   jbgrrmncb(mnb,j),jbgrtmnsb(mnb,j),jbgttmncb(mnb,j),mnb=1,mnboz)
      end do

      IF (lasym_b) THEN
      do j=1,jsize
        write(9) (rmnsb(mnb,j),zmncb(mnb,j),pmncb(mnb,j),bmnsb(mnb,j),
     1   gmnsb(mnb,j),gimnsb(mnb,j),mnb=1,mnboz)
      end do
      do j=1,jsize
        write(9) (grrmnsb(mnb,j),grtmncb(mnb,j),gttmnsb(mnb,j),
     1   grrojmnsb(mnb,j),grtojmncb(mnb,j),gttojmnsb(mnb,j),
     2   jbgrrmnsb(mnb,j),jbgrtmncb(mnb,j),jbgttmnsb(mnb,j),mnb=1,mnboz)
      end do
      END IF

      close(unit=9)

      END SUBROUTINE boozer_metric
