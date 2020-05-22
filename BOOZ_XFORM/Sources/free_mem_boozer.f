      SUBROUTINE free_mem_boozer
      USE booz_params
      USE booz_persistent
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1=0, istat2=0, istat3=0
C-----------------------------------------------
      IF (ALLOCATED(jlist)) DEALLOCATE (jlist, lsurf_boz)

      IF (ALLOCATED(bsubumnc)) DEALLOCATE(
     1    bsubumnc, bsubvmnc, bmodmnc, rmnc, zmns, lmns,
     2    xm, xn, xm_nyq, xn_nyq, hiota, phip, gpsi, 
     3    ipsi, pmns, beta_vol, pres, phi, buco, bvco,
     3    rmncb, zmnsb, pmnsb, gmncb, bmncb, bmod_b, 
     4    chip, chi, stat=istat1 )
! CRCook deallocate chi's

      IF (ALLOCATED(cosm_b)) DEALLOCATE(
     1    cosm_b, sinm_b, cosn_b, sinn_b,
     2    cosm_nyq, sinm_nyq, cosn_nyq, sinn_nyq,
     2    sfull, scl, xmb, xnb, thgrd, ztgrd, stat=istat2)

      IF (ALLOCATED(bsubumns)) DEALLOCATE(
     1    bsubumns, bsubvmns, bmodmns, rmns, zmnc, lmnc, pmnc,
     2    rmnsb, zmncb, pmncb, gmnsb, bmnsb, stat=istat3)

      IF (istat1 .ne.0 .or. istat2 .ne. 0)
     1  PRINT *,' Deallocation error in Free_mem_boozer'

      END SUBROUTINE free_mem_boozer
