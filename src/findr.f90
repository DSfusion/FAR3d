subroutine findr

! this sub sets up the r grid

! requires
! nis==number of points in the island
! ni==number of points interior to the island
! ne==number of points exterior to the island
! mj==total number of points=ni+nis+ne
! delta==width of the uniform fine grid (island)
! rc==center of the fine grid (island)
! fti==fraction of interior point in transition region
! fte==fraction of exterior points in transition region

  use mpi
  use param
  use processor
  use domain
  use findrf
  use tools
  
  implicit none

  integer :: ni0,ne0,nis0,nset,nte,nti,maxj,j,jt,jj,minj,ierror
  real(IDP) :: fe,fi,de,di,dr,ss,b

  fe=0
  fi=0
  ni0=ni
  ne0=ne
  nis0=nis
  nset=0
  if (ni < 2) then
     nset=2
     ni=2
  end if
  if (ne < 3) then
     nset=2
     ne=3
  end if
  if(nset > 1) nis=mj-ne-ni
  if(nis < 2) then
     if (myPE == 0) write(6,'(" grid: nis .lt. 2")')
     stop
  end if
  ne=ne-1
  nte=fte*ne
  nti=fti*ni
  if(nte < 1) nte=1
  if(nti < 1) nti=1
  dx=delta/(nis-1.)
  de=1.-(rc+0.5*delta)
  di=1.-delta-de
  if (de <= 0.) then
     if (myPE == 0) write(6,'(" grid: de .le. 0")')
     stop
  end if
  if (di <= 0.) then
     if (myPE == 0) write(6,'(" grid: di .le. 0")')
     stop
  end if
  n0=ne
  nt=nte
  d=de
  b=10.0_IDP
  if (nt > range(d)) b=10.0_IDP**(1.0_IDP*range(d)/nt)
  fe=zeroin(0.0_IDP,b,1.e-10_IDP)
  n0=ni
  nt=nti
  d=di
  b=10.0_IDP
  if (nt > range(d)) b=10.0_IDP**(1.0_IDP*range(d)/nt)
  fi=zeroin(0.0_IDP,b,1.e-10_IDP)
  dr=dx*fi**nti
  ss=0.
  maxj=ni-nti+1
  do j=1,maxj
     ss=ss+dr
     r(j)=ss
  end do
  do j=1,nti
     jt=nti-j+1
     dr=dx*fi**jt
     ss=ss+dr
     jj=maxj+j
     r(jj)=ss
  end do
  minj=ni+2
  maxj=ni+nis
  do j=minj,maxj
     ss=ss+dx
     r(j)=ss
  end do
  do j=1,nte
     dr=dx*fe**j
     ss=ss+dr
     jj=maxj+j
     r(jj)=ss
  end do
  minj=ni+nis+nte+1
  maxj=ni+nis+ne
  do j=minj,maxj
     ss=ss+dr
     r(j)=ss
  end do
  maxj=maxj+1
  if (mj /= maxj) then
     if (myPE == 0) write(6,'(" grid: mj .ne. ni+ne+nis")')
     CALL MPI_FINALIZE(ierror)
     stop
  end if
  r(mj)=1.0_IDP
  ne=ne0
  ni=ni0
  nis=nis0
  r(0)=0.

end subroutine findr
