    subroutine zzdisp(x,y,zzr,zzi)

      use param
      implicit none
      real(IDP) :: x,y,zzr,zzi,x1,y1,wzr1,wzi1,a,b,abr,abi

	interface
	 subroutine wzdisp(x,y,re,im)
	     use param
	     implicit none
             real(IDP) :: im,x,y,re
	 end subroutine wzdisp
	end interface

!		Subroutine required to calculate the Landau electron-ion damping effects
	
      x1=abs(x)
      y1=abs(y)
      call wzdisp(x1,y1,wzr1,wzi1)
      if(x.ge.0.0_IDP.and.y.ge.0.0_IDP) go to 1
      if(x.le.0.0_IDP.and.y.ge.0.0_IDP) go to 2
      a=2.0_IDP*x1*y1
      b=-(x1*x1-y1*y1)
      abr=2.0_IDP*exp(b)*cos(a)
      abi=-2.0_IDP*exp(b)*sin(a)
      wzr1=abr-wzr1
      wzi1=abi-wzi1
      if(x.le.0.0_IDP.and.y.le.0.0_IDP) go to 1
    2 wzi1=-wzi1
    1 zzr=-1.7724538509055_IDP*wzi1
      zzi=1.7724538509055_IDP*wzr1
      return
	  
    end subroutine zzdisp
	    
    subroutine wzdisp(x,y,re,im)

      use param
      implicit none      
      real(IDP) :: im,lambda,epsh,epsl,epsy,x,y,re,h,h2,ss,rr,ri,sr,si,tr,ti,cc,c,s
      integer :: capn, nu, nup, i, n, np1
      logical :: b

!		Subroutine required to calculate the Landau electron-ion damping effects
	  
      epsh=1.e-12_IDP
      epsl = epsh; epsy = epsh
      if(y.lt.4.29_IDP .and. x.lt.5.33_IDP) go to 10

!  (x,y) belongs to q1-r

      h=0.0_IDP
      capn=0
      nu=8
      go to 20

!  (x,y) belongs to r

   10 s=(1.0_IDP-y/4.29_IDP)*sqrt(1.0_IDP-x*x/28.41_IDP)
      h=1.6_IDP*s
      h2=2.0_IDP*h
      capn = int(6.0_IDP + 23.0_IDP*s + 0.5_IDP)
      nu = int(9.0_IDP + 21.0_IDP*s + 0.5_IDP)
      lambda=h2**capn
   20 b=(h.eq.0.0_IDP .or. lambda.lt.epsl)

!  statement (lambda.lt.epsl) covers the underflow case
!  when h(.gt.0) is very small.

      rr=0.0_IDP
      ri=0.0_IDP
      sr=0.0_IDP
      si=0.0_IDP
      nup=nu+1
      do 100 i=1,nup
      n=nup-i
      np1=n+1
      tr=y+h+np1*rr
      ti=x-np1*ri
      c=0.5_IDP/(tr*tr+ti*ti)
      rr=c*tr
      ri=c*ti
      if(.not.(h .gt. 0.0_IDP .and. n .le. capn)) go to 100
      tr=lambda+sr
      sr=rr*tr-ri*si
      si=ri*tr+rr*si
      lambda=lambda/h2
  100 continue
      cc=1.12837916709551_IDP
      if(y.lt.epsy) go to 120
      if(b) go to 110
      re=sr*cc
      go to 130
  110 re=rr*cc
      go to 130
  120 re=exp(-x*x)
  130 if(b) go to 140
      im=si*cc
      go to 150
  140 im=ri*cc
  150 return
  
    end subroutine wzdisp