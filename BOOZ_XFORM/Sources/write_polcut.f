      SUBROUTINE write_polcut
      USE booz_params
      USE booz_persistent, ONLY: sfull, xm, xn, xnb, xmb, rmnc, zmns,
     1   rmns, zmnc
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, js, kr, mn, iflag
      REAL(rprec) :: twopi, dtheta, hs, rp, rsq, delr, phik, phia, phib,
     1   phic, phit, fa, fb, fc, ft, fab, fbc, dbc, abserr, relerr
      INTEGER, DIMENSION(15) :: jp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: thetav
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: phiv, rvl, zvl
      CHARACTER(LEN=1) :: tb
      CHARACTER(LEN=16) :: confil
C-----------------------------------------------
      ALLOCATE (thetav(101))
      twopi = 8*ATAN(1.0_dp)
      dtheta=0.01*twopi
      do i=1,101
       thetav(i)=(i-1)*dtheta
      end do

      hs=one/ohs
      jp(1)=2
      jp(15)=ns
      rp=sfull(2)
      delr=(1.0_dp-sfull(2))/14.0
      do j=2,14
       rp=sfull(2)+(j-1)*delr
       rsq=rp*rp
       do js=jp(j-1),ns-1
        if (hs*js >= rsq) exit
       end do
       jp(j)=js+1
       if (js > jp(j-1) .and. (rsq-hs*(js-1)) < (hs*js-rsq)) jp(j)=js
      end do

      ALLOCATE (phiv(101,15,3), rvl(101,15,3), zvl(101,15,3))
      tb=char(9)

      phik=0.0

      do kr=1,3

       if (nboz == 0 .and. kr > 1) cycle

       do js=1,15
        j=jp(js)

        do i=1,101
         rvl(i,js,kr)=0.0
         zvl(i,js,kr)=0.0
         do mn=1,mnmax
          rvl(i,js,kr)=rvl(i,js,kr)+rmnc(mn,j)*
     1       cos(xm(mn)*thetav(i)-xn(mn)*phik)
          zvl(i,js,kr)=zvl(i,js,kr)+zmns(mn,j)*
     1       sin(xm(mn)*thetav(i)-xn(mn)*phik)
         end do
         if (lasym_b) then
          do mn=1,mnmax
           rvl(i,js,kr)=rvl(i,js,kr)+rmns(mn,j)*
     1       sin(xm(mn)*thetav(i)-xn(mn)*phik)
           zvl(i,js,kr)=zvl(i,js,kr)+zmnc(mn,j)*
     1       cos(xm(mn)*thetav(i)-xn(mn)*phik)
          end do
         end if
        end do

       end do

       write(confil,'("R_Z_",i1,".txt")') kr
       open(unit=92,file=confil,recl=512)
       write(92,'("R(",i3,")",a1,"Z(",i3,")",14(a1,"R(",i3,")",a1,"Z(",
     1       i3,")"))') jp(1),tb,jp(1),(tb,jp(j),tb,jp(j),j=2,15)
       do i=1,101
        write(92,'(1pe13.6,29(a1,1pe13.6))') rvl(i,1,kr),tb,zvl(i,1,kr),
     1       (tb,rvl(i,js,kr),tb,zvl(i,js,kr),js=2,15)
       end do
       close(unit=92)

       phik=phik+0.25*twopi/nfp
      end do

      phik=0.0
      abserr=1.e-10
      abserr=1.e-8

      do kr=1,3

       if (nboz == 0 .and. kr > 1) cycle

       do js=1,15
        j=jp(js)-1

        do i=1,101
         if (i == 1) then
          phib=phik
         else
          phib=phiv(i-1,js,kr)
         end if
         fb=phib
         do mn=1,mnboz
          fb=fb+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phib)
         end do
         if (lasym_b) then
          do mn=1,mnboz
           fb=fb+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phib)
          end do
         end if
         fb=fb-phik
         if (fb /= 0.0_dp) then 
          phic=phib+0.01
          fc=phic
          do mn=1,mnboz
           fc=fc+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phic)
          end do
          if (lasym_b) then
           do mn=1,mnboz
            fc=fc+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phic)
           end do
          end if
          fc=fc-phik
          phia=phib-0.01
          fa=phia
          do mn=1,mnboz
           fa=fa+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phia)
          end do
          if (lasym_b) then
           do mn=1,mnboz
            fa=fa+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phia)
           end do
          end if
          fa=fa-phik
          fab=fa*fb
          fbc=fb*fc
          if (fab < 0.0_dp) then
           phic=phib
           phib=phia
          else if (fbc > 0.0_dp) then
           dbc=(fc-fb)*fb
           if (dbc < 0.0_dp) then
            do
             phic=phic+0.01
             fc=phic
             do mn=1,mnboz
              fc=fc+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phic)
             end do
             if (lasym_b) then
              do mn=1,mnboz
               fc=fc+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phic)
              end do
             end if
             fc=fc-phik
             if (fb*fc < 0.0_dp) exit
            end do
            phib=phic-0.01
           else
            do
             phia=phia-0.01
             fa=phia
             do mn=1,mnboz
              fa=fa+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phia)
             end do
             if (lasym_b) then
              do mn=1,mnboz
               fa=fa+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phia)
              end do
             end if
             fa=fa-phik
             if (fb*fa < 0.0_dp) exit
            end do
            phib=phia
            phic=phia+0.01
           end if
          end if
          iflag=1
          call root(phit,ft,phib,phic,relerr,abserr,iflag)
          do while (iflag < 0)
           ft=phit
           do mn=1,mnboz
            ft=ft+pmnsb(mn,j)*sin(xmb(mn)*thetav(i)-xnb(mn)*phit)
           end do
           if (lasym_b) then
            do mn=1,mnboz
             ft=ft+pmncb(mn,j)*cos(xmb(mn)*thetav(i)-xnb(mn)*phit)
            end do
           end if
           ft=ft-phik
           call root(phit,ft,phib,phic,relerr,abserr,iflag)
          end do
         end if
         phiv(i,js,kr)=phib
        end do

        do i=1,101
         rvl(i,js,kr)=0.0
         zvl(i,js,kr)=0.0
         do mn=1,mnboz
          rvl(i,js,kr)=rvl(i,js,kr)+rmncb(mn,j)*
     1       cos(xmb(mn)*thetav(i)-xnb(mn)*phiv(i,js,kr))
          zvl(i,js,kr)=zvl(i,js,kr)+zmnsb(mn,j)*
     1       sin(xmb(mn)*thetav(i)-xnb(mn)*phiv(i,js,kr))
         end do
         if (lasym_b) then
          do mn=1,mnboz
           rvl(i,js,kr)=rvl(i,js,kr)+rmnsb(mn,j)*
     1       sin(xmb(mn)*thetav(i)-xnb(mn)*phiv(i,js,kr))
           zvl(i,js,kr)=zvl(i,js,kr)+zmncb(mn,j)*
     1       cos(xmb(mn)*thetav(i)-xnb(mn)*phiv(i,js,kr))
          end do
         end if
        end do

       end do

       write(confil,'("R_Z_",i1,"b.txt")') kr
       open(unit=92,file=confil,recl=512)
       write(92,'("R(",i3,")",a1,"Z(",i3,")",14(a1,"R(",i3,")",a1,"Z(",
     1       i3,")"))') jp(1),tb,jp(1),(tb,jp(j),tb,jp(j),j=2,15)
       do i=1,101
        write(92,'(1pe13.6,29(a1,1pe13.6))') rvl(i,1,kr),tb,zvl(i,1,kr),
     1       (tb,rvl(i,js,kr),tb,zvl(i,js,kr),js=2,15)
       end do
       close(unit=92)

       write(confil,'("phi_",i1,"b.txt")') kr
       open(unit=92,file=confil,recl=512)
       write(92,'("theta",15(a1,"phi(",i3,")"))') (tb,jp(j),j=1,15)
       do i=1,101
        write(92,'(f7.2,15(a1,f10.5))') 360.0*thetav(i)/twopi,
     1       (tb,360.0*phiv(i,js,kr)/twopi,js=1,15)
       end do
       close(unit=92)

       phik=phik+0.25*twopi/nfp
      end do

      END SUBROUTINE write_polcut
