subroutine cnvt(idir)

  use param
  use processor
  use var_para
  use cotrol
  use domain
  use equil
  use dynamo
  use transfer
  use scratch
  implicit none

  integer :: idir,i,mnum2,mnum3,l,l1,l2,l1t,lp,ind,i1,lsk,lskp,imt1,j

  select case (idir)
      case (1)
        do ind=1,nvar
           if (ind == 4) cycle
           select case (ind)
              case (1)
                 call trnsfr(psip,1,2)
              case (2)
                 call trnsfr(uztp,-1,2)
              case (3)
                 call trnsfr(prp,1,2)
              case (5)
                 call trnsfr(nfpp,1,2)
              case (6)
                 call trnsfr(vprlfp,-1,2)
              case (7)
                 call trnsfr(vthprlfp,-1,2)
              case (8)
                 call trnsfr(nalpp,1,2)
              case (9)
                 call trnsfr(vprlalpp,-1,2)
           end select
           l1=0
           do i=n_start,n_end
              mnum3=noeqn*mnumn(i)
              lsk=(ind-1)*mnumn(i)+nskpxn(i)
              do l1t=1,mnumn(i)
                 l=l1t+lnumn(i-1)
                 lp=lln(l)
                 l1=l1+1
                 lskp=l1t+lsk
                 do j=1,mj
                    imt1=lskp+mnum3*(j-1)
                    yt(imt1)=scp(l1,j)
                 end do
              end do
           end do
        end do
     case (2)
        do ind=1,nvar
           l1=0
           do i=n_start,n_end
              mnum3=noeqn*mnumn(i)
              lsk=(ind-1)*mnumn(i)+nskpxn(i)
              do l1t=1,mnumn(i)
                 l=l1t+lnumn(i-1)
                 lp=lln(l)
                 l1=l1+1
                 lskp=l1t+lsk
                 do j=1,mj
                    imt1=lskp+mnum3*(j-1)
                    scp(l1,j)=xt(imt1)
                 end do
                 if (mm(lp) == 0) then
                    scp(l1,0)=(r(2)**2*scp(l1,1)-r(1)**2*scp(l1,2))/(r(2)**2-r(1)**2)
                 else
                    scp(l1,0)=0.0_IDP
                 end if
                 if (ind > 4) scp(l1,mj)=(scp(l1,mjm1)*(r(mj)-r(mjm2))**2-scp(l1,mjm2)*(r(mj)-r(mjm1))**2)/ &
                         ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
                 ! if (((ind == 5 .or. ind == 8) .and. lp == l0) .or. ind == 6 .or. ind == 7 .or. ind == 9) &
                 !    scp(l1,mj)=(scp(l1,mjm1)*(r(mj)-r(mjm2))**2-scp(l1,mjm2)*(r(mj)-r(mjm1))**2)/ &
                 !               ((r(mj)-r(mjm2))**2-(r(mj)-r(mjm1))**2)
              end do
           end do
           select case (ind)
              case (1)
                 call trnsfr(psi,1,1)
              case (2)
                 call trnsfr(phi,-1,1)
              case (3)
                 call trnsfr(pr,1,1)
              case (4)
                 call trnsfr(uzt,-1,1)
              case (5)
                 call trnsfr(nf,1,1)
              case (6)
                 call trnsfr(vprlf,-1,1)
              case (7)
                 call trnsfr(vthprlf,-1,1)
              case (8)
                 call trnsfr(nalp,1,1)
              case (9)
                 call trnsfr(vprlalp,-1,1)
           end select
        end do
      case (3)
        do ind=1,nvar
           if (ind == 4) cycle
           l1=0
           do i=n_start,n_end
              mnum3=noeqn*mnumn(i)
              lsk=(ind-1)*mnumn(i)+nskpxn(i)
              do l1t=1,mnumn(i)
                 l=l1t+lnumn(i-1)
                 lp=lln(l)
                 l1=l1+1
                 lskp=l1t+lsk
                 do j=1,mj
                    imt1=lskp+mnum3*(j-1)
                    scp(l1,j)=yt(imt1)
                 end do
                 if (mm(lp) == 0) then
                    scp(l1,0)=(r(2)**2*scp(l1,1)-r(1)**2*scp(l1,2))/(r(2)**2-r(1)**2)
                 else
                    scp(l1,0)=0.0_IDP
                 end if
              end do
           end do
           select case (ind)
              case (1)
                 call trnsfr(psip,1,1)
              case (2)
                 call trnsfr(uztp,-1,1)
              case (3)
                 call trnsfr(prp,1,1)
              case (5)
                 call trnsfr(nfpp,1,1)
              case (6)
                 call trnsfr(vprlfp,-1,1)
              case (7)
                 call trnsfr(vthprlfp,-1,1)
              case (8)
                 call trnsfr(nalpp,1,1)
              case (9)
                 call trnsfr(vprlalpp,-1,1)
           end select
        end do
  end select

end subroutine cnvt
