!  u-zeta expression

!	WARNING: sd1 is feq-qqinv*cureq

		sd1=feq-qqinv*cureq
		do l=1,leqmax
			sceq1(:,l)=denseq*(jbgtt(:,l)-rinv*rinv*cureq*cureq*jsq(:,l)/(epsq*sd1))
		end do
		call blockj(sceq1,1,4,2,0,2,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=-denseq*(jbgrt(:,l)-cureq*bstg(:,l)/(epsq*sd1))
		end do
		call dbydreq(sceq3,sceq1,0.0_IDP,1.0_IDP,2)
		do l=1,leqmax
			sceq3(:,l)=sceq3(:,l)+rinv*sceq1(:,l)
		end do
		call dbydtheq(sceq3,sceq2,-1,1.0_IDP,1.0_IDP,2)
		call blockj(sceq3,1,4,2,0,1,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=-denseq*(feq*jbgrt(:,l)-r*r*qqinv*bsjgtt(:,l)-cureq*bstg(:,l)/epsq)/sd1
		end do
		sceq3=sceq1+sceq2
		call blockj(sceq3,-1,4,2,1,1,0,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=denseq*(feq*jbgrr(:,l)-r*r*(qqinv*bsjgrt(:,l)+bsqg(:,l)/epsq))/sd1
		end do
		call blockj(sceq2,1,4,2,2,0,0,1.0_IDP)
		call dbydreq(sceq3,sceq1,0.0_IDP,1.0_IDP,2)
		call dbydtheq(sceq3,sceq2,1,1.0_IDP,1.0_IDP,2)
		call blockj(sceq3,-1,4,2,1,0,0,1.0_IDP)
		do l=1,leqmax
			sceq1(:,l)=denseq*(rinv*cureq*jbgrt(:,l)-r*bsjgtt(:,l))/sd1
		end do
		call blockj(sceq1,-1,4,2,0,1,1,1.0_IDP)
		do l=1,leqmax
			sceq2(:,l)=-denseq*(rinv*cureq*jbgrr(:,l)-r*bsjgrt(:,l))/sd1
		end do
		call blockj(sceq2,1,4,2,1,0,1,1.0_IDP)
		call dbydreq(sceq3,sceq1,0.0_IDP,1.0_IDP,3)
		do l=1,leqmax
			sceq3(:,l)=sceq3(:,l)+rinv*sceq1(:,l)
		end do
		call dbydtheq(sceq3,sceq2,1,1.0_IDP,1.0_IDP,3)
		call blockj(sceq3,-1,4,2,0,0,1,1.0_IDP)

!	END WARNING
