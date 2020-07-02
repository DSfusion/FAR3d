	subroutine dfault

		use param
		use cotrol
		use domain
		use equil
		use dynamo
		implicit none

!		Default values of the model variables and parameters						
		
		numrun(1)="zz"
		numrun(2)="zz"
		numrun(3)="0"
		numruno(1)="fs"
		numruno(2)="##"
		numruno(3)="#"
		numruns(1)="zz"
		numruns(2)="zz"
		numruns(3)="0"
		numvac="vvvvv"
		ihist=0
		maxstp=100000
		nstep=0
		ndump=100000
		nprint=100
		lplots=0
		itime=1
		dt0=.66			
		nstep1=0
		nonlin=1	
		mj=0
		lmax=0
		leqmax=0
		mmin=0
		mmax=0
		nmin=0
		nmax=0
		mmineq=0
		mmaxeq=0 				
		nmineq=0
		nmaxeq=0
		mjm1=mj-1		
		mjm2=mj-2
		r=0.
		mm=0 	
		nn=0
		mmeq=0
		nneq=0
		l0=0
		leq0=0
		rinv=0. 	
		dc1m=0.
		dc1p=0.
		dc2m=0.
		dc2p=0.
		del2cm=0.
		del2cp=0.
		signl=0
		sgnleq=0
		lhmax=0		
		mh=0
		nh=0
		lheqmx=0
		mheq=0
		nheq=0
		m0dy=0
		nocpl=0
		ni=mj/3
		ne=mj/3
		nis=mj-ni-ne
		delta=.333
		Adens=7.0
		Bdens=0.5	
		LcA0=0.5
		LcA1=-0.88623
		LcA2=1
		LcA3=1	
        ext_prof=0	
		iflr_on=0
        ieldamp_on=0		
		epflr_on=0
        twofl_on=0	
		dpres=0.5
		DIIID_u=0
        omegar=0	
        iflr=0	
        r_epflr=0	
        r_epflralp=0		
		rc=.5
		fti=.95
		fte=.95  				
		qq=0.
		qqinv=0.
		psieq=0.
		chieq=0.
		preq=0.
		feq=0.
		cureq=0.
		denseq=1.
		teeq=1.
                                vzt_eq=0.
                                vtherm_ion=0.
		ngeneq=0.
		ndevice(1)="        "
		ndevice(2)="        "
		eps=0.		
		bet0=1.
		rs=0.
		leq=1
		nfp=1
		eta=0.
		etann=0.
		dt=0.
		time=0.
		etascl=1.
		reta=0.
		eta0=0.
		etalmb=0.
		ietaeq=0
		s=1.e+5
		gamma=0.
		ipert=0
		widthi=0.
		gammai=0.
		pertscl=1.
		stdifp=0.
		stdifu=0.
		stdifv=0.
		stdifnf=0.
		stdifvf=0.
		cnep=0.
		ctep=0.
		cvep=0.
        LcA0alp=0.
		LcA1alp=0.
		LcA2alp=0.
		LcA3alp=0.
		alpha_on=0.
		stdifnalp=0.
		stdifvalp=0.
		Adensalp=0.
		Bdensalp=0.
		bet0_alp=0.		
        cnfpalp=0.
		cvfpalp=0.		
		omcyalp=0.
		betath_factor=1.
		EP_dens_on=0.
		Alpha_dens_on=0.
		EP_vel_on=0.
		Alpha_vel_on=0.
		ext_prof_len=100
                vtk_on=0
        deltaq=0.
		deltaiota=0.	
                                Eq_vel_on=0.
                                q_prof_on=0.
                                omcyb=0.
	                rbound=0.
                                Trapped_on=0.

	end subroutine dfault