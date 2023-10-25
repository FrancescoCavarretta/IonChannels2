:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn 
        RANGE msh, mk1, mk2, mmin
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gIhbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)
        msh = 0
        mk1 = 0
        mk2 = 0
        mmin = 0
}

ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	gIh	(S/cm2)
	mInf
	mTau    (ms)
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIh = gIhbar*m
	ihcn = gIh*(v-ehcn)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

FUNCTION efun(z) {
	 if (fabs(z) < 1e-4) {
	    efun = 1 - z/2
	 }else{
	    efun = z/(exp(z) - 1)
         }
}

PROCEDURE rates(){
	UNITSOFF

		mAlpha =  0.001*6.43*11.9* efun((v+154.9+msh)/(11.9* (1+mk1)))
		mBeta  =  0.001*193*exp((v+msh)/(33.1* (1+mk2)))
		mInf = mmin + (1-mmin) * mAlpha/(mAlpha + mBeta)

        	mAlpha =  0.001*6.43*11.9* efun((v+154.9)/11.9)
		mBeta  =  0.001*193*exp(v/33.1)
		mTau = 1/(mAlpha + mBeta)
	UNITSON
}
