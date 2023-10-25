:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica 
        RANGE msh, mk, mmin, hsh, hk, hmin
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_HVAbar = 0.00001 (S/cm2)
        msh = 0
        mk1 = 0
        mk2 = 0
        mmin = 0
        hsh = 0
        hk1 = 0
        hk2 = 0
        hmin = 0         
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_HVA	(S/cm2)
	mInf    (1)
	mTau    (ms)
	mAlpha  (1)
	mBeta   (1)
	hInf    (1)
	hTau    (ms)
	hAlpha  (1)
	hBeta   (1)
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_HVA = gCa_HVAbar*m*m*h
	ica = gCa_HVA*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
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
		mAlpha =  (0.055*3.8)*efun(-(27+v+msh)/(3.8*(1+mk1)))
		mBeta  =  (0.94*exp(-(75+v+msh)/(17*(1+mk2))))
		mInf = mmin + (1-mmin) * mAlpha/(mAlpha + mBeta)
        
        	mAlpha =  (0.055*3.8)*efun(-(27+v)/(3.8))
		mBeta  =  (0.94*exp(-(75+v)/(17)))
		mTau = 1/(mAlpha + mBeta)
        
		hAlpha =  (0.000457*exp(-(13+v+hsh)/(50*(1+hk1))))
		hBeta  =  (0.0065/(exp(-(v+15+hsh)/(28*(1+hk2)))+1))
		hInf = hmin + (1-hmin) * hAlpha/(hAlpha + hBeta)


		hAlpha =  (0.000457*exp(-(13+v)/(50)))
		hBeta  =  (0.0065/(exp(-(v+15)/(28))+1))
                hTau = 1/(hAlpha + hBeta)
	UNITSON
}
