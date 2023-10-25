:Reference :Colbert and Pan 2002

NEURON	{
	SUFFIX NaTa_t
	USEION na READ ena WRITE ina
	RANGE gNaTa_tbar, gNaTa_t, ina
        RANGE msh, mk, mmin, hsh, hk, hmin
        
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNaTa_tbar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0
        hsh = 0
        hk = 0
        hmin = 0  
}

ASSIGNED	{
          celsius (degC)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa_t	(S/cm2)
	mInf
	mTau    (ms)
	mAlpha
	mBeta
	hInf
	hTau    (ms)
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNaTa_t = gNaTa_tbar*m*m*m*h
	ina = gNaTa_t*(v-ena)
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
  LOCAL qt
  qt = 2.3^((celsius-21)/10)
	
  UNITSOFF
		mAlpha = (0.182 * 6)* efun(-(v- -38)/6)
		mBeta  = (0.124 * 6)* efun(-(-v -38)/6)
		mTau = (1/(mAlpha + mBeta))/qt
		mInf = mmin + (1-mmin) * 1/(1+exp(-(v+40.302+msh)/(6* (1+mk))))


		hAlpha = (0.015 *6)* efun((v- -66)/6)
		hBeta  = (0.015 *6)* efun((-v -66)/6)
		hTau = (1/(hAlpha + hBeta))/qt
		hInf = hmin + (1-hmin) * 1/(1+exp((v+66+hsh)/(6* (1+hk))))
	UNITSON
}
