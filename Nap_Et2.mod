:Comment : mtau deduced from text (said to be 6 times faster than for NaTa)
:Comment : so I used the equations from NaT and multiplied by 6
:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Nap_Et2
	USEION na READ ena WRITE ina
	RANGE gNap_Et2bar, gNap_Et2, ina
        RANGE msh, mk, mmin, hsh, hk, hmin
        RANGE mtmin, mtmax, mtsh, mtk1, mtk2, htmin, htmax, htsh, htk1, htk2        
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNap_Et2bar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0
        hsh = 0
        hk = 0
        hmin = 0

        mtmin = 0
        mtmax = 0
        mtsh = 0
        mtk1 = 0
        mtk2 = 0

        htmin = 0
        htmax = 0
        htsh = 0
        htk1 = 0
        htk2 = 0        
}

ASSIGNED	{
          celsius (degC)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNap_Et2	(S/cm2)
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
	gNap_Et2 = gNap_Et2bar*m*m*m*h
	ina = gNap_Et2*(v-ena)
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
		mInf = mmin + (1-mmin) * 1.0/(1+exp((v+52.6+msh)/(-4.6* (1+mk))))

		mAlpha = (0.182 *6) * efun(-(v + 38 + mtsh)/(6 * (1 + mtk1)))
		mBeta  = (0.124 *6) * efun( (v + 38 + mtsh)/(6 * (1 + mtk2)))
		mTau = ( mtmin + (1 + mtmax) * 6 / (mAlpha + mBeta) )/qt



		hInf = hmin + (1-hmin) * 1.0/(1+exp((v+48.8+hsh)/(10* (1+hk))))
  
                hAlpha = -2.88e-6*4.63 * efun( (v + 17 + htsh)/(4.63  * (1 + htk1)))
                hBeta  = -6.94e-6*2.63 * efun(-(v + 64.4 + htsh)/(2.63 * (1 + htk2)))
		hTau = ( htmin + (1 + htmax) * 1 / (hAlpha + hBeta) )/qt
	UNITSON
}
