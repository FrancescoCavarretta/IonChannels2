:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Ca_LVAst
	USEION ca READ eca WRITE ica
	RANGE gCa_LVAstbar, gCa_LVAst, ica
        RANGE msh, mk, mmin, hsh, hk, hmin
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_LVAstbar = 0.00001 (S/cm2)
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
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAst	(S/cm2)
	mInf
	mTau    (ms)
	hInf
	hTau    (ms)
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_LVAst = gCa_LVAstbar*m*m*h
	ica = gCa_LVAst*(v-eca)
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

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((celsius-21)/10)

	UNITSOFF
		v = v + 10
		mInf = mmin + (1-mmin) * 1.0000/(1+ exp((v - -30.000 + msh)/(-6*(1+mk))))
		mTau = (5.0000 + 20.0000/(1+exp((v - -25.000)/5)))/qt
		hInf = hmin + (1-hmin) * 1.0000/(1+ exp((v - -80.000 + hsh)/(6.4*(1+hk))))
		hTau = (20.0000 + 50.0000/(1+exp((v - -40.000)/7)))/qt
		v = v - 10
	UNITSON
}
