:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Im
	USEION k READ ek WRITE ik
	RANGE gImbar, gIm, ik
        RANGE msh, mk, mmin
        RANGE mtmin, mtmax, mtsh, mtk
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gImbar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0

        mtmin = 0
        mtmax = 0
        mtsh = 0
        mtk = 0
}

ASSIGNED	{
          celsius (degC)
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm	(S/cm2)
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
	gIm = gImbar*m
	ik = gIm*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((celsius-21)/10)

	UNITSOFF
		mAlpha = 3.3e-3*exp(2.5*0.04*(v + 35 + mtsh) / (1 + mtk))
		mBeta = 3.3e-3*exp(-2.5*0.04*(v + 35 + mtsh) / (1 + mtk))
		mTau = mtmin + (1 + mtmax) * (1/(mAlpha + mBeta))/qt
		mInf = mmin + (1-mmin) * 1/(1 + exp(-(v+35+msh)/ ( 5 * (1+mk) ) ))
	UNITSON
}
