:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX K_Tst
	USEION k READ ek WRITE ik
	RANGE gK_Tstbar, gK_Tst, ik
        RANGE msh, mk, mmin, hsh, hk, hmin
        RANGE mtmin, mtmax, mtsh, mtk, htmin, htmax, htsh, htk  
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gK_Tstbar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0
        hsh = 0
        hk = 0
        hmin = 0
        
        mtmin = 0
        mtmax = 0
        mtsh = 0
        mtk = 0

        htmin = 0
        htmax = 0
        htsh = 0
        htk = 0 
        
}

ASSIGNED	{
          celsius (degC)
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Tst	(S/cm2)
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
	gK_Tst = gK_Tstbar*(m^4)*h
	ik = gK_Tst*(v-ek)
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
		mInf =  mmin + (1-mmin) *1/(1 + exp(-(v+0+msh)/(19*(1+mk))))
		mTau =  ( mtmin+0.34+(1+htmax)*0.92*exp(-((v+71+mtsh)/(59*(1+mtk)))^2) )/qt
		hInf =  hmin + (1-hmin) *1/(1 + exp(-(v+66+hsh)/(-10*(1+hk))))
		hTau =  ( htmin+8+(1+htmax)*49*exp(-((v+73+htsh)/(23*(1+htk)))^2) )/qt
		v = v - 10
	UNITSON
}
