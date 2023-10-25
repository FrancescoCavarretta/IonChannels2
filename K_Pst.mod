:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21


NEURON	{
	SUFFIX K_Pst
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, ik
        RANGE msh, mk, mmin, hsh, hk, hmin
        RANGE mtmin, mtmax1, mtmax2, mtsh, mtk1, mtk2, htmin, htmax1, htmax2, htsh, htk1, htk2 
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gK_Pstbar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0
        hsh = 0
        hk = 0
        hmin = 0
        
        mtmin = 0
        mtmax1 = 0
        mtmax2 = 0
        mtsh = 0
        mtk1 = 0
        mtk2 = 0

        htmin = 0
        htmax1 = 0
        htmax2 = 0
        htsh = 0
        htk1 = 0     
        htk2 = 0         
}

ASSIGNED	{
          celsius (degC)
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Pst	(S/cm2)
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
	gK_Pst = gK_Pstbar*m*m*h
	ik = gK_Pst*(v-ek)
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
		mInf =  mmin + (1-mmin) * (1/(1 + exp(-(v+1+msh)/(12*(1+mk)))))
  
        if (v < (-50-mtsh)) {
            mTau = ( mtmin + 1.25 + (1+mtmax1)*175.03*exp( (v + mtsh) * 0.026 / (1 + mtk1) ) )/qt
        } else {
            mTau = ( mtmin + 1.25 + (1+mtmax2)*13    *exp(-(v + mtsh) * 0.026 / (1 + mtk2) ) )/qt
        }
		hInf =  hmin + (1-hmin) * 1/(1 + exp(-(v+54+hsh)/(-11*(1+hk))))

		hTau =  ( htmin + 360 + htmax1 * 1010  * exp(-((v+75+htsh)/(48 * (1 + htk1)))^2) + htmax2 * 24 * (v + 55 + htsh)  * exp(-((v+75+htsh)/(48 * (1 + htk2)))^2))/qt
  
		v = v - 10
	UNITSON
}
