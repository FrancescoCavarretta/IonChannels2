:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)

NEURON	{
	SUFFIX SKv3_1
	USEION k READ ek WRITE ik
	RANGE gSKv3_1bar, gSKv3_1, ik
        RANGE msh, mk, mmin

        RANGE mtmin, mtmax, mtsh, mtk
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gSKv3_1bar = 0.00001 (S/cm2)
        msh = 0
        mk = 0
        mmin = 0

        mtmin = 0
        mtmax = 0
        mtsh = 0
        mtk = 0
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gSKv3_1	(S/cm2)
	mInf	(1)
	mTau	(ms)
}

STATE	{ 
	m      (1)
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gSKv3_1 = gSKv3_1bar*m
	ik = gSKv3_1*(v-ek)
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
	UNITSOFF
		mInf =  mmin + (1-mmin) * 1/(1+exp(((v - 18.700 + msh)/(-9.700 * (1+mk) ))))
        
		mTau =  mtmin + (1 + mtmax) * (0.2*20.000) / (1 + exp(((v + 46.560 + mtsh)/(-44.140 * (1 + mtk)))))
	UNITSON
}
