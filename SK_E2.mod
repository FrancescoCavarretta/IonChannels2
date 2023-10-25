: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK_E2
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gSK_E2bar, gSK_E2, ik
       RANGE msh, mk, mmin
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gSK_E2bar = .000001 (mho/cm2)
          zTau = 1              (ms)

        msh = 0
        mk = 0
        mmin = 0
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_E2	       (mho/cm2)
         ek           (mV)
         cai          (mM)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           gSK_E2  = gSK_E2bar * z
           ik   =  gSK_E2 * (v - ek)
}

DERIVATIVE states {
        rates(cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
          :if(ca < 1e-7 (mM)){
	  :            ca = ca + 1e-07 (mM)
          :}
          :zInf = 1/(1 + (0.00043 (mM)/ ca)^4.8)

          if(ca < 1e-300) {
            zInf = mmin
          } else {
            zInf = mmin + (1-mmin) * 1 / (1 + exp(-(log(ca) + 7.752 + msh) / (0.208 * (1+mk))))
          }
}

INITIAL {
        rates(cai)
        z = zInf
}
