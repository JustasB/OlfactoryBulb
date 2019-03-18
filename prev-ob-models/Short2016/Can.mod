TITLE Calcium dependent non-specific cation channel 
 
NEURON { 
 SUFFIX Ican 
 NONSPECIFIC_CURRENT i 
 USEION ca READ cai
 RANGE gbar, i, minf, tauinf, m, cai, g, MCa, mm, Shm 
 GLOBAL tauca
} 
 
UNITS { 
 (mA) = (milliamp) 
 (mV) = (millivolt) 
 (molar) = (1/liter) 
 (mM) = (millimolar)
 (uM) = (micromolar) 
} 
 
PARAMETER { 
 v (mV) 
 dt (ms) 
 gbar= 0.003 (mho/cm2) <0,1e9> 
 e = 10 (mV) 
 tauca = 800 (ms) 
 kdcan = .200 (mM)
 Shm = -2
} 
 
STATE { 
 m 
 g (mho/cm2)
 MCa 
 mm
} 
 
ASSIGNED { 
 cai (mM)
 i (mA/cm2)
 minf   
 tauinf (ms)    
} 
 
INITIAL { 
 rate(v) 
 m = minf
} 
 
BREAKPOINT { 
 SOLVE state METHOD cnexp
 MCa =  (cai/(kdcan + cai)) : different ways to read out all the different gating variables
 mm = m*MCa
 g = gbar*mm
 i = g*(v - e) 
} 
 
DERIVATIVE state { 
 rate(v)    
 m' = (minf -m)/tauinf
} 
 
PROCEDURE rate(v(mV)) { 
UNITSOFF
minf   = 1/(1+exp(-(43+v-Shm)/5.2))
tauinf = 1.6 + 2.7 /(exp((-55-v)/15) + exp((-55-v)/(-15))) 
UNITSON
} 
 
