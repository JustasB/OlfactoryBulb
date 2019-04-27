TITLE HH P/N calcium channel

NEURON { 
 SUFFIX Icapn
 USEION ca WRITE ica 
 RANGE gbar, ica, h, m, g 
 GLOBAL minf, hinf, mtau, htau 
} 
 
UNITS { 
 (mA) = (milliamp) 
 (mV) = (millivolt) 
} 
 
PARAMETER { 
 v (mV) 
 gbar = 0.1 (mho/cm2) <0,1e9> 
 e = 100 (mV) 
} 
 
STATE { 
 m
 h
 g (mho/cm2)
} 
 
ASSIGNED { 
 ica (mA/cm2) 
 minf 
 hinf 
 mtau (ms) 
 htau (ms) 
} 
 
INITIAL { 
 rates(v) 
 m = minf 
 h = hinf 
} 
 
BREAKPOINT { 
 SOLVE states METHOD cnexp 
 g = gbar*m*m*h
 ica = g*(v - e) 
} 
 
DERIVATIVE states { 
 : computes state variables m and h at present v, t 
 rates(v) 
 m' = (minf - m)/mtau 
 h' = (hinf - h)/htau 
} 
 

 PROCEDURE rates(v(mV)) { 
 UNITSOFF
 mtau = 0.4 + 0.7/(exp((-5-v)/15) + exp((-5-v)/(-15))) 
 minf = 1/(1+exp(-10-v)/4)  
 htau = 300 + 100/(exp((-40-v)/9.5) + exp((-40-v)/(-9.5))) 
 hinf = 1/(1+exp((-25-v)/(-2))) 
 UNITSON
} 



