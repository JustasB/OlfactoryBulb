TITLE HH T-type calcium channel 
 
NEURON { 
 SUFFIX Icat 
 USEION ca WRITE ica 
 RANGE gbar, ica, g, h, m, sha, shi, k_tauH
 GLOBAL minf, hinf, mtau, htau 
} 
 
UNITS { 
 (mA) = (milliamp) 
 (mV) = (millivolt) 
} 
 
PARAMETER { 
 v (mV) 
 gbar = 0.036 (mho/cm2) <0,1e9> 
 e = 100 (mV) 
 sha = 0
 shi = 0
 k_tauH = 1
 
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
 mtau = 1.5 + 3.5/(exp(-(v+30-sha)/15) + exp((v+30-sha)/(15))) 
 minf = 1/(1+exp(-(v+44-sha)/5.5))  
 
 htau = 10 + 40/(exp(-(v+50-shi)/15) + exp((v+50-shi)/(15)))
 hinf = k_tauH*1/(1+exp((v+70-shi)/(4))) 
 UNITSON
} 
