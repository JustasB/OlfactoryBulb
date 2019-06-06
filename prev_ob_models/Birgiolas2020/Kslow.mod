TITLE Slow inactivating current Iks

COMMENT
     From: Wang (1993) Ionic basis for intrinsic 40 Hz neuronal oscillations, Neuroreport
ENDCOMMENT


UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX Kslow
	USEION k READ ek WRITE ik
	RANGE gbar, ik
}

PARAMETER { 
	gbar =  0.0 	(mho/cm2)
	ek   = -70	(mV)
	k    =  1.5       : 1.5 !!!
	sha  =  0   (mV)
	shi  =  -3   (mV) : -3 !!!
	tauM = 10   (ms) : 10
	celsius
	q10 = 3
} 

ASSIGNED { 
    v  (mV)
	ik  		(mA/cm2) 
	minf 		
	mtau 		(ms)
	hinf
	htau		(ms)
	qt
} 

STATE {
	m h
}

INITIAL {
	qt=q10^((celsius-35)/10)
	Rates(v) 
	m = minf
	h = hinf
} 

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gbar * m * h * ( v - ek )
} 


DERIVATIVE states { 
	Rates(v) 
	m' = ( minf - m ) / mtau
	h' = ( hinf - h ) / htau 
}

UNITSOFF
 
PROCEDURE Rates(v) { 
	minf  = 1 / ( 1 + exp( -(v + 34 - sha) / 6.5 ) )
	mtau  = tauM / qt
	hinf = 1 / ( 1 + exp((v + 65 - shi) / 6.6 ))
	htau = 200 + k*220/( 1 + exp( -(v + 71.6) / 6.85 ))
	htau = htau / qt
}
UNITSON

