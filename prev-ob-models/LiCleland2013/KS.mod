TITLE Slow inactivating current Iks

COMMENT
     Slow inactivating current: Wang 1993 Neuroreport
ENDCOMMENT

:INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX IKs
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
} 

ASSIGNED { 
    v  (mV)
	ik  		(mA/cm2) 
	minf 		
	mtau 		(ms)
	hinf
	htau		(ms)
} 

STATE {
	m h
}

INITIAL { 
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
	:TABLE minf, mtau, hinf, htau FROM -120 TO 40 WITH 641
	minf  = 1 / ( 1 + exp( -(v + 34 - sha) / 6.5 ) ) 
	mtau  = tauM 
	hinf = 1 / ( 1 + exp((v + 65 - shi) / 6.6 ))
	htau = 200 + k*220/( 1 + exp( -(v + 71.6) / 6.85 ))
}
UNITSON

