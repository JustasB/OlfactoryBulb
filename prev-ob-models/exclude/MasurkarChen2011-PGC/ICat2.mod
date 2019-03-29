TITLE Sustained calcium current from plateau potential-firing rat olfactory glomerular neurons

: originally implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
: original title - Calcium high-threshold L type current for RD Traub, J Neurophysiol 89:909-921, 2003
: Modified by Arjun Masurkar to resemble HVA current in plateau-firing
: juxtaglomerular neurons

: changes: carev=35
: to change: IV curve should match (peak I at 0mv, activation similar to IT)


INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
	SUFFIX iCat2
	USEION ca WRITE ica
	RANGE  gbar, ica
}


UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}
 

PARAMETER { 
	v  		(mV)  
	gbar = 0.0 	(mho/cm2)
	cai	= 8e-5 (mM)		: eca=? mV
	cao	= 2    (mM)
}
 

STATE {
	m
}

ASSIGNED { 
	ica 		(mA/cm2) 
	alpha beta	(/ms)
	carev	(mV)
	i	(mA/cm2)
}
 


BREAKPOINT { 
	SOLVE states METHOD cnexp
	

:	Will assume m^2 model with no inactivation
	carev = 35

	ica = gbar * m * m * ( v - carev) 
        i = ica		: diagnostic i added to display the current
}
 
INITIAL { 
	settables(v) 
	m = alpha / ( alpha + beta )
	m = 0
}
 
DERIVATIVE states { 
	settables(v) 
	m' = alpha * ( 1 - m ) - beta * m 
}

UNITSOFF 

PROCEDURE settables(v) { LOCAL tmp
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	alpha = 1.11 / ( 1 + exp( - 0.058 * ( v +10 ) ) )
	tmp = v +23.9
	if ( fabs( tmp ) < 1e-6 ) {
		beta  = 0.1 * exp( - tmp / 5 ) 
	}else{
		beta  = 0.02 * tmp / ( exp( tmp / 5 ) - 1 )
	}
}

UNITSON