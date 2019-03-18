
COMMENT

kd.mod

Potassium channel, Hodgkin-Huxley style kinetics
Kinetic rates based on Sah et al. and Hamill et al. (1991)

Use with na.mod

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu
	
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar
	GLOBAL ninf, ntau
	GLOBAL Ra, Rb, tha, qa
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 500   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
								
	tha  = 1.037727	(mV)		: v 1/2 for inf
	qa   = 13.422609	(mV)		: inf slope		
	
	Ra   = 0.17832937	(/ms)		: max act rate
	Rb   = 0.0052501573	(/ms)		: max deact rate	

	dt		(ms)
	celsius		(degC)
	temp = 16	(degC)		: original temp 	
	q10  = 2.3			: temperature sensitivity

	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	rates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp

DERIVATIVE states {   :Computes state variable n 
        rates(v)      :             at the current v and dt.
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	a = trap0(v, tha, Ra, qa)
	b = trap0(-v, -tha, Rb, qa)
        ntau = 1/(a+b)
	ninf = a*ntau
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	


