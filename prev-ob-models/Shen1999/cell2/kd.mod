
COMMENT

kd.mod

Use with na.mod
	
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar
	GLOBAL ninf, ntau, ik
	GLOBAL nna, nnc, nnq1, nnq2, ncen, nslp
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 74.338   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)

	ncen = -46.104	(mV)
	nslp = 1.4503	(mV)
	nna   = 5.8073	(/ms)		: max act rate
	nnc   = -68.210	(mV)		: max deact rate
	nnq1  = 10.625	(mV)		: v 1/2 for inf
	nnq2  = 71.411	(mV)		: inf slope		

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
	trates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n
	ik = (1e-4) * gk * (v - ek)
} 

LOCAL nexp

DERIVATIVE states {   :Computes state variable n 
        trates(v)      :             at the current v and dt.
	n' = (ninf - n)/ntau
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        TABLE ninf, ntau
	DEPEND celsius, temp, ncen, nslp, nna, nnc, nnq1, nnq2
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable_hh == 1

}


PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

        ntau = nna / ( exp((v-nnc)/nnq1) + exp(-(v-nnc)/nnq2) )
	ninf = 1 / ( 1 + exp(-(v-ncen)/nslp) )
}

