
COMMENT

na.mod

Sodium channel, Hodgkin-Huxley style kinetics.  

use with kd.mod

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar, vshift
	GLOBAL mslp, mcen, ma, mc, mq1, mq2
	GLOBAL hslp, hcen, ha, hc, hq1, hq2
	GLOBAL minf, hinf, mtau, htau, ina
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gbar = 258.272   	(pS/um2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)
								
	mslp = 6.7513		(mV)
	mcen = -31.235		(mV)
	ma = 0.23118		(ms)
	mc = 19.098		(mV)
	mq1 = 1		(mV)
	mq2 = 54.927		(mV)

	hslp  = 2.631		(mV)		: v 1/2 for inact 	
	hcen  = -47.953		(mV)		: v 1/2 for inact 	
	ha   = 14.042		(ms)		: inact (v)	
	hc = -73.517		(mV)		: inact inf slope	
	hq1  = 24.053		(mV)		: inact inf slope
	hq2   = 19.627		(mV)	        : inact tau slope

	temp = 23	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states  METHOD cnexp
        gna = gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 

DERIVATIVE states {   :Computes state variables m, h, and n 
        trates(v+vshift)      :             at the current v and dt.
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}


PROCEDURE trates(v (mV)) {  
                      
        TABLE minf, mtau , hinf, htau
	DEPEND dt, mcen, mslp, ma, mc, mq1, mq2, hcen, hslp, ha, hc, hq1, hq2
	
	FROM vmin TO vmax WITH 199

UNITSOFF
	rates(v): not consistently executed from here if usetable == 1
UNITSON

}

UNITSOFF

PROCEDURE rates(vm) {  
        LOCAL  a, b

	mtau = 	xtau(vm, ma, mc, mq1, mq2)
	minf =  xinf(vm, mcen, mslp)

		:"h" inactivation 

	htau = xtau(vm, ha, hc, hq1, hq2)
	hinf = xinf(-vm, -hcen, hslp)
}


FUNCTION xinf(v, xcen, xslp) {
	xinf = 1/( 1 + exp(-(v-xcen)/xslp) )
}

FUNCTION xtau(v, xa, xc, xq1, xq2) {
	xtau = xa / (exp(-(v-xc)/xq2) + exp((v-xc)/xq1))
}

UNITSON
