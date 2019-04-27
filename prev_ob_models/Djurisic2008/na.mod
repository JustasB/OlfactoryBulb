
COMMENT

na.mod

Sodium channel, Hodgkin-Huxley style kinetics.  


qi is not well constrained by the data, since there are no points
between -80 and -55.  So this was fixed at 5 while the thi1,thi2,Rg,Rd
were optimized using a simplex least square proc

voltage dependencies are shifted approximately +5mV from the best
fit to give higher threshold

use with kd.mod

Author: Upinder S. Bhalla, California Institute of Technology
J. of Neurophysiology, V69, N6, 1993

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar, vshift
	GLOBAL thm1, thm2, qm1, qm2, thi1, thi2, qi, qinf, thinf
	GLOBAL minf, hinf, mtau, htau, ina
	GLOBAL Am1, Am2, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gbar = 258.272   	(pS/um2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)
								
	thm1  = -70.3833	(mV)		: v 1/2 for act		(-42)
	thm2  = -21.8432	(mV)		: v 1/2 for act		(-15)
	Am1   = 0.242621	(/ms)		: open (v)		
	Am2   = 0.819229	(/ms)		: close (v)		
	qm1   = 3.51809		(mV)		: act slope		
	qm2   = 3.9834		(mV)		: act slope		

	thi1  = -39.1689	(mV)		: v 1/2 for inact 	
	thi2  = -38.4483	(mV)		: v 1/2 for inact 	
	qi   = 5.63879		(mV)	        : inact tau slope
	thinf = -48.1801	(mV)		: inact inf slope	
	qinf  = 3.73406		(mV)		: inact inf slope
	Rg   = 0.00422366	(/ms)		: inact (v)	
	Rd   = 0.0802232	(/ms)		: inact recov (v) 

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
	DEPEND dt, celsius, temp, Am1, Am2, Rd, Rg, thm1, thm2, thi1, thi2, qm1, qm2, qi, qinf, thinf
	
	FROM vmin TO vmax WITH 199

UNITSOFF
	rates(v): not consistently executed from here if usetable == 1
UNITSON

}

UNITSOFF

PROCEDURE rates(vm) {  
        LOCAL  a, b

	a = trap0(vm,thm1,Am1,qm1)
	b = trap0(-vm,-thm2,Am2,qm2)
	mtau = 1/(a+b)
	minf = a*mtau

		:"h" inactivation 

	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))
}


FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	

UNITSON
