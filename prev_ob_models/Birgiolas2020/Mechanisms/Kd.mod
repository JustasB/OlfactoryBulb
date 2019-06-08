TITLE K-DR
: K-DR current for Mitral Cells from Wang et al (1996)
: M.Migliore Jan. 2002
: ik = gbar*m*m*(v - ek) !!!
: Original: :ik = gbar*m*(v - ek)  !!!


NEURON {
	SUFFIX Kd
	USEION k READ ek WRITE ik
	RANGE  gbar
	GLOBAL minf, mtau
}

PARAMETER {
	gbar = 0.002   	(mho/cm2)	
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0m=0.0035
	vhalfm=-50
	zetam=0.055
	gmm=0.5
	celsius (degC)
	q10=3
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf
    mtau (ms)
    qt
}
 

STATE { m}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*m*(v - ek)
} 

INITIAL {
	qt=q10^((celsius-24)/10)
	trates(v)
	m=minf
}

DERIVATIVE states {   
	trates(v)
	m' = (minf-m)/mtau
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp(-(v-21)/10))
	mtau = betm(v)/(qt*a0m*(1+alpm(v)))
	mtau = mtau / qt
}

FUNCTION alpm(v(mV)) {
    alpm = exp(zetam*(v-vhalfm))
}

FUNCTION betm(v(mV)) {
    betm = exp(zetam*gmm*(v-vhalfm))
}
