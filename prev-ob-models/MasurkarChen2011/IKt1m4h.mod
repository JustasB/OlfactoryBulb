TITLE K-A
: K-A current for Mitral Cells from Wang et al (1996)
: M.Migliore Jan. 2002

NEURON {
	SUFFIX Ikt1m4h
	USEION k READ ek WRITE ik
	RANGE  gbar, ik
	GLOBAL minf, mtau, hinf, htau
}

PARAMETER {
	gbar = 0.002   	(mho/cm2)	
								
	celsius
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0m=0.04
	vhalfm=-31.5
	zetam=0.17
	gmm=0.35

	a0h=0.018
	vhalfh=-16
	zetah=0.6
	gmh=0.8

	sha=9.9
	shi=5.7
	
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
	minf 		mtau (ms)	 	
	hinf 		htau (ms)	 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*m*m*m*h*(v + 98)
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(v) {  
	LOCAL qt
        qt=q10^((celsius-24)/10)
        minf = 1/(1+exp(-(v+39.1)/18.56))
	mtau=betm(v)/(0.54*(3+alpm(v)))+0.23
        hinf = 1/(1+exp((v+39.57)/3.68))
	htau = beth(v)/(0.2*(0.000001+alph(v)))+5
}

FUNCTION alpm(v(mV)) {
  alpm = exp(zetam*(v-vhalfm)) 
}

FUNCTION betm(v(mV)) {
  betm = exp(zetam*gmm*(v-vhalfm)) 
}

FUNCTION alph(v(mV)) {
  alph = exp(zetah*(v-vhalfh)) 
}

FUNCTION beth(v(mV)) {
  beth = exp(zetah*gmh*(v-vhalfh)) 
}
