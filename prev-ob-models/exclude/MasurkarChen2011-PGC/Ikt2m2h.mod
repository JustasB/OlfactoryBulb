TITLE K-A
: K-A current for Mitral Cells from Wang et al (1996)
: M.Migliore Jan. 2002

NEURON {
	SUFFIX Ikt2m2h
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
	vhalfm=-17.91
	zetam=0.562
	gmm=0.978

	a0h=0.018
	vhalfh=21.04
	zetah=0.695
	gmh=0.966

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
	ik = gbar*m*m*h*(v + 98)
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
        minf = 1/(1+exp(-(v+9.83)/13.18))
	mtau=betm(v)/(0.2*(23+alpm(v)))+2.45
        hinf = 1/(1+exp((v+40.3)/12.3))
	htau = beth(v)/(0.00818*(0.00000005945+alph(v)))+115.9
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
