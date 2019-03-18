:A-Type Potassium current for intrinsic excitability

NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
	RANGE gk, gkbar, gcomp, glearn, i
	RANGE Tau, Te, Tp, eps, float, r, K
	RANGE z, e, p, thresh, delay
}

UNITS {
	(S) = (siemens)
	(uS) = (microsiemens)
	(mV)= (millivolt)
	(mA)= (milliamp)
}

PARAMETER {
	Tau = 20.0 (ms) <1e-9,1e9>
	Te = 200.0 (ms) <1e-9,1e9>
	Tp = 1000.0 (ms) <1e-9,1e9>
	eps = 1e-6
	r = 0.8
	gkbar = 54.8 (uS/cm2)
	thresh = -20 (mV)
	delay = 7
	glearn = 0
	K = 1 <1e-9,1e9>
}

ASSIGNED {
	v  (mV)
	ek (mV)
	ik (mA)
	i  (mA)
	gk  (S/cm2)
	gcomp
	firing
	up
	time
	counter
	ready
}

STATE {
	z
	e
	p
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcomp = g_comp(p)
	if (K<eps) {
		gk = gkbar * glearn
	} else {
		gk = gkbar * gcomp
	}
	i = 1e-6*gk*(v-ek)
	ik = i
}

DERIVATIVE states {
	detect(v)
	if (firing == 1) {z =  z + r*(1-z)}
	z' = -z/Tau
	e' = (z-e)/Te
	p' = (e-p)/Tp
}

INITIAL {
	z = 0.01
	e = 0.01
	p = 0.01
	up = 0
	firing = 0
	time = 0
	counter = 0
	ready = 0	
}

FUNCTION g_comp(p) {
	if (p < eps) {gcomp = 1}
	else {g_comp = log(p)/log(eps)}
}

PROCEDURE detect(v (mV)) {
	if ( v>thresh && up==0 ) {
		counter = delay
		up = 1
		ready = 1
	}
	if( ready==1 && counter>0) {counter = counter-1}
	if( ready==1 && counter<=0) {
		firing = 1
		ready = 0
		time = t
		}
	if ( t>time ) { firing = 0 }	
	if ( v < thresh) { up = 0 }
}
