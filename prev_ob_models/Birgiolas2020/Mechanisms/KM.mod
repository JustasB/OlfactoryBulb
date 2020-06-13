TITLE HH KM channel
: Hodgkin - Huxley KM channel with parameters from
: US Bhalla and JM Bower, J. Neurophysiol. 69:1948-1983 (1993)
: Andrew Davison, The Babraham Institute, 1998.

NEURON {
	SUFFIX KM
	USEION k READ ek WRITE ik
	RANGE gbar, ik
	GLOBAL xinf, xtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gbar=.036 (mho/cm2) <0,1e9>
	ek = -70 (mV)
	celsius (degC)
	q10=3
}
STATE {
	x
}
ASSIGNED {
	ik (mA/cm2)
	xinf
	xtau (ms)
	qt
}

INITIAL {
	qt=q10^((celsius-35)/10)
	rates(v)
	x = xinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gbar*x*(v - ek)
}

DERIVATIVE states {
	rates(v)
	x' = (xinf - x)/xtau
}

PROCEDURE rates(v(mV)) {
	TABLE xinf, xtau FROM -100 TO 100 WITH 200
	xinf = 1/(1 + exp(-(v*1(/mV) + 35)/5))
	xtau = 1000(ms)/(3.3*exp((v*1(/mV) + 35)/40) + exp(-(v*1(/mV) + 35)/20))
	xtau = xtau / qt
}

