TITLE HH slow potassium channel with FUNCTION_TABLEs  (.inf and .tau files)
: Implemented in Rubin and Cleland (2006) J Neurophysiology
: Parameters from Bhalla and Bower (1993) J Neurophysiology
: Adapted from /usr/local/neuron/demo/release/nachan.mod - squid
:   by Andrew Davison, The Babraham Institute  [Brain Res Bulletin, 2000]

NEURON {
	SUFFIX kslowtab
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL ninf, kinf, ntau, ktau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar= 0.120 (mho/cm2) <0,1e9>
	ek = -70 (mV)
}
STATE {
	n k
}
ASSIGNED {
	ik (mA/cm2)
	ninf
	kinf
	ntau (ms)
	ktau (ms)
}

INITIAL {
	rates(v)
	n = ninf
	k = kinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*n*n*k*(v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
	k' = (kinf - k)/ktau
}

FUNCTION_TABLE tabninf(v(mV))
FUNCTION_TABLE tabntau(v(mV)) (ms)
FUNCTION_TABLE tabkinf(v(mV))
FUNCTION_TABLE tabktau(v(mV)) (ms)

PROCEDURE rates(v(mV)) {
	ninf = tabninf(v)
	ntau = tabntau(v) 
	kinf = tabkinf(v)
	ktau = tabktau(v)
}
