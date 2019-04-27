TITLE Calcium dependent potassium channel
: Implemented in Rubin and Cleland (2006) J Neurophysiology
: Parameters from Bhalla and Bower (1993) J Neurophysiology
: Adapted from /usr/local/neuron/demo/release/nachan.mod - squid
:   by Andrew Davison, The Babraham Institute  [Brain Res Bulletin, 2000]

NEURON {
	SUFFIX kca3
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkbar, ik, Yconcdep, Yvdep
	GLOBAL Yalpha, Ybeta
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gkbar= 0.120 (mho/cm2) <0,1e9>
	ek = -80 (mV)
	Ybeta = 0.05 (/ms)
	cai (mM) := 1e-5 (mM)
}


STATE {
	Y
}

ASSIGNED {
	ik (mA/cm2)
	Yalpha   (/ms)
	Yvdep    
	Yconcdep (/ms)
}

INITIAL {
	rate(v,cai)
	Y = Yalpha/(Yalpha + Ybeta)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gkbar*Y*(v - ek)
}

DERIVATIVE state {
	rate(v,cai)
	Y' = Yalpha*(1-Y) - Ybeta*Y
}

PROCEDURE rate(v(mV),cai(mM)) {
	vdep(v)
	concdep(cai)
	Yalpha = Yvdep*Yconcdep
}

PROCEDURE vdep(v(mV)) {
	TABLE Yvdep FROM -100 TO 100 WITH 100
	Yvdep = exp((v*1(/mV)+70)/27)
}

PROCEDURE concdep(cai(mM)) {
	TABLE Yconcdep FROM 0 TO 0.01 WITH 1000
	if (cai < 0.01) {
		Yconcdep = 500(/ms)*( 0.015-cai*1(/mM) )/( exp((0.015-cai*1(/mM))/0.0013) -1 )
	} else {
		Yconcdep = 500(/ms)*0.005/( exp(0.005/0.0013) -1 )
	}
}
