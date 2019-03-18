TITLE LCa calcium channel with fixed reversal potential
: LCa channel with parameters from US Bhalla and JM Bower,
: J. Neurophysiol. 69:1948-1983 (1993)
: Adapted from /usr/local/neuron/demo/release/nachan.mod - squid
: by Andrew Davison, The Babraham Institute.
: 25-08-98

NEURON {
	SUFFIX LCa3_mit_usb
	USEION ca WRITE ica
	RANGE gcabar, ica
	GLOBAL sinf, rinf, stau, rtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

CONSTANT { eca = 70 (mV) }

PARAMETER {
	v (mV)
	dt (ms)
	gcabar	= 0.120 (mho/cm2) <0,1e9>
:	eca = 70 (mV)
}

STATE {
	r s
}

ASSIGNED {
	ica (mA/cm2)
	sinf
	rinf
	stau (ms)
	rtau (ms)
}

INITIAL {
	rates(v)
	s = sinf
	r = rinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar*s*r*(v - eca)
}

DERIVATIVE states {
	rates(v)
	s' = (sinf - s)/stau
	r' = (rinf - r)/rtau
}

FUNCTION alp(v(mV),i) (/ms) {
	if (i==0) {
		alp = 7.5(/ms)/(1 + exp((-v *1(/mV) + 13)/7))
	}else if (i==1){
		alp = 0.0068(/ms)/(1 + exp((v *1(/mV) + 30)/12))
	}
}

FUNCTION bet(v(mV),i)(/ms) {
	if (i==0) {
		bet = 1.65(/ms)/(1 + exp((v *1(/mV) - 14)/4))
	}else if (i==1){
		bet = 0.06(/ms)/(1 + exp(-v* 1(/mV)/11))
	}
}

PROCEDURE rates(v(mV)) {LOCAL a, b
	TABLE sinf, rinf, stau, rtau FROM -100 TO 100 WITH 200
	a = alp(v,0)  b=bet(v,0)
	stau = 1/(a + b)
	sinf = a/(a + b)
	a = alp(v,1)  b=bet(v,1)
	rtau = 1/(a + b)
	rinf = a/(a + b)
}

