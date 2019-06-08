TITLE Calcium decay
: as described in Bhalla and Bower, J. Neurophysiol. 69:1948-1983 (1993)
: Andrew Davison, The Babraham Institute, 1998
: partially based on cadecay.mod by Alain Destexhe, Salk Institute 1995.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms) }

NEURON{
	SUFFIX CaPool
	USEION ca READ ica, cai WRITE cai
	RANGE channel_flow, depth, B
	GLOBAL cai, cainf, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
}

CONSTANT {
	FARADAY = 96154 (coul)
}

PARAMETER {
	tau = 10 	    (ms)    : cai decay constant
	cainf = 1e-5	(mM)	: baseline calcium concentration
	ica 	        (mA/cm2): calcium current from other channels
    diam            (um)
}

ASSIGNED {
	channel_flow	(mM/ms)
	B		        (mM cm2/ms/mA)
    depth  	        (um)           : shell within which cai is calculated
}

STATE {
	cai		(mM)
}

INITIAL {
    depth = diam/4.0        : to match Bhalla and Bower 1993 set
					        : depth = diam/4 for each compartment

	B = -(1e4)/(2*FARADAY*depth)
	cai = cainf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	channel_flow = B*ica

    : one way flow in channel
	if (channel_flow < 0.0 ) {
		channel_flow = 0.0
	}

	cai' = channel_flow  - (cai - cainf)/tau
}
	





