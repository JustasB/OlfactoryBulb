TITLE Calcium decay
: Implemented in Rubin and Cleland (2006) J Neurophysiology
: written by Andrew Davison et al (2000) Brain Res Bulletin
: as described in Bhalla and Bower (1993) J Neurophysiology 69:1948-1983 (1993)
: partially based on cadecay.mod by Alain Destexhe, Salk Institute 1995.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms) }

NEURON{
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ica, channel_flow, depth, B
	GLOBAL cai, tau, cainf
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
	:FARADAY = 93149 (coul)		: moles do not appear in units
						: note this value is chosen to fit with
						: Genesis
}

PARAMETER {
	dt (ms)
	depth = 1 	(um)		: shell within which cai is calculated
					: to match Bhalla and Bower 1993 set
					: depth = diam/4 for each compartment
	tau = 1e-5 	(ms)		: cai decay constant
	cainf = 1e-5	(mM)		: baseline calcium concentration
	ica		(mA/cm2)
}

STATE {
	cai		(mM)
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	channel_flow	(mM/ms)
	B		(mM cm2/ms/mA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	B = -(1e4)/(2*FARADAY*depth)
	channel_flow = B*ica
	if (channel_flow <= 0.0 ) { channel_flow = 0.0 }	: one way flow in channel
	cai' = channel_flow  - (cai - cainf)/tau
}
	





