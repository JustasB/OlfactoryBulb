
: A passive non-specific cationic current to simulate nicotinic effect

NEURON {
	SUFFIX Inic
	NONSPECIFIC_CURRENT i
	RANGE i, enic, gbar
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gbar = 1e-3 (siemens/cm2) < 0, 1e9 >
	enic = 3.2  (mV)
}

ASSIGNED {
	v (mV)
	i (mA/cm2)
}

BREAKPOINT { 
	i = gbar*(v - enic)
}
