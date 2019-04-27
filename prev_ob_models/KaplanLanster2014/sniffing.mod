NEURON {
	SUFFIX sniffinginput 
	NONSPECIFIC_CURRENT i
	RANGE i, e, gor, tstart, or, tau, sniffperiod, tshift
}

PARAMETER {
	or = 0.0 <0,1>
	gor = 0.001 (siemens/cm2) <0, 1e9>
	e = 0 (millivolt)
	tstart = 0  (ms) :will be overwritten by cell constructor parameters
	tau = 20 (ms) 
	tstop = 500 (ms) :will be overwritten by cell constructor parameters
	sniffperiod = 80 (ms)
	tshift = 40 (ms)
}

ASSIGNED {
	i (milliamp/cm2)
	v (millivolt)
}

BREAKPOINT {
	i = 0
		
	if (t>=tstart && t<=tstop){ 	: sniffing input
		: i = (sin(t / sniffperiod - tshift))^2 * or * gor
		i = (sin(t / sniffperiod - tshift))^2 * or * gor * (v - e)
	}
	: else i = 0
}

