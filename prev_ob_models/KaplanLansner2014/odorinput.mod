NEURON {
	SUFFIX odorinput
	NONSPECIFIC_CURRENT i
	RANGE i, e, gor, tstart,or
}

PARAMETER {
	or = 0.0 <0,1>
	gor = 0.001 (siemens/cm2) <0, 1e9>
	e = 0 (millivolt)
	tstart = 0  (ms) :will be overwritten by cell constructor parameters
	tau = 20 (ms) 
	tstop = 500 (ms) :will be overwritten by cell constructor parameters
	: good value
	:tstop = tstart + 25 * tau(ms) :will be overwritten by cell constructor parameters
	: these values give t_rise = 87 ms, t_stim =372 ms
}

ASSIGNED {
	i (milliamp/cm2)
	v (millivolt)
}

BREAKPOINT {
	i = 0
		
	if (t >= tstart){ 	
			i = (1 / (1 + exp(-(t - 12*tau) / tau))) * or * gor * (v - e)
	}
	if (t > tstop){ : stimulus decrease
        i = (1 / (1 + exp(-(tstop + 10 * tau - t) / tau))) * or * gor * (v - e)
	}
	: else i = 0
}

