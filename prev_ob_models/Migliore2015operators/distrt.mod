TITLE ...just to signal thr crossing
: M.Migliore Mar 2012

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
	thr=-10
}


NEURON {
	SUFFIX th
        GLOBAL thr
		RANGE flag
}

ASSIGNED {
	flag
	vpre
}

INITIAL {
	vpre=v
	flag=0
}


BREAKPOINT {
	if (v>thr && vpre<thr) {flag=1} else {flag=0}
	vpre=v	
}
