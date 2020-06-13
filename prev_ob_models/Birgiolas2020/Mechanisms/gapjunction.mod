

NEURON {

    POINT_PROCESS GapJunction
    RANGE v_other
    RANGE g, i
    NONSPECIFIC_CURRENT i
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
    g = 0 (uS)
}

ASSIGNED {

    v       (mV)
    v_other (mV)
    i       (nA)
}

BREAKPOINT {

	if (g>0) {
	    i = g * (v-v_other)
	}

}