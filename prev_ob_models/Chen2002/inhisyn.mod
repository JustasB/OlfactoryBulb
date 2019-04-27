COMMENT
This is a model for a single input pulse of NMDA type. The conductance of 
synaptic input is affected by extracellular Magnesium concentration.

ENDCOMMENT

NEURON {
        POINT_PROCESS InhiSyn
        RANGE  erev, g, i, onset
        RANGE tp, gmaxampa
        RANGE tau1, tau2, gmaxnmda
        NONSPECIFIC_CURRENT i
}

UNITS {
	(mV) = (millivolt)
        (mA) = (milliampere)
        (nA) = (nanoampere)
        (uS) = (microsiemens)
        (umho) = (micromho)
        (mM) = (milli / litre)
}

PARAMETER {
        v               (mV)

        onset = 1      (ms)

	tp = 5	(ms)
        gmaxampa = 0.5  (uS)

	tau1 = 80	(ms)
	tau2 = 10	(ms)
        gmaxnmda = 0.2  (uS)

        erev = -70        (mV)
}

ASSIGNED {
        i  (nA)
        g  (uS)
}

INITIAL {
UNITSOFF
	syncon(v)
UNITSON
}

BREAKPOINT {
UNITSOFF
        syncon(v)
UNITSON
        i = g * (v - erev)
}

UNITSOFF
PROCEDURE syncon(v) { LOCAL gampa, gnmda
	at_time(onset)
        if ( t < onset ) {
                g = 0
        } else {
	gampa = gmaxampa * (exp(1)/tp) * (t - onset) * exp(- (t - onset) / tp)
        gnmda = gmaxnmda * ( exp(- (t - onset) / tau1) - exp(- (t - onset) / tau2))
        g = gampa + gnmda
        }
}
UNITSON

