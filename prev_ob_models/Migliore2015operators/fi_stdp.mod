: Weight adjuster portion based on stdwa_songabbott.mod in ModeDB 64261
: Conductance portion based on Exp2Syn.
: This model is intended for use as the mitral side of the reciprocal
: synapse and as such the pre events come from the granule side ThreshDetect
: instance of the Mitral Granule Reciprocal Synapse (MGRS) with non-negative
: weight (positive delay required),
: and the post events come from the mitral side ThreshDetect instance
: of the MGRS with negative weight (0 delay allowed).

COMMENT
Spike Timing Dependent Weight Adjuster
based on Song and Abbott, 2001.
Andrew Davison, UNIC, CNRS, 2003-2004
ENDCOMMENT

NEURON {
	POINT_PROCESS FastInhibSTDP

	: conductance
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE gmax
	RANGE mgid, ggid, srcgid

	: weight adjuster
	RANGE interval, tlast_pre, tlast_post, M, P
	RANGE deltaw, wmax, aLTP, aLTD
	RANGE wsyn
	GLOBAL tauLTP, tauLTD, on
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	: conductance
	tau1=1 (ms) <1e-9,1e9>
	tau2 = 200 (ms) <1e-9,1e9>
	gmax = .003 (uS) 
	e = -80	(mV)

	: weight adjuster
	tauLTP  = 20	(ms)    : decay time for LTP part ( values from           )
	tauLTD  = 20	(ms)    : decay time for LTD part ( Song and Abbott, 2001 )
	wmax    = 1		: min and max values of synaptic weight
	aLTP    = 0.001		: amplitude of LTP steps
	aLTD    = 0.00106	: amplitude of LTD steps
	on	= 1		: allows learning to be turned on and off globally

	: administrative
	mgid = -1 : associated mitral gid
	ggid = -1 : associated granule gid
	srcgid = -1 : the gid of the granule detector
}

ASSIGNED {
	: conductance
	v (mV)
	i (nA)
	g (uS)
	factor

	: weight adjuster
	interval	(ms)	: since last spike of the other kind
	tlast_pre	(ms)	: time of last presynaptic spike
	tlast_post	(ms)	: time of last postsynaptic spike
	M			: LTD function
	P			: LTP function
	deltaw			: change in weight
	wsyn			: weight of the synapse

}

STATE {
	A
	B
}

INITIAL {
	: conductance
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	: weight adjuster
	interval = 0
	tlast_pre = 0
	tlast_post = 0
	M = 0
	P = 0
	deltaw = 0
	wsyn = 0
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	g = (B - A)*gmax
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE (w) {
	if (w >= 0) {				: this is a pre-synaptic spike
		P = P*exp((tlast_pre-t)/tauLTP) + aLTP
		interval = tlast_post - t	: interval is negative
		tlast_pre = t
		deltaw = wmax * M * exp(interval/tauLTD)
	} else {				: this is a post-synaptic spike
		M = M*exp((tlast_post-t)/tauLTD) - aLTD
		interval = t - tlast_pre	: interval is positive
		tlast_post = t
		deltaw = wmax * P * exp(-interval/tauLTP)
	}
	if (on) {
		wsyn = wsyn + deltaw
		if (wsyn > wmax) {
			wsyn = wmax
		}
		if (wsyn < 0) {
			wsyn = 0
		}
	}
	A = A + wsyn*factor
	B = B + wsyn*factor
}
