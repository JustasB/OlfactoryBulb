: copied by Hines from Exp2syn and added spike dependent plasticity
COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS FastInhib
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE gmax
	RANGE x, mgid, ggid, srcgid
	GLOBAL ltdinvl, ltpinvl, sighalf, sigslope

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1= 1 (ms) <1e-9,1e9> : was 1
	tau2 = 200 (ms) <1e-9,1e9> : was 200
	gmax = .003 (uS) 
	e = -80	(mV)
	ltdinvl = 250 (ms)		: longer intervals, no change
	ltpinvl = 33.33 (ms)		: shorter interval, LTP
	sighalf = 25 (1)
	sigslope = 3 (1)
	x = 0 (um) : cartesian synapse location
	mgid = -1 : associated mitral gid
	ggid = -1 : associated granule gid
	srcgid = -1 : the gid of the granule detector
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	w (uS)
	total (uS)
}

STATE {
	A
	B
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
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

FUNCTION plast(step(1))(1) {
	plast = 1 - 1/(1 + exp((step - sighalf)/sigslope))
}

FUNCTION norm_weight_to_sig(w) {
    norm_weight_to_sig = floor(0.4999 + log(((-1/(w-1))-1)/exp(-sighalf/sigslope))*sigslope)
}

NET_RECEIVE(weight, s, w, tlast (ms)) {
	INITIAL {
		s = 0
		w = 0
		tlast = -1e9(ms)
	}
	if (t - tlast < ltpinvl) { : LTP
		s = s + 1
		if (s > 2*sighalf) { s = 2*sighalf }
	}else if (t - tlast > ltdinvl) { : no change
	}else{ : LTD
		s = s - 1
		if (s < 0) { s = 0 }
	}
	tlast = t
	w = weight*plast(s)
	A = A + w*factor
	B = B + w*factor
}
