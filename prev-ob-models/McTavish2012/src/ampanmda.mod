TITLE simple NMDA receptors

: Hines combined AMPA and NMDA and spike dependent plasticity

: Modified from the original AMPA.mod, M.Migliore Jan 2003
: A weight of 0.0035 gives a peak conductance of 1nS in 0Mg

COMMENT
-----------------------------------------------------------------------------

	Simple model for glutamate AMPA receptors
	=========================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate
    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used
    to estimate the parameters of the present model; the fit was performed
    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.
    1: 195-230, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that 
    included the release of transmitter in adjacent terminals, its lateral 
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.  

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.


-----------------------------------------------------------------------------
ENDCOMMENT



NEURON {
	POINT_PROCESS AmpaNmda
	RANGE R, g, mg, inmda, iampa, gnmda, gampa
	RANGE x, mgid, ggid, srcgid, gmax
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, Alpha, Beta, E, Rinf, Rtau, ampatau
	GLOBAL gampafactor, nmdafactor
	GLOBAL ltdinvl, ltpinvl, sighalf, sigslope
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur	= 1		(ms)	: transmitter duration (rising phase)
	Alpha	= 0.35		(/ms)	: forward (binding) rate
	Beta	= 0.035		(/ms)	: backward (unbinding) rate
	E	= 0	(mV)		: reversal potential
	mg	= 1    (mM)		: external magnesium concentration
	gmax = 2 (umho)
	gampafactor = 0.001 (1)
	nmdafactor = 0.0035 (1)
	ltdinvl = 250 (ms)		: longer intervals, no change
	ltpinvl = 33.33 (ms)		: shorter interval, LTP
	sighalf = 25 (1)
	sigslope = 3 (1)
	ampatau = 3 (ms)
	x = 0 (um) : cartesian synapse location
	mgid = -1 : associated mitral gid
	ggid = -1 : associated granule gid
	srcgid = -1 : gid of the mitral detector
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: total current = iampa+inmda
	inmda 		(nA)		: current = gnmda*(v - E)
	iampa 		(nA)		: current = gampa*(v - E)
	gnmda 		(umho)		: 
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
}

STATE {Ron Roff
	gampa 		(umho)
}

INITIAL {
	PROTECT Rinf = Alpha / (Alpha + Beta)
	PROTECT Rtau = 1 / (Alpha + Beta)
	synon = 0
	gampa = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	gnmda = mgblock(v)*(Ron + Roff)*gmax*nmdafactor
	inmda = gnmda*(v - E)
	iampa = gampa*(v - E)
	i = iampa + inmda
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
	gampa' = -gampa/ampatau
}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first


FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

FUNCTION plast(step(1))(1) {
	plast = 1 - 1/(1 + exp((step - sighalf)/sigslope))
}

FUNCTION norm_weight_to_sig(w) {
    norm_weight_to_sig = floor(0.4999 + log(((-1/(w-1))-1)/exp(-sighalf/sigslope))*sigslope)
}

NET_RECEIVE(weight, s, w, tlast (ms), r0, t0 (ms)) {
	INITIAL {
		s = 0
		w = 0
		tlast = -1e9 (ms)
		r0 = 0
		t0 = -1e9 (ms)
	}
	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		: plasticity affects this spike. If desired to affect
		: the next spike then put following group after
		: net_send
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
		gampa = gampa + w*gmax*gampafactor
		r0 = r0*exp(-Beta*(t - t0))
		t0 = t
		synon = synon + w
		Ron = Ron + r0
		Roff = Roff - r0
		: come again in Cdur with flag = current value of w+1
		net_send(Cdur, w + 1)
        }else{ : turn off what was added Cdur ago
		r0 = (flag-1)*Rinf + (r0 - (flag-1)*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - (flag-1)
		Ron = Ron - r0
		Roff = Roff + r0
	}
}

