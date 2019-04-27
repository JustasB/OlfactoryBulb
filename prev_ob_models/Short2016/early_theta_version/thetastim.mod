: thetastim.mod derived from
: $Id: netstim.mod 2212 2008-09-08 14:32:26Z hines $
: comments at end

NEURON	{ 
  ARTIFICIAL_CELL ThetaStim
  RANGE interval, number, start, actual_start
  RANGE noise
  RANGE outer_interval, outer_number, outer_start
  RANGE outer_noise
  THREADSAFE : only true if every instance has its own distinct Random
  POINTER donotuse
}

PARAMETER {
	interval	= 25 (ms) <1e-9,1e9>: time between spikes (msec)
	number	= 4 <0,1e9>	: number of spikes (independent of noise)
	start		= 25 (ms)	: start of first spike relative to group
	noise		= 0 <0,1>	: amount of randomness (0.0 - 1.0)
	outer_interval	= 200 (ms) <1e-9,1e9>: time between spikes (msec)
	outer_number	= 5 <0,1e9>	: number of grous of spikes (independent of noise)
	outer_start		= 25 (ms)	: start of first spike group
	outer_noise		= 0 <0,1>	: amount of randomness (0.0 - 1.0)
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
	outer_event (ms)
	outer_on
	outer_ispike
	actual_start : start time for each train of spikes: changed as they occur
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0 : off
	ispike = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (outer_noise < 0) {
		outer_noise = 0
	}
	if (outer_noise > 1) {
		outer_noise = 1
	}
	if ((outer_start >= 0 && outer_number > 0) && (number > 0)) {
		outer_on = 1
		: randomize the first spike group so on average it occurs at
		: outer_start + outer_noise*outer_interval
		outer_event = outer_start + invl(outer_interval) - outer_interval*(1. - outer_noise)
		: but not earlier than 0
		if (outer_event < 0) {
			outer_event = 0
		}
		net_send(outer_event, 13) : flag=13 starts a net_send to send flag=3 to start the spikes
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}
PROCEDURE init_outer_sequence(t(ms)) {
	if (outer_number > 0) {
		outer_on = 1
		outer_event = 0
		outer_ispike = 0
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	if (noise == 0) {
		invl = mean
	}else{
		invl = (1. - noise)*mean + noise*mean*erand()
	}
}
FUNCTION outer_invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	if (outer_noise == 0) {
		outer_invl = mean
	}else{
		outer_invl = (1. - outer_noise)*mean + outer_noise*mean*erand()
	}
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		/* only can be used in main thread */
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		: the old standby. Cannot use if reproducible parallel sim
		: independent of nhost or which host this instance is on
		: is desired, since each instance on this cpu draws from
		: the same stream
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval)
	}
	if (ispike >= number) {
		on = 0
	}
}
PROCEDURE next_outer_invl() {
	if (outer_number > 0) {
		outer_event = outer_invl(outer_interval)
	}
	if (outer_ispike >= outer_number) {
		outer_on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			: but not if a netsend is on the queue
			init_sequence(t)
			: randomize the first spike so on average it occurs at
			: noise*interval (most likely interval is always 0)
			next_invl()
			event = event - interval*(1. - noise)
			net_send(event, 1)
		}else if (w < 0) { : turn off spiking definitively
			on = 0
		}
	}
	if (flag == 3) { : from flag 13 event (start of a group of spikes)
		if (on == 1) { : but ignore if turned off by external event
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		ispike = ispike + 1
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event, 1)
		}
	}
	if (flag == 13) { : from INITIAL
		if (outer_on == 1) { : but ignore if turned off by external event
			init_outer_sequence(t)
			on = 1
			net_send(0, 11)           : setup for group of spikes right away
		}
	}
	if (flag == 11 && outer_on == 1) {
		on = 1
		net_send(start, 3) : starts group of spikes
		outer_ispike = outer_ispike + 1
		: no spikes for outer group so no net_event(t)
		next_outer_invl()
		if (outer_on == 1) {
			net_send(outer_event, 11) : setup for next group
		}
	}
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed)

Parameters;
   noise: 	between 0 (no noise-periodic) and 1 (fully noisy)
   interval: 	mean time between spikes (ms)
   number: 	number of spikes (independent of noise)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode
The intrinsic bursting parameters have been removed since
generators can stimulate other generators to create complicated bursting
patterns with independent statistics (see below)

Modified by Michael Hines to use logical event style with NET_RECEIVE
This stimulator can also be triggered by an input event.
If the stimulator is in the on==0 state (no net_send events on queue)
 and receives a positive weight
event, then the stimulator changes to the on=1 state and goes through
its entire spike sequence before changing to the on=0 state. During
that time it ignores any positive weight events. If, in an on!=0 state,
the stimulator receives a negative weight event, the stimulator will
change to the on==0 state. In the on==0 state, it will ignore any ariving
net_send events. A change to the on==1 state immediately fires the first spike of
its sequence.

Modified by Tom Morse to provide a ThetaStim protocol:
Added modified duplicate "outer" functions and testing to the
NET_RECEIVE section that were allegorical to a NetStim's control of one group of
spikes to add control for groups of spike bursts, i.e. a Theta Stimulation.
The "outer" functions and variables control the start time and number of groups,
and the original NetStim variable names control the spikes within each group.

ENDCOMMENT

