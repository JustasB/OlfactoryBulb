: $Id: monx.mod,v 1.3 2007/03/09 03:45:56 ted Exp ted $

COMMENT
Does 3 things:
1. Has storage that can be used for various things, 
   e.g. distance from soma.
2. At initialization, stores local area in a variable
3. Monitors peak v during a simulation and stores in another variable.


CAVEATS
=======
Beware of one drawback of this approach to storing values of range variables:
attempts to access values at 0 or 1 will instead return or set the values at 
the first or last internal nodes!  For example, in a section with nseg = 1,
  for (x) dist_monx(x) = distance(x)
will result in dist_monx(0.5) actually storing the value for distance(1), 
and subsequently invoking
  for (x) print dist_monx(x)
will print the same number 3 times!


If monx is to track currents, a custom initialization is needed 
so that the initialization reaches a stable point.

According to Michael Hines, 
"The problem is that ina is not set until just before return from 
finitialize when the equivalent of an fcurrent() call is executed 
(cvode will also end up calling the statements in the BREAKPOINT's 
SOLVE block)."

If one is using the GUI, the simple solution would be to just do 
an Init followed by Init & Run, since that is how many times 
finitialize() must be called to reach a fixed point.

A solution that is less prone to user error is to replace the 
run time system's init() with a custom init() that executes a 
pair of finitialize() calls.

For example, this will do:

proc init() {
	finitialize(v_init)
	finitialize(v_init)
	if (cvode.active()) {
		cvode.re_init()
	} else {
		fcurrent()
	}
	frecord_init()
}

ENDCOMMENT

NEURON {
	SUFFIX monx
:	USEION na READ ina
:	RANGE inamin, tinamin, inamax, tinamax, qna
	RANGE vmax, dist, surf, temp, zin, gesyn, gisyn
}

ASSIGNED {
:	ina (milliamp/cm2)
:	inamin (milliamp/cm2)
:	inamax (milliamp/cm2)
:	tinamin (ms)
:	tinamax (ms)
	area (micron2)
	surf (micron2)
	v (millivolt)
	vmax (millivolt)
	tvmax (ms)
	dist (micron)
	gesyn (microsiemens)
	gisyn (microsiemens)
	zin (megohm)
	temp (1) : for temporary/intermediate variables!
}

: STATE { qna (picocoulomb) }

INITIAL {
:	inamin = ina
:	inamax = ina
: printf("ina = %f\n", ina) 
:	qna = 0
	vmax = 0
	tvmax = 0
	surf = area
	: dist must be initialized at the hoc level
}

: BREAKPOINT {
:	SOLVE integrate METHOD cnexp
:	
: }

: DERIVATIVE integrate {
: 	qna' = (0.01)*ina*area
: }

AFTER SOLVE {
:	if (ina<inamin) {
:		inamin = ina
:		tinamin = t
:	}
:	if (ina>inamax) {
:		inamax = ina
:		tinamax = t
:	}
	if (v>vmax) {
		vmax = v
		tvmax = t
	}
}
