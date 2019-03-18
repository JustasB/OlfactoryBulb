COMMENT
-----------------------------------------------------------------------------

GradGABAa.mod
Graded ionotropic GABA-A synaptic mechanism

Simple synaptic mechanism for the thresholded, graded release of GABA
(e.g., in dendodendritic synapses).  

Thomas A. Cleland (tac29@cornell.edu) and Praveen Sethupathy
Cornell University
Summer 2003, January 2004
adapted by Francois David 2005
-----------------------------------------------------------------------------

There is a pointer called "PreActiv" which must be set to the (presynaptic) variable 
that is supposed to trigger synaptic release.  This variable is usually the
presynaptic voltage (mV) but it could be another variable, such as the presynaptic 
calcium concentration.

Once pre has crossed the threshold value given by Prethresh, a quantity of
C (i.e., GABA) is released which increases monotonically with (pre-Prethresh)
until Cmax (the maximum permissible concentration of GABA in the cleft) is
reached.  Postsynaptic conductance (g) is currently scaled directly to C such that 
g=gmax when C=Cmax.  

This mechanism is perfectly adequate for most network simulations; graded
release schema tend to be forgiving of kinetic imprecision.  It will, however,
be found wanting in studies of or critically depending upon GABA receptor kinetics.  

-----------------------------------------------------------------------------
ENDCOMMENT


NEURON {				
	POINT_PROCESS GradGABAa		
	POINTER PreActiv
	RANGE C, g, gmax, Erev, timestep, tau, g_inf, vref, thres, slop
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Erev	 = -70	 (mV)			: reversal potential 
	gmax	 = 0	 (uS)			: maximum conductance
	tau 	 = 3 	 (ms)

	g_inf 	 = 0 	 (umho)
	vref     = 1    (mV)                    : sert juste a l'homogeneite des calculs
	thres    = -45  (mV) 
	slop	 = 0.2			 	: no unit	
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C				: analog of transmitter concentration or synapse activity
	PreActiv			: pointer to presynaptic variable
}

INITIAL {
	C = 0
	g = 0
}

BREAKPOINT { LOCAL delta_g
	  
	C = 1/(1+exp(4*slop*(thres-PreActiv)/vref))  : sigmoid 0.6/4 is the max slope at -45mV
					   	     : values fitted to get the wanted effect

	g_inf = gmax* C	
	delta_g = (dt/2) * (g_inf-g)/tau      : exponential function
					      : breakpoint called twice per time step then dt/2

	g = g + delta_g	

	:i=0

	i = (g*(v-Erev))
	:if ( i > 5 ){
	:	i=5	:to avoid bug
	:} 
	:optional to make the GABA conductnace  only inhibiting
	:if ( i < 0 ) {
	:	i=0	:to avoid bug 
	:}
}



