: $Id: fvpre.mod,v 1.10 2003/07/29 23:37:22 billl Exp $
COMMENT
synapse taken from Wang, X.-J. and Buzsaki G. (1996) Gamma oscillations by
synaptic inhibition in a hippocampal interneuronal network.  
J. Neurosci. 16, 6402-6413.
ENDCOMMENT
					       
NEURON {
  POINT_PROCESS gradAMPA
  RANGE gmax, g, i, alpha, beta, thetasyn,e, sigma
  :GLOBAL thetasyn,e
  NONSPECIFIC_CURRENT i
  POINTER vpre
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  gmax  = 1e-3 (uS)
  alpha = 1 (/ms)      : 1 ms  
  beta  = 0.1818 (/ms) :  5.5 ms
  e =  0	  (mV)     :  Reserval potential 
  thetasyn =  0 (mV)   :  Activation threshold
  sigma    = 2  : !!!
}

ASSIGNED { vpre (mV) v (mV) i (nA)  g (uS)}

STATE { s }

INITIAL {
  s =  alpha*F(vpre)/(alpha*F(vpre)+beta)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  g = gmax * s
  i = g*(v - e)
}

DERIVATIVE state {
  s' = alpha*F(vpre)*(1-s) - beta*s
}

FUNCTION F (v1 (mV)) {
  F = 1/(1 + exp(-(v1-thetasyn)/sigma))
}  
