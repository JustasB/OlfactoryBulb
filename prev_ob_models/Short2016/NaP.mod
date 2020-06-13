TITLE Sodium persistent current

COMMENT
      From Wang 1993 NeuroReport  
ENDCOMMENT

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 

NEURON { 
	SUFFIX INaP
	USEION na READ ena WRITE ina
	RANGE gbar, ina
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	ena  = 45	(mV)
	shm  = 1   (mV)  : 0!!!
} 

ASSIGNED { 
    v   (mV)
	ina	(mA/cm2) 
} 


BREAKPOINT { 
	ina = gbar * minf(v) * ( v - ena )
} 


UNITSOFF
 
FUNCTION minf (v) {
  minf = 1/(1+exp(-(v+51 - shm)/5))

}
UNITSON

