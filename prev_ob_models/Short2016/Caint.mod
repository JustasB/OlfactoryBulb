TITLE Mechanism  internal calcium concentration cai

NEURON {
	SUFFIX Cacon
	USEION ca READ ica, cai WRITE cai
	GLOBAL tauca, A, camin
}

UNITS { 
 (mA) = (milliamp) 
 (mV) = (millivolt) 
 (molar) = (1/liter) 
 (mM) = (millimolar)
 (uM) = (micromolar) 
} 

PARAMETER {  
 dt (ms) 
 tauca = 800 (ms) 
: A = 1.03e-7 (mM cm2 / ms mA) : same as in paper 
 A = 0.2  :0.103
 camin = 1e-8 (mM)  : arbitrary

} 

STATE {
	cai		(mM) 
}

INITIAL {
	cai = camin
}

ASSIGNED { 
 ica (mA/cm2)    
} 

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	 cai' = -A*ica - (cai-camin)/tauca 
}
