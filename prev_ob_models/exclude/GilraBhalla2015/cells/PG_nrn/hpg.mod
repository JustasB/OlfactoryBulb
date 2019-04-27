TITLE I-h channel for periglomerular cells from Cadetti and Belluzzi (2001)

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
        eh  		(mV)        
	celsius 	(degC)
	ghbar=.00005 	(mho/cm2)
        vhalft=-65   	(mV)
        a0t=0.00085      	(/ms)
        zetat=0.085    	(1)
        gmt=.5   	(1)
	q10=4.5
}


NEURON {
	SUFFIX Ih_cb : Aditya changed this to match with his MOOSE model
	NONSPECIFIC_CURRENT i
        RANGE ghbar, eh
        GLOBAL linf,ltau
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
        linf      
        ltau
}

INITIAL {
    eh = -30 (mV) : Aditya added this else one has to set it from an external .hoc script
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*l*(v-eh)

}


FUNCTION alpt(v(mV)) {
  alpt = exp(zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/ltau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL qt
        qt=q10^((celsius-30)/10)
        linf = 1/(1+ exp((v+80)/10))
        ltau = bett(v)/(qt*a0t*(1+alpt(v)))
}














