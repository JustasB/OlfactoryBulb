COMMENT
-----------------------------------------------------------------------------

hpg.mod	
Mechanism for hyperpolarization-activated cation channel ("H-channel")

Based on a mechanism written by Cadetti & Belluzzi (2001) for H-current
in rat olfactory bulb juxtaglomerular cells (i.e., primarily periglomerular
and external tufted cells), as described in:
  
	Cadetti L, Belluzzi O (2001) Hyperpolarisation-activated current in
	glomerular cells of the rat olfactory bulb.  Neuroreport 12(14):3117-20.

Slightly modified for redistribution by

	Thomas A. Cleland (tac29@cornell.edu) and Praveen Sethupathy
	Cornell University
	Summer 2003, January 2004

-----------------------------------------------------------------------------
ENDCOMMENT

TITLE I-h channel for juxtaglomerular cells 

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v 			(mV)
  	eh  			(mV)        
	celsius 		(degC)
	ghbar  = 5e-05 	(mho/cm2)
      vhalft=-65   	(mV)
      a0t=0.00085      	(/ms)
      zetat=0.085    	(1)
      gmt=.5   		(1)
	q10=4.5
}

NEURON {
	SUFFIX hpg
	NONSPECIFIC_CURRENT i
      RANGE ghbar, eh, g
      GLOBAL linf,taul
}

STATE { l }

ASSIGNED {
	i (mA/cm2)
      linf      
      taul
	g
}

INITIAL {
	rate(v)
	l=linf
	g=10
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = ghbar*l
	i = g*(v-eh)
}

FUNCTION alpt(v(mV)) {
  	alpt = exp(zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  	bett = exp(zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
      rate(v)
      l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
      LOCAL qt
      qt=q10^((celsius-30)/10)
      linf = 1/(1+ exp((v+80)/10))
      taul = bett(v)/(qt*a0t*(1+alpt(v)))
}














