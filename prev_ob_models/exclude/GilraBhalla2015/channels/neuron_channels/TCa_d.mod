TITLE Low threshold calcium current
:
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by standard equations (NOT GHK)
:   using a m2h format, according to the voltage-clamp data
:   (whole cell patch clamp) of Huguenard & Prince, J Neurosci.
:   12: 3804-3817, 1992.
:   See also: http://www.cnl.salk.edu/~alain ,  http://cns.fmed.ulaval.ca
:
:    - Kinetics adapted to fit the T-channel of reticular neuron
:    - Q10 changed to 5 and 3
:    - Time constant htau fitted from experimental data
:    - shift parameter for screening charge
:
:   ACTIVATION FUNCTIONS FROM EXPERIMENTS (NO CORRECTION)
:
:   Reversal potential taken from Nernst Equation
:
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX TCa_d
	USEION ca READ cai, cao WRITE ica
	RANGE gmax, shift : aditya changed gcabar to gmax
	GLOBAL minf, mtau, hinf, htau : aditya shifted these from RANGE to GLOBAL for ease in plotting them
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC) : aditya - this will be overridden by external global celsius anyway
:	eca	= 120	(mV)
	gmax	= .0008	(mho/cm2)
	shift	= 0 	(mV)
	cai	: aditya put the default values in INITIAL block
	cao	: aditya put the default values in INITIAL block
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	minf
	mtau	(ms)
	hinf
	htau	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD euler
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gmax * m*m*h * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

UNITSOFF
INITIAL {
:
:   Activation functions and kinetics were obtained from
:   Huguenard & Prince, and were at 23-25 deg.
:   Transformation to 36 deg assuming Q10 of 5 and 3 for m and h
:   (as in Coulter et al., J Physiol 414: 587, 1989)
:
	cai	= 2.4e-4 (mM)		: adjusted for eca=120 mV : aditya put them in INITIAL block so that they don't need to be supplied externally
	cao	= 2	(mM)    : aditya put them in INITIAL block so that they don't need to be supplied externally

	phi_m = 5.0 ^ ((celsius-24)/10)
	phi_h = 3.0 ^ ((celsius-24)/10)

	evaluate_fct(v)

	m = minf
	h = hinf
}

PROCEDURE evaluate_fct(v(mV)) { 
:
:   Time constants were obtained from J. Huguenard
:

	minf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
	hinf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

	mtau = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
	htau = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h
}
UNITSON
