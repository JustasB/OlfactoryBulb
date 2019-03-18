TITLE Transient calcium current from plateau potential-firing rat olfactory glomerular neurons
:   based on recordings by Arjun Masurkar; code based on -
:
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:   
:   code adapted for Vitko et al 2005 (Perez-Reyes group)


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iCat1m3h
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, i
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
	celsius	= 35	(degC)
	gcabar	= .0008	(mho/cm2)
	cai	= 8e-5 (mM)		: eca=? mV
	cao	= 2	(mM)
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	i	(mA/cm2)
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = 35

:	Will assume m3h model
	ica = gcabar * m*m*m*h * (v-carev)
	i = ica		: diagnostic i added to display the current
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {


	evaluate_fct(v)

	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 

	m_inf = 1.0/(1+exp((-51.44-v)/7.23))
	h_inf = 1.0 / ( 1 + exp((v+73.43)/6.04) )

	tau_m = 2.29 + 1.0 / ( exp((v+68.03)/-27.68) + exp((v+39.08)/2.74) ) 
	if (v < -50)
	{tau_h = exp((v+770)/162.5)}
	else
	{tau_h = 37 + exp((v+27)/-6)}
}

UNITSON
