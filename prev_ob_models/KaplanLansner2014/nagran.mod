TITLE HH fast sodium channel with FUNCTION_TABLEs
: Hodgkin - Huxley granule cell sodium channel using the data given in
: US Bhalla and JM Bower, J. Neurophysiol. 69:1948-1983 (1993)
: Needs the files tabchannels.dat and tabchannels.hoc
: Andrew Davison, The Babraham Institute, 1998.

NEURON {
	SUFFIX nagrantab
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
	GLOBAL minf, hinf, mtau, htau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gnabar= 0.120 (mho/cm2) <0,1e9>
	ena = 45 (mV)
}
STATE {
	m h
}
ASSIGNED {
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

FUNCTION_TABLE tabminf(v(mV))
FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhinf(v(mV))
FUNCTION_TABLE tabhtau(v(mV)) (ms)

PROCEDURE rates(v(mV)) {
	minf = tabminf(v)
	mtau = tabmtau(v) 
	hinf = tabhinf(v)
	htau = tabhtau(v)
}
