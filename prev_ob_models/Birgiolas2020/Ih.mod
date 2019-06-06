TITLE Hyperpolarization-activated cation current
: Implemented in Rubin and Cleland (2006) J Neurophysiology
: Adapted from Saraga et al (2003) J. Physiology

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX Ih
    USEION h READ eh WRITE ih VALENCE 1
    RANGE gbar,ih
    GLOBAL rinf, rexp, tau_r
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    p = 5 (degC)
    dt (ms)
    gbar = 0.00 (mho/cm2)
    eh = -32.9 (mV)
    q10=4.5
    celsius
}
 
STATE {
    r
}
 
ASSIGNED {
    ih (mA/cm2)
	rinf rexp
	tau_r
	qt
}
 
BREAKPOINT {
    SOLVE deriv METHOD derivimplicit
    ih = gbar*r*(v - eh)
}
 
INITIAL {
    qt=q10^((celsius-35)/10)
	rates(v)
	r = rinf
}

DERIVATIVE deriv { :Computes state variable h at current v and dt.
	rates(v)
	r' = (rinf - r)/tau_r
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    TABLE rinf, rexp, tau_r DEPEND dt, p FROM -200 TO 100 WITH 300

	rinf = 1/(1 + exp((v+84.1)/10.2))
	rexp = 1 - exp(-dt/(tau_r))
	tau_r = 100 + 1/(exp(-17.9-0.116*v)+exp(-1.84+0.09*v))
	tau_r = tau_r / qt
}
 
UNITSON

