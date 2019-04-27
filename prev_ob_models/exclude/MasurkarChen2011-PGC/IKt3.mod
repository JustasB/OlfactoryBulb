
COMMENT

kd.mod

Potassium channel, Hodgkin-Huxley style kinetics
Kinetic rates based on Sah et al. and Hamill et al. (1991)

Use with na.mod

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu
	
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Ikt3	USEION k READ ek WRITE ik
	RANGE n, gk, gbar, actshift, i
	GLOBAL ninf, ntau, ik
	GLOBAL Ra, Rb, tha, qa
	GLOBAL tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 41.8448   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
								
	tha  = 40	(mV)		: v 1/2 for inf
	qa   = 16.7	(mV)		: inf slope		
	
	Ra   = 0.6	(/ms)		: max act rate
	Rb   = 0.012	(/ms)		: max deact rate	

	dt		(ms)
	celsius		(degC)
	

	vmin = -120	(mV)
	vmax = 60	(mV)
  actshift = 0 (mV)} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
        i               (mA/cm2)
}
 

STATE { n }

INITIAL { 
	trates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = gbar*n*n*n*n
	ik = (1e-4) * gk * (v +  98)
	i = ik

} 

LOCAL nexp

DERIVATIVE states {   :Computes state variable n 
        trates(v)      :             at the current v and dt.
	n' = (ninf - n)/ntau
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        TABLE ninf, ntau
	DEPEND celsius, Ra, Rb, tha, qa
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable_hh == 1

}


PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)
	ninf = a*ntau
}


