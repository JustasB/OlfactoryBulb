
TITLE Burst of IPSC

: Burst of IPSC

: written by Francois David 2006


COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
Each burst is composed of single IPSCs.
Burst parameters and IPSC parameters can be chosen in a large variety of ranges
ENDCOMMENT
 
NEURON {
	POINT_PROCESS ipsc_gauss
	NONSPECIFIC_CURRENT  ip
	RANGE ip, gp
	RANGE dep, dur, delm, sigi, amp, tau, nipsc
	RANGE noise_seed
}
UNITS {
	(nA) = (nanoamp)
:	(uS) = (microsiemens)
}

PARAMETER {
	dep  = 5 (ms)		 : start      }
	dur  = 80 (ms) <0,1e9>   : duree       }  parameters of the IPSC burst
	delm = 10 (ms)		 : delay mean }
	sigi = 50 (ms)		 : delay std }
	dt   = 0.01 (ms) 
	v (mV)			 : fix by Tom McTavish added 20110504

	amp  = 0.005    (nS)     : conductance 			}
	tau  = 5 	(ms)     : synaptic time constant	} parameter of individual IPSCs
	Erev = -70	(mV)	 : GABAA reversal potential	}


	nipsc = 5		 : nb total of ipsc
	noise_seed = 1
}

PROCEDURE seed1(x) {
	set_seed(x)
}

ASSIGNED { ip (nA)
	   gp (uS)
	   indic1
	   events[20000]
	   amplitude[20000]
	   tevents[200]
	   count
}

LOCAL j,k 

INITIAL {LOCAL indic3
	ip = 0
	gp = 0
	count = 1
	indic1 = dur/dt-1
	seed1(noise_seed)
	FROM j = 0 TO indic1 {
		events[j]    = 0  : by default, no event
		amplitude[j] = 0  : by default all the amplitudes are zero
	}

	FROM k = 0 TO nipsc-1 {
		indic3 = -1
		WHILE (indic3<0 || indic3>dur/(dt)) {
		      tevents[k] = normrand(delm,sigi)
		      indic3 = floor((tevents[k]-dep)/(dt)) : indice of the 
		}
		amplitude[indic3]   = scop_random() : pick a random amplitude factor for the conductance
		
	        events[indic3]      = events[indic3] + amplitude[indic3]  : nb of events per time step
	}
}

BREAKPOINT { LOCAL indic2
	if (t<dep+dur && t >= dep) {
		indic2 = floor((t-dep)/(dt)) : indic qui varie avec le temps
		gp = (-gp/tau)*(dt/2) + amp * (amplitude[indic2]/2) * events[indic2] + gp : breakpoint called twice per dt
	
		ip = gp * (v - Erev)
		
		:printf("duree%g\n",ip)
		if (ip<0 || ip>100) {
			ip = 0			:printf("toto\n")
		}
		if (ip>10) {
			ip = 10
		}
	}
	else {
		gp = 0
	}
}
