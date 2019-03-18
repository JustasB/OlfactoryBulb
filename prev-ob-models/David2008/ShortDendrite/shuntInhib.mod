COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS shuntI
	RANGE del, dur, amp, tau, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
	Erev = -65.4  (mV)
}
ASSIGNED { i (nA) }

INITIAL {
	i = 0
}

BREAKPOINT { LOCAL g
	if (t<del+dur && t >= del) {
		g = amp
		i = -g * (v - Erev)
	}
	else {
		i = 0
	}
}
