
NEURON {
  POINT_PROCESS OdorInput
  RANGE i,del,dur,f0,f1,r,torn,bias
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  del = 50    (ms)
  dur = 200   (ms)
  torn= 500   (ms)
  f0  = 0.2   (nA)
  f1  = 0.8   (nA)
  r   = 60
  bias = 0    (nA)
  pi = 3.1415926
}

ASSIGNED {
  ival (nA)
  i (nA)
  amp (nA)
  on (1)
}

INITIAL {
  i  = 0
  on = 0
  net_send(del, 1)
}

BEFORE BREAKPOINT {
  if  (on) {
    amp = f0 + 0.5*(f1-f0)*(tanh((t-torn)/(r/3)/(1(ms))-3)+1)
	ival = amp
  } else {
    ival = 0
  }
}

BREAKPOINT {
  i = ival
}

NET_RECEIVE (w) {
  if (flag == 1) {
    if (on == 0) {
      : turn it on
      on = 1
      : prepare to turn it off
      net_send(dur, 1)
    } else {
      : turn it off
      on = 0
    }
  }
}


