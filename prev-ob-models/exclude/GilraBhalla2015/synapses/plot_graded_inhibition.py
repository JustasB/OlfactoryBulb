from pylab import *

figure()
title("Graded inhibitory synapse: granule-|mitral")
# adjust activation curve so that all the dynamics happens between -52 and -28mV.
act = [0.0] # at -52mV
act.extend( [1/(1+math.exp(-(vm + 40.5e-3)/1.0e-3)) for vm in arange(-50e-3,-30.00001e-3,2e-3)] )
act.extend([1.0]) # for -28mV
plot(arange(-52e-3,-28.00001e-3,2e-3),act,'r,-')
show()
