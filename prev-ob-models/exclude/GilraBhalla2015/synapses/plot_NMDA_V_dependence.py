from pylab import *
from synapseConstants import *

figure()
title("NMDA synapse - Mg dependence")
Vrange = arange(-80e-3,0e-3,1e-3)
eta = 1.0 / mitral_granule_NMDA_KMg_A
gamma = 1.0 / mitral_granule_NMDA_KMg_B
Vdep = [ 1.0/(1+eta*MG_CONC*exp(-gamma*V)) for V in Vrange]
plot(Vrange,Vdep,'r,-')
show()
