# -*- coding: utf-8 -*-
"""
Simple quadratic integrate-and-fire (QIF) model for the granule cells.
Parameters were fitted to match model f-I curve with granule f-I curve of Davison (2001, PhD Thesis)
"""

from brian import *

print "Granule equation initialization"

# Granule membrane parameters
V_T=-60.*mV
I_T=0.02*nA # threshold of the Davison f-I curve

# All next 3 variables are linked to reproduce the Davison f-I curve (# example with a multiplicative factor)
Delta_T=0.1*mV # Delta_T*factor
tau_m=60*ms # tau_m/factor
gL_G=16.666*nS # gL_G*factor => not a real leak constant (because such a constant has no meaning in a QIF model)


# Synaptic current parameters
tauE=3.*ms
Ee = 0.*mV
gE_max=4.*nS # gE_max*factor <-- if excitation must remain the same when a factor is applied above
 
QIF_eqs="""
dV/dt =  ((1/(2*Delta_T))*(V-V_T)**2 - I_T/gL_G + Iinj/gL_G - (sE*gE_max + gEinj)/gL_G*(V-Ee))/tau_m  : volt
Iinj : amp
gEinj : siemens
dsE/dt = -sE/tauE : 1
"""

def QIF_reset(P, spikes):
    P.V[spikes]=-70.*mV # -70*mV in UB4.c