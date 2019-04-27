# -*- coding: utf-8 -*-
"""
This script provides equations for the mitral model initially
published in David et al. (2009, PLOS Comp Biology)
"""

from brian import *

print "Mitral equation initialization"

# Mitral membrane and ionic current params
Cm     = 0.01*farad*meter**-2 
gL     = 0.1*siemens*meter**-2 
EL     = -66.5*mV 
gNa    = 500*siemens*meter**-2
gNap   = 1.1*siemens*meter**-2
gKA    = 100.*siemens*meter**-2
gKS    = 310.*siemens*meter**-2
gKF    = 100.*siemens*meter**-2
mhKA   = 0.004
ENa    = 45.*mV
EK     = -75.*mV
taumKF = 2.6*ms

# Synaptic current params
Ee = 0.*mV
Ei = -70.*mV
tauIr = 2.*ms # Weak inhibition synapse rise time
gI_max=0.18*siemens*meter**-2 # Weak inhibition max conductance (conductance time integral depends on rise time)
gI_G=3.0*siemens*meter**-2 # Strong inhibition max conductance 

Mitral_eqs=Equations('''
dV/dt=(-gL*(V-EL)-gNa*mNa*mNa*mNa*(V-ENa)-gKA*mhKA*(V-EK)-gKS*W*X*(V-EK)-gNap/(1.+exp(-(V+51.*mV)/(5.*mV)))*(V-ENa)-gKF*Y*(V-EK) + Iinj - gEinj * (V-Ee) - (gI_cst + sI * gI_max + sI_G * gI_G)  * (V-Ei) )/Cm : volt
dW/dt = ( 1./( exp( -(V+34.*mV)/(6.5*mV) ) + 1.) - W ) / taumKS : 1
dX/dt = 2.*( 1./( exp(  (V+65.*mV)/(6.6*mV) ) + 1.) - X ) / ((200. + 220./( exp(-(V+71.6*mV)/(6.85*mV)) + 1.))*ms) : 1
dY/dt = -Y/taumKF : 1
aNa=0.32 * (-(V+50.*mV))/( exp(-(V+50.*mV)/(4.*mV)) - 1. ) : mV
bNa=0.28 * (V+23.*mV) / (exp((V+23.*mV)/(5.*mV)) - 1. )  : mV
mNa=aNa/(aNa+bNa) : 1
dsI/dt = (rI -sI) /tauI: 1
drI/dt= -rI/tauIr : 1
dsI_G/dt = -sI_G /tauI_G : 1
Iinj : amp*meter**-2
gEinj : siemens*meter**-2
taumKS : ms
tauI : ms
tauI_G : ms
gI_cst : siemens*meter**-2

dsI2/dt = (rI2 -sI2) /tauI: 1
drI2/dt= -rI2/tauIr : 1
LFP=sI2 * gI_max :siemens*meter**-2

Isyn=- (sI * gI_max + sI_G * gI_G)  * (V-Ei) :amp*meter**-2
Isyn_all=(- gEinj * (V-Ee) - (gI_cst + sI * gI_max + sI_G * gI_G)  * (V-Ei)) :amp*meter**-2
Gsyn=(sI * gI_max + sI_G * gI_G):siemens*meter**-2
''')

def Mitral_reset(P, spikes):
    P.V[spikes]=-65.*mV # -70*mV in UB4.c
    P.W[spikes]+=0.03 # mKs, 0.03 in UB4.c
    P.X[spikes]+=0.002 # hKs
    P.Y[spikes]+=0.4 # mKf