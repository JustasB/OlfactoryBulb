# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
## USAGE: python2.6 fit_odor_pulses.py ../results/odor_pulses/2011_02_10_10_38_odorpulses_SINGLES_JOINTS_PGS_numgloms10.pickle

from scipy import optimize
from pylab import *

Mgconc = 1 # mM
vals = [0.44,0.78,0.96]

def Mgblock(V,eta,gamma):
    """
    Mgblock is (1 - (the V factor in g_NMDA))
    """
    return 1 - ( 1/(1+eta*Mgconc*exp(-gamma*V)) )

def chisqfunc(params):
    chisqsum=0
    chisqarray = []
    for i,V in enumerate([-28e-3,-48e-3,-78e-3]):
        chisq = ( vals[i]-Mgblock(V,params[0],params[1]) )**2
        chisqsum += chisq
        chisqarray.append(chisq)
    print chisq
    return array(chisqarray)
    
if __name__ == "__main__":
    eta0 = 0.33 # /mM
    gamma0 = 60 # /V
    params0 = [eta0,gamma0]
    params = optimize.leastsq( chisqfunc, params0, \
        full_output=1, maxfev=100000, ftol=1e-15, xtol=1e-15 )
    print params[3] # print the status message
    params = params[0]
    print params
