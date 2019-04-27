import sys

sys.path.extend(["..","../simulations","../networks"])

from simset_odor import *
from poisson_utils import *
from stimuliConstants import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

## USAGE: python2.6 generate_firefiles_constantrate.py

if varied == 'mainrate':
    firingrates = varied_mainrate # from stimuliConstants.py
    filebase = '../firefiles/firefiles_whitenoise/'
else:
    firingrates = arange(0.0,20.0,0.5)
    #firingrates = arange(0.0,60.0,3.0)
    filebase = '../firefiles/firefiles_constrate/'

def generate_constrate_files(numtrials):
    seed([100.0]) ##### Seed numpy's random number generator.
    ## Just a Poisson spike generator which outputs spike times
    ## for a series of constant firing rates into individual files.
    for firingrate in firingrates: # firingrate in Hz.
        for trialnum in range(numtrials):
            fn = filebase+'firetimes_constrate'+str(firingrate)+'_trial'+str(trialnum)+'.txt'
            firefile = open(fn,'w')
            for i in range(NUM_ORN_FILES_PER_GLOM):
                ornstimvector = poissonTrain(PULSE_RUNTIME,firingrate,REFRACTORY)
                firefile.write(' '.join([str(t) for t in ornstimvector])+'\n')
            firefile.close()
            print "wrote ",fn

if __name__ == "__main__":
    #generate_constrate_files(NUMWHITETRAINS)
    generate_constrate_files(40)
