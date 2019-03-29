import sys, pickle

from generate_firefiles_odors import *

## generate firefiles i.e. list of spike times
## from the firerates computed by generate_firerates_....py
## USAGE:
## from node000:
## mpiexec -machinefile ~/hostfile -n <numtrains*numtrials+1> ~/Python-2.6.4/bin/python2.6 generate_firefiles_whitenoise.py
## mpiexec -machinefile ~/hostfile -n 251 ~/Python-2.6.4/bin/python2.6 generate_firefiles_whitenoise.py
## typical value for numtrains = 250
## typical value for num of trials = 1
## (depends on number of available processing nodes and number of odorfiles generated)

# time points for the firing rate which is read from a pickled file
pulsetsteps = arange(0,PULSE_RUNTIME,NOISEDT)
numtpulse = len(pulsetsteps)

f = open('firerates/firerates_whitenoise_seed'+\
    str(stim_rate_seednum)+'_dt'+str(NOISEDT)+'_trains'+str(NUMWHITETRAINS)+'.pickle','r')
frateResponseList = pickle.load(f)
f.close()

def whitenoise_stimuli(trainnum,trialnum):
    ## mitral and PG odor ORNs firing files
    for glomnum in range(NUM_GLOMS):
        frate = frateResponseList[glomnum][trainnum]
        mitralfirefilename = '../firefiles/firefiles_whitenoise/firetimes_whitenoise_glom'\
            +str(glomnum)+'_train'+str(trainnum)+'_trial'+str(trialnum)
        ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
            mitralfirefilename, PULSE_RUNTIME, pulsetsteps)
    return True
    
if __name__ == "__main__":
    if mpirank == boss: # boss collates
        numavgs = (mpisize-1)/NUMWHITETRAINS
        for avgnum in range(numavgs):
            for trainnum in range(NUMWHITETRAINS):
                procnum = avgnum*NUMWHITETRAINS + trainnum + 1
                print 'waiting for process '+str(procnum)+'.'
                ## below: you get a numpy array of 
                ## rows=NUM_GLOMS*MIT_SISTERS and cols=spike times
                ## mitral responses has spike times
                ## we calculate STA to get kernel from spike times.
                ok = mpicomm.recv(source=procnum, tag=0)
                print 'received from process '+str(procnum)+'.'
    else:
        seed([500.0*mpirank]) ##### Seed numpy's random number generator.
        trialnum = (mpirank-1)/NUMWHITETRAINS
        trainnum = (mpirank-1)%NUMWHITETRAINS
        ok = whitenoise_stimuli(trainnum,trialnum)
        mpicomm.send( ok, dest=boss, tag=0 )
