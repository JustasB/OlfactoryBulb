import sys, pickle

from generate_firefiles_odors import *

## generate firefiles i.e. list of spike times
## from the firerates computed by generate_firerates_....py
## USAGE:
## from node000:
## mpiexec -machinefile ~/hostfile -n <num_sines*numtrials+1> ~/Python-2.6.4/bin/python2.6 generate_firefiles_sinusoids.py
## mpiexec -machinefile ~/hostfile -n 201 ~/Python-2.6.4/bin/python2.6 generate_firefiles_sinusoids.py
## typical value for num_sines = 5
## typical value for num of trials = 40
## (depends on number of available processing nodes and number of odorfiles generated)

# time points for the firing rate which is read from a pickled file
sinepulsetime = arange(0,SIN_RUNTIME,FIRINGFILLDT)
len_pulsetime = len(sinepulsetime)
num_sins = len(sine_frequencies)

fn = 'firerates/firerates_sinusoids_seed'+str(stim_rate_seednum)+\
    '_ampl'+str(sine_amplitude)+'_mean'+str(sine_ORN_mean)+'.pickle'
print "Loading rate file:",fn
f = open(fn,'r')
frateResponseList = pickle.load(f)
f.close()

def sinusoid_stimuli(fnum,trialnum):
    ## mitral and PG odor ORNs firing files
    for glomnum in range(NUM_GLOMS):
        frate = frateResponseList[glomnum][fnum]
        mitralfirefilename = '../firefiles/firefiles_sins/firetimes_sin_glom'\
            +str(glomnum)+'_fnum'+str(fnum)+'_trial'+str(trialnum)
        ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
            mitralfirefilename, SIN_RUNTIME, sinepulsetime)
    return True
    
if __name__ == "__main__":
    if mpirank == boss: # boss collates
        numavgs = (mpisize-1)/num_sins
        for avgnum in range(numavgs):
            for fnum in range(num_sins):
                procnum = avgnum*num_sins + fnum + 1
                print 'waiting for process '+str(procnum)+'.'
                ## below: you get a numpy array of 
                ## rows=NUM_GLOMS*MIT_SISTERS and cols=spike times
                ## mitral responses has spike times
                ## we calculate STA to get kernel from spike times.
                ok = mpicomm.recv(source=procnum, tag=0)
                print 'received from process '+str(procnum)+'.'
    else:
        seed([500.0*mpirank]) ##### Seed numpy's random number generator.
        trialnum = (mpirank-1)/num_sins
        fnum = (mpirank-1)%num_sins
        ok = sinusoid_stimuli(fnum,trialnum)
        mpicomm.send( ok, dest=boss, tag=0 )
