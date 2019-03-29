import sys, pickle

from generate_firefiles_odors import *

## generate firefiles i.e. list of spike times
## from the firerates computed by generate_firerates_....py
## USAGE: python2.6 generate_firefiles_variedinhibition.py

RUNTIME = REALRUNTIME + SETTLETIME

# binning to plot odor responses
bindt = (RUNTIME-SETTLETIME)/respbins
tlist = arange(SETTLETIME+bindt/2.0,RUNTIME,bindt)
pulsebindt = PULSE_RUNTIME/pulsebins
pulsetlist = arange(pulsebindt/2.0,PULSE_RUNTIME,pulsebindt)

# time points for the firing rate which is read from a pickled file
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
numt = len(firingtsteps)
extratime = arange(0,2*RESPIRATION,FIRINGFILLDT)
pulsetsteps = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
numtpulse = len(pulsetsteps)

f = open('firerates/firerates_2sgm_variedinh.pickle','r')
frateOdorList,kernels = pickle.load(f)
f.close()

def variedinh_stimuli():
    ## mitral and PG odor ORNs firing files
    ## Seed numpy's random number generator.
    ## If no parameter is given, it uses current system time
    seed([700.0])
    for glomnum in range(NUM_GLOMS):
        figure()
        title(str(glomnum))
        ylabel('Hz')
        xlabel('time (s)')
        totinh = float(len(frateOdorList[glomnum]))
        for inhidx,frate in enumerate(frateOdorList[glomnum]):
            for avgnum in range(MAXNUMAVG):
                mitralfirefilename = '../firefiles/firefiles_variedinh/firetimes_2sgm_glom_'\
                    +str(glomnum)+'_inhnum'+str(inhidx)+'_avgnum'+str(avgnum)
                ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
                    mitralfirefilename, RUNTIME, firingtsteps)
            ## plotBins here returns firing rate of NUM_ORN_FILES_PER_GLOM combined.
            ## So divide by NUM_ORN_FILES_PER_GLOM.
            ## We just plot ornstimvector_merged for the last avgnum 
            ratebins = [rate/NUM_ORN_FILES_PER_GLOM for rate in plotBins(ornstimvector_merged,\
                respbins, RUNTIME, SETTLETIME)]
            plot(tlist, ratebins, color=(inhidx/totinh,1-inhidx/totinh,0), marker=',')

if __name__ == "__main__":
    ## seeds are set in each of the below functions:
    ## so you may comment / uncomment below functions depending on what you want to generate.
    variedinh_stimuli()
    show()

