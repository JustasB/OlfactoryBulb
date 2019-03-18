import sys, pickle

from pylab import *
import numpy

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *

#### Generate random trains of firing rates
#### to feed to gloms individually to obtain their kernels.
#### At each time point, the firing rate is chosen from a Gaussian distribution
#### cut off at zero below, and mean as mean firing rate for odors,
#### and variance same as mean.
#### The firing rate as a function of time is fed
#### to a Poisson spike generator in generate_firefiles_randompulses.py .

frateResponseList = []

len_pulsetime = int(PULSE_RUNTIME/NOISEDT)
noisepulsetime = arange(0,PULSE_RUNTIME,NOISEDT)

def firingRateWhiteNoise():
    ## in Hz
    ## array of Gaussian distributed firing rates at each time point
    ## mean = FIRINGMEANA, standard deviation = sqrt(FIRINGMEANA)
    #pulse_steps = normal(loc=FIRINGMEANA,scale=sqrt(FIRINGMEANA),size=len_pulsetime)
    
    ## invert flat real FFT to get the time series:
    ## gaussian white noise with mean FIRINGMEANA and variance 0.1*FIRINGMEANA/NOISEDT
    ## See Dayan & Abbot Chap 1 eqn 1.25 for the 1/dt factor
    pulse_rFFT = [ 0.1*FIRINGMEANA/NOISEDT*exp(1j*uniform(0,2*pi)) \
        for idx in range(len_pulsetime/2+1) ]
    pulse_steps = irfft(pulse_rFFT) + FIRINGMEANA
    
    ## clip firing rates below zero; in-place hence pulse_steps is also the output
    clip(pulse_steps,0,1e6,pulse_steps)
    return array(pulse_steps)

def noise_stimuli():
    ## firing rates to generate Poisson input to mitrals and PGs
    for glomnum in range(NUM_GLOMS):
        frateResponseList.append([])
        for trainnum in range(NUMWHITETRAINS):
            # firing rates
            frate = firingRateWhiteNoise()
            ## important to put within [] or (...,) for extend
            frateResponseList[-1].extend([frate])

def fleshout_frate(frate, xtimes):
    fleshed_frate = []
    for f in frate:
        fleshed_frate.extend([f]*xtimes)
    return array(fleshed_frate)

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    seed([123.0])#[stim_rate_seednum]) ##### Seed numpy's random number generator.

    noise_stimuli()

    filename = 'firerates/firerates_whitenoise_seed'\
        +str(stim_rate_seednum)+'_dt'+str(NOISEDT)+'_trains'+str(NUMWHITETRAINS)+'.pickle'
    fireratefile = open(filename,'w')
    pickle.dump( frateResponseList, fireratefile)
    fireratefile.close()
    print "wrote",filename
    
    figure(facecolor='w')
    title('clipped noisetrain')
    plot(frateResponseList[0][0])

    figure(facecolor='w')
    title('avg psd of (noisetrain-mean)')
    flesh_factor = 1
    avg_fftsq = zeros(len_pulsetime*flesh_factor)
    for trainnum in range(NUMWHITETRAINS):
        fleshed_frate = fleshout_frate(frateResponseList[0][0],flesh_factor)
        avg_fftsq += abs(fft(array(fleshed_frate)-fleshed_frate.mean()))**2.0
    avg_fftsq /= float(NUMWHITETRAINS)
    plot(avg_fftsq**0.5)

    # glom0 & glom1
    figure(facecolor='w')
    title('Glomerulus 0 & 1')
    xlabel('time (s)', fontsize='large')
    ylabel('firing rate (Hz)', fontsize='large')
    bar(noisepulsetime, frateResponseList[0][0], width=NOISEDT, color=(1,0,0))
    bar(noisepulsetime, frateResponseList[1][0], width=NOISEDT, color=(0,1,0))
    
    show()
