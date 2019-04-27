import sys, pickle

from pylab import *
import numpy

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *

#### Generate sinusoids of different frequencies and amplitudes
#### to feed to gloms individually to obtain frequency response.
#### Cut off at zero below, and stimuliConstants has mean firing rate for odors,
#### keep amplitude less than mean firing rate to avoid cutoff,
#### and variance in firing rate of different ORNs is same as mean.
#### The firing rate as a function of time is fed
#### to a Poisson spike generator in generate_firefiles_sinusoids.py .

frateResponseList = []

sinepulsetime = arange(0,SIN_RUNTIME,FIRINGFILLDT)
len_pulsetime = len(sinepulsetime)

def firingRateSinusoid(DC,ampl,f):
    ## in Hz
    ## array of Gaussian distributed firing rates at each time point
    ## mean = FIRINGMEANA, standard deviation = sqrt(FIRINGMEANA)
    #pulse_steps = normal(loc=FIRINGMEANA,scale=sqrt(FIRINGMEANA),size=len_pulsetime)
    
    pulse_steps = array( [ DC + ampl*sin(2*pi*t*f) for t in sinepulsetime ] )
    
    ## clip firing rates below zero; in-place hence pulse_steps is also the output
    clip(pulse_steps,0,1e6,pulse_steps)
    return array(pulse_steps)

def sinusoid_stimuli():
    ## firing rates to generate Poisson input to mitrals and PGs
    for glomnum in range(NUM_GLOMS):
        frateResponseList.append([])
        for sine_f in sine_frequencies:
            ## firing rates
            frate = firingRateSinusoid(sine_ORN_mean,sine_amplitude,sine_f)
            ## important to put within [] or (...,) for extend
            frateResponseList[-1].extend([frate])

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    seed([123.0])#[stim_rate_seednum]) ##### Seed numpy's random number generator.

    sinusoid_stimuli()

    filename = 'firerates/firerates_sinusoids_seed'+str(stim_rate_seednum)+\
        '_ampl'+str(sine_amplitude)+'_mean'+str(sine_ORN_mean)+'.pickle'
    fireratefile = open(filename,'w')
    pickle.dump( frateResponseList, fireratefile)
    fireratefile.close()
    print "wrote",filename
    
    figure(facecolor='w')
    title('psd of sinusoid')
    frate = frateResponseList[0][0]
    fftsq = abs(fft(array(frate)-frate.mean()))**2.0
    plot(fftsq**0.5)

    # glom0 & glom1
    figure(facecolor='w')
    title('Glomerulus 0 & 1')
    xlabel('time (s)', fontsize='large')
    ylabel('firing rate (Hz)', fontsize='large')
    plot(sinepulsetime, frateResponseList[0][0], color=(1,0,0))
    plot(sinepulsetime, frateResponseList[1][0], color=(0,1,0))
    
    show()
