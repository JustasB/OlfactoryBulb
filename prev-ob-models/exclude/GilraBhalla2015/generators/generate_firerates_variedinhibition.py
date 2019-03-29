import sys, pickle

sys.path.extend(["..","../networks","../simulations"])

from generate_firerates_odors import *

#### We have dual sigmoid functions as responses for
#### respiration, odor A and odor B (each is different for every glomerulus).
#### These responses are deconvolved with respiration to get kernels.
#### For Adil type odor morph experiments, these kernels are again convolved with respiration.
#### For Priyanka's data, these kernels are convolved with odor and air pulses.
#### This only gives firing rates.
#### The firing rate as a function of time is fed
#### to a Poisson spike generator in generate_firefiles.py .

frateOdorList = []
kernels = []

def compute_doublesigmoid_params(delay, risetime, duration):
    """ return the double sigmoid params for the given delay, risetime and duration."""
    ######### get double sigmoid params
    spread2_factor = 4.0 # not used presently
    spread1, center2, spread2 = \
        compute_dblsigmoid_params(risetime, duration, spread2_factor)
    #spread2 = spread2_factor*spread1
    ######### check if the solved params are reasonable
    ######### risetime and duration match and spread2 is at least twice spread1.
    actual_risetime = invert_dblsigmoid(0.0, spread1, center2, spread2, 0.9) - \
        invert_dblsigmoid(0.0, spread1, center2, spread2, 0.1)
    actual_duration = invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5, risephase=False) - \
        invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5)
    if abs(actual_risetime-risetime)>1e-4 or abs(actual_duration-duration)>1e-4:
        print "The actual risetime and duration are", actual_risetime, actual_duration
        print "But the required risetime and duration are", risetime, duration
        sys.exit(1)
    if spread2 < 2*spread1:
        print "spread2 < 2*spread1, spread1 =",spread1,", spread2 =",spread2
        sys.exit(1)
    ########## shift curve to include latency
    ## I had taken center1 = 0.0, now shift the curve
    ## so that t=0 is actually where first sigmoid is 0.05*peak
    ## and then add the delay / latency to it!
    offset = - invert_dblsigmoid(0.0, spread1, center2, spread2, 0.05) + delay
    return [offset, spread1, center2+offset, spread2]

def variedinh_stimuli():
    ## mitral and PG odor ORNs firing files

    ## glom0 - odor kernel
    odorparams0_e = [FIRINGMEANA]
    odorparams0_e.extend( compute_doublesigmoid_params(2*delay_mean,1.5*risetime_mean,1.5*duration_mean) )
    odorparams0_i = [0.0]
    odorparams0_i.extend( compute_doublesigmoid_params(delay_mean,risetime_mean,duration_mean) )
    kernel0odor = getkernel(odorparams0_e,odorparams0_i)
    ## glom0 - air kernel
    odorparams0_e = [0.0]#[FIRINGRATEAIR_MEAN]
    odorparams0_e.extend( compute_doublesigmoid_params(delay_mean,risetime_mean,duration_mean) )
    odorparams0_i = [0.0]
    odorparams0_i.extend( compute_doublesigmoid_params(delay_mean,risetime_mean,duration_mean) )
    kernel0air = getkernel(odorparams0_e,odorparams0_i)
    kernels.append((kernel0air,kernel0odor))
    ## glom0 - firing rates
    ## odorparams are relics of non-kernel era,
    ## actually not used by receptorFiringRate(), only kernels are used.
    frate = receptorFiringRate(1.0, 1.0, 0.0,\
        odorparams0_e, odorparams0_i,\
        odorparams0_e, odorparams0_i,\
        odorparams0_e, kernel0odor, kernel0air, kernel0air)
    frateOdorList.append([frate])

    risetime = risetime_mean
    duration = duration_mean
    for glomnum in range(1,NUM_GLOMS):
        print "Computing kernels and firing rates for glomerulus", glomnum

        ### For each glomerulus, set its responses to odor and air
        ### as difference of excitatory and inhibitory random dual-exponentials
        ### for each different respiration phase
        
        kernels.append([])        
        frateOdorList.append([])
        delta_delay = (delay_mean+3*delay_sd)/float(NUMINHS)
        for delay in arange(0.0,delay_mean+3*delay_sd,delta_delay):
            print "Computing kernels and firing rates for delay =", delay

            odorparamsA_e = [FIRINGMEANA]
            odorparamsA_e.extend( compute_doublesigmoid_params(delay,risetime,duration) )
            odorparamsA_i = [0.0]
            odorparamsA_i.extend( compute_doublesigmoid_params(delay,risetime,duration) )

            ## kernel for Priyanka can be obtained
            ## by numerically deconvolving with the respiration pulse.
            kernelA = getkernel(odorparamsA_e,odorparamsA_i)

            ### respiration also has both excitatory and inhibitory dual exponentials.
            odorparamsR_e = [FIRINGRATEAIR_MEAN]
            odorparamsR_e.extend( compute_doublesigmoid_params(delay,risetime,duration) )
            odorparamsR_i = [0.0]
            odorparamsR_i.extend( compute_doublesigmoid_params(delay,risetime,duration) )

            kernelR = getkernel(odorparamsR_e,odorparamsR_i)

            kernels[-1].extend( (kernelR,kernelA) )

            # firing rates
            frate = receptorFiringRate(1.0, 1.0, 0.0,\
                odorparamsA_e, odorparamsA_i,\
                odorparamsR_e, odorparamsR_i,\
                odorparamsR_e, kernelA, kernelR, kernelR)
            ## important to put within [] or (...,) for extend
            frateOdorList[-1].extend([frate])

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    seed([stim_rate_seednum]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time

    variedinh_stimuli()

    filename = 'firerates/firerates_2sgm_variedinh.pickle'
    fireratefile = open(filename,'w')
    pickle.dump( (frateOdorList,kernels), fireratefile)
    fireratefile.close()
    print "wrote",filename

    # glom0 & glom1
    figure(facecolor='w')
    title('Glomerulus 0 & 1')
    xlabel('time (s)', fontsize='large')
    ylabel('firing rate (Hz)', fontsize='large')
    plot(firingtsteps, frateOdorList[0][0], color=(0,0,0), marker=',')
    totinh = float(len(frateOdorList[1]))
    for inhidx,frate in enumerate(frateOdorList[1]):
        plot(firingtsteps, frate, color=(inhidx/totinh,1-inhidx/totinh,0), marker=',')
    
    show()
