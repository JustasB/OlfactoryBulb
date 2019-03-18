
## USAGE:
## python2.6 generate_firefiles_gran_baseline.py <invitro|noresp>
## no extra arg for invivo, use extra arg for invitro and noresp options

import sys
sys.path.extend(['..','../networks'])

from networkConstants import * # for gran_spines
from stimuliConstants import *
from neuro_utils import *
    
from pylab import * # part of matplotlib that depends on numpy but not scipy

from generate_firerates_odors import *
from generate_firefiles_odors import *

## granule baseline firing - write files
def gran_files(mit_base_rate, filebase, invitro_str):
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    seed([100.0]) ##### Seed numpy's random number generator.

    for i in range(MAXNUMAVG_GRANS):
        firefilename = filebase+\
            '/firetimes_gran_baseline'+invitro_str+'_'+str(i)+'.txt'
        firefile = open(firefilename,'w')
        for num in range(num_gran_baseline_files):            
            #### for each spine on the granule,
            #### make a mitral fire at a baseline rate (slightly varied randomly).
            #ornstimvector_merged = []
            ## there will be at least one non-baseline synapse connected already.
            #for i in range(gran_spines-1):
            #    mit_base_f_rand = uniform(mit_base_rate*0.75,mit_base_rate*1.25)
            #    ornstimvector = poissonTrain(MAXRUNTIME+SETTLETIME,\
            #        mit_base_f_rand,REFRACTORY) # from moose_utils.py
            #    ornstimvector_merged.extend(ornstimvector)
            #ornstimvector_merged.sort()
            #firefile.write(' '.join([str(t) for t in ornstimvector_merged])+'\n')
            
            ### instead of merging and sorting large number of granules,
            ### (smaller variance due to averaging);
            ### just generate the full baseline excitation to granule
            ### with frate as gaussian/normal with variance ~= mean
            ### roughly matches Carleton (in vitro) & Cang&Isaacson (in vivo)
            gran_base_rate = mit_base_rate*(gran_spines-1)
            gran_base_f_rand = normal(gran_base_rate,sqrt(gran_base_rate))
            if gran_base_f_rand<0: gran_base_f_rand = 0
            ornstimvector = poissonTrain(MAXRUNTIME+SETTLETIME,\
                gran_base_f_rand,REFRACTORY) # from moose_utils.py
            firefile.write(' '.join([str(t) for t in ornstimvector])+'\n')
            
        firefile.close()
        print "wrote ", firefilename

def gran_files_resp(filebase, extrastr, weight, showfig):
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    seed([100.0]) ##### Seed numpy's random number generator.

    delay = delay_mean
    risetime = risetime_mean
    duration = duration_mean
    ######### get double sigmoid params
    spread2_factor = 4.0 # not used presently
    spread1, center2, spread2 = \
        compute_dblsigmoid_params(risetime, duration, spread2_factor)
    #spread2 = spread2_factor*spread1
    ########## shift curve to include latency
    ## I had taken center1 = 0.0, now shift the curve
    ## so that t=0 is actually where first sigmoid is 0.05*peak
    ## and then add the delay / latency to it!
    offset = - invert_dblsigmoid(0.0, spread1, center2, spread2, 0.05) + delay
    odorparamsR_e = [ gran_baseline_for_resp_tuning ]
    odorparamsR_e.extend( [offset, spread1, center2+offset, spread2] )

    ## inhibition kicks in 200ms later, delay_mean above is 154ms
    delay = 200e-3+delay_mean # ms
    risetime = risetime_mean
    duration = duration_mean/1.5
    ######### get double sigmoid params
    spread2_factor = 4.0 # not used presently
    spread1, center2, spread2 = \
        compute_dblsigmoid_params(risetime, duration, spread2_factor)
    #spread2 = spread2_factor*spread1
    ########## shift curve to include latency
    ## I had taken center1 = 0.0, now shift the curve
    ## so that t=0 is actually where first sigmoid is 0.05*peak
    ## and then add the delay / latency to it!
    offset = - invert_dblsigmoid(0.0, spread1, center2, spread2, 0.05) + delay
    odorparamsR_i = [ gran_baseline_for_resp_tuning ]
    odorparamsR_i.extend( [offset, spread1, center2+offset, spread2] )

    kernelR = getkernel(odorparamsR_e,odorparamsR_i)
    frate = float(weight) * (
        receptorFiringRate(0, 1.0, 0, \
            kernelR, kernelR, kernelR) \
            + gran_baseline_apart_from_resp_tuning
        )
    figure()
    plot(firingtsteps, frate, color=(1,0,0), marker=',')
    ylim(0,60)
    ## Take the last respiration period from the end and integrate
    lastfrate = frate[-int(RESPIRATION/FIRINGFILLDT):]
    frateavg = sum([fratei*FIRINGFILLDT for fratei in lastfrate]) / RESPIRATION
    print "Average firing rate for respiratory tuned response is",frateavg
    firefilename = filebase+'/firetimes_gran_baseline'+extrastr+'_'
    ## RUNTIME is defined in simset_odor.py imported via generate_firefiles_odors.py
    for i in range(MAXNUMAVG_GRANS):
        ## write_odor_files() is in generate_firefiles.py
        ornstimvector_merged = write_odor_files(num_gran_baseline_files,
            frate, firefilename+str(i), RUNTIME, firingtsteps, 
            vary=(frateavg,sqrt(frateavg)) )
    if showfig:
        figure()
        ratebins = [rate/float(num_gran_baseline_files)\
            for rate in plotBins(ornstimvector_merged, respbins, RUNTIME, SETTLETIME)]
        plot(tlist, ratebins, marker=',')
        show()

if __name__ == "__main__":
    ## seed for every function separately

    filebase = '../firefiles/firefiles_baseline'
    #filebase = '../firefiles/firefiles_whitenoise'
    #filebase = '../firefiles/firefiles_variedinh'
    
    if len(sys.argv)>1:
        arg1 = sys.argv[1]
        if arg1=='invitro':
            ### granule baseline firing in vitro
            ### gran_files has a loop of gran_spines iterations
            ### that effectively makes firing rate = mit_base_f_invitro*gran_spines
            gran_files(mit_base_f_invitro, filebase, '_invitro')
        elif arg1=='noresp':
            ### granule baseline firing in vivo :
            ### for activity-dep inhibition & random pulses for tracheotomized rat
            ### constant rate
            ### gran_files has a loop of gran_spines iterations
            ### that effectively makes firing rate = mit_base_f*gran_spines
            gran_files(mit_base_f, filebase, '_noresp')
            gran_files(extraexc_factor*mit_base_f, filebase, '_noresp_extra')
        else:
            print "Unrecognized param",arg1
    else:
        ### In vivo, we have a constant baseline firing and a respiratory tuned firing on top of it,
        ### since average mitral response is respiration tuned - see network_constants.py
        ### mean is roughly 34-35Hz when the weight (2nd param below) is 1.0
        gran_files_resp(filebase, '', 1.0, showfig=False)
        gran_files_resp(filebase, '_extra', extraexc_factor, showfig=False)
