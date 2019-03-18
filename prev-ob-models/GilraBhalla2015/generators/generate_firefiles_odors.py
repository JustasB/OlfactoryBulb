import sys, pickle

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *

from moose_utils import *
from poisson_utils import *
from data_utils import * # has mpi import and variables also

from pylab import * # part of matplotlib that depends on numpy but not scipy

### generate firefiles i.e. list of spike times from the firerates computed by generate_firerates.py
############## For generating random pulses -- takes maximum time, hence parallelized
## from node000 (not from gj, else master process on node000 cannot connect to display for graphs):
## mpiexec -machinefile ~/hostfile -n <numpulses*NUM_GLOMS+1> ~/Python-2.6.4/bin/python2.6 generate_firefiles_odors.py
## typical value for numpulses = RANDOM_PULSE_NUMS*2
## typically for 10 gloms,
## mpiexec -machinefile ~/hostfile -n 61 ~/Python-2.6.4/bin/python2.6 generate_firefiles_odors.py
############## For generating odor morphs, etc. not parallelized
## python2.6 generate_firefiles_odors.py

## REALRUNTIME is from simset_odor.py
RUNTIME = REALRUNTIME + SETTLETIME

## binning to plot odor responses
bindt = (RUNTIME-SETTLETIME)/respbins
## to make bins causal, keep bin_t at right edge of bin
tlist = arange(SETTLETIME+bindt,RUNTIME+bindt/2.0,bindt)
pulsetlist = arange(pulsebindt,PULSE_RUNTIME+pulsebindt/2.0,pulsebindt)

# time points for the firing rate which is read from a pickled file
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
numt = len(firingtsteps)
extratime = arange(0,2*RESPIRATION,FIRINGFILLDT)
pulsetsteps = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
numtpulse = len(pulsetsteps)

stim_firefiles_dirname = '../firefiles/firefiles'+str(stim_rate_seednum)
if NONLINEAR_ORNS: stim_firefiles_dirname += '_NL'+NONLINEAR_TYPE

def logisticfn(x,offset,width):
    return 1.0/(1.0+exp(-(x-offset)/width))

def thresholded_sigmoid(response):
    ## take a numpy array and threshold and sigmoid() it
    response[ where(response<=0)[0] ] = 0 # clip negative values to 0

    ## I want a sigmoid for only non-negative Reals
    ## Rather than a truncated sigmoid,
    ## I use the cumulative distribution function of the gamma distribution (non-negative domain).
    ## Gamma distribution with scale = 9.0 and shape = 0.5 is close to Gaussian.
    ## Its integral i.e. CDF is a 'sigmoid' over the non-negative Reals.
    #### No I can't use it as it requires scipy.stats.gamma and it is not available on cluster nodes
    #response = scipy.stats.gamma.cdf(response,[9.0]*scipy.stats.gamma.numargs,loc=0,scale=0.5)
    
    ## erf(0.5)~=0.5; set NONLINEAR_CENTER as the center of the sigmoid.
    ## Width should be normal NONLINEAR_WIDTH (spans 2-3 orders of mag in concentration).
    #### No I can't use it as it requires scipy.stats.gamma and it is not available on cluster nodes
    #response = scipy.special.erf((response-NONLINEAR_CENTER)/NONLINEAR_WIDTH)
    ## set the response to zero where negative
    #response[ where(response<=0)[0] ] = 0

    ## finally since gammaCDF and erf can't be used, we go for logistic function.
    ## 1.5*FIRINGMEANA is the center of the sigmoid, FIRINGMEANA/2.0 is width
    response = logisticfn(response,NONLINEAR_CENTER,NONLINEAR_WIDTH)

    ## peak ORN firing set as 3*FIRINGMEANA
    return 3*FIRINGMEANA*response

def NLcondition(frate,glomnum):
    if not NONLINEAR_ORNS: return frate
    if NONLINEAR_TYPE == 'P':
        if glomnum == 0: return thresholded_sigmoid(frate)
        else: return frate
    else:
        if glomnum == 0: return frate
        else: return thresholded_sigmoid(frate)

def odor_and_air_stimuli():
    ## mitral and PG odor ORNs firing files
    seed([stim_rate_seednum]) # Seed numpy's random number generator.
    for glomnum in range(NUM_GLOMS):
        if SHOWFIGS:
            figure()
            title(str(glomnum))
            ylabel('Hz')
            xlabel('time (s)')
        frateperglomList = []
        for odoridx,(odorA, odorB) in enumerate(inputList):
            frate = frateOdorList[glomnum][odoridx]
            frate = NLcondition(frate,glomnum) # output non-linearity if switch is set
            for avgnum in range(MAXNUMAVG):
                mitralfirefilename = stim_firefiles_dirname+\
                    '/firetimes_2sgm_glom_'+str(glomnum)+\
                    '_odor_'+str(odorA)+'_'+str(odorB)+'_avgnum'+str(avgnum)
                ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
                    mitralfirefilename, RUNTIME, firingtsteps)
            if SHOWFIGS:
                ## plotBins here returns firing rate of NUM_ORN_FILES_PER_GLOM combined.
                ## So divide by NUM_ORN_FILES_PER_GLOM.
                ## We just plot ornstimvector_merged for the last avgnum 
                ratebins = [rate/NUM_ORN_FILES_PER_GLOM \
                    for rate in plotBins(ornstimvector_merged, respbins, RUNTIME, SETTLETIME)]
                plot(tlist, ratebins, color=(odorA,odorB,0), marker=',')
            ## extra unmodelled sister cell excitation to granules.
            ## No need of separate files, above files will do for sister excitation
            #mitralfirefilename = '../firefiles/firetimes_2exp_sisters_glom_'+str(glomnum)+\
            #    '_odor_'+str(odorA)+'_'+str(odorB)+'_avgnum'+str(avgnum)
            #write_odor_files(NUM_ORN_FILES_PER_GLOM, frate, mitralfirefilename)

def random_pulse_stimuli(glomnum, pulse_i):
    numpulses = 2*RANDOM_PULSE_NUMS
    ## mitral and PG odor ORNs firing files
    frate_air = randomResponseList[glomnum][0]
    ## the first 'frate' in randomPulseList[glomnum] is for air (rectangle pulse),
    ## the second is a random pulsed air pulse.
    ## then alternately odorA and odorB.
    ## the last odorA and odorB frates are combined and given simultaneously.
    frate = randomResponseList[glomnum][pulse_i+1]
    ## randomPulseList[glomnum][[pulses,frate],...]
    ## pulse_i goes from 0 to numpulses-1 as the
    ## mpi slave allocator skips the first (rectangle air)
    ## and last (combined with 2nd last) pulse
    if pulse_i == 0:
        ## this a just a random air pulse sequence,
        ## so rectangular air pulse 'frate_air' is not added
        frate = frate + ORN_BGND_MEAN
    elif pulse_i < (numpulses-1):
        ## the ORN background is added here, and not used in
        ## fitting the kernel in generate_firerates_odors.py,
        ## as dF/F is not sensistive to background.
        frate = frate + frate_air + ORN_BGND_MEAN
    else:
        ## the last odorA and odorB frates are combined and given simultaneously.
        ## the ORN background is added here, and not used to fit the kernel,
        ## as dF/F is not sensistive to background.
        frate = frate + randomResponseList[glomnum][pulse_i+2] + frate_air + ORN_BGND_MEAN
    frate = NLcondition(frate,glomnum) # output non-linearity if switch is set
    for avgnum in range(MAXNUMAVG):
        mitralfirefilename = stim_firefiles_dirname+\
            '/firetimes_rndpulse_glom_'+str(glomnum)+\
            '_pulse_'+str(pulse_i)+'_avgnum'+str(avgnum)
        ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
            mitralfirefilename, PULSE_RUNTIME, pulsetsteps)
    ## plotBins here returns firing rate of NUM_ORN_FILES_PER_GLOM combined.
    ## So divide by NUM_ORN_FILES_PER_GLOM.
    ## We just plot ornstimvector_merged for the last avgnum 
    ratebins = [rate/NUM_ORN_FILES_PER_GLOM \
        for rate in plotBins(ornstimvector_merged, pulsebins, PULSE_RUNTIME, 0.0)]
    ## extra unmodelled sister cell excitation to granules.
    ## No need of separate files, above files will do for sister excitation
    #mitralfirefilename = '../firefiles/firetimes_2exp_sisters_glom_'+str(glomnum)+\
    #    '_odor_'+str(odorA)+'_'+str(odorB)+'_avgnum'+str(avgnum)
    #write_odor_files(NUM_ORN_FILES_PER_GLOM, frate, mitralfirefilename)
    return ratebins

def interglom_stimuli():
    """
    This reproduces SA small-world excitatory network to PGs via ETs.
    I HAVE NOT YET INCLUDED THE BACKGROUND DUE TO MEAN RESPIRATORY ACTIVITY
    """
    seed([600.0]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    # before taking mean activity over all glomeruli, we should clip negative values
    # because SA cells only receive _excitatory_ input from ET cells.
    frateOdorListClipped = clip(frateOdorList,0.0,1e6)
    # mean activity over gloms (axis=0) as a function of time
    fratemean = mean(frateOdorListClipped,axis=0)
    # Short axon network delays and neuronal integration is modelled as 25ms integrated excitation to PGs
    # the delay should be put in as part of the synapse
    integrate_window_length = int(SA_integrate_time/FIRINGFILLDT)
    # Make a normalized integration window
    # ones() makes array of floats by default, so no integer division problems
    rect = ones((integrate_window_length,))/integrate_window_length
    figure()
    title('SAs')
    ylabel('Hz')
    xlabel('time (s)')
    for odoridx,(odorA,odorB) in enumerate(inputList):
        # convolution takes moving average,
        # I must use mode='full' and take [0:numt], else the SA->(ET->)PG excitation is acausal.
        # of course, the initial part from t = 0 till 'SA_integrate_time' is invalid, but that is within SETTLETIME.
        interglom_rate = convolve(fratemean[odoridx], rect, mode='full')[0:numt]
        SAfirefilename = stim_firefiles_dirname+\
            +'/firetimes_SA_odor_'+str(odorA)+'_'+str(odorB)
        ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM,
            interglom_rate, SAfirefilename, RUNTIME, firingtsteps)
        ratebins = [rate/NUM_ORN_FILES_PER_GLOM \
        for rate in plotBins(ornstimvector_merged, respbins, RUNTIME, SETTLETIME)]
        plot(tlist, ratebins, color=(odorA,odorB,0), marker=',')

def pulseinresp_stimuli(_showfigs):
    seed([stim_rate_seednum]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    stim_frate_filename = '../generators/firerates/firerates_pulseinresp_'+str(stim_rate_seednum)
    stim_frate_filename += '.pickle'
    f = open(stim_frate_filename,'r')
    pulseInRespList,kernels = pickle.load(f)
    f.close()

    pulseinresp_RUNTIME = (NUM_RESPS+1)*RESPIRATION
    pulseinresp_tsteps = arange(0,pulseinresp_RUNTIME+1e-10,FIRINGFILLDT)
    pulseinresp_bintsteps = arange(0,pulseinresp_RUNTIME,pulsebindt)

    for glomnum in range(NUM_GLOMS):
        (_,frates) = pulseInRespList[glomnum]
        if glomnum==0 and _showfigs: figure()
        for ratenum,frate in enumerate(frates):
            frate = NLcondition(frate,glomnum) # output non-linearity if switch is set
            for avgnum in range(MAXNUMAVG):
                mitralfirefilename = stim_firefiles_dirname+\
                    '/firetimes_pulseinresp_glom_'+str(glomnum)+\
                    '_odor'+['R','A','B'][ratenum]+'_avgnum'+str(avgnum)
                ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
                    mitralfirefilename, pulseinresp_RUNTIME, pulseinresp_tsteps)
            ## plotBins here returns firing rate of NUM_ORN_FILES_PER_GLOM combined.
            ## So divide by NUM_ORN_FILES_PER_GLOM.
            ## We just plot ornstimvector_merged for the last avgnum 
            pulseinrespbins = int((pulseinresp_RUNTIME/pulsebindt))
            ratebins = [rate/NUM_ORN_FILES_PER_GLOM \
                for rate in plotBins(ornstimvector_merged, pulseinrespbins, pulseinresp_RUNTIME, 0.0)]
            if glomnum==0 and _showfigs:
                plot(pulseinresp_bintsteps,ratebins,['k','r','b'][ratenum])

def scaledpulses_stimuli(_showfigs):
    seed([stim_rate_seednum]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    stim_frate_filename = '../generators/firerates/firerates_scaledpulses_width'+\
        str(scaledWidth)+'_'+str(stim_rate_seednum)
    stim_frate_filename += '.pickle'
    f = open(stim_frate_filename,'r')
    scaledPulses,kernels = pickle.load(f)
    f.close()

    scaledpulses_tsteps = arange(0,SCALED_RUNTIME+1e-10,FIRINGFILLDT)
    scaledpulses_bintsteps = arange(0,SCALED_RUNTIME,pulsebindt)

    for glomnum in range(NUM_GLOMS):
        frates = scaledPulses[glomnum]
        if glomnum==0 and _showfigs: figure()
        for scalenum,scale in enumerate(scaledList):
            for fratenum,frate in enumerate(frates[scalenum]):
                frate = NLcondition(frate,glomnum) # output non-linearity if switch is set
                for avgnum in range(MAXNUMAVG):
                    mitralfirefilename = stim_firefiles_dirname+\
                        '/firetimes_scaledpulses_width'+str(scaledWidth)+'_glom_'+str(glomnum)+\
                        '_odor'+['A','B'][fratenum]+'_scale'+str(scalenum)+'_avgnum'+str(avgnum)
                    ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM, frate,\
                        mitralfirefilename, SCALED_RUNTIME, scaledpulses_tsteps)
                if glomnum==0 and _showfigs:
                    ## plotBins here returns firing rate of NUM_ORN_FILES_PER_GLOM combined.
                    ## So divide by NUM_ORN_FILES_PER_GLOM.
                    ## We just plot ornstimvector_merged for the last avgnum 
                    scaledpulsebins = int((SCALED_RUNTIME/pulsebindt))
                    ratebins = [rate/NUM_ORN_FILES_PER_GLOM \
                        for rate in plotBins(ornstimvector_merged, scaledpulsebins, SCALED_RUNTIME, 0.0)]
                    plot(scaledpulses_bintsteps[:len(ratebins)],ratebins)

def interglom_pulsestimuli():
    """
    This reproduces SA small-world excitatory network to PGs via ETs.
    I HAVE NOT YET INCLUDED THE BACKGROUND DUE TO MEAN RESPIRATORY ACTIVITY
    """
    seed([800.0]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    # before taking mean activity over all glomeruli, we should clip negative values
    # because SA cells only receive _excitatory_ input from ET cells.
    # randomPulseList[glomnum][[pulses,frate],...]
    numpulses = len(randomPulseList[0])
    # the first pulse is air, then A and B alternate, the last has to be combined with the second last.
    frateList = empty([NUM_GLOMS,numpulses-2,len(randomPulseList[0][0][1])])
    for glomnum in range(NUM_GLOMS):
        frate_air = randomPulseList[glomnum][0][1]
        for pulsenum in range(1,numpulses-1):
            if pulsenum == 1:
                frateList[glomnum][pulsenum-1] =\
                    array( randomPulseList[glomnum][pulsenum][1] + ORN_BGND_MEAN )
            elif pulsenum < numpulses-2:
                frateList[glomnum][pulsenum-1] =\
                    array( randomPulseList[glomnum][pulsenum][1] + frate_air + ORN_BGND_MEAN)
            else:
                frateList[glomnum][pulsenum-1] = array( randomPulseList[glomnum][pulsenum+1][1] + \
                    randomPulseList[glomnum][pulsenum+1][1] + frate_air + ORN_BGND_MEAN)
    fratePulseListClipped = clip(frateList,0.0,1e6)
    # mean activity over gloms (axis=0) as a function of time
    fratemean = mean(fratePulseListClipped,axis=0)
    # Short axon network delays and neuronal integration is modelled as 25ms integrated excitation to PGs
    # the delay should be put in as part of the synapse
    integrate_window_length = int(SA_integrate_time/FIRINGFILLDT)
    # Make a normalized integration window
    # ones() makes array of floats by default, so no integer division problems
    rect = ones((integrate_window_length,))/integrate_window_length
    figure()
    title('SAs Pulses')
    ylabel('Hz')
    xlabel('time (s)')
    for pulse_i in range(numpulses-2):
        ## convolution takes moving average,
        ## I must use mode='full' and take [0:numt], else the SA->(ET->)PG excitation is acausal.
        ## of course, the initial part from t = 0 till 'SA_integrate_time' is invalid,
        ## but that is within SETTLETIME.
        interglom_rate = convolve(fratemean[pulse_i], rect, mode='full')[0:numtpulse]
        SAfirefilename = stim_firefiles_dirname+\
            +'/firetimes_SA_rndpulse'+'_pulse_'+str(pulse_i)
        ornstimvector_merged = write_odor_files(NUM_ORN_FILES_PER_GLOM,\
            interglom_rate, SAfirefilename, PULSE_RUNTIME, pulsetsteps)
        ratebins = [rate/NUM_ORN_FILES_PER_GLOM\
            for rate in plotBins(ornstimvector_merged, pulsebins, PULSE_RUNTIME, 0.0)]
        plot(pulsetlist, ratebins, color=(pulse_i/float(numpulses),1-pulse_i/float(numpulses),0), marker=',')

def write_odor_files(numfiles, frate, firefilename, runtime, tsteps, vary=None):
    frate = clip(frate,0.0,1e6)
    PEAK_FIRING_RATE = max(frate)*2
    if PEAK_FIRING_RATE < 1.0: PEAK_FIRING_RATE = 1.0
    ornstimvector_merged = []
    firefile = open(firefilename+'.txt','w')
    for i in range(numfiles):
        ## if variation is called for, then vary the firing rate by a factor for each file
        if not vary is None: frate_rand = frate * normal(vary[0],vary[1])/vary[0]
        else: frate_rand = frate
        ornstimvector = poissonTrainVaryingRate(runtime,PEAK_FIRING_RATE,REFRACTORY,tsteps,frate_rand)
        firefile.write(' '.join([str(t) for t in ornstimvector])+'\n')
        ornstimvector_merged.extend(ornstimvector)
    firefile.close()
    ornstimvector_merged.sort()
    print "wrote ", firefilename
    return ornstimvector_merged

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    #seed([600.0]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time

    if 'NOSHOW' not in sys.argv: SHOWFIGS=True
    else: SHOWFIGS=False

    stim_frate_filename = '../generators/firerates/firerates_2sgm_'+str(stim_rate_seednum)
    stim_frate_filename += '.pickle'
    f = open(stim_frate_filename,'r')
    frateOdorList,fratePulseList,randomPulseList,\
    randomPulseStepsList,randomResponseList,kernels\
        = pickle.load(f)
    f.close()

    if mpisize == 1:
        ## seeds are set in each of the below functions:
        ## so you may comment / uncomment below functions
        ## depending on what you want to generate.
        #interglom_stimuli()
        if 'PULSEINRESP' in sys.argv: pulseinresp_stimuli(SHOWFIGS)
        elif 'SCALEDPULSES' in sys.argv: scaledpulses_stimuli(SHOWFIGS)
        else: odor_and_air_stimuli()
        #interglom_pulsestimuli()
        if SHOWFIGS: show()
    else:
        ## Random pulses, not doing two pulse stimulus
        ## since random pulses convey all the info.
        numpulses = RANDOM_PULSE_NUMS*2
        numpulses_fl = float(numpulses)
        if mpirank == boss: # boss collates
            for glomnum in range(NUM_GLOMS):
                if SHOWFIGS:
                    figure()
                    title("Pulses for glom "+str(glomnum))
                    ylabel('Hz')
                    xlabel('time (s)')
                for pulsenum in range(numpulses):
                    procnum = glomnum*numpulses + pulsenum + 1
                    print 'waiting for process '+str(procnum)+'.'
                    ratebins = mpicomm.recv(source=procnum, tag=0)
                    if SHOWFIGS:
                        plot(pulsetlist, ratebins, \
                            color=(pulsenum/numpulses_fl,1-pulsenum/numpulses_fl,0), marker=',')
                    print 'received from process '+str(procnum)+'.'
            if SHOWFIGS: show()
        else:
            seed([500.0*mpirank]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
            glomnum = (mpirank-1)/numpulses
            pulsenum = (mpirank-1)%numpulses
            ratebins = random_pulse_stimuli(glomnum,pulsenum)
            mpicomm.send( ratebins, dest=boss, tag=0 )
