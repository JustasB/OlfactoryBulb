import sys, pickle

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *

from moose_utils import *
from neuro_utils import * # has the dual exponential, DoG (subtractedGaussian) and gamma fns
from data_utils import *

from pylab import * # part of matplotlib that depends on numpy but not scipy
import scipy.optimize
import scipy.interpolate
import scipy.special # special function erf()
import scipy.stats
## for scipy.signal.deconvolve but it doesn't work
## if any of the respiration pulse coeff's are zero - see below.
#import scipy.signal

##### USAGE:
## python2.6 generate_firerates_odors.py


"""	1) respirationPulse1 (one respiration pulse) is peak normalized to the absolute peak!
        Not to the positive peak!!! And then rectified.
        respirationPulses (multiple (2) respiration pulses) is also peak normalized.
	2) odor_and_air_stimuli(): For each of odorA, odorB and respiration, an excitatory & inhibitory component is generated:
		a. compute_odor_params(): For each exc / inh component first draw latency, and duration from gamma functions.
            Set risetime from latency. Then fit a double sigmoid (4 params) to this to get a continuous curve.
		b. get_kernel(): Combine the exc and inh components and for each of odorA and odorB,
            fit a kernel so that convolving with the peak normalized respirationPulse1 gives the combined exc+inh response.
            [We have used same params for exc and inh components,
            Carey et al characterized only exc components.
            They measure dF/F, so air response is already subtracted.]
		c. Adil's & Ashesh's ORN responses: receptorFiringRate(): convolve respirationPulses with each kernel,
            and weight by odour conc (percent conc maps to weight, 1% maps to 1.0 factor).
            Add for each odour (conc% as weight) and for air (weight=1),  and add const ORN_BGND_MEAN.
		d. Priyanka's ORN responses:
			i. receptorFiringRateRespPulse(): Create an air pedestal of weight 0.25 (avg of peak=1 respirationPulse1?).
			ii. receptorFiringRateRandomPulse(): Create m-sequence pulses
                with on-off times as olfactometer_tau and saturation/peak weight = 1. (thus 1% odor at peak).
                [But if my air pedestal is 0.25, I should keep this peak also at 0.25!
                Suction is constant 0.25 of air peak, volume of air is integral of this.
                Number of molecules will also be integral of 0.25 const, instead of integral of respiration pulse with peak=1.]
			iii. pulseConvolve(): convolve the above resp pedestal and random pulses with respective air and odour kernels.
			iv. In generate_firefiles_odor.py, combine the resp pedestal and odor random pulse and ORN_BGND_MEAN
                to get the ORN frate for an actual m-sequence delivery.
"""

#### We have dual sigmoid functions as responses for
#### respiration, odor A and odor B (each is different for every glomerulus).
#### These responses are deconvolved with respiration to get kernels.
#### For Adil type odor morph experiments, these kernels are again convolved with respiration.
#### For Priyanka's data, these kernels are convolved with odor and air pulses.
#### This only gives firing rates.
#### The firing rate as a function of time is fed
#### to a Poisson spike generator in generate_firefiles.py .

#### A design philosophy has been to ensure that random number call is made
#### even if that number is not used, due to some switch / if statement.
#### This ensures that the same response gets generated
#### irrespective of switch / number of glomeruli.

## salient responses
if stim_rate_seednum<0:
    stim_nums = len(salient_seed_glom_mapA)
    stim_seed = -int(stim_rate_seednum)-1 ## -1 to start indices from 0 instead of 1.
    ampl_scale = stim_seed / stim_nums # integer division
    stim_seed = stim_seed % stim_nums # modulo
    ampl_scale = ampl_scale % 2 + 1 # only two amplitude scales: 1.0 and 2.0
    responseA_indices = salient_seed_glom_mapA[stim_seed] # odorA response indices for gloms0,1,2
    responseB_indices = salient_seed_glom_mapB[stim_seed] # odorB response indices for gloms0,1,2    

frateOdorList = []
fratePulseList = []
randomPulseList = []
randomPulseStepsList = []
randomResponseList = []
pulseInRespList = []
scaledPulses = []
kernels = []

def periodic_respiration_fn(t,negclip=True):
    """
    Respiration function but periodic.
    """
    tmodrespn = ((t/RESPIRATION)-floor(t/RESPIRATION))*RESPIRATION
    return respiration_function(tmodrespn,negclip)

def respiration_function(t,negclip=True):
    """
    Normalized (integral~=1) respiration function (approx norm as fn does not go to inf).
    Inspiration is modelled as a fast positive dual exponential (integral ~= +0.6)
    and a slow dual exponential function (integral ~= +0.5).
    Expiration as a dip is modelled as a negative flipped dual exponential (integral ~= -0.1)
    But the response is then rectified,
     else the negative dip will act like 'negative' odor in convolution.
    After rectification, the already approximate normalization does not hold!!
    Actually, there is no meaning to above integral normalization here.
    There is a peak normalization performed later, but that is arbitrary anyways.
    What matter are the ratios of 'suction' air and odor pulse heights to this 'respiration' peak.
    """
    value = 0.6*dualexpfunction(t,0.0,tau1rins,tau2rins) +\
        0.5*dualexpfunction(t,0.0,tau1rinsslow,tau2rinsslow) -\
        0.1*dualexpfunction(-(t-RESPIRATION),0.0,tau1rexp,tau2rexp)
    if negclip: return clip(value,0.0,1e6)
    else: return value

def periodic_impulse(t):
    if abs((t/RESPIRATION)-floor(t/RESPIRATION)) < FIRINGFILLDT/2: return 1.0
    else: return 0.0

def normresponse(wavefm):
    """
    Peak-normalized waveform.
    """
    wavemax = absolute(wavefm).max()
    if wavemax!=0:
        return wavefm/wavemax
    else:
        return wavefm

firingtsteps = arange(0,ODORRUNTIME+1e-10,FIRINGFILLDT) # include the last RUNTIME point also.
numt = len(firingtsteps)
extratime = arange(0,(NUM_RESPS+1)*RESPIRATION+1e-10,FIRINGFILLDT)
respirationtime = arange(0,RESPIRATION-1e-10,FIRINGFILLDT) # don't include last point
randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
binsize = int(fitting_dt/FIRINGFILLDT) # ensure these are multiples
extratime_decimated = extratime[::binsize]

## NUM_RESPS respiration pulses.
#respirationPulses = array([math.sin(t/RESPIRATION * 2*math.pi) for t in extratime])
respirationPulses = array([periodic_respiration_fn(t) for t in extratime])
## Peak normalize the respiration pulses
respirationPulses = normresponse(respirationPulses)
## A single respiration pulse extended with zeroes over one more respiration cycle.
respirationPulse1 = [respiration_function(t) for t in respirationtime]
respirationPulse1.extend([0]*(len(respirationtime)+1))
## Peak normalize the respiration pulse. But this is arbitrary.
## What matter are the ratios of air and odor 'suction' pulse heights to this 'respiration' peak.
respirationPulse1 = normresponse(respirationPulse1)
# NUM_RESP impulses to convolve with (instead of respiration) for Adil's experiment
impulses = array([periodic_impulse(t) for t in extratime])
startidx = int((RESPIRATION-SETTLETIME)/FIRINGFILLDT)

def respiration_integrate():
    """ This fn is just used once to find out integral and mean and peak of the respiration signal.
    Respiration is peak normalized later, so take ratio of integral or mean to the peak. """
    resp1cycle = array([respiration_function(t) for t in respirationtime])
    print "Integral of the rectified respiration pulse over one cycle is",sum(resp1cycle*FIRINGFILLDT)
    print "Peak of the respiration pulse is",max(resp1cycle)
    print "Average of the rectified respiration pulse over one cycle is",mean(resp1cycle)

#def odorResponse(odorparams_e, odorparams_i):
def odorResponse(kernel):
    """
    For Adil's experiments, just calculate response as a double sigmoid.
    The parameters (Carey et al) are for the respiration tuned response so just convolve with impulses.
    Just in case I use params for kernel, I have a RESP_CONVOLVE switch.
    Assumes resp has two or more cycles.
    Discard most of the first cycle except the last 'settletime'.
    Take the second cycle fully i.e. 'settletime' to 'runtime' contains a full cycle.
    The first cycle is put in only to generate a valid convolved second cycle.
    """
    #if RESP_CONVOLVE:
        #resp = respirationPulses2
        #### I must actually deconvolve the response and use that as kernel
    #else: resp = impulses2
    ## This is not the true kernel (using Carey et al params), it is already respiration tuned.
    ## I just convolve with impulses, rather than with respiration waveform.
    #ekernel = [dblsigmoid(t,odorparams_e[1],odorparams_e[2],odorparams_e[3],odorparams_e[4]) for t in extratime]
    ## The peak response for this odor is set, after convolution with respiration
    #eresponse = odorparams_e[0] * normconvolve( resp, ekernel )
    #ikernel = [dblsigmoid(t,odorparams_i[1],odorparams_i[2],odorparams_i[3],odorparams_i[4]) for t in extratime]
    ## The peak response for this odor is set, after convolution with respiration
    #iresponse = odorparams_i[0] * normconvolve( resp, ikernel )
    ## take only the end SETTLETIME amount of the first respiratory cycle and the full 2nd & 3rd cycle.
    ## odor response is zero for negative time, then SETTLETIME of end-resp cycle, then 2 full resp cycles.
    ## for long kernels, the 1st full and 2nd full resp odor responses will be different.
    #return array(eresponse-iresponse)[startidx:startidx+numt]

    ## pre-chop off the respirationPulses till startidx (only end SETTLETIME of 1st cycle remains),
    ## and then convolve.
    ## always multiply convolution_result by dt to leave result invariant on decimation.
    response = convolve(respirationPulses[startidx:],kernel,mode='full')*FIRINGFILLDT
    ## take only the end SETTLETIME amount of the first respiratory cycle and NUM_RESP cycles fully.
    ## THE FIRST CYCLE CONVOLUTION IS ALSO NOT VALID!
    return response[0:numt]

def normconvolve(resp,kernel):
    """
    Peak-normalized convolve.
    """
    ## always multiply convolution_result by dt to leave result invariant on decimation.
    conv = convolve(resp,kernel,mode='full')*FIRINGFILLDT
    convmax = absolute(conv).max()
    if convmax!=0:
        return conv/convmax
    else:
        return conv

def receptorFiringRate(respwt, concA, concB, kernelA, kernelB, kernelR, glomnum = None):
    combinedresponse = concA * CONC_SCALE * odorResponse(kernelA) + \
        concB * CONC_SCALE * odorResponse(kernelB) + \
        respwt * odorResponse(kernelR)
    ## constant response
    if stim_rate_seednum<=-19 and CONST_GLOM0 and glomnum==0:
        return [ 0.5*FIRINGMEANA*ampl_scale ] * len(combinedresponse)
        #return [ 0.5*FIRINGMEANA*ampl_scale + ORN_BGND_MEAN ] * len(combinedresponse)
    ## the ORN background is added here, and not used to fit the kernel,
    ## as dF/F is not sensistive to background.
    return combinedresponse + ORN_BGND_MEAN

def receptorFiringRatePulse(pulsedelayA, pulsedurationA, pulsedelayB, pulsedurationB,
        kernelA, kernelB, kernelR):
    ##### in Hz
    len_pulsetime = int(PULSE_DURATION/FIRINGFILLDT)+1
    pulseA = zeros([len_pulsetime])
    pulseA[int(pulsedelayA/FIRINGFILLDT):int((pulsedelayA+pulsedurationA)/FIRINGFILLDT)] = \
                                    1*resp_pulsewt*CONC_SCALE
    pulseB = zeros([len_pulsetime])
    pulseB[int(pulsedelayB/FIRINGFILLDT):int((pulsedelayB+pulsedurationB)/FIRINGFILLDT)] = \
                                    1*resp_pulsewt*CONC_SCALE
    pulseR = zeros([len_pulsetime])
    # respiration pulse starts 50 ms before any of the odors. If negative, clip to zero.
    pulsedelayR = clip( min(pulsedelayA,pulsedelayB)-50e-3 , 0.0, 1e6 )
    # respiration pulse stops 50ms after both the odors.
    pulsedurationR = max(pulsedelayA+pulsedurationA,pulsedelayB+pulsedurationB)+50e-3
    pulseR[int(pulsedelayR/FIRINGFILLDT):int(pulsedurationR/FIRINGFILLDT)]=1*resp_pulsewt
    combinedresponse = convolve(pulseA, kernelA, mode='full') + \
        convolve(pulseB, kernelB, mode='full') + \
        convolve(pulseR, kernelR, mode='full')
    ## always multiply convolution_result by dt to leave result invariant on decimation.
    combinedresponse = combinedresponse * FIRINGFILLDT
    combinedresponse = combinedresponse[0:len_pulsetime]
    # the ORN background is added here, and not used to fit the kernel,
    # as dF/F (from which kernel is computed) is not sensistive to background.
    return combinedresponse + ORN_BGND_MEAN

def receptorFiringRateScaledPulse(scaleA, scaleB, kernelA, kernelB, kernelR):
    ##### in Hz
    len_pulsetime = int(SCALED_RUNTIME/FIRINGFILLDT)+1
    pulseA = zeros([len_pulsetime])
    pulseA[int(PULSE_START/FIRINGFILLDT):int((PULSE_START+scaledWidth)/FIRINGFILLDT)] = \
                                    1*scaleA*resp_pulsewt
    pulseB = zeros([len_pulsetime])
    pulseB[int(PULSE_START/FIRINGFILLDT):int((PULSE_START+scaledWidth)/FIRINGFILLDT)] = \
                                    1*scaleB*resp_pulsewt
    pulseR = zeros([len_pulsetime])
    # respiration pulse starts after settle time
    pulsedelayR = SETTLETIME
    # respiration pulse stops at the end of full runtime
    pulsedurationR = SCALED_RUNTIME
    pulseR[int(pulsedelayR/FIRINGFILLDT):int(pulsedurationR/FIRINGFILLDT)]=1*resp_pulsewt
    combinedresponse = convolve(pulseA, kernelA, mode='full') + \
        convolve(pulseB, kernelB, mode='full') + \
        convolve(pulseR, kernelR, mode='full')
    ## always multiply convolution_result by dt to leave result invariant on decimation.
    combinedresponse = combinedresponse * FIRINGFILLDT
    combinedresponse = combinedresponse[0:len_pulsetime]
    # the ORN background is added here, and not used to fit the kernel,
    # as dF/F (from which kernel is computed) is not sensistive to background.
    return combinedresponse + ORN_BGND_MEAN

def receptorFiringRateRespPulse():
    """ Make an analog and a digital pedestal for air,
    on which the air random pulse will ride, (added up in generate_firefiles) .
    Since respiration periodic function is peak normalized,
    resp_pulsewt is typically set low at 1/3.0,
    to get (250 to 300ml/min) / (900ml/min) ratio and
    similar physiological firing rate on convolving with the same kernel:
    the avg respiration < 0.25 times peak.
    Pedestal mean=peak, but binary pulses mean will be half their peak."""
    ##### in Hz
    len_pulsetime = int(PULSE_RUNTIME/FIRINGFILLDT)+1
    # create the respiratory constant pulse around the random odor pulse
    # with AIR_BUFFER period on both sides.
    pulseR = zeros([len_pulsetime])
    pulseRsteps = zeros([len_pulsetime])
    startidx = int((PULSE_START-PULSE_AIR_BUFFER)/FIRINGFILLDT)
    stopidx = int((PULSE_START+PULSE_DURATION+PULSE_AIR_BUFFER)/FIRINGFILLDT)
    pulseR[startidx:stopidx] = 1*resp_pulsewt
    pulseRsteps[startidx:stopidx] = 1
    return pulseR,pulseRsteps

len_pulsetime = int(PULSE_RUNTIME/FIRINGFILLDT)+1
pulsetsteps = arange(0,PULSE_RUNTIME+1e-6,FIRINGFILLDT)

def receptorFiringRateRandomPulse():
    """ returns a digital and analog (time constant olfactometer_tau) m-sequence.
    for analog m-sequence, y_infinity = 1 when pulse is on, 0 when off.
    Since respiration periodic function is peak normalized,
    resp_pulsewt is typically set low at 1/3.0,
    to get (250 to 300ml/min) / (900ml/min) ratio and
    similar physiological firing rate on convolving with the same kernel:
    the avg respiration < 0.25 times peak.
    Pedestal mean=peak, but binary pulses mean will be half their peak. """
    ##### in Hz
    ## pulse is analog, pulse_steps is digital
    pulse = zeros([len_pulsetime])
    pulse_steps = zeros([len_pulsetime])
    t = PULSE_START
    ## In the NUMBITS-sized binbuffer, set the bits to be randomly True or False
    ## But binbuffer should not have all False. Else shift register will cycle only False-s!
    binbuffer = [False]*NUMBITS
    while binbuffer == [False]*NUMBITS:
        for i in range(NUMBITS):
            if uniform(0,1)>=0.5: binbuffer[i]=True
    ### I use exponential euler as described on pg 334 of Book of Genesis
    D = exp(-FIRINGFILLDT/olfactometer_tau)
    y = 0.0
    PULSE_END = PULSE_START+PULSE_DURATION
    ## create pulse_on and pulse_off pulses of BIT_INTERVAL width
    ## as m-sequence using shift register with xor feedback
    ## as in fig 3.2.3 and table 3.2.1 in Norton's book
    ## odor pulse duration is set at (2^NUMBITS-1)*BIT_INTERVAL
    while t<PULSE_END-BIT_INTERVAL+1e-9: # add a small number to not screw up floating point comparison
        toinsert = binbuffer[3] ^ binbuffer[6] # 4th and 7th bits are xor-ed
        binbuffer.insert(0,toinsert)
        pulse_on = binbuffer.pop() # take the last element off
        ## increase or decay the rate during a pulse
        for subt in range(int(t/FIRINGFILLDT),int((t+BIT_INTERVAL)/FIRINGFILLDT)):
            ## dy/dt = A - B*y
            ## B = 1/olfactometer_tau
            ## A/B = y_infinity = resp_pulsewt
            ## y(t+dt) = y(t)*D + A/B*(1-D)
            if pulse_on:
                y = y*D + (1-D)*resp_pulsewt
                ystep = 1
            else:
                y = y*D
                ystep = 0
            pulse[subt] = y
            pulse_steps[subt] = ystep
        t += BIT_INTERVAL
    ## decay the rate for the rest of the AIR_BUFFER period.
    for subt in range(int(t/FIRINGFILLDT),int((t+PULSE_AIR_BUFFER)/FIRINGFILLDT)):
        y = y*D
        pulse[subt] = y
        pulse_steps[subt] = 0
    return array(pulse),array(pulse_steps)

pulseinresp_width = 50e-3 # ms
pulseinresp_delay = 50e-3 # ms

def receptorFiringRatePulseinResp(kernelR,kernelodor):
    """ returns an analog pulse (time constant olfactometer_tau) shaped by respiration.
    y_infinity = 1 when pulse is on, 0 when off.
    Here, this pulse is in freely breathing animal, not in tracheotomized.
    So, pulse height is 1.0, i.e. not scaled by resp_pulsewt.
    respiration periodic function is peak normalized.
    """
    ## pulse with olfactometer dynamics is created first
    len_resptime = len(extratime)
    pulse = zeros([len_resptime])
    t = 0.0
    tnum = 0
    ### I use exponential euler as described on pg 334 of Book of Genesis
    D = exp(-FIRINGFILLDT/olfactometer_tau)
    y = 0.0
    while t<extratime[-1]+1e-9: # add a small number to not screw up floating point comparison
        tmodresp = (t/RESPIRATION)-floor(t/RESPIRATION)
        if pulseinresp_delay < tmodresp < pulseinresp_delay+pulseinresp_width: pulse_on = True
        else: pulse_on = False
        ## increase or decay the rate during a pulse
        ### I use exponential euler as described on pg 334 of Book of Genesis
        ## dy/dt = A - B*y
        ## B = 1/olfactometer_tau
        ## A/B = y_infinity = 1.0
        ## y(t+dt) = y(t)*D + A/B*(1-D)
        if pulse_on:
            y = y*D + (1-D)
        else:
            y = y*D
        pulse[tnum] = y
        t += FIRINGFILLDT
        tnum += 1

    ##### in Hz
    resp_air = convolve(respirationPulses,kernelR,mode='full')[0:len_resptime]*FIRINGFILLDT
    ## pulses is the odor time course of olfactometer,
    ## rat takes the instantaneous odor in with respiration waveform
    ## Hence just an element-wise multiplication
    pulseinresp = pulse*respirationPulses
    frate = resp_air + pulseConvolve(pulseinresp,kernelodor) + ORN_BGND_MEAN

    ## both are numpy arrays
    return pulse,resp_air+ORN_BGND_MEAN,frate

def pulseConvolve(pulse, kernel):
    ## finally the firing rate by convolution of pulse with kernel
    ## since kernel was calculated with decimated pulse, do the same here
    ## always multiply convolution_result by dt to leave result invariant on decimation.
    frate = convolve(pulse,kernel,mode='full')[:len(pulse)] * FIRINGFILLDT
    return array(frate)

def dblsigmoid_eqns((spread1, center2, spread2), spread2_factor, risetime, duration):
    #spread2 = spread2_factor*spread1
    if SHARP_KERNELS: spread2 = spread1
    eqns = [
    invert_dblsigmoid(0.0, spread1, center2, spread2, 0.9) - \
        invert_dblsigmoid(0.0, spread1, center2, spread2, 0.1) - risetime,
    invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5, risephase=False) - \
        invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5) - duration,
    0.0
    ]
    return eqns

def compute_dblsigmoid_params(risetime, duration, spread2_factor):
    """
    Compute double sigmoid spread1, center2, spread2 for given risetime and duration.
    center1 is taken as 0.0. Initial values of spread1 and center2 for fsolve are risetime, duration.
    """
    return fsolve(dblsigmoid_eqns, [risetime, duration, risetime],
        args=(spread2_factor,risetime,duration))

def compute_odor_params(glomnum, componenttype = ''):
    """ first compute delay, risetime and duration. risetime and delay are correlated.
    then get the double sigmoid params for the generated delay, risetime and duration."""

    while True: # keep generating params, until reasonable ones get generated
        ######### generate delay, risetime and duration

        ### trunc_stdnormal() not to be used since it takes only positive side of the gaussian!
        #x1 = trunc_stdnormal()
        #x2 = trunc_stdnormal()
        #delay = x1*delay_sd + delay_mean
        #### delay and rise_time are correlated normal variables
        ### but since I am generating them as gamma distributed variables,
        ### I cannot use below correlation based generator,
        ### which is valid only for normally distributed variables.
        #risetime = ( delay_risetime_correlation*x1 + \
        #    sqrt(1-delay_risetime_correlation**2)*x2 ) * risetime_sd + risetime_mean
        ## trunc_stdnormal() not to be used since it takes only positive side of the gaussian!
        #duration = trunc_stdnormal()*duration_sd + duration_mean

        ## use Gamma functions for approximating non-negative distributions
        ## correlation between delay and risetime is ignored
        ## mean is k*theta, variance is k*theta^2,
        ## thus: shape = k = mean/theta, scale = theta = sd^2/mean
        shape,scale = gamma_shape_scale(delay_mean,delay_sd)
        delay = gamma(shape,scale)
        shape,scale = gamma_shape_scale(risetime_mean,risetime_sd)
        risetime = delay*risetime_mean/delay_mean#gamma(shape,scale)
        shape,scale = gamma_shape_scale(duration_mean,duration_sd)
        duration = gamma(shape,scale)
        
        ## salient responses for negative stim_rate_seednums
        if stim_rate_seednum<0 and glomnum<3:
            ## salient response params for this glomnum
            ## response_indices selected using seednum at the very top
            if 'B' in componenttype: # for odor B
                response_params = salient_responsesB[responseB_indices[glomnum]]
            else: # for odorA and respiration
                response_params = salient_responsesA[responseA_indices[glomnum]]
            
            if 'e' in componenttype:
                delay = response_params[0]
                risetime = response_params[1]
                duration = response_params[2]
            else: # 'i' in componenttype
                delay = response_params[4]
                risetime = response_params[5]
                duration = response_params[6]

        ######## Forced responses for 3 gloms
        if FORCE_RESPONSES:
            if 'A_e' in componenttype:
                if glomnum == 0:
                    delay = 2.0*delay_mean
                    risetime = 1.5*risetime_mean
                    duration = 1.5*duration_mean
                if glomnum == 1:
                    delay = 0e-3#delay_mean - delay_sd
                    risetime = risetime_mean# - risetime_sd
                    duration = duration_mean - duration_sd
                if glomnum == 2:
                    delay = delay_mean + 3*delay_sd
                    risetime = risetime_mean# + risetime_sd
                    duration = duration_mean# + duration_sd
            elif 'A_i' in componenttype:
                if glomnum == 1:
                    delay = delay_mean + delay_sd
                    risetime = risetime_mean + risetime_sd
                    duration = duration_mean + duration_sd
                if glomnum == 2:
                    delay = delay_mean - delay_sd
                    risetime = risetime_mean - risetime_sd
                    duration = duration_mean - duration_sd
            elif 'B_e' in componenttype:
                if glomnum == 0:
                    delay = delay_mean + delay_sd
                    risetime = risetime_mean + risetime_sd
                    duration = duration_mean + 2*duration_sd
                if glomnum == 1:
                    delay = delay_mean - delay_sd
                    risetime = risetime_mean - risetime_sd
                    duration = duration_mean - duration_sd
                if glomnum == 2:
                    delay = delay_mean + delay_sd
                    risetime = risetime_mean - risetime_sd
                    duration = duration_mean - duration_sd
            elif 'B_i' in componenttype:
                if glomnum == 0:
                    delay = delay_mean - delay_sd
                    risetime = risetime_mean - risetime_sd
                    duration = duration_mean + duration_sd
                if glomnum == 2:
                    delay = delay_mean + delay_sd
                    risetime = risetime_mean + risetime_sd
                    duration = duration_mean + duration_sd

        ######### fit to double sigmoid, get params
        spread2_factor = 4.0 # not used presently
        spread1, center2, spread2 = \
            compute_dblsigmoid_params(risetime, duration, spread2_factor)
        if SHARP_KERNELS: spread2=spread1
        #spread2 = spread2_factor*spread1
        ######### check if the solved params are reasonable
        ######### risetime and duration match and spread2 is at least 1.5*spread1.
        actual_risetime = invert_dblsigmoid(0.0, spread1, center2, spread2, 0.9) - \
            invert_dblsigmoid(0.0, spread1, center2, spread2, 0.1)
        actual_duration = invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5, risephase=False) - \
            invert_dblsigmoid(0.0, spread1, center2, spread2, 0.5)

        ## if fit is not reasonable, complain and repeat loop, else break out
        repeatloop_condn = abs(actual_risetime-risetime)>1e-3 or abs(actual_duration-duration)>1e-3
        if not SHARP_KERNELS: repeatloop_condn = repeatloop_condn or spread2 < 1.5*spread1
        if repeatloop_condn:
            print "The params do not fit well ..."
            print "The actual risetime and duration are", actual_risetime, actual_duration
            print "The required risetime and duration are", risetime, duration
            print "spread2 should be > 1.5*spread1, spread1 =",spread1,", spread2 =",spread2
        else: break
    
    ########## shift curve to include latency
    ## I had taken center1 = 0.0, now shift the curve
    ## so that t=0 is actually where first sigmoid is 0.05*peak
    ## and then add the delay / latency to it!
    offset = - invert_dblsigmoid(0.0, spread1, center2, spread2, 0.05) + delay
    return [offset, spread1, center2+offset, spread2]

def kernelleastsqfunc(kernel, response, respiration):
    ## don't use mode='same': it takes only the valid part of the convolution.
    ## always multiply convolution_result by dt to leave the response invariant on decimation
    responsecalc = convolve( respiration, kernel, mode='full' )[0:len(response)] * fitting_dt
    return [ (val-response[i]) for i,val in enumerate(responsecalc) ]

def getkernel(odorparams_e, odorparams_i):
    """ To get the kernel deconvolve the odor response obtained using Carey's params.
    """
    eresponse = odorparams_e[0] * normresponse(
        [dblsigmoid(t,odorparams_e[1],odorparams_e[2],odorparams_e[3],odorparams_e[4]) for t in extratime] )
    iresponse = odorparams_i[0] * normresponse( 
        [dblsigmoid(t,odorparams_i[1],odorparams_i[2],odorparams_i[3],odorparams_i[4]) for t in extratime] )
    response = eresponse - iresponse
    ############ deconvolve does not work, gives:
    ## return sigtools._linear_filter(b, a, x, axis)
    ## ValueError: BUG: filter coefficient a[0] == 0 not supported yet
    # kernel = scipy.signal.deconvolve(response,respirationPulse1)
    
    ############ Doing rfft -> irfft does not work either :
    ## kernel is slightly non-zero at the start and does not follow response there.
    #response_fft = rfft(response)
    #respiration = respirationPulse1[0:len(response)/2]
    #respirationPulse1_fft = rfft(respirationPulse1,len(response))
    #kernel_fft = response_fft/respirationPulse1_fft
    ## see numpy.fft.irfft documentation for why len(response) is required.
    #kernel = irfft(kernel_fft,len(response))

    ############# FINALLY, obtain the kernel by fitting!
    ## initial value of kernel is the response itself.
    ## Do not fit for all the points - too slow.
    ## only fit for number of bins in kernel i.e. decimate
    response = response[::binsize]
    ## The respiration is peak normalized!
    respiration = respirationPulse1[::binsize]
    kernel_info = scipy.optimize.leastsq( kernelleastsqfunc,\
        response, args=(response,respiration), full_output=1 )
    kernel = kernel_info[0]
    #figure()
    #plot(kernel,'rx-')
    #plot(respiration,'go-')
    #plot(convolve( respiration, kernel, mode='full' )[0:len(response)]*fitting_dt,'b+-')
    #plot(response,'kx',linestyle='dashed')
    print 'message :',kernel_info[3]
    ## kernel was a decimated one, now interpolate to have dt = FIRINGFILLDT
    kernel_fn = scipy.interpolate.interp1d(extratime_decimated, kernel, kind='cubic')
    kernel_full = [kernel_fn(t) for t in extratime]
    return array(kernel_full)

def get_high_frate(cutoffrate,meanfrate,glomnum):
    #if glomnum == 0: None
    #else: cutoffrate = meanfrate
    fraterand = 0
    while fraterand<cutoffrate:
        fraterand = exponential(meanfrate)
    return fraterand

def compute_Guassian_params(odor,glomnum):
    while True:
        center = normal(DoG_center_mean,DoG_center_sd)
        if center>0: break ## allow only positive values for Gaussian's center
    if glomnum==0: thinfactor = 1.0
    else: thinfactor = 1.0#0.75
    if odor=='R':
        width = uniform(DoG_airwidth_min*thinfactor,DoG_airwidth_max*thinfactor)
    else:
        width = uniform(DoG_width_min*thinfactor,DoG_width_max*thinfactor)
    return center, width

def get_DoG_kernel(odor,factor,glomnum):
    if odor=='A': fratemean = FIRINGMEANA
    elif odor=='B': fratemean = FIRINGMEANB
    else: fratemean = FIRINGRATEAIR_MEAN

    if glomnum==0: # central glom
        ecenter, ewidth = compute_Guassian_params(odor,glomnum)
        icenter, iwidth = compute_Guassian_params(odor,glomnum)
        if DoG_INHON:
            ### older exp distribution, with peak2peak and upper cutoffs after returning from this function
            #epeak = factor * get_high_frate(fratemean,fratemean,glomnum) * kernel2ORNfrate_ratio
            ## exc+inh kernels uniform distrib
            epeak = factor * uniform(fratemean,2.5*fratemean) * kernel2ORNfrate_ratio
            ### older exp distribution, with peak2peak and upper cutoffs after returning from this function
            #ipeak = factor * get_high_frate(fratemean/2.0,fratemean,glomnum) * kernel2ORNfrate_ratio
            ## exc+inh kernels uniform distrib
            ipeak = factor * uniform(0.0*fratemean,0.75*fratemean) * kernel2ORNfrate_ratio
        else:
            ## purely exc kernels
            #epeak = factor * uniform(0.2*fratemean,1.5*fratemean) * kernel2ORNfrate_ratio
            epeak = factor * uniform(0.4*fratemean,1.5*fratemean) * kernel2ORNfrate_ratio # new pg firing regime
            ## purely exc kernels
            ipeak = 0.0
        DoGkernel = subtractedAlphaGaussians(\
            extratime, epeak, ipeak, ecenter, icenter, ewidth/2.0, iwidth/2.0)
    else: # lateral glom
        while True:
            ecenter, ewidth = compute_Guassian_params(odor,glomnum)
            icenter, iwidth = compute_Guassian_params(odor,glomnum)
            #epeak = factor * uniform(0.2*fratemean,1.5*fratemean) * kernel2ORNfrate_ratio
            epeak = factor * uniform(0.4*fratemean,1.5*fratemean) * kernel2ORNfrate_ratio # new pg firing regime
            if DoG_INHON:
                ipeak = factor * uniform(0.2*fratemean,1.5*fratemean) * kernel2ORNfrate_ratio
            else:
                ipeak = 0.0 # purely exc input
            DoGkernel = subtractedAlphaGaussians(\
                extratime, epeak, ipeak, ecenter, icenter, ewidth/2.0, iwidth/2.0)
            if DoG_INHON:
                ## exc+inh lat i/p -- ensure reasonable undulations/excursions amplitude to get 2nd order effects
                if (max(DoGkernel)-min(DoGkernel))>1.0*fratemean*kernel2ORNfrate_ratio: break
            else: break
            ### below is obsolete
            #if max(DoGkernel)>1.0*fratemean*kernel2ORNfrate_ratio and \
            #    min(DoGkernel)<-0.5*fratemean*kernel2ORNfrate_ratio: break
    return DoGkernel

def odor_and_air_stimuli():
    # mitral and PG odor ORNs firing files

    ## Priyanka's random pulse sequences: should be the same for all glomeruli.
    ## The kernels for each glom are different, but ORN response is obtained by
    ## convolving different kernels with the same random pulse sequence.
    ## randomPulseList[pulsenum]
    ## if you redefine randomPulseList = [...] here, it will exist only locally.
    ## first a rectangular air pulse
    ## receptorFiringRateRespPulse() uses resp_pulsewt.
    randomPulse,randomPulseSteps = receptorFiringRateRespPulse()
    randomPulseList.append( randomPulse )
    randomPulseStepsList.append( randomPulseSteps )
    ## an air random pulse, then RANDOM_PULSE_NUMS number of odor A and odor B random pulses.
    for i in range(2*RANDOM_PULSE_NUMS+1):
        randomPulse,randomPulseSteps = receptorFiringRateRandomPulse()
        randomPulseList.append( randomPulse )
        randomPulseStepsList.append( randomPulseSteps )

    for glomnum in range(NUM_GLOMS_MAX):
        print "Computing firing rates for odor morphs and pulses for glomerulus", glomnum
        Afactor,Bfactor = GLOMS_ODOR_weights[glomnum]

        if not DoG_KERNELS:
            ### For each glomerulus, set its responses to odor A and B
            ### as difference of excitatory and inhibitory random dual-exponentials
            ### the inhibitory component should be activated only for 50% of odors
            ### thus, 50+50/2=75% of responses will be excitatory and 50/2=25% inhibitory
            ### thus, inh:exc = 25:75 = 20:60 as in Savigner et al 2009.
            ### peak firing rate for exc and inh is distributed as a falling exponential
            ### with 10 spikes/s mean (figs of Duchamp-Viret et al 2000)

            if SHARP_KERNELS: # higher exc as large inh, also better kernel fit
                odorparamsA_e = [Afactor * get_high_frate(FIRINGMEANA,FIRINGMEANA,glomnum)]
            else:
                odorparamsA_e = [Afactor * get_high_frate(FIRINGMEANA/2.0,FIRINGMEANA,glomnum)]
            if FORCE_RESPONSES and glomnum in [0,1,2]: odorparamsA_e = [FIRINGMEANA]
            elif stim_rate_seednum<0 and glomnum in [0,1,2]: # salient odorA
                odorparamsA_e = [ FIRINGMEANA*ampl_scale*salient_responsesA[responseA_indices[glomnum]][3] ]
            odorparamsA_e.extend( compute_odor_params(glomnum, 'A_e') )
            if SHARP_KERNELS: ifactor = 1.0
            else:
                if uniform(0.0,1.0)<0.5: ifactor = 1.0
                else: ifactor = 0.0
            if SHARP_KERNELS: # all of Priyanka's mitral kernels have inh components, so force a minimum
                odorparamsA_i = [Afactor * ifactor * get_high_frate(FIRINGMEANA/4.0,FIRINGMEANA,glomnum)]
            else: odorparamsA_i = [Afactor * ifactor * exponential(FIRINGMEANA)]
            if FORCE_RESPONSES and glomnum in [0,1,2]: odorparamsA_i = [0.0]
            elif stim_rate_seednum<0 and glomnum in [0,1,2]:
                odorparamsA_i = [ FIRINGMEANA*ampl_scale*salient_responsesA[responseA_indices[glomnum]][7] ]
            odorparamsA_i.extend( compute_odor_params(glomnum, 'A_i') )

            ### kernel for Priyanka can be obtained by numerically deconvolving with the respiration pulse.
            kernelA = getkernel(odorparamsA_e,odorparamsA_i)

            if SHARP_KERNELS: # higher exc as large inh, also better kernel fit
                odorparamsB_e = [Bfactor * get_high_frate(FIRINGMEANB,FIRINGMEANB,glomnum)]
            else:
                odorparamsB_e = [Bfactor * get_high_frate(FIRINGMEANB/2.0,FIRINGMEANB,glomnum)]
            if FORCE_RESPONSES and glomnum in [0,1,2]: odorparamsB_e = [FIRINGMEANB]
            elif stim_rate_seednum<0 and glomnum<3: # salient odorB
                odorparamsB_e = [ FIRINGMEANB*ampl_scale*salient_responsesB[responseB_indices[glomnum]][3] ]
            if stim_rate_seednum<=-19 and CONST_GLOM0 and glomnum<3: odorparamsB_e = [ 0.0 ] # const odorA, rest zero
            odorparamsB_e.extend( compute_odor_params(glomnum, 'B_e') )
            if SHARP_KERNELS: ifactor = 1.0
            else:
                if uniform(0.0,1.0)<0.5: ifactor = 1.0
                else: ifactor = 0.0
            if FORCE_RESPONSES and glomnum in [1]: ifactor = 0.0
            if SHARP_KERNELS: # all of Priyanka's mitral kernels have inh components, so force a minimum
                odorparamsB_i = [Bfactor * ifactor * get_high_frate(FIRINGMEANB/4.0,FIRINGMEANB,glomnum)]
            else: odorparamsB_i = [Bfactor * ifactor * exponential(FIRINGMEANB)]
            if FORCE_RESPONSES and glomnum in [0,2]: odorparamsB_i = [FIRINGMEANB]
            elif stim_rate_seednum<0 and glomnum<3: # salient odorB
                odorparamsB_i = [ FIRINGMEANB*ampl_scale*salient_responsesB[responseB_indices[glomnum]][7] ]
            if stim_rate_seednum<=-19 and CONST_GLOM0 and glomnum<3: odorparamsB_i = [ 0.0 ] # const odorA, rest zero
            odorparamsB_i.extend( compute_odor_params(glomnum, 'B_i') )

            kernelB = getkernel(odorparamsB_e,odorparamsB_i)

            ### respiration also has both excitatory and inhibitory dual exponentials.
            #if FORCE_RESPONSES and glomnum in [0]: odorparamsR_e = [0.0] # zero exc resp
            odorparamsR_e = [get_high_frate(FIRINGRATEAIR_MEAN/2.0,FIRINGRATEAIR_MEAN,glomnum)]
            if FORCE_RESPONSES and glomnum in [0,1,2]: odorparamsR_e = [FIRINGRATEAIR_MEAN]
            elif stim_rate_seednum<0 and glomnum<3:
                odorparamsR_e = [ 2*FIRINGRATEAIR_MEAN*salient_responsesA[responseA_indices[glomnum]][3] ]
            odorparamsR_e.extend( compute_odor_params(glomnum,'R_e') )
            if SHARP_KERNELS: ifactor = 1.0
            else:
                if uniform(0.0,1.0)<0.5: ifactor = 1.0
                else: ifactor = 0.0
            if FORCE_RESPONSES and glomnum in [0,1,2]: odorparamsR_i = [0.0] # zero inh resp
            elif stim_rate_seednum<0 and glomnum<3:
                odorparamsR_i = [ 2*FIRINGRATEAIR_MEAN*salient_responsesA[responseA_indices[glomnum]][7] ]
            elif SHARP_KERNELS: # all of Priyanka's mitral kernels have inh components, so force a minimum
                odorparamsR_i = [ifactor * get_high_frate(FIRINGRATEAIR_MEAN/8.0,FIRINGRATEAIR_MEAN,glomnum)]
            else: odorparamsR_i = [ifactor * exponential(FIRINGRATEAIR_MEAN)]
            odorparamsR_i.extend( compute_odor_params(glomnum,'R_i') )

            kernelR = getkernel(odorparamsR_e,odorparamsR_i)

            if stim_rate_seednum<=-19 and not CONST_GLOM0:
                if glomnum==0:
                    kernelR = kernelR/1.5
                    kernelA = kernelA/1.5
                    kernelB = kernelB/1.5
                else:
                    kernelR = kernelR
                    kernelA = kernelA*1.5
                    kernelB = kernelB*1.5

        else: # DoG (difference of Gaussians) kernels
            ####################### Generate kernels first, not retro-fit them
            
            ## air and air+odor must have minimum positive mean for central glom.
            ## min 1Hz mean mitral response rate, so 0.2Hz ORN frate,
            ## so 0.2*kernel2ORNfrate_ratio mean kernel value;
            ## also factor 0.5 since long blank post-padded kernel.
            ## All side glomeruli should also have some response, since I simulate
            ## only odor-responsive glomeruli that will cause lateral inhibition.
            ## Note that mean is over full kernel_time which is 1.5s or 2s,
            ## so don't set too high lowcutoff for mean.
            lowcutoff_odorkernelmean = 0.0#0.45*kernel2ORNfrate_ratio # 0.5* is for exc+inh kernels
            lowcutoff_airkernelmean = 0.0
            #lowcutoff_odorkernelampl = 5*kernel2ORNfrate_ratio # min 5Hz ORN odor amplitude (max-min)
            #highcutoff_odorkernelampl = 9*kernel2ORNfrate_ratio # max 9Hz ORN odor amplitude (max-min)
            #highcutoff_airkernelampl = 4*kernel2ORNfrate_ratio # max 4Hz ORN air amplitude (max-min)
            #highcutoff_odorkernelmean = 6*kernel2ORNfrate_ratio # mean < 6Hz ORN odor amplitude
            #highcutoff_airkernelmean = 3*kernel2ORNfrate_ratio # mean < 3Hz ORN air amplitude
            while True:
                kernelA = get_DoG_kernel('A',Afactor,glomnum)
                kernelB = get_DoG_kernel('B',Bfactor,glomnum)
                kernelR = get_DoG_kernel('R',1.0,glomnum)
                meanAR = mean(kernelA+kernelR)
                meanBR = mean(kernelB+kernelR)
                meanR = mean(kernelR)
                #peak2peakA = (max(kernelA)-min(kernelA))
                #peak2peakB = (max(kernelB)-min(kernelB))
                #peak2peakR = (max(kernelR)-min(kernelR))
                if ( meanAR < lowcutoff_odorkernelmean or meanBR < lowcutoff_odorkernelmean \
                    or meanR < lowcutoff_airkernelmean ):
                    continue
                else: break

        kernels.extend( (kernelR,kernelA,kernelB) )

        ## Adil's odor morph responses
        frateperglomList = []
        for (odorA, odorB) in inputList:
            #if glomnum == 0:
            #    odorA = 1.0/(1+exp(-(odorA-0.5)/0.1))
            #    odorB = 1.0/(1+exp(-(odorB-0.5)/0.1))
            #print glomnum,odorA, odorB
            frate = receptorFiringRate(respwt, odorA, odorB, kernelA, kernelB, kernelR, glomnum)
            frateperglomList.append(frate)
        ## frateOdorList[glomnum][inputnum]
        frateOdorList.append(frateperglomList)

        ## Priyanka's paired pulse responses
        fratepulseperglomList = []
        for (pulsedelayA,pulsedurationA,pulsedelayB,pulsedurationB) in pulseList:
            frate = receptorFiringRatePulse(\
                pulsedelayA, pulsedurationA, pulsedelayB, pulsedurationB,\
                kernelA, kernelB, kernelR)
            fratepulseperglomList.append(frate)
        ## fratePulseList[glomnum][0-4][timeidx]
        fratePulseList.append((fratepulseperglomList,kernelA, kernelB, kernelR))

        ## Priyanka's random pulse responses
        ## randomResponseList[glomnum][frate,...]
        ## for the stimuli below -18, if CONST_GLOM0, I choose glom0's response as constant odorA
        if stim_rate_seednum<=-19 and CONST_GLOM0 and glomnum==0:
            ## glom0 is constant odorA, zero air responses!
            per_glom_random_pulse_list = [randomPulseList[0]*0.0]
            per_glom_random_pulse_list.append( randomPulseList[1]*0.0 )            
            for i in range(RANDOM_PULSE_NUMS):
                per_glom_random_pulse_list.append( randomPulseList[2*(i+1)]*(0.5*FIRINGMEANA*ampl_scale) )
                per_glom_random_pulse_list.append( randomPulseList[2*(i+1)+1]*0.0 )
        else:
            per_glom_random_pulse_list = [pulseConvolve(randomPulseList[0], kernelR)]
            per_glom_random_pulse_list.append( pulseConvolve(randomPulseList[1], kernelR) )
            for i in range(RANDOM_PULSE_NUMS):
                per_glom_random_pulse_list.append( \
                                pulseConvolve(randomPulseList[2*(i+1)], kernelA)*CONC_SCALE )
                per_glom_random_pulse_list.append( \
                                pulseConvolve(randomPulseList[2*(i+1)+1], kernelB)*CONC_SCALE )
        randomResponseList.append(per_glom_random_pulse_list)
        
        ## pulseinresp: Predict kernel from pulse in respiration of freely-breathing animal?
        pulses,frateR,frate = receptorFiringRatePulseinResp(kernelR,kernelA)
        per_glom_pulseinresp = [frateR]
        per_glom_pulseinresp.append(frate)
        per_glom_pulseinresp.append(receptorFiringRatePulseinResp(kernelR,kernelB)[2])
        ## Saving pulses for every glom is redundant, yet programmatically easier
        ## pulseInRespList[glomnum][0|1][pulsesbins | frateR,frateA,frateB][nodimpresent | ratebinnum]
        pulseInRespList.append((pulses,per_glom_pulseinresp))

        ## scaled pulses for linearity vs concentration
        fratescaledpulseperglomList = []
        for scale in scaledList:
            frateA = receptorFiringRateScaledPulse(scale, 0.0, kernelA, kernelB, kernelR)
            frateB = receptorFiringRateScaledPulse(0.0, scale, kernelA, kernelB, kernelR)
            fratescaledpulseperglomList.append((frateA,frateB))
        scaledPulses.append(fratescaledpulseperglomList)

def generate_firerates(rateseed):
    ## Seed numpy's random number generator.
    ## If no parameter is given, it uses current system time
    seed([rateseed])

    odor_and_air_stimuli()

    if 'PULSEINRESP' in sys.argv:
        ## write pulse in resp separately
        filename = '../generators/firerates/firerates_pulseinresp_'+str(rateseed)
        filename += '.pickle'
        fireratefile = open(filename,'w')
        pickle.dump(  \
            (pulseInRespList,kernels), \
            fireratefile)
        fireratefile.close()
        print "wrote",filename
    elif 'SCALEDPULSES' in sys.argv:
        ## write pulse in resp separately
        filename = '../generators/firerates/firerates_scaledpulses_width'+\
            str(scaledWidth)+'_'+str(rateseed)
        filename += '.pickle'
        fireratefile = open(filename,'w')
        pickle.dump(  \
            (scaledPulses,kernels), \
            fireratefile)
        fireratefile.close()
        print "wrote",filename
    else:
        filename = '../generators/firerates/firerates_2sgm_'+str(rateseed)
        filename += '.pickle'
        fireratefile = open(filename,'w')
        pickle.dump(  \
            (frateOdorList,fratePulseList,randomPulseList, \
            randomPulseStepsList,randomResponseList,kernels), \
            fireratefile)
        fireratefile.close()
        print "wrote",filename

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    generate_firerates(stim_rate_seednum)

    ################# Only plotting below this ##################
    if 'NOSHOW' not in sys.argv:
        ## Only plot gloms 0,1,2 which should be interesting enough to simulate with
        for glomi in range(min(NUM_GLOMS_MAX,3)):
            figure(facecolor='w')
            title('Glomerulus '+str(glomi))
            xlabel('time (s)', fontsize='large')
            ylabel('firing rate (Hz)', fontsize='large')
            for odoridx,(odorA, odorB) in enumerate(inputList):
                frate = frateOdorList[glomi][odoridx]
                plot(firingtsteps, frate, color=(odorA,odorB,0), marker=',')
            #twinx()
            plot(extratime[:numt], respirationPulses[startidx:startidx+numt], color=(0,0,1),marker=',')

        fig1 = figure(facecolor='w')
        title('Glomerulus 0 : Pulse Responses')
        xlabel('time (s)', fontsize='large')
        ylabel('firing rate (Hz)', fontsize='large')
        ax1 = fig1.add_subplot(111)
        lenpulselist = len(pulseList)
        pulsetime = arange(0.0,PULSE_DURATION+1e-10,FIRINGFILLDT)
        for pulseidx,(pulsedelayA,pulsedurationA,pulsedelayB,pulsedurationB) in enumerate(pulseList):
            frate = fratePulseList[0][0][pulseidx]
            colorfrac = pulseidx/float(lenpulselist)
            pulseA = zeros([len(pulsetime)])+pulseidx*10
            pulseA[int(pulsedelayA/FIRINGFILLDT):int((pulsedelayA+pulsedurationA)/FIRINGFILLDT)]=pulseidx*10+9
            ax1.plot(pulsetime, pulseA, color=(colorfrac,1-colorfrac,1), marker=',')
            ax1.plot(pulsetime, frate, color=(colorfrac,1-colorfrac,0), marker=',')

        fig2 = figure(facecolor='w')
        ax2 = fig2.add_subplot(111)
        title('Glomerulus 0 : Kernels')
        xlabel('time (s)', fontsize='large')
        ylabel('firing rate (Hz)', fontsize='large')
        ax2.plot(extratime, fratePulseList[0][1], color=(1,0,0),
            marker=',',linestyle='solid',label='kernelA')
        ax2.plot(extratime, fratePulseList[0][2], color=(0,1,0),
            marker=',',linestyle='solid',label='kernelB')
        ax2.plot(extratime, fratePulseList[0][3], color=(0,0,0),
            marker=',',linestyle='solid',label='kernelR')
        legend()

        fig3 = figure(facecolor='w')
        ax3 = fig3.add_subplot(111)
        title('Respiration waveform')
        xlabel('time (s)', fontsize='large')
        ylabel('firing rate (Hz)', fontsize='large')
        ax3.plot(extratime, respirationPulses, color=(1,0,0),
            marker=',',linestyle='solid',label='multiple respirations')
        ax3.plot(arange(0,2*RESPIRATION+1e-10,FIRINGFILLDT), respirationPulse1, color=(0,0,1),
            marker=',',linestyle='solid',label='single respirations')
        legend()
        
        for glomnum in range(3):#NUM_GLOMS_MAX):
            figure(facecolor='w')
            title('Glomerulus '+str(glomnum)+' : Random Pulse Responses - odorA')
            xlabel('time (s)', fontsize='large')
            ylabel('firing rate (Hz)', fontsize='large')
            ### odor pulse + air pulse + background (adding different units! - just schematic)
            plot(randompulsetime, randomPulseList[2]+randomPulseList[0]+ORN_BGND_MEAN,\
                color=(1,0,0), marker=',',linestyle='solid',label='pulseA')
            plot(randompulsetime, randomPulseList[3]+randomPulseList[0]+ORN_BGND_MEAN,\
                color=(1,0,0.5), marker=',',linestyle='solid',label='pulseB')
            #twinx()
            ### odor response + air pedestal response + background
            plot(randompulsetime, randomResponseList[glomnum][2]+randomResponseList[glomnum][0]+ORN_BGND_MEAN,\
                color=(0,1,0), marker=',',linestyle='solid',label='ORN responseA')
            plot(randompulsetime, randomResponseList[glomnum][3]+randomResponseList[glomnum][0]+ORN_BGND_MEAN,\
                color=(0,1,0.5), marker=',',linestyle='solid',label='ORN responseB')
            legend()

            fig2 = figure(facecolor='w')
            ax2 = fig2.add_subplot(111)
            title('Glomerulus '+str(glomnum)+' : Kernels')
            xlabel('time (s)', fontsize='large')
            ylabel('arb units', fontsize='large')
            ax2.plot(extratime, fratePulseList[glomnum][1], color=(1,0,0),
                marker=',',linestyle='solid',label='kernelA')
            ax2.plot(extratime, fratePulseList[glomnum][2], color=(0,1,0),
                marker=',',linestyle='solid',label='kernelB')
            ax2.plot(extratime, fratePulseList[glomnum][3], color=(0,0,0),
                marker=',',linestyle='solid',label='kernelR')
            legend()

            fig3 = figure(facecolor='w')
            ax3 = fig3.add_subplot(111)
            title('Glomerulus '+str(glomnum)+' : PulseInResp')
            xlabel('time (s)', fontsize='large')
            ylabel('arb units', fontsize='large')
            ax3.plot(extratime, respirationPulses, color=(1,0,1),
                marker=',',linestyle='solid',label='Respiration')
            ax3.plot(extratime, pulseInRespList[glomnum][0], color=(0,1,0),
                marker=',',linestyle='solid',label='Pulse')
            ax3.plot(extratime, pulseInRespList[glomnum][1][0], color=(0,0,0),
                marker=',',linestyle='solid',label='OdorA')
            ax3.plot(extratime, pulseInRespList[glomnum][1][1], color=(1,0,0),
                marker=',',linestyle='solid',label='OdorA')
            ax3.plot(extratime, pulseInRespList[glomnum][1][2], color=(0,0,1),
                marker=',',linestyle='solid',label='OdorB')
            legend()

        show()
