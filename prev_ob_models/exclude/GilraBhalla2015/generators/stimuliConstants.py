# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

from stimuliConstantsMinimal import *

## Priyanka's kernels are sharper with all of them having inh components.
## Both below should not be True, only none or one should be set to True.
## (SHARP_KERNELS is used at multiple places, never combine it with DoG_KERNELS.)
SHARP_KERNELS = False # sharper kernels but using retro-fitting of Carey et al type responses.
DoG_KERNELS = True # Use difference of Gaussians kernels, no reference to Carey et al.
DoG_INHON = False # whether the inh component is on/off in the DoG kernels, default off/False.

## generate stimuli for NUM_GLOMS_MAX glomeruli. If you increase NUM_GLOMS_MAX,
## the stimuli will change even for the glomnums of old NUM_GLOMS_MAX (random number generation order).
## though in generate_firerates_odors.py, I have tried to ensure that this doesn't happen
## by generating random numbers irrespective of switches. Still to be on the safe side...
NUM_GLOMS_MAX = 10#20
## which glom responds to which odor: for each glom weights (odorA, odorB)
## be careful when you set NUM_GLOMS = 3. You then need to permute below for all possible permutations.
## Glom0 responds to both odors as that is how Adil selected his cells.
## Should not give both odors to the others, as they may have mixture interactions.
## Actually gloms 2 and 3 can each be considered a proxy
## for two different gloms with different odorA and odorB responses.
## There could be mixture interactions though for the intermediate morphs!
## Am allowing each glom to respond to both odours, Each glom can be considered to be
## different gloms responding to one odor each, but at roughly the same distance!
## Basically, I only want to see how decorr and linearity changes with number of glomeruli.
GLOMS_ODOR_weights = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),
    (1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]

# time for the cells to settle
SETTLETIME = 0.25 # seconds
RESPIRATION = 0.5 # s respiration cycle of anaesthetized rat is roughly 1 s.
NUM_RESPS = 2
ODORRUNTIME = NUM_RESPS*RESPIRATION + SETTLETIME # seconds

## Scale odor responses by this amount. 1.0 corresponds to 1% SV.
CONC_SCALE = 1.0 # default is 1.0 for 1% saturated vapour.
# inputList = [ (odorA,odorB) , ... ]
odorAon = True
odorBon = True
if odorAon and odorBon:
    inputList = [ (0.0,1.0), (0.2,0.8), (0.4,0.6), (0.6,0.4), (0.8,0.2), (1.0,0.0), (0.0,0.0) ]
elif odorAon:
    inputList = [ (0.0,0.0), (0.2,0.0), (0.4,0.0), (0.6,0.0), (0.8,0.0), (1.0,0.0), (0.0,0.0) ]
elif odorBon:
    inputList = [ (0.0,1.0), (0.0,0.8), (0.0,0.6), (0.0,0.4), (0.0,0.2), (0.0,0.0), (0.0,0.0) ]
# pulseList = [ (delayA,durationA,delayB,durationB) , ... ]
pulseList = [(100e-3,50e-3,500e-3,0), (100e-3,100e-3,500e-3,0), (100e-3,200e-3,500e-3,0),\
    (100e-3,500e-3,500e-3,0), (100e-3,1000e-3,500e-3,0), (100e-3,2000e-3,500e-3,0)]
# for pulses of different concentrations, 0.0 is needed as the air response.
scaledList = [0.0,1/3.0,2/3.0,1.0,2.0,5.0] # concentration scalings
#numneighbours = 1 # number of mitral cells connected to recorded mitral
NUMBITS = 7 ## Can't change this! The bits to be xor-ed are hard-coded in generate_firerates_odors.py
BIT_INTERVAL = 50e-3
PULSE_DURATION = BIT_INTERVAL*(2.0**NUMBITS-1) # 6.35 seconds # 7-bit m-sequence
SCALED_DURATION = 3 # seconds, duration of scaled pulses to test linearity
PULSE_AIR_BUFFER = 1.0 # seconds # on both sides of the random pulse
PULSE_START = SETTLETIME+PULSE_AIR_BUFFER # seconds
PULSE_RUNTIME = PULSE_START+PULSE_DURATION+PULSE_AIR_BUFFER # seconds
SCALED_RUNTIME = PULSE_START+SCALED_DURATION
if DoG_KERNELS: olfactometer_tau = 40e-3 # seconds
else: olfactometer_tau = 50e-3 # seconds
#olfactometer_tau_off = 70e-3 # seconds # taking same tau-s for the olfactometer
# number of random pulse stimuli
RANDOM_PULSE_NUMS = 3
pulsebindt = 0.1 # seconds
pulsebins = int((PULSE_RUNTIME/pulsebindt)) ## pulsebins scale with runtime

## max number of trials for morphs or pulses
## used to generate ORN firefiles
## Do average over more than these many trials as these many stimuli only.
MAXNUMAVG = 10
## Each (odornum,trialnum) pair has a different set of granule baseline firefiles.
MAXNUMAVG_GRANS = 65

## granule non-respiratory tuned background files are present only up to MAXRUNTIME.
MAXRUNTIME = 30.0 #s
## respiratory tuned are present upto 1 resp cycle.

############## Try to keep FIRINGFILLDT a multiple of clock 0 i.e. SIMDT
## the time step for which one must evaluate firing rate to fill the kernel, odor waveform, rate vector, etc.
FIRINGFILLDT = 1e-3 # s
## 35Hz EPSPs in granule cells at baseline Cang & Isaacson 2003
## translates roughly to average 0.35Hz for each of the 100 mitrals into a granule.
mit_base_f = 0.35 # Hz
## 3.45Hz EPSPs in granule cells in vitro (in mice Carleton et al)
## translates roughly to average 0.0345Hz for each of the 100 mitrals into a granule.
mit_base_f_invitro = 0.0345 # Hz
num_gran_baseline_files = 1000
## exc factor instead of 1.0 to ALL granules during odor period if CLUB_MITRALS = False in simset_odor.py
## based on Ashesh's and Adil's air vs odor response averages (though they sampled high firing cells)
## But these are created once only -- just switch between firefiles/firefiles_baseline_1.5x and _1.0x
extraexc_factor = 1.0#1.5

######## ORN input parameters
NUM_ORN_FILES_PER_GLOM = 1000
respwt=1
## resp_pulsewt is the scaling in peak airflowrate in pulses compared to respiration.
## It is constrained by experiment, and should be kept at 1/3.0.
## larger exc and inh components and sharper kernels, keep at 1/2.0 -- OBSOLETE.
if SHARP_KERNELS: resp_pulsewt=1.0/2.0#0.25
else: resp_pulsewt=1/3.0#0.25
respbins=100
REFRACTORY = 1.0e-3 # s refractory period
################ Constant baseline:
# arbitrarily put in a baseline to ensure that excitatory responses are favoured over inhibitory responses.
#BASELINEORNFIRING = 1.25 # spikes/s # not used for the dual exp form of virtual ORNs

################# Odor rates:
# Duchamp Viret et al 2000, mean firing rate for odors (from figures) is <10 Hz. Peaks are from 15 to 30 Hz.
# Double sigmoid can be approximated as a rectangle of half the peak height.
# So take mean of the distribution as twice the mean firing rate.
# Using FIRINGRATEA = mean*4 = 10*4 = 40Hz, I got very high 70Hz pegged firing of mitral in vivo.
# Reduce the mean firing rates by half, since for Priyanka's stimuli, the peak rates go double that for Adil's!
FIRINGMEANA = 6.0#10/2.0 # spikes/s
FIRINGMEANB = 6.0#10/2.0 # spikes/s

################# Air rates:
# Duchamp Viret et al 2000: peak = mean*2 (double sigmoid - mean~=peak/2); mean = 1.9Hz
# but Duchamp-Viret 2005 reported 1.25Hz, so I take 1.5Hz
#FIRINGRATERESP = 1.5*2 # spikes/s typical firing rate for respiration
# Exponential distribution with below mean is chosen to correspond to fig 1 of Duchamp Viret et al 2000.
## Now using uniform distribution as we are actually using mean of 400 ORNs,
## not the exp distributed frates of individual ORNs given a mean
## With new _mod_withspikeinit mitral use 1.5Hz air mean; for default _mod mitral, use 2Hz.
FIRINGRATEAIR_MEAN = 2.0#1.5 # spikes/s

################# Constant background rate to have inhibition on mitrals from beginning
ORN_BGND_MEAN = 0.5 # spikes/s

## obsolete - i now use 2*maxfiringrate to set the max rate before spike thinning.
#MAXFIRINGRATE = 4.0*FIRINGMEANA # spikes/s - maximum attainable firing rate for spike thinning

NONLINEAR_CENTER = 1.2*FIRINGMEANA # 7.2Hz
NONLINEAR_WIDTH = FIRINGMEANA/FIRINGMEANA # 1Hz

########## ORN waveforms and respiration function parameters:
#respiration: Values set by me to be similar to Carey et al 2009 Fig 1:
# RESPIRATION is respiration period (1s).
tau1rins = 0.06*RESPIRATION
tau2rins = 0.05*RESPIRATION
tau1rinsslow = 0.3*RESPIRATION
tau2rinsslow = 0.1*RESPIRATION
tau1rexp = 0.03*RESPIRATION
tau2rexp = 0.025*RESPIRATION

### Values of spreads taken from Carey et al (Wachowiak lab), 2009
### I should not take the values from the paper as values for the kernels (for Priyanka's experiments),
### rather, I should deconvolve numerically with respiration to get the kernel.

## delay - table 1 of Carey et al 2009: -- this is only for one rat, but table 2 has similar values
## 154 +- 59 ms from 'inspiration start' to 'max of slope times curvature of response'
delay_mean = 154.0e-3
delay_sd = 59.0e-3
## rise-time - table 2 (anaethetized case): 122 +- 32 ms for 10% to 90% rise time
## Cannot use Gaussian for this - it has a long tail!
## but it is correlated with delay; and the easy formula
## for correlated random variables is for Gaussian random variables.
risetime_mean = 121e-3 # should have been 122e-3 from table 2 - typo - but doesn't matter much
risetime_sd = 32e-3
delay_risetime_correlation = 0.3 # delay and rise time are positively correlated.
## duration - table 2 (anaesthetized case): 443 +- 119 ms for response above 50% of peak
if SHARP_KERNELS:
    duration_mean = 225e-3
else: #default Carey et al param
    duration_mean = 443e-3
duration_sd = 119e-3
## shape and scale of gamma functions are decided by the mean and SD.
#delay_shape = 10 # shape of gamma function estimated by eye from fig 4 of Carey et al 2009
#risetime_shape = 2 # shape of gamma function estimated by eye from fig 4 of Carey et al 2009
#duration_shape = 5 # shape of gamma function estimated by eye from fig 4 of Carey et al 2009

## force shape of first 3 gloms' responses rather than generate randomly
FORCE_RESPONSES = False
## FORCE_RESPONSES above is different from the salient responses below
## broad, early narrow and late narrow are the salient responses
## salient_responses = [ (exc_delay,exc_risetime,exc_duration,exc_ampl_factor,\
##    inh_delay,inh_risetime,inh_duration,inh_ampl_factor), ... ]
## even if ampl_factor is zero for the first response's inh component,
## can't have zero delay,risetime,duration, else fitting routine pukes.
salient_responsesA = [\
    (delay_mean-delay_sd,risetime_mean-risetime_sd,duration_mean+2*duration_sd,0.15,\
        delay_mean,risetime_mean,duration_mean,0.0),\
    (delay_mean-delay_sd,risetime_mean-risetime_sd,duration_mean,0.6,\
        delay_mean+2*delay_sd,risetime_mean+2*risetime_sd,duration_mean,0.3),\
    (delay_mean+2*delay_sd,risetime_mean+2*risetime_sd,duration_mean+2*duration_sd,0.8,\
        delay_mean,risetime_mean,duration_mean+duration_sd,0.8)\
    ]
salient_responsesB = [\
    (delay_mean+delay_sd,risetime_mean+risetime_sd,duration_mean+duration_sd,0.25,\
        delay_mean,risetime_mean,duration_mean,0.0),\
    (delay_mean,risetime_mean,duration_mean,0.8,\
        delay_mean+2*delay_sd,risetime_mean+2*risetime_sd,duration_mean,0.3),\
    (delay_mean+2*delay_sd,risetime_mean+2*risetime_sd,duration_mean+2*duration_sd,1.0,\
        delay_mean,risetime_mean,duration_mean+duration_sd,0.8)\
    ]
## if stim_rate_seednum<-19, this sets whether to give const odorA and zero air and zero odorB to glom0
CONST_GLOM0 = False
## stimseed / len(salient_seed_glom_map) maps to 1% or 2% peak concentration
## stimseed % len(salient_seed_glom_map) maps to the salient_response indices for first 3 glomeruli
## salient_seed_glom_map = [ (response_index_glom0,response_index_glom1,response_index_glom2), ... ] 
salient_seed_glom_mapA = [ (0,1,2),(0,0,1),(0,0,2), (1,1,2),(1,0,1),(1,0,2), (2,1,2),(2,0,1),(2,0,2) ]
salient_seed_glom_mapB = [ (0,1,2),(1,2,0),(2,0,1), (0,2,1),(1,0,2),(2,1,0), (0,1,2),(1,0,1),(2,0,2) ]

## DoG (difference of Gaussians) kernel params
kernel2ORNfrate_ratio = 10.0 # rough factor from retro-fit kernels
DoG_center_mean = 250e-3 # s
DoG_center_sd = 100e-3 # s
DoG_width_min = 250e-3#300e-3#200e-3#100e-3 # s
DoG_width_max = 450e-3#450e-3#400e-3#250e-3 # s
DoG_airwidth_min = 250e-3#300e-3#200e-3#100e-3 # s
DoG_airwidth_max = 450e-3#450e-3#400e-3#250e-3 # s

## These two params are adjusted to obtain 35Hz
## average background (calculated in the generator program)
## Respiration tuned granule background baseline apart from respiratory tuning
gran_baseline_apart_from_resp_tuning = 27
gran_baseline_for_resp_tuning = 23

## Fitting the kernel cannot work on len(extratime) number of params.
## Decimate input, response and kernel to have the
## same large fitting_dt which must be a multiple of FIRINGFILLDT
kernel_time = 1.5 # seconds # length of kernel # 1.5s=(NUM_RESPS+1)*RESPIRATION, but diff from Priyanka
fitting_dt = 50e-3 # 50 ms same as Priyanka uses

## introduce random inhibition in mitrals from file. not used.
random_inh = False

## Number of different variations of stimuli to non-central glomeruli
## to measure effect of inhibition 
NUMINHS = 6

## Vary the main mitral firing rate or distance between the two or granule baseline frequency
## depending on the 'varied' variable in simset_odor.py
varied_mainrate = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0] # vary mitral A from 0 to 10 Hz in 6 steps
varied_distance = []
varied_gran_baseline = [1.0,1.5,2.0]

## White noise time constant, this causes the bandlimiting.
## should be a multiple of FIRINGFILLDT
NOISEDT = 5e-3
## Number of white noise trains to calculate a kernel
NUMWHITETRAINS = 250#60*4#6

## amount of slice on each side of the mitral that is healthy.
healthy_slice_width = 100e-6 # metres

sine_ORN_mean = 3#3 # Hz
sine_amplitude = 3#1 # Hz
sine_frequencies = [1.0,2.0,5.0,10.0,15.0] # Hz
#sine_frequencies = [20.0,30.0,40.0,50.0,60.0] # Hz
num_sins = len(sine_frequencies)
SIN_RUNTIME = 1.0+1.0/sine_frequencies[0] * 2 # 1s + two periods of the lowest frequency
