# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3


from pylab import *
rc('path',simplify= False) # To ensure matplotlib connects points at expense of speed.

# assuming sys.path already has '../networks'
from networkConstants import * # has SETTLETIME, PULSE_RUNTIME, pulsebins

from simset_odor_minimal import *

IN_VIVO = True

## set to 1 for 2MITS / 2 for 2GLOMS set during neuroml model generation by generate_neuroML.py
mitralsidekickidx = 2

## set varied to one of (None,'variedinh', 'mainrate', 'distance', 'gran_baseline')
## varied = None is for usual odor morphs and pulses, keep all mits.
## rest are for 2GLOMS:
##      'variedinh' is for changing delay etc to the odor response of glomB (variedinh sim)
##      others (mainrate, distance, gran_baseline) change conditions
##           under which kernel of glomB is calculated using whitenoise sim.
varied = None#'variedinh'#None#'mainrate'#'variedinh'
if varied == 'variedinh':
    ORNpathINHstr = "../firefiles/firefiles_variedinh/"
    ORNpathseedstr = "../firefiles/firefiles_variedinh/"
elif varied == 'mainrate':
    ORNpathINHstr = "../firefiles/firefiles_whitenoise/"
    ORNpathseedstr = "../firefiles/firefiles_whitenoise/"
granfilebase = "../firefiles/firefiles_baseline/firetimes_gran_baseline"

## for sinusoid inputs and related odor_sins.py, etc. only
LAT_SINS = True # whether sinusoids are given laterally or directly
SIN_MAIN_CONSTRATE = 3.0 # Hz # constant rate to give to central glom, with sin at lateral.

SIMDT = 5e-5 # seconds
## Correction factor for event-based vs graded synchans
## In event-based case, activation = weight/simdt
## while in graded case activation is used directly.
## Further, due to every sim point in graded case acting like an event,
## moose gives a convolved Gk which is much larger
## (say 1000 times) than the event based one.
## So I should correct for simdt and the activation weighted multiple events at every simdt.
#synchan_activation_correction = 1/SIMDT / 1000 # for SIMDT=1e-5
synchan_activation_correction = 1/SIMDT / 250 # for SIMDT=5e-5

PLOTDT = 5e-5 # seconds
## add SETTLETIME from stimuliConstants.py to the REALRUNTIME to get the simulation RUNTIME
REALRUNTIME = 1.0 # seconds # 2 respiration cycles

VOLTAGE_CLAMP = False

NUMBINS = 10
###### Somehow mpiexec doesn't like number of processes > 1000 and 
## gives errors somewhere in the middle of the run while trying to send to the boss,
## after which boss decides the below:
## rank 0 in job 35  gulabjamun.ncbs.res.in_56625 caused collective abort of all ranks
############ So ensure that NUMAVG+1 < 1000
####### For a network rather than isolated mitral cell, even 500 processes make the cluster hang.
##### Has worked only for 32 processes, didn't test inbetween yet.
##### Basically memory overflow problems - use less memory by using less data tables

## If unmodeled odor-responsive mitrals connect to same granule as their sister modelled mitrals,
## provide extra excitation to granule from modeled sister mitral
## as compensation for unmodelled sister odor-responsive mitrals, if CLUB_MITRALS=True.
## If CLUB_MITRALS=False, provide extra generic background to ALL granules during odor.
CLUB_MITRALS = False#True
## Only two mitrals for STA/whitenoise
## i.e. varied == 'mainrate'|'distance'|'gran_baseline'
if varied is None or varied == 'variedinh': ONLY_TWO_MITS = False
else: ONLY_TWO_MITS = True

ODOR_GIVEN = True
