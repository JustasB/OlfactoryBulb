# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

# These are the simulation parameters for testing recurrent and lateral inhibition on mitrals via granules
# Only 2 gloms, 2 mitral sisters per glom.

from pylab import *
rc('path',simplify= False) # To ensure matplotlib connects points at expense of speed.

## assuming sys.path already has '../networks'
from networkConstants import * # has PROXIMAL_CONNECTION, DIRECTED and FRAC_DIRECTED

## set to 1 for 2MITS / 2 for 2GLOMS option set in generate_neuroML.py
## note: mit 2 is randomly connected to mit 0; mit 3 is super-inh onto mit 0.
mitralsidekickidx = 2

SIMDT = 5e-5 # seconds
## Correction factor for event-based vs graded synchans
## In event-based case, activation = weight/simdt while in graded case activation is used directly.
## Further, due to every sim point in graded case acting like an event,
## moose gives a convolved Gk which is much larger (say 1000 times) than the event based one.
## So I should correct for simdt and the activation weighted multiple events at every simdt.
#synchan_activation_correction = 1/SIMDT / 1000 # for SIMDT=1e-5
synchan_activation_correction = 1/SIMDT / 250 # for SIMDT=5e-5

PLOTDT = 5e-5 # seconds
# add SETTLETIME from stimuliConstants.py to the REALRUNTIME to get the simulation RUNTIME
REALRUNTIME = 1.0 # seconds - for active inhibition - stimulus only for 500ms

VOLTAGE_CLAMP = False

NUMBINS = 10
###### Somehow mpiexec doesn't like number of processes > 1000 (perhaps memory overflow)
###### and gives errors somewhere in the middle of the run while trying to send to the boss,
###### after which boss decides the below:
## rank 0 in job 35  gulabjamun.ncbs.res.in_56625 caused collective abort of all ranks
############ So ensure that NUMAVG+1 < 1000
####### For a network rather than isolated mitral cell, even 500 processes make the cluster hang.
####### Has worked only for 32 processes, didn't test inbetween yet.
##### Basically memory overflow problems - use less memory by using less data tables
NUMAVG = 50 # don't make NUMAVG higher than 200.

netseedstr = "300.0"#"500.0"
mitdistancestr = "_mitdist50.0" #microns
IN_VIVO = False

## Better to compare MAXSYNS when comparing PGvsGranule and synaptic strengths,
## else there will be +-0.25x random variation.
OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1_numgloms2_seed'\
    +netseedstr+mitdistancestr
granfilebase = "../firefiles/firefiles_baseline/firetimes_gran_baseline_invitro"
if DIRECTED:
    OBNet_file += '_directed'+str(FRAC_DIRECTED)
    if PROXIMAL_CONNECTION:
        OBNet_file += '_proximal'
    else:
        OBNet_file += '_distal'
if IN_VIVO: OBNet_file += '_2GLOMS.xml'
else: OBNet_file += '_2GLOMS_INVITRO.xml'

## SPIKEBLOCK below: TTX and TEA: set Na and Kfast to zero
## ONLY implemented for mitral channels!!!
SPIKEBLOCK = False # only for mitral cells presently

NO_SPINE_INH = False
NO_SINGLES = False
NO_JOINTS = False
NO_MULTIS = False
NO_PGS = False

CLUB_MITRALS = False # If unmodeled mitrals connect to same granule as modeled mitrals, provide extra excitation to granule from modeled mitral
ONLY_TWO_MITS = True#False # Only the two central mitrals with indices (neuroml id-s) given in networkConstants.py to test lateral inhibition.

ODOR_GIVEN = False # No odor given only current injection as below
offInject = 750e-12#200e-12#450e-12 # current in A, roughly 24Hz (8-9 in 400ms) firing w/o lateral inh -- Arevian et al
#offInject = 1000e-12 # current in A, roughly 40Hz firing w/o lateral inh -- Friedman & Strowbridge 2000
onInject = 2000e-12#1500e-12#1750e-12 # current in B, roughly 80Hz firing w/ lateral - Arevian et al 2008 use 60 to 100Hz in mitral B.
