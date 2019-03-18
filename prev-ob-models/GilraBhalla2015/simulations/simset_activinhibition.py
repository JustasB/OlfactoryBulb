# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## These are the simulation parameters for testing activity dependent inhibition between two mitrals via granules
## Only 2 gloms, 2 mitral sisters per glom.

from pylab import *
rc('path',simplify= False) # To ensure matplotlib connects points at expense of speed.

## assuming sys.path already has '../networks'
from networkConstants import * # has SETTLETIME, PULSE_RUNTIME, pulsebins
## below file is programmatically generated for repeats with different params.
from simset_activinhibition_minimal import *

SIMDT = 5e-5 # seconds
## Correction factor for event-based vs graded synchans
## In event-based case, activation = weight/simdt while in graded case activation is used directly.
## Further, due to every sim point in graded case acting like an event,
## moose gives a convolved Gk which is much larger (say 1000 times) than the event based one.
## So I should correct for simdt and the activation weighted multiple events at every simdt.
#synchan_activation_correction = 1/SIMDT / 1000 # for SIMDT=1e-5
synchan_activation_correction = 1/SIMDT / 250 # for SIMDT=5e-5

PLOTDT = 5e-5 # seconds
## add SETTLETIME from stimuliConstants.py to the REALRUNTIME to get the simulation RUNTIME
## very important to use 400ms stimulus (a la Arevian),
## as the granule starts firing after 400ms for ~100Hz
REALRUNTIME = 0.4 # seconds - for active inhibition - stimulus only for 400ms

VOLTAGE_CLAMP = False

NUMBINS = 10
###### Somehow mpiexec doesn't like number of processes > 1000 and
###### gives errors somewhere in the middle of the run
###### while trying to send to the boss, after which boss decides the below:
## rank 0 in job 35  gulabjamun.ncbs.res.in_56625 caused collective abort of all ranks
############ So ensure that NUMAVG+1 < 1000
####### For a network rather than isolated mitral cell, even 500 processes make the cluster hang.
####### Has worked only for 32 processes, didn't test inbetween yet.
##### Basically memory overflow problems - use less memory by using less data tables
NUMAVG = 50 # don't make NUMAVG higher than 200.

if IN_VIVO:
    #OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed10_jointsclubbed1_numgloms1_seed'+netseedstr
    OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1'\
        '_numgloms2_seed'+netseedstr+mitdistancestr
    granfilebase = "../firefiles/firefiles_baseline/firetimes_gran_baseline_noresp"
else:
    ## half the granules are lost in slicing :
    ## look at my old pre-CNS simulations with various r-dependent connectivity.
    ## even for columns, I assume half of the granules on the prim dend will die.
    ## Instead of half the number of syns, grans more than healthy_slice_width away in y are pruned.
    OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1'\
        '_numgloms2_seed'+netseedstr+mitdistancestr
    granfilebase = "../firefiles/firefiles_baseline/firetimes_gran_baseline_invitro"
if DIRECTED:
    OBNet_file += '_directed'+str(FRAC_DIRECTED)
    if PROXIMAL_CONNECTION:
        OBNet_file += '_proximal'
    else:
        OBNet_file += '_distal'
if IN_VIVO: OBNet_file += '_2GLOMS.xml'
else: OBNet_file += '_2GLOMS_INVITRO.xml'

PLOT_EXTRAS = False # whether to plot the first few singles and joints for each data point.
if PLOT_EXTRAS: SPIKETABLE = False #make the mitral soma tables and SOME interneurons record Vm-s and not spiketimes.
else: SPIKETABLE = True #make the mitral soma tables and ALL interneurons record spiketimes and not Vm-s.

PG_raiseRMP = 0.0#7.5e-3 # 5mV
granule_raiseRMP = 0.0#7.5e-3 # 5mV
PG_inject = 0.0 #1e-12 # 1pA
granule_inject = 0.0 #1e-12 # 1pA

## unmodelled sister mitrals are not clubbed for activity dep inh 
## into modelled sister mitrals to multiply excitation to granules,
## as there is no common input via the glomerulus.
CLUB_MITRALS = False # If unmodeled sister mitrals connect to same granule as modeled sister mitrals, provide extra excitation to granule from modeled mitral
ONLY_TWO_MITS = True#False # Only two mitrals with indices (neuroml id-s) given in sim_utils.py to test lateral inhibition.

ODOR_GIVEN = False # No odor given only current injection as below
if IN_VIVO:
    onInject = 2000e-12 # current in B, roughly 70Hz firing w/o lateral - Actually Arevian et al 2008 use 60 to 100Hz in mitral B.
else:
    onInject = oninject_ext#1500e-12 # current in B, roughly 80Hz firing w/o lateral - Actually Arevian et al 2008 use 60 to 100Hz in mitral B.
if ASYM_TEST: Imax = 1001.0e-12#A # endpoint of activDDI plot
else: Imax = 2501e-12#3501.0e-12#3001.0e-12# A # endpoint of activDDI plot
ORNmax = 20.0
