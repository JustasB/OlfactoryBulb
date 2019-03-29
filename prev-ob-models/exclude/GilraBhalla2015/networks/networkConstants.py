import sys
sys.path.extend(["../synapses"])
from synapseConstants import * # for mitral_granule_NMDA_tau2
from networkConstantsMinimal import *

#### z coordinates of each cell layer:
MITRAL_Z = 0.0003 # meters # all mitral somas are in a layer at this z
GRANULE_Z = 0.0 # meters # all granule somas are in a layer at this z
PG_Z = 0.0003 + 0.00075 # meters # all PG somas are in a layer at this z - place them slightly higher than mitrals's gloms for easy visualization

#### parameters for rat olfactory bulb
#### - using Egger and Urban's supplementary table - (obtained via email)
gran_spines = 100 # number of granule spines
# had to increase mit_dend_r to accomodate BBmit1993 - required for anisotropic granule connectivity.
mit_dend_r = 1700e-6#850e-6 # meters - mitral dendritic radius
gran_dend_r = 50e-6#25e-6#50e-6 # meters - granule dendritic radius # using smaller radius to make thinner column
allow_multi_gran_conn = True # if gran_dend_r is small, cannot enforce only one mit-gran connection, so allow.
bulb_r = 2e-3 # meters - olfactory bulb radius
mits_per_glom = 50 # actual number of mitrals per glom # take two tufted as one mitral 25+50/2=50
num_gloms = 2400 # actual number of glomeruli in the olfactory bulb - not modelled or even used presently.
# whether granule connectivity is columnar or uniform along dendrite.
# if True this prunes away granules not in a column
# No longer pruning, columnar due to reciprocal synapses over the primary dendrite.
#COLUMNAR = True
mit_syns = 10000 # number of mitral synapses, In a slice, there will be less!
# Number of MIT_SISTERS is extremely hard-coded everywhere. So don't change it from 2.
# num of modeled sister mitrals at each glomerulus
# - give same ORN inputs - as if gap junctions in glomerulus.
MIT_SISTERS = 2
# These are the indices of the central mitral and central glomerulus
mitralidx = 0
central_glom = mitralidx/2 # integer division - presently glom num 0.

grans_per_mit = 50 # number of granules per mitral/tufted.
glom_r = 50e-6 # meters radius of glomerulus including PG cells. I use this to calculate density of granule cells.
soma_r = 10e-6 # meters radius of mitral soma - taken from Gordon Shepherd

gran_density = mits_per_glom*grans_per_mit/(2*glom_r)**2 ## taking a square of side 2*glom_radius
# assume the same constant density throughout the bulb i.e. glomeruli are closely packed.
gran_spacing = 2*glom_r/(mits_per_glom*grans_per_mit)**0.5

## COLUMNAR specifies if extra granules outside glomeruli are to be pruned?
## For 2 mits, essentially granules near mit0 will remain.
## If mit1 is placed farther away, it will not have many granules near its soma.
#COLUMNAR = True
## Whether FRAC_DIRECTED of mits_per_syns will be
## connected between pairs listed in DIRECTED_CONNS.
## Keep DIRECTED True for simulating odors, not for simulating activity dependent inhibition.
## For ADI, we choose any two mitrals, not just connected ones.
DIRECTED = directed
## ensure that FRAC_DIRECTED * num of mitrals directed < 1.
## For NUM_GLOMS=10, 20mits all connected to mit0, FRAC_DIRECTED < 0.05.
## Can set FRAC_DIRECTED to 0.0 keeping DIRECTED=True. This will ensure that
## other mits lat dends are over directed centralmit's soma, if PROXIMAL_CONNECTION = True
FRAC_DIRECTED = frac_directed # I think you need to set this to 0.05 to get reasonable phase separation?
## STRONG_SYNAPSES sets whether directed synapses are increased by directed_synfactor
if DIRECTED and FRAC_DIRECTED>0.0: STRONG_SYNAPSES = True
else: STRONG_SYNAPSES = False
STRONG_IHNSYNS_CENTRAL_ONLY = False # Strongify the inh gran-|mit only on the central 'receiver' mits
## directed synapse exc and inh are increased by strongsynfactor-s below
## gran-|mit is strengthened uniformly, but mit->gran is enhanced more distally than proximally 
strongsynfactorinh = 4.0
strongsynfactorexclateral = 3.0#2.0#6.0
## WARNING! below strongsynfactorexccentral only affects strongified synapses, rest remain at 1.0
## So setting it below 1.0 counteracts strongification!
## Actually, setting it different from 1.0 is difficult to achieve by learning,
## since the proximal synapses are probably determined by the mean activity mostly,
## i.e. self-inhibition is set by the long-term mean activity, not a rare odor that gets decorrelated.
## Distal synapses have more scope, since they do not contribute much to self-inhibition.
strongsynfactorexccentral = 1.0
## some synapses are 1x, some at 2x, but baseline is set at 1x of mitral_granule_AMPA_Gbar
exc_avg_factor = 1.0#2.5 # avg exc mit->gran weight for baseline poisson input to granules
## exccentral is given from below proximal mitral compartments, rest provide exclateral to granules
## 0 is soma, 15-20 are prim dend, rest are nearby secondary dends upto level 4 (directed granules reach level 4)
#proximal_mitral_segids = [0,15,16,17,18,19,20,115,158,201,254]
proximal_mitral_segids = [0, 15,16,17,18,19,20,  115,116,117,153,118,  158,159,160,196,161,  201,202,203,249,204,  254,255,256,282,257]
## mitnum and centralmitnum are connected if DIRECTED_CONNS[mitnum] = centralmitnum
## DIRECTED_CONNS[mitnum] = None means that it is not connected to any of the two mits of centralglom
DIRECTED_CONNS = [None,None]*10 ## initialize
DIRECTED_CONNS[central_glom]=None ## central_glom
DIRECTED_CONNS[central_glom+1]=None ## central_glom
###### ALWAYS directed connect the 0th mitral, and keep the 1st mitral None or connect that also.
###### Gloms 1,2,3,4 get odorA, 5,6,7,8 get odorB and 9 only air.
###### So connect these gloms mits alternately to central glom's mits 0 and 1.
DIRECTED_CONNS[2] = 0
DIRECTED_CONNS[3] = None
DIRECTED_CONNS[4] = 1
DIRECTED_CONNS[5] = None
DIRECTED_CONNS[6] = 0
DIRECTED_CONNS[7] = None
DIRECTED_CONNS[8] = 1
DIRECTED_CONNS[9] = None
DIRECTED_CONNS[10] = 0
DIRECTED_CONNS[11] = None
DIRECTED_CONNS[12] = 1
DIRECTED_CONNS[13] = None
DIRECTED_CONNS[14] = 0
DIRECTED_CONNS[15] = None
DIRECTED_CONNS[16] = 1
DIRECTED_CONNS[17] = None
DIRECTED_CONNS[18] = 0
DIRECTED_CONNS[19] = None
## if PROXIMAL_CONNECTION, then
##   reasonable length of dendrite of mitnum must fall on centralmitnum soma.
##   if DIRECTED, FRAC_DIRECTED of granules to centralmitnum are connected near soma to mitnum.
## if not PROXIMAL_CONNECTION,
##   then roles of centralmitnum and mitnum are reversed.
##   but centralmitnum's dendrites are not forced to lie over the mitnum-s it connects to.
##   since it gets progressively more difficult to align all.
##   rather, mitnum's glomerulus is placed on to the centralmitnum's dendrite.
##     if a glom's two mitnums are connected to different centralmitnums,
##     then the glom is placed on dendrite of centralmitnum for the first mitnum.
PROXIMAL_CONNECTION = True
## Force all directed joints of mit0 to connect only to primary dendrite of mit0.
PRIMDEND_FORCE = False
## Force the first two gloms to be near the central glom
NEAR_2GLOMS = False#True

GRANS_CLUB_JOINTS = 1 # club GRANS_CLUB_JOINTS number of joint granules together
GRANS_CLUB_SINGLES = 100 #10 # club GRANS_CLUB_SINGLES number of single granules together
## WARNING: Despite setting 2GLOMS club factor for joints and singles below, it is used only for generation.
##          During simulation, the filebasename specified in simset_inhibition and simset_activinhibition
##          contains the club factor for joints and singles for 2GLOMS.
##          Also the club factor  for testing existence of netfiles (& running generator) for tuftADI and ADI,
##          is hardcoded in inhibition_tuftinput_repeats.py and activdep_inhibition_repeats.py
GRANS_CLUB_JOINTS_2GLOMS = 1#2 # as above but for ADI, STA and varied sims with 2 gloms.
GRANS_CLUB_SINGLES_2GLOMS = 100#20#100 # as above but for ADI, STA and varied sims with 2 gloms.
## Even though singles are clubbed, they can activate a realistic number of synapses; however
## computationally/neuroml-ly expensive to have all GRANS_CLUB_SINGLES number of gran-|mit synapses.
## So, per clubbed single, there will be SYNS_PER_CLUBBED_SINGLE number of inh syns,
## each of weight GRANS_CLUB_SINGLES/float(SYNS_PER_CLUBBED_SINGLE) times unit weight.
SYNS_PER_CLUBBED_SINGLE = 10

##### 12000 ORNs per glomerulus, 25 mitral cells per glomerulus
##### each ORN typically makes 8 connections within the glomerulus to interneurons and mitral: say half to mitral.
##### 12000*4 ~= 48K synapses. If an ORN does not connect to the same mitral twice,
##### a mitral has 48K/75 = 640 ORN synapses (25mitrals+50 tufted)
##### Since I am modelling all ORN connections via a single glomerulus compartment,
##### I can afford to connect 640 ORNs - odors are not generated at runtime.
##### I actually use 200 ORNs * 4 synapses per mitral = 800 NUM_ORN_MITRAL_SYNS
##### Set NUM_ORN to 1 for 'synchronized' inputs, vary only firing rate.
##### Or set this higher and vary firing rate over smaller range to distribute firings.
NUM_ORN = 200
NUM_SYN_PER_ORN = 4
## Without below dark ORN factor, mitral firing saturates for merely 5Hz ORN firing.
## Both Meredith's and Duchamp-Viret's experiments are extracellular.
## I predict that there are many dark ORNs a la dark mitrals,
## reducing the mean firing rate by an order of magnitude.
DARKORNFACTOR = 2 # should be an int, that divides NUM_ORN*NUM_SYN_PER_ORN exactly
NUM_ORN_MITRAL_SYNS = NUM_ORN*NUM_SYN_PER_ORN / DARKORNFACTOR # 800 ORN synapses per mitral, 50 spines per PG.
NUM_PG_SPINES = 50 # Shao et al 2009 says 25 for mice, 2 spines per 10microns for 10*(radius=25microns) dendritic length gives ~50 spines for rats.
NUM_ORN_PG_SYNS = NUM_PG_SPINES #int(0.6*NUM_PG_SPINES) # ORN->ET->PG is folded in, not just ORN->PG.
NUM_PG_to_M_T_ET_SPINES = NUM_PG_SPINES/2 # integer division, 50% of PG spines are reciprocal with ET or mitral. Rest to only ORNs.
#NUM_PG_PER_MITRAL = 20 # 20 PGs per mitral (Shepherd in Syn Org of Brain)
NUM_PG_PER_GLOM = 1000 # 1500-2000 JG cells per glom. 50% of these are PGs i.e. 1000 (Shipley 1996).
ET_PER_GLOM = 200 # Shipley private email. 10% of JG cells are ET. & Parrish-Aungst et al 2007 in mice.
NUM_SA_SYNS_PER_PG = 25 # 50 spines per PG. But actually SA->PG is actually SA->ET->PG. Also ORN->ET->SA. Unknown transfer function ORN->ET->SA->ET'->PG.
PG_CLUB = 1#5 # Thus the simulated num of PGs per mitral will be 20/PG_CLUB! Make it a divisor of NUM_PG_PER_MITRAL else integer division.
PGSYNS_PER_MITRAL = int( NUM_PG_to_M_T_ET_SPINES*NUM_PG_PER_GLOM/float(mits_per_glom+ET_PER_GLOM)/PG_CLUB ) # 25*1000/250/1 = 100 clubbed PG syns
## I spread out mitrals from (-glom_r/2,-glom_r/2) to (glom_r/2,glom_r/2)
## since mitral somas are quite spread out.
## Hence I need to spread the PGs also.
PG_halfextent = 2*glom_r
PG_spacing = 2*PG_halfextent/(NUM_PG_PER_GLOM/PG_CLUB)**0.5
PGindexsize = int( 2*PG_halfextent / PG_spacing ) + 1 # make a PGindexsize*PGindexsize array of PGs at every glomerulus.

# Presently I do not put in any centrifugal inhibition
# either as an EPSP factor or as a current injection - hence below values
GRAN_INH_FACTOR = 1.0#1.0/6.0 # Due to GABA-ergic inhibitory input in vivo, heights of EPSPs drop - not sure if it is centrifugal or present in slice also - Wellis and Kauer 1994, JPhysiol.
GRAN_INH_INJECT = 0.0#-20e-12 # Amp # Rather than an inhibition factor on the synapses as above, put in a current injection into the granules.

if DIRECTED:
    AMPA_factor=GRAN_INH_FACTOR # was 1.0
    NMDA_factor=GRAN_INH_FACTOR # was 1.25
    ## for FRAC_DIRECTED=0.0, use GABA_factor=1.5,
    ## for FRAC_DIRECTED=0.05, use GABA_factor=0.25
    GABA_factor=1.0 #0.15
else:
    #AMPA_factor=0.26*GRAN_INH_FACTOR
    #NMDA_factor=0.21*GRAN_INH_FACTOR
    AMPA_factor=GRAN_INH_FACTOR
    NMDA_factor=GRAN_INH_FACTOR
    GABA_factor=1.0 #0.15    

SATURATING_SINGLE_INHIBITION = False
SYN_VARIATION = 0.25 # +-25% variation in synaptic weights
## whether synaptic weights must be lognormal distributed vs uniformly distributed.
lognormal_weights = True#False
## whether the granule resting membrane potentials should be varied or not
VARY_GRANS_RMP = True
## SD of the granule RMP: varying Eleak by 2.25mV varies RMP by 1.5mV
gran_RMP_SD = 2.25e-3#1.5e-3 # 1.5 mV from Cang and Isaacson 2003
## whether the PG resting membrane potentials should be varied or not
VARY_PGS_RMP = True
## SD of the PG RMP: varying Eleak by 3mV varies RMP by 2mV for PG 2010
## PG 2013 RMP hardly varies on changing Eleak
PG_RMP_SD = 3e-3#2.0e-3 # 2.0 mV from Hayar et al 2004b

USE_SECDEND_DECAY = True
## gram--|mit inhibitory synapse conductance density (per area)
## decays exponentially with dendritic distance
mit_prim_decay_l = 100e-6 # meters approx from Lowe 2002
mit_sec_decay_l = 150e-6 # meters approx from Lowe 2002

## synaptic delays are for the first cellspike respresented by a clubbed cell
## spikes for the rest of the cells, that the clubbed cell represents, are 'spread' by values below...
exc_synaptic_delay = 1.8e-3 # seconds
inh_synaptic_delay = 0.6e-3 # seconds
ORN_synaptic_delay = 0.0 # seconds # fixed delay doesn't matter much, as I input Poisson spikes.

## unmodelled mitrals are put in via modelled mitral with random delays.
mitral_PG_delay_spread = 40e-3 # seconds
## unmodelled PGs are put in via modelled PGs with random delays.
## since exc is delayed by mitral_PG_delay spread, advance the inh delay below
PG_mitral_delay_spread = 160e-3#mitral_granule_NMDA_tau2 - mitral_PG_delay_spread
## unmodelled mitrals are put in via modelled mitral with random delays.
mitral_granule_delay_spread = 40e-3 # seconds
## unmodelled granules are put in via modelled granules with random delays.
## since exc is delayed by mitral_granule_delay spread, advance the inh delay below
granule_mitral_delay_spread = 160e-3#mitral_granule_NMDA_tau2 - mitral_granule_delay_spread
## delay_distrib: 0=none, 1=uniform, 2=exponential, 3=gamma
exc_delay_distrib = 1 # uniform
inh_delay_distrib = 2 # exponential

SA_integrate_time = 25e-3 # seconds
## integration - convolve taking the invalid initial period
## also automatically puts in a delay of SA_integrate_time
## so don't put in any extra delay in the synapse.
SA_delay = 0.0#25e-3 # seconds

## PGs like to stay pegged at -20mV or so, hence -15mV threshold.
## tuft and lateral dendrite do not reach 0mV and barely reach -10mV, hence this threshold.
THRESHOLD = -15.0e-3 #-35.0e-3 # V

## Reduce the extra proxy excitation from mitrals to granules by this factor
## because mitrals are experimentally spread out over 200-300 microns [Sosulski et al 2011],
## but, I only spread them over glom_r for their tufts to remain near the glom.
## As per my onenote labnote of 10 may 2011, the extra excitation reduces by half
## with 6*glom_r spread of mitrals. (set mitspread_extraexc_redux = 0.5)
## I'm now spreading out the mitrals by 2*glom_r (not glom_r), but still let the factor be 0.5.
mitspread_extraexc_redux = 0.5
