## This file used to be programmatically generated for converging to best fit Activity Dependent Inhibition curve.
## But that doesn't give decent result, so set by hand.

import sys
sys.path.extend(["../networks"])
## do not import networkConstants as that imports this file, and it's circular then!!!
from networkConstantsMinimal import *
## STRONG_SYNAPSES is defined in networkConstants, but can't import it due to reason above,
## so duplicating the directed and frac_directed check below again.

## For STRONG_SYNAPSES i.e differential connectivity set mitral -> granule base excitation to 0.2nS
## else, for random / uniform connectivity, set the base value to 0.3nS
## This is to get the same amount of activity dependent inhibition (Arevian et al)
## for the different network connectivities...
if directed and frac_directed>0.0:
    mitral_granule_AMPA_Gbar = 0.2e-9 # Siemens
    granule_mitral_GABA_Gbar = 1.0e-9#12.0e-09 # Siemens
else: #### confirm ADI for 0% frac_directed setting below
    ## 0.3e-9 for 3% frac_directed, _mod mitral,
    ## but 0.2e-9 for 1% frac_directed, _mod_spikeinit mitral
    mitral_granule_AMPA_Gbar = 0.2e-9#0.3e-9 # Siemens
    granule_mitral_GABA_Gbar = 1.5e-9#12.0e-09 # Siemens
## For the _mod mitral with _spikeinit,
## self Gbar below must be reduced to 5 pS, else huge self-inhibition
## For the _mod mitral, 50 pS is fine, it doesn't get affected much by inhibition!
self_mitral_GABA_Gbar = 5e-12#5e-12#50e-12 # Siemens
