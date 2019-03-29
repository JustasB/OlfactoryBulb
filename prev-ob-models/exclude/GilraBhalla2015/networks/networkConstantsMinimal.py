## actual number of modelled gloms could be 10 (for odor testing)
## or 2 (for inhibition testing) decided during neuroml generation.
## can set number of modelled glom to whatever you like.
## Randomly half of them will lie on central glom's mit0 or mit1.
## First half will receive odor A. Rest will receive odor B.
NUM_GLOMS = 2

## Whether FRAC_DIRECTED of mits_per_syns will be
## connected between pairs listed in DIRECTED_CONNS.
## Keep directed True for simulating odors,
## Even for ADI, choose two connected mitrals.
directed = True

## ensure that FRAC_DIRECTED * num of mitrals directed < 1.
## For NUM_GLOMS=10, 20mits all connected to mit0, FRAC_DIRECTED < 0.05.
## Can set FRAC_DIRECTED to 0.0 keeping DIRECTED=True. This will ensure that
## other mits lat dends are over directed centralmit's soma, if PROXIMAL_CONNECTION = True
frac_directed = 0.01 # I think you need to set this to 0.05 to get reasonable phase separation?
