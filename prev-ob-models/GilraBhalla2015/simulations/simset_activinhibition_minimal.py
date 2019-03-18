## This file is programmatically generated.

netseedstr = "100.0"
mitdistance = 400.0 # microns
mitdistancestr = "_mitdist400.0" # microns

NO_SINGLES = False
## spine inhibition and singles are self-inh
## toggle them on/off together
NO_SPINE_INH = NO_SINGLES
NO_JOINTS = False
NO_MULTIS = NO_JOINTS
NO_PGS = False

## When testing ADI (ASYM_TEST = False), fixed current in mitB to generate 80Hz. 1mM Mg++.
## When testing asymmetry in inhibition (ASYM_TEST=True), same currents in mitA and mitB, and 0.2mM Mg++.
ASYM_TEST = False
## reverse roles of mitA and mitB in activity dependent inhibition
REVERSED_ADI = True
IN_VIVO = True
## tuftinput: if ODORINH, use higher inputs to ORNs, higher gran bgnd not used.
ODORINH = True
oninject_ext = 15.0 # Hz 
