# -*- coding: cp1252 -*-
from custom_params import *

try:
    
    from neuron import h
    h.celsius = 35
    
except:
    pass

# load the custom parameters
try:
  module = __import__(filename)
  globals().update(vars(module))
except:
  print 'error during params import'

from copy import copy
from math import pi,exp


#plasticity
fi_tau1 = 1.0
fi_tau2 = 100.0


Ngloms = 127 # glomerolous
Nmitral_per_glom = 5 # mitral per glomerolous
Nmitral = Ngloms * Nmitral_per_glom


# Random123 secondary stream identifiers
# note: the primary stream index is the "gid" which is ordered as
# Nmitral, Ngranule, synapses
# Not all secondary identifiers are used for all types
stream_last=1
stream_soma = stream_last; stream_last += 1
stream_dend = stream_last; stream_last += 1
stream_apic = stream_last; stream_last += 1
stream_tuft = stream_last; stream_last += 1
stream_glom = stream_last; stream_last += 1
stream_granule_pos = stream_last; stream_last += 1
stream_mitraldendconnect = stream_last; stream_last += 1000 #allows per dendrite stream
stream_dummysyn = stream_last; stream_last +=  Nmitral * 5 #needs to be greater than max_dummies in dymmysyns which is currently 75.
stream_dummygen = stream_last; stream_last += Nmitral * 5

# for the odorstim
stream_last += stream_ods_shift
stream_ods_w = stream_last; stream_last += 6
stream_ods_act = stream_last; stream_last += 1
stream_ods_bg = stream_last; stream_last += 1

def ranstream(id1, id2):
    #print 'ranstream ', id1, id2
    r = h.Random()
    r.Random123(id1, id2)
    return r


# glomerulus radius
GLOM_RADIUS = 40.
GLOM_DIST = 100.

# glom_axis, define axis size for glomerulus
bulbCenter = [ 50., 1275., 0. ]
bulbAxis = [ 2100., 2800., 2100.]
glomAxis = copy(bulbAxis)
somaAxis = [ copy(bulbAxis), copy(bulbAxis) ]

for i in range(3):
    somaAxis[0][i] -= 600 + 100
    somaAxis[1][i] -= 600
    glomAxis[i] += 300

granule_max_depth = 400.
granAxisInf = [ somaAxis[0][0] - 2 * granule_max_depth,
             somaAxis[0][1] - 2 * granule_max_depth,
             somaAxis[0][2] - 2 * granule_max_depth ]

from glom import *
loadGloms() # put glomerulus coords in glomRealCoords

# soma max x and max y
SOMA_MAX_X = 1000
SOMA_MAX_Y = 2000





granule_field_radius = 50. #microns

gid_mitral_begin = 0
gid_granule_begin = gid_mitral_begin + Nmitral

granule_diam = 10. 

granule_priden2_len = 250.


Nx_granule = int(somaAxis[0][0] / grid_dim) + 1
Ny_granule = int(somaAxis[0][1] / grid_dim) + 1
Nz_granule = int(somaAxis[0][2] / grid_dim) + 1
granule_origin = [ int(bulbCenter[0] - somaAxis[0][0] * .5),
                   int(bulbCenter[1] - somaAxis[0][1] * .5),
                   int(bulbCenter[2] - somaAxis[0][2] * .5) ]

mg_on_granule_std = 25.
