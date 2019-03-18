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
from math import pi

use_fi_stdp = False # FastInhibSTDP vs FastInhib
m2g_multiple = False # multiple connections within glomerulus

#plasticity
fi_tau1 = 1.0
fi_tau2 = 100.0


# specific grow parameters

GROW_MAX_ITERATIONS = 2000
GROW_MAX_ATTEMPTS = 100

##### ATTENTION: parameters may change in the future

## apical params
APIC_DIAM = 5. # modified from Francesco's value
APIC_LEN_MAX = 550.

## tuft params
TUFT_DIAM = .8 # modified from Francesco's value
N_MIN_TUFT = 4 
N_MAX_TUFT = 6
TUFT_MIN_LEN = 40
TUFT_MAX_LEN = 80

# dendrites paramaters
N_MIN_DEND = 4
N_MAX_DEND = 6
DEND_DIAM = 2

# dendrites definition
N_MIN_DEND = 4
N_MAX_DEND = 9
N_MEAN_DEND = 5

# dendrites max length, normal distribution
DEND_LEN_MU = 837.97
DEND_LEN_VAR = 283.04 ** 2
DEND_LEN_MIN = 50.
DEND_LEN_MAX = 1800.

# dendrites bifurcation parameters, exponential distribution
BIFURC_LEN_MU_1 = 85.32
BIFURC_LEN_MIN_1 = 2.95
BIFURC_LEN_MAX_1 = 325.
BIFURC_LEN_MU_2 = 226.61
BIFURC_LEN_MIN_2 = .5
BIFURC_LEN_MAX_2 = 825.
BIFURC_PROB = [ 0.75, 0.63, 0.42, 0.28, 0.06 ]

# neg. exp
BIFURC_PHI_MU = .5

def gen_dend_diam(dist): return 0.9 + 2.6 * exp(-0.0013 * dist) #value's estimated with a fitting

#class DiamDend:
    # cambiare con cautela!
  #  __DIAM_MIN = 0.2
  #  __DIAM_MAX = 3.

    # NOISE
 #   __ALPHA_NOISE = 0.3 # peso del noise
    
 #   @staticmethod
 #   def __f_diam(dist): return 0.9 + 2.6 * exp(-0.0013 * dist)

 #   @staticmethod
 #   def __ns_diam(dist): return 0.6 + 1.5 * exp(-0.003 * dist)

 #   @staticmethod
 #   def genFirstDiam(r): return DiamDend.genDiam(0, r)

 #   @staticmethod
 #   def genDiam(dist, r):
 #       d = DiamDend.__f_diam(dist) + r.normal(0, DiamDend.__ns_diam(dist) ** 2) * DiamDend.__ALPHA_NOISE
 #       while d < DiamDend.__DIAM_MIN or d > DiamDend.__DIAM_MAX:
 #           d = DiamDend.__f_diam(dist) + r.normal(0, DiamDend.__ns_diam(dist)) * DiamDend.__ALPHA_NOISE
 #       return d

## random walk, noise
WALK_RHO = 20.
BIFURC_PHI_MIN = pi / 24.
BIFURC_PHI_MAX = pi / 5.
NS_PHI_B = 0.16
NS_PHI_MIN = -6.26
NS_PHI_MAX = 6.26
NS_THETA_B = 0.1407
NS_THETA_MIN = -1.56
NS_THETA_MAX = 1.18

GROW_RESISTANCE = 1.5


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
stream_m2g = stream_last; stream_last += 1
stream_dummysyn = stream_last; stream_last +=  Nmitral #needs to be greater than max_dummies in dymmysyns which is currently 75.
stream_dummygen = stream_last; stream_last += Nmitral 

# for the odorstim
stream_last += stream_ods_shift
stream_ods_w = stream_last; stream_last += N_MAX_TUFT
stream_ods_act = stream_last; stream_last += 1
stream_ods_bg = stream_last; stream_last += 1

def ranstream(id1, id2):
    #print 'ranstream ', id1, id2
    r = h.Random()
    r.Random123(id1, id2)
    return r


# glomerulus radius
GLOM_RADIUS = 50.
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

from Glom import *
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
