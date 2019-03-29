#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math

# The PYTHONPATH should contain the location of moose.py and _moose.so
# files.  Putting ".." with the assumption that moose.py and _moose.so
# has been generated in ${MOOSE_SOURCE_DIRECTORY}/pymoose/ (as default
# pymoose build does) and this file is located in
# ${MOOSE_SOURCE_DIRECTORY}/pymoose/examples
# sys.path.append('..\..')
try:
    import moose
except ImportError:
    print "ERROR: Could not import moose. Please add the directory containing moose.py in your PYTHONPATH"
    import sys
    sys.exit(1)

from synapseConstants import *
    
class mitral_granule_saturatingAMPA(moose.KinSynChan):
    """Saturating AMPA synapse from mitral to granule cell."""
    def __init__(self, *args):
        moose.KinSynChan.__init__(self,*args)
        self.Ek = mitral_granule_AMPA_Ek
        self.Gbar = mitral_granule_saturatingAMPA_Gbar
        # KinSynChan is implemented from Destexhe, Mainen and Sejnowski, 1994
        # pulseWidth is the time for which the neurotransmitter is on
        self.pulseWidth = mitral_granule_saturatingAMPA_pulseWidth
        # decay time after neurotransmitter is switched off 1/beta in Destexhe, Mainen and Sejnowski, 1994
        self.tau1 = mitral_granule_saturatingAMPA_tau1
        # the fraction of bound/open receptors in infinite time with one synaptic event.
        # rise time tau2 or tau_r in the paper is calculated as tau_1*(1-rInf) and cannot be set.
        self.rInf = mitral_granule_saturatingAMPA_rInf
        self.addField('graded')
        self.setField('graded','False')
        self.addField('mgblock')
        self.setField('mgblock','False')
