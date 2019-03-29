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
    
class mitral_granule_AMPA(moose.SynChan):
    """Non-saturating AMPA synapse from mitral to granule cell."""
    def __init__(self, *args):
        moose.SynChan.__init__(self,*args)
        self.Ek = mitral_granule_NMDA_Ek
        self.Gbar = mitral_granule_AMPA_Gbar
        self.tau1 = mitral_granule_AMPA_tau1
        self.tau2 = mitral_granule_AMPA_tau2
        self.addField('graded')
        self.setField('graded','False')
        self.addField('mgblock')
        self.setField('mgblock','False')
