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

## if plastic synapse, choose the short term plasticity syn chan
if GABA_plastic: GABAbaseSynChan = moose.STPSynChan
else: GABAbaseSynChan = moose.SynChan

class mitral_self_GABA(GABAbaseSynChan):
    """Non-saturating GABA synapse from granule to mitral cell."""
    def __init__(self, *args):
        GABAbaseSynChan.__init__(self,*args)
        self.Ek = granule_mitral_GABA_Ek
        self.Gbar = self_mitral_GABA_Gbar
        self.tau1 = granule_mitral_GABA_tau1
        self.tau2 = granule_mitral_GABA_tau2
        self.addField('graded')
        self.setField('graded','False')
        self.addField('mgblock')
        self.setField('mgblock','False')
        
        ## Only depression, no facilitation:
        ## Venki Murthy 2005 paper
        self.tauD1 = GABA_recovery_time
        self.deltaF = 0.0
        self.d1 = GABA_depression_factor
        self.d2 = 1.0
