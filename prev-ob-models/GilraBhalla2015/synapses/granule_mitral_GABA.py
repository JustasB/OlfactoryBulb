#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math

from pylab import *

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
from moose.utils import * # for BSplineFill

## if plastic synapse, choose the short term plasticity syn chan
if GABA_plastic: GABAbaseSynChan = moose.STPSynChan
else: GABAbaseSynChan = moose.SynChan

class granule_mitral_GABA(GABAbaseSynChan):
    """Non-saturating GABA synapse from granule to mitral cell."""
    def __init__(self, synchan_activation_correction, *args):
        GABAbaseSynChan.__init__(self,*args)
        self.Ek = granule_mitral_GABA_Ek
        self.Gbar = granule_mitral_GABA_Gbar
        self.tau1 = granule_mitral_GABA_tau1
        self.tau2 = granule_mitral_GABA_tau2
        self.addField('graded')
        if GRANULE_INH_GRADED:
            self.setField('graded','True')
            ######## Graded synapse
            inhsyntable = moose.Table(self.path+"/graded_table")
            inhsyntable.xmin = -43e-3#-52e-3 # V
            inhsyntable.xmax = -19e-3#-28e-3 # V
            inhsyntable.xdivs = 12
            # adjust activation curve so that all the dynamics happens between -52 and -28mV.
            act = [0.0] # at -52mV
            act.extend( [1/(1+math.exp(-(vm + 31.5e-3)/1.0e-3)) for vm in arange(-41e-3,-21.00001e-3,2e-3)] )
            act.extend([1.0]) # for -28mV
            for i,activation in enumerate(act):
                inhsyntable[i] = activation*synchan_activation_correction # synchan_activation_correction depends on SIMDT!
            inhsyntable.tabFill(1000,BSplineFill)
            inhsyntable.connect("outputSrc",self,"activation")
        else:
            self.setField('graded','False')
        self.addField('mgblock')
        self.setField('mgblock','False')

        ## Only depression, no facilitation:
        ## Venki Murthy 2005 paper
        self.tauD1 = GABA_recovery_time
        self.deltaF = 0.0
        self.d1 = GABA_depression_factor
        self.d2 = 1.0

        ########## Remember that the Vm graded synapse is actually a proxy for Ca dependent synapse:
        ########## Graded Ca dependent synapse
        #inhsyntable = moose.Table(granule._gran.path+"/InhSynTable_"+mcstr+"_"+gcstr+"_"+nsstr)
        #inhsyntable.xmin = 0 # Ca concentration - mMolar = mol/m^3 SI
        #inhsyntable.xmax = 4.0e-4 # Ca concentration - mMolar = mol/m^3 SI
        #inhsyntable.xdivs = 10
        ##for i,activation in enumerate([0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,1.0,1.0]):
        #for i,activation in enumerate([0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,1.0,1.0]):
        #    inhsyntable[i] = activation*synchan_activation_correction
        #inhsyntable.tabFill(1000,BSplineFill)
        #inhsyntable.connect("outputSrc",inhsyn,"activation")
        #granule._granPeriCaPool.connect('concSrc',inhsyntable,'msgInput')
