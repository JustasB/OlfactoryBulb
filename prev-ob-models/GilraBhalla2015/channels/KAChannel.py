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

from channelConstants import *
    
VKA = -70.0e-3 # Volts

GKA = 58.7*sarea # Siemens, from mit4.hoc

taup=1.38e-3 # seconds

tauq=150.0e-3 # seconds

def calc_KA_pinf(v):
    return 1/(1+math.exp(-(42e-3+v)/13e-3))

def calc_KA_qinf(v):
    return 1/(1+math.exp((110e-3+v)/18e-3))

class KAChannel(moose.HHChannel):
    """Na channel inherits from HHChannel."""
    def __init__(self, *args):
        """Setup the Na channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VKA
        self.Gbar = GKA
        self.addField('ion')
        self.setField('ion','K')
        self.Xpower = 1 # This will create HHGate instance xGate inside the Na channel
        self.Ypower = 1 # This will create HHGate instance yGate inside the Na channel
        ## Below gates get created after Xpower or Ypower are set to nonzero values
        ## I don't anymore have to explicitly create these attributes in the class
        #self.xGate = moose.HHGate(self.path + "/xGate")
        #self.yGate = moose.HHGate(self.path + "/yGate")
        self.xGate.A.xmin = VMIN
        self.xGate.A.xmax = VMAX
        self.xGate.A.xdivs = NDIVS
        self.xGate.B.xmin = VMIN
        self.xGate.B.xmax = VMAX
        self.xGate.B.xdivs = NDIVS
        self.yGate.A.xmin = VMIN
        self.yGate.A.xmax = VMAX
        self.yGate.A.xdivs = NDIVS
        self.yGate.B.xmin = VMIN
        self.yGate.B.xmax = VMAX
        self.yGate.B.xdivs = NDIVS
        
        v = VMIN

        for i in range(NDIVS+1):
            self.xGate.A[i] = calc_KA_pinf(v)/taup
            self.xGate.B[i] = 1.0/taup
            self.yGate.A[i] = calc_KA_qinf(v)/tauq
            self.yGate.B[i] = 1.0/tauq
            v = v + dv
