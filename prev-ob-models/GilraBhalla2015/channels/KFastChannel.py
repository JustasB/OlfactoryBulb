#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math
import os

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

VK = -70.0e-3 # Volts

GK = 1956*sarea #Siemens, from mit4.hoc


class KFastChannel(moose.HHChannel):
    """K fast channel inherits from HHChannel."""
    def __init__(self, *args):
        """Setup the Na channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VK
        self.Gbar = GK
        self.addField('ion')
        self.setField('ion','K')
        self.Xpower = 2 # This will create HHGate instance xGate inside the Na channel
        self.Ypower = 1 # This will create HHGate instance yGate inside the Na channel
        ## Below gates get created after Xpower or Ypower are set to nonzero values
        ## I don't anymore have to explicitly create these attributes in the class
        #self.xGate = moose.HHGate(self.path + "/xGate")
        #self.yGate = moose.HHGate(self.path + "/yGate")
        selfdir = os.path.dirname(__file__)+os.sep
        fkfast_ninf=open(selfdir+"kfast_n.inf")
        fkfast_ntau=open(selfdir+"kfast_n.tau")
        fkfast_kinf=open(selfdir+"kfast_k.inf")
        fkfast_ktau=open(selfdir+"kfast_k.tau")
        kfast_NDIVS=1000 # ensure that there are 1001 lines in these above data files
        # VMIN and VMAX for kfast should be hard set to -100mV and 50mV, as per the data files, irrespective of other channels.
        VMIN_k = -0.1 # V
        VMAX_k = 0.050 # V
        self.xGate.A.xmin = VMIN_k
        self.xGate.A.xmax = VMAX_k
        self.xGate.A.xdivs = kfast_NDIVS
        self.xGate.B.xmin = VMIN_k
        self.xGate.B.xmax = VMAX_k
        self.xGate.B.xdivs = kfast_NDIVS
        self.yGate.A.xmin = VMIN_k
        self.yGate.A.xmax = VMAX_k
        self.yGate.A.xdivs = kfast_NDIVS
        self.yGate.B.xmin = VMIN_k
        self.yGate.B.xmax = VMAX_k
        self.yGate.B.xdivs = kfast_NDIVS

        for i in range(kfast_NDIVS+1):
            ## The files are from Davison's model, tau in physiological units ms, so convert
            ninf=float(fkfast_ninf.readline().split()[1]) # split each line in the file on whitespace and take the second value (first value is the voltage).
            ntau=1.0e-3*float(fkfast_ntau.readline().split()[1])
            self.xGate.A[i] = ninf/ntau
            self.xGate.B[i] = 1.0/ntau
            kinf=float(fkfast_kinf.readline().split()[1])
            ktau=1.0e-3*float(fkfast_ktau.readline().split()[1])
            self.yGate.A[i] =kinf/ktau
            self.yGate.B[i] =1.0/ktau

        fkfast_ninf.close()
        fkfast_ntau.close()
        fkfast_kinf.close()
        fkfast_ktau.close()
