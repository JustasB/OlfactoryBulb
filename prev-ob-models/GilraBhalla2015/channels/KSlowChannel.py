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

GK = 28*sarea # Siemens, from mit4.hoc


class KSlowChannel(moose.HHChannel):
    """K slow channel inherits from HHChannel."""
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
        fkslow_ninf=open(selfdir+"kslow_n.inf")
        fkslow_ntau=open(selfdir+"kslow_n.tau")
        fkslow_kinf=open(selfdir+"kslow_k.inf")
        fkslow_ktau=open(selfdir+"kslow_k.tau")
        kslow_NDIVS=1000 # ensure that there are 1001 lines in these above data files
        # VMIN and VMAX for kfast should be hard set to -100mV and 50mV, as per the data files, irrespective of other channels.
        VMIN_k = -0.1 # V
        VMAX_k = 0.050 # V
        self.xGate.A.xmin = VMIN_k
        self.xGate.A.xmax = VMAX_k
        self.xGate.A.xdivs = kslow_NDIVS
        self.xGate.B.xmin = VMIN_k
        self.xGate.B.xmax = VMAX_k
        self.xGate.B.xdivs = kslow_NDIVS
        self.yGate.A.xmin = VMIN_k
        self.yGate.A.xmax = VMAX_k
        self.yGate.A.xdivs = kslow_NDIVS
        self.yGate.B.xmin = VMIN_k
        self.yGate.B.xmax = VMAX_k
        self.yGate.B.xdivs = kslow_NDIVS

        for i in range(kslow_NDIVS+1):
            ## The files are from Davison's model, physiological units ms^-1, so convert
            ninf=float(fkslow_ninf.readline().split()[1]) # split each line in the file on whitespace and take the second value (first value is the voltage).
            ntau=1.0e-3*float(fkslow_ntau.readline().split()[1])
            self.xGate.A[i] = ninf/ntau
            self.xGate.B[i] = 1.0/ntau
            kinf=float(fkslow_kinf.readline().split()[1])
            ktau=1.0e-3*float(fkslow_ktau.readline().split()[1])
            self.yGate.A[i] =kinf/ktau
            self.yGate.B[i] =1.0/ktau

        fkslow_ninf.close()
        fkslow_ntau.close()
        fkslow_kinf.close()
        fkslow_ktau.close()
