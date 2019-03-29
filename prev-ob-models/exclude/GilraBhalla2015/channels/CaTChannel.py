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

# ALL SI UNITS

VCa = 0.070 # Volts

GCa = 40*sarea # Siemens, from mit4.hoc

def calc_Ca_alpha_s(v):
    return 7.5e3/(1+math.exp((-30e-3-v)/7e-3))

def calc_Ca_beta_s(v):
    return 1.65e3/(1+math.exp((v+30e-3)/4e-3))

def calc_Ca_alpha_r(v):
    return 6.8/(1+math.exp((v+60e-3)/12e-3))

def calc_Ca_beta_r(v):
    return 60/(1+math.exp(-v-30e-3/11e-3))


class CaTChannel(moose.HHChannel):
    """Ca channel inherits from HHChannel."""
    def __init__(self, *args):
        """Setup the Ca channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VCa
        self.Gbar = GCa
        self.addField('ion')
        self.setField('ion','Ca')
        self.Xpower = 1 # This will create HHGate instance xGate inside the Ca channel
        self.Ypower = 1 # This will create HHGate instance yGate inside the Ca channel
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
            self.xGate.A[i] = calc_Ca_alpha_s(v)
            self.xGate.B[i] = calc_Ca_alpha_s(v) + calc_Ca_beta_s(v)
            self.yGate.A[i] = calc_Ca_alpha_r(v)
            self.yGate.B[i] = calc_Ca_alpha_r(v) + calc_Ca_beta_r(v)
            v = v + dv
