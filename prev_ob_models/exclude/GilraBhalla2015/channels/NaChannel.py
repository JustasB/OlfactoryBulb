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

VNa = 45e-3 # Volts

GNa = 1532*sarea # Siemens, from mit4.hoc and then multiplying by area of mitral soma

#### V IMP: 10^6 factor in alpha_m and beta_m, as there is V in Nr, so the constant in front had units ms-1 mV^-1, so the 10^6 factor.
def calc_Na_alpha_m(v):
    return 0.32e6*(v+42e-3)/(1-math.exp(-(42e-3+v)/4e-3)) # From Bhalla and Bower 1993 paper

def calc_Na_beta_m(v):
    return 0.28e6*(v+15e-3)/(-1+math.exp((15e-3+v)/5e-3)) # From Bhalla and Bower 1993 paper

def calc_Na_alpha_h(v):
    return 0.128e3/math.exp((v+38e-3)/18e-3) # From Bhalla and Bower 1993 paper

def calc_Na_beta_h(v):
    return 4.0e3/(1+math.exp(-(v+15e-3)/5.0e-3)) # From Bhalla and Bower 1993 paper


class NaChannel(moose.HHChannel):
    """Na channel inherits from HHChannel."""
    def __init__(self, *args):
        """Setup the Na channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VNa
        self.Gbar = GNa
        self.addField('ion')
        self.setField('ion','Na')
        self.Xpower = 3 # This will create HHGate instance xGate inside the Na channel
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
            self.xGate.A[i] = calc_Na_alpha_m(v)
            self.xGate.B[i] = calc_Na_alpha_m(v) + calc_Na_beta_m(v)
            self.yGate.A[i] = calc_Na_alpha_h(v)
            self.yGate.B[i] = calc_Na_alpha_h(v) + calc_Na_beta_h(v)
            v = v + dv
