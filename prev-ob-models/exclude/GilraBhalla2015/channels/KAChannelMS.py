#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math

# KA channel of Migliore and Shepherd 2005 to 2008 - translated from kamt.mod

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
    
VKA = -90.0e-3 # Volts

GKA = 20.0*sarea # Siemens

a0m=0.04
vhalfm=-45e-3 # V
zetam=0.1e3 # /V
gmm=0.75 # dimensionless

a0h=0.018
vhalfh=-70e-3 # V
zetah=0.2e3 # /V
gmh=0.99 # dimensionless

sha=9.9e-3 # V
shi=5.7e-3 # V

q10=3.0
qt=q10**((CELSIUS-24.0)/10.0) # CELSIUS is a global constant
# (CELSIUS-24)/10 - integer division!!!! Ensure floating point division
# FOR INTEGERS: be careful ^ is bitwise xor in python - used here qt=2, but above qt=3!

def calc_KA_malp(v):
    return math.exp(zetam*(v-vhalfm))

def calc_KA_mbet(v):
    return math.exp(zetam*gmm*(v-vhalfm))

def calc_KA_halp(v):
    return math.exp(zetah*(v-vhalfh))

def calc_KA_hbet(v):
    return math.exp(zetah*gmh*(v-vhalfh))

def calc_KA_mtau(v):
    return calc_KA_mbet(v)/(qt*a0m*(1+calc_KA_malp(v))) * 1e-3 # convert to s

def calc_KA_htau(v):
    return calc_KA_hbet(v)/(qt*a0h*(1+calc_KA_halp(v))) * 1e-3 # convert to s

def calc_KA_minf(v):
    return 1/(1 + math.exp(-(v-sha-7.6e-3)/14e-3))

def calc_KA_hinf(v):
    return 1/(1 + math.exp((v-shi+47.4e-3)/6e-3))

class KAChannelMS(moose.HHChannel):
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
            mtau = calc_KA_mtau(v)
            self.xGate.A[i] = calc_KA_minf(v)/mtau
            self.xGate.B[i] = 1.0/mtau
            htau = calc_KA_htau(v)
            self.yGate.A[i] = calc_KA_hinf(v)/htau
            self.yGate.B[i] = 1.0/htau
            v = v + dv
