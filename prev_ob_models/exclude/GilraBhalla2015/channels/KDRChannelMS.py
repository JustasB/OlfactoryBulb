#!/usr/bin/env python
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

VKDR = -90.0e-3 # Volts

gransarea = 5e-9 # m^2 default surface area of soma from granule.tem
GKDR = 20*gransarea # Siemens

a0m=0.0035
vhalfm=-50e-3 # V
zetam=0.055e3 # /V
gmm=0.5 # dimensionless

q10=3.0
qt=q10**((CELSIUS-24.0)/10.0) # CELSIUS is a global constant
# (CELSIUS-24)/10 - integer division!!!! Ensure floating point division
# FOR INTEGERS: be careful ^ is bitwise xor in python - used here qt=2, but above qt=3!

def calc_KA_malp(v):
    return math.exp(zetam*(v-vhalfm))

def calc_KA_mbet(v):
    return math.exp(zetam*gmm*(v-vhalfm))

def calc_KA_mtau(v):
    return calc_KA_mbet(v)/(qt*a0m*(1+calc_KA_malp(v))) * 1e-3 # convert to seconds

def calc_KA_minf(v):
    return 1/(1 + math.exp(-(v-21e-3)/10e-3))


class KDRChannelMS(moose.HHChannel):
    """K Delayed Rectifier channel translated from Migliore and Shepherd 2007."""
    def __init__(self, *args):
        """Setup the KDR channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VKDR
        self.Gbar = GKDR
        self.addField('ion')
        self.setField('ion','K')
        self.Xpower = 1 # This will create HHGate instance xGate inside the Na channel
        #self.Ypower = 0 # This will create HHGate instance yGate inside the Na channel
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
        
        v = VMIN

        for i in range(NDIVS+1):
            mtau = calc_KA_mtau(v)
            self.xGate.A[i] = calc_KA_minf(v)/mtau
            self.xGate.B[i] = 1.0/mtau
            v = v + dv
