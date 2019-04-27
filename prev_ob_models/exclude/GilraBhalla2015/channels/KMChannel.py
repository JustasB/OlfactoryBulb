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

VKM = -70.0e-3 # Volts

gransarea = 5e-9 # m^2 default surface area of soma from granule.tem
GKM = 1334*gransarea # Siemens

def calc_KM_taux(v):
    return 1/( 3.3*math.exp((v+35.0e-3)/40.0e-3) + math.exp(-(v+35.0e-3)/20.0e-3) ) # seconds

def calc_KM_xinf(v):
    return 1/(1+math.exp(-(35.0e-3+v)/5.0e-3))

class KMChannel(moose.HHChannel):
    """K muscarinic channel inherits from HHChannel."""
    def __init__(self, *args):
        """Setup the KM channel with defaults"""
        moose.HHChannel.__init__(self,*args)
        self.Ek = VKM
        self.Gbar = GKM
        self.addField('ion')
        self.setField('ion','K')
        self.Xpower = 1 # This will create HHGate instance xGate inside the Na channel
        #self.Ypower = 0 # This will create HHGate instance yGate inside the Na channel
        ## Below gate gets created after Xpower or Ypower are set to nonzero values
        ## I don't anymore have to explicitly create this attribute in the class
        #self.xGate = moose.HHGate(self.path + "/xGate")
        self.xGate.A.xmin = VMIN
        self.xGate.A.xmax = VMAX
        self.xGate.A.xdivs = NDIVS
        self.xGate.B.xmin = VMIN
        self.xGate.B.xmax = VMAX
        self.xGate.B.xdivs = NDIVS
        
        v = VMIN

        for i in range(NDIVS+1):
            taux = calc_KM_taux(v)
            self.xGate.A[i] = calc_KM_xinf(v)/taux
            self.xGate.B[i] = 1/taux
            v = v + dv
