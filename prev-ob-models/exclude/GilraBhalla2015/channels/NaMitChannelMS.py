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

# This Na channel is the same for mitral and granule 'naxn.mod' in Migliore and Shepherd 2008,
# except for VNa and sh.
# sh and VNa are passed to the constructor

#VNa = 50e-3 # Volts # Migliore and Shepherd 2008 mitral - different from granule 60mV

GNa = 200*sarea # Siemens, from mit4.hoc and then multiplying by area of mitral soma

## Lifted directly from Migliore and Shepherd 2008 Na gran channel
#sh = 10e-3 # V # (5mV in nmdol but 10mV in hoc file - mitral.hoc) - 15mV for granule
mmin=0.02e-3 # s
hmin=0.5e-3 # s
qt = 2.0**((CELSIUS-24.0)/10.0) # CELSIUS is defined in globalConstants.py
#Na_alpha_m:
tha = -30e-3 #V
Ra=0.4e6 # /s/V,
qa=7.2e-3 #V
#Na_beta_m:
tha = -30e-3 #V
Rb=0.124e6 # /s/V
qa=7.2e-3 # V
#Na_alpha_h:
thi1 = -45e-3 #V
Rd=0.03e6 # /s/V
qd=1.5e-3 #V 
#Na_beta_h:
thi2 = -45e-3 #V
Rg=0.01e6 # /s/V
qg=1.5e-3 #V 
### They use the above to calculate minf, mtau, htau, BUT finally have a different formula for hinf!!!
thinf = -50e-3 #V
qinf = 4e-3 #V

def trap0(vminusth,a,q):
	if (abs(vminusth) > 1e-9):
	    return a * (vminusth) / (1 - math.exp(-(vminusth)/q))
	else:
	    return a * q

#### V IMP: 10^6 factor in alpha_m and beta_m, as there is V in Nr, so the constant in front had units ms-1 mV^-1, so the 10^6 factor.
#def calc_Na_alpha_m(v):
#    #return 0.32e6*(v+42e-3)/(1-math.exp(-(42e-3+v)/4e-3)) # From Bhalla and Bower 1993 paper
#    return 400.0*(v+15e-3)/(1-math.exp(-(15e-3+v)/7.2e-3)) # tha+sh=-15mV. Use 15mV instead of 42mV
#
#def calc_Na_beta_m(v):
#    #return 0.28e6*(v+15e-3)/(-1+math.exp((15e-3+v)/5e-3)) # From Bhalla and Bower 1993 paper
#    return 124.0*(v-15e-3)/(-1+math.exp((-15e-3+v)/7.2e-3)) # -tha-sh=15mV. Use -15mV instead of 15mV
#
#def calc_Na_alpha_h(v):
#    #return 0.128e3/math.exp((v+38e-3)/18e-3) # From Bhalla and Bower 1993 paper
#    return 30.0*(v+30e-3)/math.exp((v+30e-3)/1.5e-3) # thi1+sh = -30mV, Use 30mV instead of 38mV.
#
#def calc_Na_beta_h(v):
#    #return 4.0e3/(1+math.exp(-(v+15e-3)/5.0e-3)) # From Bhalla and Bower 1993 paper
#    return 10.0*(v-30e-3)/(1+math.exp(-(v-30e-3)/1.5e-3)) # -thi2-sh = 30mV, Use -30mV instead of 15mV.
#####

class NaMitChannelMS(moose.HHChannel):
    """Na channel inherits from HHChannel."""
    def __init__(self, sh, VNa, *args):
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
            a = trap0(v-(tha+sh),Ra,qa)
            b = trap0(-v+(tha+sh),Rb,qa)
            mtau = 1/(a+b)/qt
            if mtau<mmin: mtau=mmin
            minf = a/(a+b)
            self.xGate.A[i] = minf/mtau
            self.xGate.B[i] = 1/mtau
            #self.xGate.tweakTau() # convert mtau and minf to A and B tables
            
            a = trap0(v-(thi1+sh),Rd,qd)
            b = trap0(-v+(thi2+sh),Rg,qg)
            htau = 1/(a+b)/qt
            if htau<hmin: htau=hmin
            hinf = 1/(1+math.exp((v-thinf-sh)/qinf))
            self.yGate.A[i] = hinf/htau
            self.yGate.B[i] = 1/htau
            #self.yGate.tweakTau() # convert htau and hinf tables to A and B tables.
            v = v + dv

