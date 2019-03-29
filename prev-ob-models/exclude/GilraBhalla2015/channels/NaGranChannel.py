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

VNa = 45e-3 # Volts

gransarea = 5e-9 # m^2 default surface area of soma from granule.tem
GNa = 1611*gransarea # Siemens, from granule.tem and then multiplying by area of mitral soma

class NaGranChannel(moose.HHChannel):
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
        # These gates get created only after Xpower or Ypower are set to nonzero values
        # So we have to explicitly insert these fields in the class
        self.xGate = moose.HHGate(self.path + "/xGate")
        self.yGate = moose.HHGate(self.path + "/yGate")
        tabchan_file = open("channels/tabchannels.dat")
        nagran_NDIVS = 1000 # ensure that there are 1001 lines in these above data files
        self.xGate.A.xmin = VMIN
        self.xGate.A.xmax = VMAX
        self.xGate.A.xdivs = nagran_NDIVS
        self.xGate.B.xmin = VMIN
        self.xGate.B.xmax = VMAX
        self.xGate.B.xdivs = nagran_NDIVS
        self.yGate.A.xmin = VMIN
        self.yGate.A.xmax = VMAX
        self.yGate.A.xdivs = nagran_NDIVS
        self.yGate.B.xmin = VMIN
        self.yGate.B.xmax = VMAX
        self.yGate.B.xdivs = nagran_NDIVS

        v = VMIN

        for i in range(nagran_NDIVS+1):
            ## The files are from Davison's model, physiological units ms^-1, so convert
            (minf_s,mtau_ms_s,hinf_s,htau_ms_s) = tabchan_file.readline().split()[9:13] # split each line in the file on whitespace and take the 10th to 13th columns (see tabchannels.hoc for parsing).
            minf = float(minf_s)
            mtau = float(mtau_ms_s) * 1.0e-3
            hinf = float(hinf_s)
            htau = float(htau_ms_s) * 1.0e-3
            self.xGate.A[i] = minf/mtau
            self.xGate.B[i] = 1.0/mtau ### Python version below 3.0 will treat 1/x as integer division! Use 1.0/x !!!!
            self.yGate.A[i] = hinf/htau
            self.yGate.B[i] = 1.0/htau
            v = v + dv

        tabchan_file.close()
