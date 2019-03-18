#!/usr/bin/env python
# -*- coding: utf-8 -*-

### This program plots a channel's state variables / hinf, htau etc. as a function of voltage.

mechanism_names = [
'Na_rat_ms','KDR_ms','KA_ms', # migliore and shepherd
'TCa_d','Ih_cb', # PG cell of Cleland and Sethupathy uses channels from various places
'Na_mit_usb','K2_mit_usb','K_mit_usb','LCa3_mit_usb','KA_bsg_yka']
#'Gran_CaHVA_98','cal'] # extra channels from neuroConstruct examples
mechanism_vars = [
['minf','mtau','hinf','htau'],
['minf','mtau'],
['minf','mtau','hinf','htau'],
['minf','mtau','hinf','htau'],
['linf','ltau'],
['minf','mtau','hinf','htau'],
['ninf','ntau','kinf','ktau'],
['ninf','ntau','kinf','ktau'],
['sinf','stau','rinf','rtau'],
['pinf','ptau','qinf','qtau']
#['minf','mtau','hinf','htau'],
#['minf','mtau']
]

# index of the mechanism to test in above arrays.
test_mechanism_index = 8

import sys
import math

# The PYTHONPATH should contain the location of moose.py and _moose.so
# files.  Putting ".." with the assumption that moose.py and _moose.so
# has been generated in ${MOOSE_SOURCE_DIRECTORY}/pymoose/ (as default
# pymoose build does) and this file is located in
# ${MOOSE_SOURCE_DIRECTORY}/pymoose/examples
try:
    import moose
except ImportError:
    print "ERROR: Could not import moose. Please add the directory containing moose.py in your PYTHONPATH"
    import sys
    sys.exit(1)

sys.path.append('..')
from mooseConstants import *

from load_channels import *

from pylab import *
                        
if __name__ == "__main__":
    load_channels()
    idx = test_mechanism_index
    for varidx in range(len(mechanism_vars[idx])/2): # loop over each inf and tau i.e. two elements of mechanism_vars[idx] at a time
        var = ['x','y','z'][varidx]
        gate = moose.HHGate('/library/'+mechanism_names[idx]+'/'+var+'Gate')
        VMIN = gate.A.xmin
        VMAX = gate.A.xmax
        NDIVS = gate.A.xdivs # will use same VMIN, VMAX and NDIVS for A and B tables.
        dv = (VMAX-VMIN)/NDIVS
        vrange = [VMIN+i*dv for i in range(NDIVS+1)]
        figure()
        plot(vrange,[gate.A[i]/gate.B[i] for i in range(NDIVS+1)],'b-,') # Assume A and B have corresponding number of entries
        title('state variable '+mechanism_vars[idx][2*varidx]+' of '+mechanism_names[idx]+' vs Voltage (V)')
        figure()
        plot(vrange,[1.0/gate.B[i] for i in range(NDIVS+1)],'b-,')
        title('state variable '+mechanism_vars[idx][2*varidx+1]+' of '+mechanism_names[idx]+' vs Voltage (V)')
    show()

