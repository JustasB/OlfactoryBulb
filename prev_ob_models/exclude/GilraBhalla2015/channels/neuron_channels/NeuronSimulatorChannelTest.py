#!/usr/bin/env python

### This program plots a channel's state variables / hinf, htau etc. as a function of voltage.

mechanism_names = [
'Na_rat_ms','KDR_ms','KA_ms', # migliore and shepherd
'TCa_d','Ih_cb', # PG cell of Cleland and Sethupathy uses channels from various places
'Na_mit_usb','K2_mit_usb','K_mit_usb','LCa3_mit_usb','KA_bsg_yka'] # Bhalla and Bower 1993 mitral cell
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
]

mechanism_vars_global_or_ranges = ['global','global','global','global','global','global','global','global','global','global']

# index of the mechanism to test in above arrays.
test_mechanism_index = 2

import sys
import math

try:
    from neuron import *
except ImportError:
    print "ERROR: Could not import neuron. Please add the directory containing neuron.py in your PYTHONPATH"
    import sys
    sys.exit(1)

from pylab import *

sys.path.append('../..')
from globalConstants import * # for CELSIUS

h.load_file("stdlib.hoc")
h.load_file("stdrun.hoc") # for h.stop ,etc.
h.load_file("mit_memb.hoc")

VMIN = -100.0 #mV
VMAX = 100.0 # mV
VSTEPS = 400
DT = 0.001 # ms
# USING PHYSIOLOGICAL UNITS AS NEURON USES THEM
TSTOP =  DT*VSTEPS # ms

class ChannelTest:

    def __init__(self, mechanism_name, plotvars):
        h('create soma')
        self.soma = h.soma
        self.soma.insert(mechanism_name)
        print "Inserted mechanism ",mechanism_name,"in soma"
        
        self.vclamp = h.SEClamp(0.5, sec=self.soma) # single electrode clamp in soma
        # Ensure that electrode impedance is 0, thus a perfect voltage clamp!
        # But setting to exactly zero seems not to work (Perhaps divide by zero internally returning NaN-s)
        self.vclamp.rs = 1e-10
        self.vclamp.dur1 = TSTOP
        
        self.vclampV = h.Vector(arange(VMIN,VMAX,(VMAX-VMIN)/VSTEPS)) # mV
        self.vclampV.play(self.vclamp._ref_amp1,DT)
        self.somavecV = h.Vector()
        self.vect = h.Vector()
        self.somavecV.record(self.soma(0.5)._ref_v)
        #self.somavecV.record(self.vclamp._ref_vc)
        self.vect.record(h._ref_t)
        self.vararray = []
        for i,var in enumerate(plotvars):
            #self.vararray.append(h.Vector())
            #self.vararray[-1].record(h.ref(self.soma(0.5).Na_rat_ms.gmax)) # works for gmax but not for minf, etc.
            h('objref vec'+str(i))
            h('vec'+str(i)+' = new Vector()')
            self.vararray.append(eval('h.vec'+str(i)))
            
            # USE ONE OR THE OTHER OPTIONS BELOW:
            if mechanism_vars_global_or_ranges[test_mechanism_index] == 'global':
                # if the var say minf in naxn.mod has been declared global / parameter, it is not accessed as a section / assigned variable but as &minf_Na_rat_ms
                # As long as there is only one compartment having Na_rat_ms, I suppose minf_Na_rat_ms will have the value returned for that compartment's voltage.
                h('vec'+str(i)+'.record(&'+var+'_'+mechanism_name+')')
            else:
                # if the var say minf in naxn.mod has been declared as a range, it is accessed as a section variable i.e. as &soma.minf_Na_rat_ms(0.5)
                h('vec'+str(i)+'.record(&soma.'+var+'_'+mechanism_name+'(0.5))')
            
    def run(self):
        h.celsius = CELSIUS # nrnivmodl ignored my settings in .mod file! default is 6.3! Also ena, etc are also ignored in the .mod file!
        h.finitialize(VMIN) # Set all section voltages to VMIN, so that state variables do not start from say v=-65mV.
        h.frecord_init() # initialize the recording vectors.
        h.tstop = TSTOP
        h.dt = DT
        h.run()
                
if __name__ == "__main__":
    idx = test_mechanism_index
    ct = ChannelTest(mechanism_names[idx],mechanism_vars[idx])
    ct.run()
    plot(ct.vect,ct.somavecV,'+-r')
    title('soma V (mV) vs time (ms)')
    for i,var in enumerate(ct.vararray):
        figure()
        plot(ct.somavecV,var,',-b')
        title('state variable '+mechanism_vars[idx][i]+' of '+mechanism_names[idx]+' vs Voltage (mV)')
    show()
