#!/usr/bin/env python
# -*- coding: utf-8 -*-

### This program plots a synapse's Ik and Gk on an event

synapse_names = ['mitral_granule_NMDA']

# index of the synapse to test in above arrays.
test_synapse_index = 0

import sys,os
sys.path.extend(["..","../synapses"])
from load_synapses import *
from moose_utils import *

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

from pylab import *

SETTLETIME = 0.25 # s
RUNTIME = 1.0+SETTLETIME # s
SIMDT = 5e-6 #s
PLOTDT = 5e-6 #s
                        
if __name__ == "__main__":
    load_synapses(1.0)
    idx = test_synapse_index
    syn_name = synapse_names[idx]
    context = moose.PyMooseBase.getContext()

    soma = moose.Compartment('/soma')
    soma.length = 1e-6 # m
    soma.diameter = 1e-6 # m
    soma.Cm = 0.01 * 3.14159 * soma.length * soma.diameter

    synid = context.deepCopy(context.pathToId('/library/'+syn_name),soma.id,syn_name)
    syn = moose.SynChan(synid)
    #### connect the soma to the synapse
    if syn.getField('mgblock')=='True': # If NMDA synapse based on mgblock, connect to mgblock
        mgblock = moose.Mg_block(syn.path+'/mgblock')
        compartment_connection = mgblock
    else: # if SynChan or even NMDAChan, connect normally
        compartment_connection = syn

    soma.connect("channel", compartment_connection, "channel")
    tt = moose.TimeTable('/soma/tt')
    tt.connect("event", syn,"synapse")
    syn.setWeight(syn.numSynapses-1, 1.0)
    syn.setDelay(syn.numSynapses-1, 0)
    fn = 'temp_spikefile.txt'
    f = open(fn,'w')
    f.write(str(SETTLETIME))
    f.close()
    tt.filename = fn
    os.remove(fn)
    
    somaVm = setupTable('/soma/somaVm', soma, 'Vm')
    chanIk = setupTable('/soma/chanIk', compartment_connection, 'Ik')
    chanGk = setupTable('/soma/chanGk', compartment_connection, 'Gk')

    resetSim(context, SIMDT, PLOTDT)
    context.step(RUNTIME)

    tvec = arange(0.0,RUNTIME+1e-12,PLOTDT)
    figure()
    title('somaVm')
    plot(tvec,somaVm, 'r,')
    figure()
    title('Synapse Ik')
    plot(tvec,chanIk, 'g,')
    figure()
    title('Synapse Gk')
    plot(tvec,chanGk, 'b,')
    
    show()
