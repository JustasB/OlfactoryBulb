#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program creates a Migliore & Shepherd 2008 gran cell model along with tables to pull data.
# Only the two biggest compartments are modelled.
import sys
import math
# The PYTHONPATH should contain the location of moose.py and _moose.so
# files.  Putting ".." with the assumption that moose.py and _moose.so
# has been generated in ${MOOSE_SOURCE_DIRECTORY}/pymoose/ (as default
# pymoose build does) and this file is located in
# ${MOOSE_SOURCE_DIRECTORY}/pymoose/examples
# sys.path.append('../..')

import moose

from mooseConstants import *
from globalConstants import *
    
from pylab import * # part of matplotlib that depends on numpy but not scipy
import random
import pickle

class BBGranule:

    def __init__(self, moosename='/granule', table=False):
        self.context = moose.PyMooseBase.getContext()
        self.path = moosename
        self.setupClocks()
        self.loadCell()
        if table:
            self.setupTables()
        # print self._gran.method ######## By default Cell object uses hsolve method
        print("<< "+moosename+" fully loaded >>")

    def loadCell(self):
        self.context.readCell('gran_aditya_migliore.p',self.path)
        self._gran = moose.Cell(self.path)
        self._granSoma = moose.Compartment(self.path+'/soma')
        self._granSomaKA = moose.Compartment(self.path+'/soma/KA_ms')
        self._granPeri = moose.Compartment(self.path+'/periphery')
        self._granPeriNa = moose.HHChannel(self.path+'/periphery/Na_rat_ms')
        self._granPeriKA = moose.HHChannel(self.path+'/periphery/KA_ms')
                
    def setupTables(self):
        self._data = moose.Neutral(self.path+"/data")
        # Setup the tables to pull data
        self._vmTableSoma = moose.Table("vmTableSoma", self._data)
        self._vmTableSoma.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        self._vmTableSoma.connect("inputRequest", self._granSoma, "Vm")
        self._vmTableSoma.useClock(PLOTCLOCK)
        
    def setupClocks(self):
        ###### I suppose by default all elements use clock 0 the global clock. Not sure, Niraj checked that hsolve uses clock 1
        self.context.setClock(0, SIMDT, 0)
        self.context.setClock(1, SIMDT, 0) #### The hsolve and ee methods use clock 1
        self.context.setClock(2, SIMDT, 0) #### hsolve uses clock 2 for mg_block, nmdachan and others.
        self.context.setClock(PLOTCLOCK, PLOTDT, 0)

def setup_iclamp(compartment, name, delay1, width1, level1):
    moose.Neutral('/elec') # If /elec doesn't exists it creates /elec and returns a reference to it. If it does, it just returns its reference.
    pulsegen = moose.PulseGen('/elec/pulsegen'+name)
    iclamp = moose.DiffAmp('/elec/iclamp'+name)
    iclamp.saturation = 1e6
    iclamp.gain = 1.0
    pulsegen.trigMode = 0 # free run
    pulsegen.baseLevel = 0.0
    pulsegen.firstDelay = delay1
    pulsegen.firstWidth = width1
    pulsegen.firstLevel = level1
    pulsegen.secondDelay = 1e6
    pulsegen.secondLevel = 0.0
    pulsegen.secondWidth = 0.0
    pulsegen.connect('outputSrc',iclamp,'plusDest')
    iclamp.connect('outputSrc',compartment,'injectMsg')
    return pulsegen


if __name__ == "__main__":
    sys.path.append("channels/")
    from NaGranChannelMS import *
    from KAChannelMS import *
    from KDRChannelMS import *
    NaGranChannelMS("/library/Na_rat_ms") # maintaining the old name to avoid changing .p file even though the Na channel here is from Migliore and Shepherd 2008
    KAChannelMS("/library/KA_ms")
    KDRChannelMS("/library/KDR_ms")
    
    seed() ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    gran = BBGranule(table=True)
    iclamp = setup_iclamp(gran._granSoma, '_gran', 50e-3, 300e-3, 50e-12) # 500pA for 100ms to get a single AP
    gran.context.reset()
    gran.context.step(500e-3)
    plot(gran._vmTableSoma,'r,-')
    #gran._granPeriKA.Gbar = 0
    #gran._granSomaKA.Gbar = 0
    #gran.context.reset()
    #gran.context.step(500e-3)
    #plot(gran._vmTableSoma,'g,-')
    show()
    
    
