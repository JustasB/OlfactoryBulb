#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program creates a Bhalla and Bower 1993 (upgraded to genesis 2 by Beeman) mitral cell model along with tables to pull data
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

EREST = -65.0e-3 # V
############# Set these such that 100 ORNs at approx 50Hz make the mitral cell fire in the middle of its linear range.
SYN_EXC_G = 1 * 8.6516e-9 # Siemens
SYN_INH_G = 1 * 2.2126e-9 # Siemens


class BBMitral:

    def __init__(self, moosename='/mitral', table=False):
        self.context = moose.PyMooseBase.getContext()
        self.path = moosename
        self.setupClocks()
        self.loadCell()
        if table:
            self.setupTables()
        self.addSynapses()
        self.connectORNs()
        # print self._mitral.method ######## By default Cell object uses hsolve method
        print("<< "+moosename+" fully loaded >>")
                                
    def loadCell(self):
        self.context.readCell('mit_aditya_davison_reduced.p',self.path)
        self._mitral = moose.Cell(self.path)
        self._mitralSoma = moose.Compartment(self.path+'/soma')
        self._mitralGlom = moose.Compartment(self.path+'/glom')
        self._mitralDend = moose.Compartment(self.path+'/dend')
        self._mitralDendNa = moose.HHChannel(self.path+'/dend/Na_mit_usb')
        self._mitralSomaNa = moose.HHChannel(self.path+'/soma/Na_mit_usb')
        self._mitralSomaCaPool = moose.CaConc(self.path+'/soma/Ca_mit_conc')
        self._mitralSomaLCa = moose.HHChannel(self.path+'/soma/LCa3_mit_usb')
        self._mitralSomaKCa = moose.HHChannel(self.path+'/soma/Kca_mit_usb')
        # Connect the LCa current to the Ca Pool
        self._mitralSomaLCa.connect('IkSrc',self._mitralSomaCaPool,'current')
        # Connect the KCa channel to the Ca concentration of Ca Pool
        self._mitralSomaCaPool.connect('concSrc',self._mitralSomaKCa,'concen')

    def addSynapses(self):
        ##### Excitatory + Inhibitory Synpase combo taken from the paper Djurisic etal 2008 JNeurosci.
        ##### Actually it's only an excitatory synapse, but they have used the inhibitory one to model the later time course.
        ##### Though this might be needed to account for PG cell inhibition?
        ##### Gbar-s have been set to ensure 16mV EPSP at the glom tuft as done in the paper.
        ##### Paper's defaults give only 8mV EPSP peak at glom tuft. So here multiplied by two. Cannot use supplemental table values as no cytoplasmic resistivity Ri in this model. Only axial resistance from glom to prim which doesn't matter much [maybe it does - loading rather than input resistance?]. Makes sense only if dendritic tree with many compartments.
        ##### Also no idea of how much to change the inhibitory part without Ri, so multiplied that also by 2
        self._ExcORNSyn = moose.SynChan(self.path+"/ExcSyn")
        self._ExcORNSyn.Ek = EREST + 0.07 # Volts, 70mV above EREST
        self._ExcORNSyn.Gbar = 0.3 * SYN_EXC_G # Siemens
        ### The delay and weight can be set only after connecting a spike event generator.
        ### delay and weight are arrays: multiple event messages can be connected to a single synapse
        self._ExcORNSyn.tau1 = 1.6094e-3 # seconds
        self._ExcORNSyn.tau2 = 1.6094e-3 # seconds

        self._InhORNSyn = moose.SynChan(self.path+"/InhSyn")
        self._InhORNSyn.Ek = EREST - 0.01 # Volts, 10mV below EREST
        self._InhORNSyn.Gbar = 0.3 * SYN_INH_G # Siemens
        ### The delay and weight can be set only after connecting a spike event generator.
        ### delay and weight are arrays: multiple event messages can be connected to a single synapse
        self._InhORNSyn.tau1 = 5.6631e-3 # seconds
        self._InhORNSyn.tau2 = 5.6631e-3 # seconds

        self._mitralGlom.connect("channel", self._ExcORNSyn, "channel")
        self._mitralGlom.connect("channel", self._InhORNSyn, "channel")

    def connectORNs(self):
        self.spiketable = moose.TimeTable(self.path+'/tt')
        #### SynChan's synapse MsgDest takes time as its argument. Thus spiketable should contain a list of spike times.
        self.spiketable.connect("event", self._ExcORNSyn,"synapse")
        self.spiketable.connect("event", self._InhORNSyn,"synapse")
        self._ExcORNSyn.setWeight(0, 1) # 0th element in synaptic array set to weight 1
        self._InhORNSyn.setDelay(0, 7.6337e-3) # seconds
        self._InhORNSyn.setWeight(0, 1) # 0th element in synaptic array set to weight 1
                
    def setupTables(self):
        self._data = moose.Neutral(self.path+"/data")
        # Setup the tables to pull data
        self._vmTableSoma = moose.Table("vmTableSoma", self._data)
        self._vmTableSoma.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        self._vmTableSoma.connect("inputRequest", self._mitralSoma, "Vm")
        self._vmTableSoma.useClock(PLOTCLOCK)
        self._vmTableGlom = moose.Table("vmTableGlom", self._data)
        self._vmTableGlom.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        self._vmTableGlom.connect("inputRequest", self._mitralGlom, "Vm")
        self._vmTableGlom.useClock(PLOTCLOCK)
        self._vmTableDend = moose.Table("vmTableDend", self._data)
        self._vmTableDend.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        self._vmTableDend.connect("inputRequest", self._mitralDend, "Vm")
        self._vmTableDend.useClock(PLOTCLOCK)
        
    def setupClocks(self):
        ###### I suppose by default all elements use clock 0 the global clock. Not sure, Niraj checked that hsolve uses clock 1
        self.context.setClock(0, SIMDT, 0)
        self.context.setClock(1, SIMDT, 0) #### The hsolve and ee methods use clock 1
        self.context.setClock(2, SIMDT, 0) #### hsolve uses clock 2 for mg_block, nmdachan and others.
        self.context.setClock(PLOTCLOCK, PLOTDT, 0)

def testMitral():
    mitral = BBMitral() ## Use the default randomly generated synConns
    mitral.context.reset()
    mitral.context.step(RUNTIME)
    plot(mitral._vmTableSoma,label='Soma Vm')    
    plot(mitral._vmTableGlom,label='Glom Vm')  
    show()  
    

if __name__ == "__main__":
    random.seed() ##### If no parameter is given, it uses current system time
    testMitral()
