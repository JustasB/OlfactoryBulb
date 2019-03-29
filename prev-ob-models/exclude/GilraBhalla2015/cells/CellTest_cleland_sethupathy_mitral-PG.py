#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../channels","../synapses","../networks","../simulations"])

from networkConstants import *
from synapseConstants import *
from moose_utils import * # has printCellTree(cell)
from moose.neuroml import *

from load_channels import *
from load_synapses import *
from simset_activinhibition import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

######## NOT WORKING -- BROKEN ------- USE inhibition_tuftinput.py INSTEAD.
## See attachInhibition() method -- basically:
## MorphML reader no longer sets up spikegen-s in potential locations suggested in neuroml cell file.
## So in attachInhibition() we have an empty spikegen_list and error later.

########## Need to run inside this directory:
##### From node000, run
## mpiexec -machinefile ~/hostfile -n <len(ORNfiringrates)+1> ~/Python2.6.4/bin/python2.6 CellTest_cleland_sethupathy_mitral-PG.py <mitral|PG|both>
## mpiexec -machinefile ~/hostfile -n 11 ~/Python-2.6.4/bin/python2.6 CellTest_cleland_sethupathy_mitral-PG.py both
##### 0 rank process is for collating all jobs. (rank starts from 0)
##### I assume rank 0 process always runs on the machine whose X window system has a Display connected and can show the graphs!
##### The rank 0 stdout is always directed to the terminal from which mpiexec was run.
##### X output is to the display attached to the terminal use ssh -Y node...
##### For long simulations save results in a text file for replotting later and avoid above ambiguity.
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size() # Total number of processes
mpirank = mpicomm.Get_rank() # Number of my process
mpiname = MPI.Get_processor_name() # Name of my node
# The 0th process is the boss who collates/receives all data from workers
boss = 0
print 'Process '+str(mpirank)+' on '+mpiname+'.'

ORNfiringrates = arange(0.0,20.0,2.0)

class CellTest:

    def __init__(self, outputtype):
        load_channels()
        ## synchan_activation_correction depends on SIMDT,
        ## hence passed all the way down to
        ## granule_mitral_GABA synapse from the top level
        load_synapses(synchan_activation_correction)
        MML = MorphML({'temperature':CELSIUS})
        self.outputType = outputtype
        self.cellList = []
        ## order matters, as self.cellname etc remain that in last iteration
        for cellname in ['PG','mitral']:
            if cellname == 'PG':
                #filename = 'PG_aditya2010_neuroML_L1_L2_L3.xml'
                filename = 'PG_aditya2012_neuroML_L1_L2_L3.xml'
                self.ORNCellSynapse = 'ORN_PG'
                self.somaName = 'soma_0'
                self.numORNSyns = NUM_ORN_PG_SYNS
            elif cellname == 'mitral':
                #filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3_mod.xml'
                filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3_mod_withspikeinit.xml'
                #filename = 'mitral_migliore_etal_2007.xml'
                self.ORNCellSynapse = 'ORN_mitral'
                self.somaName = 'Seg0_soma_0'
                self.numORNSyns = NUM_ORN_MITRAL_SYNS
            self.cellname = cellname
            ## returns {cellname:segDict}
            ## where segDict = { segid1 : [ segname,(proximalx,proximaly,proximalz),
            ##  (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] , ... }
            cellDict = MML.readMorphMLFromFile(filename,{})
            segDict = cellDict[cellname]
            self.context = moose.PyMooseBase.getContext()
            self.cell = moose.Cell(self.context.deepCopy(\
                self.context.pathToId('/library/'+cellname),self.context.pathToId('/'),cellname))
            self.soma = moose.Compartment('/'+self.cellname+'/'+self.somaName)
            self.soma_vm = setupTable('soma_vm', self.soma,'Vm')
            self.soma_vm.stepMode = TAB_SPIKE # store only spikes
            self.soma_vm.stepSize = -0.020 # that cross -20 mV
            self.spikeTableList = []
            self.makeSynapses(segDict,self.cell)
            self.attachORNs()
            #printCellTree(self.cell)
            #self.soma.inject = 1e-9 # Ampere
            ## 0th cell is PG, 1st cell is mitral
            self.cellList.append((self.cell,self.soma,self.soma_vm,self.spikeTableList))
        if self.outputType == 'both':
            self.attachInhibition()
        ### testing - channels
        # setup only one of the Tables below
        #self.var = setupTable('PG_mitral', moose.SynChan('/mitral/Seg0_glom_0_21/PG_mitral'), 'Gk')
        #self.var = setupTable('PG_mitral', moose.SpikeGen('/PG/soma_0/PG_mitral_spikegen'), 'state')
        #self.var = setupTable('soma_vm', moose.Compartment('/mitral/Seg0_soma_0'), 'Vm')
        #self.var = setupTable('soma_vm', moose.Compartment('/PG/soma_0'), 'Vm')
        #self.var = setupTable('mitral_glom', moose.Compartment('/mitral/Seg0_glom_0_21'), 'Vm')
        #self.var = setupTable('soma_ca', moose.CaConc(self.soma.path+'/Ca_mit_conc'), 'Ca')

    def makeSynapses(self,segDict,cell):
        ## segDict = { segid1 : [ segname,(proximalx,proximaly,proximalz),
        ##  (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] , ... }
        for seg in segDict.values():
            segname = seg[0]
            segpath = cell.path+'/'+segname
            if not self.context.exists(segpath):
                print "Missing compartment",segpath
                sys.exit(1)
            compartment = moose.Compartment(segpath)
            for syn_name in seg[5]:
                ## no need to make granule_mitral synapses, also these don't exist directly in /library
                ## rather they exist as granule_mitral_GABA: historical relic, shouldn't have been so.
                if ('ORN' not in syn_name) and ('PG' not in syn_name): continue
                synapse = self.context.deepCopy(self.context.pathToId('/library/'+syn_name),\
                    self.context.pathToId(compartment.path),syn_name)
                compartment.connect("channel", synapse, "channel")

    def attachInhibition(self): # 600 PG cells inhibit each mitral
        PG_data = self.cellList[0]
        mitral_data = self.cellList[1]
        synapse_list = []
        for compartment in mitral_data[0].getChildren(mitral_data[0].id): # compartments
            for child in moose.Neutral(compartment).getChildren(compartment): # channels and synapses
                if moose.Neutral(child).name in ['PG_mitral']:
                    synapse_list.append(moose.SynChan(child)) # mitral has the synapse being the 'post' cell
        spikegen_list = []
        ## MorphML reader no longer sets up spikegen-s in potential locations suggested in neuroml cell file.
        ## So below returns an empty spikegen_list and error later.
        ## ALSO MITRAL->PG IS NOT BEING SET UP!!!
        for compartment in PG_data[0].getChildren(PG_data[0].id): # compartments
            for child in moose.Neutral(compartment).getChildren(compartment): # channels and synapses
                if moose.Neutral(child).name in ['PG_mitral_spikegen']:
                    spikegen_list.append(moose.SpikeGen(child))
        for i in range(PGSYNS_PER_MITRAL):
            synapse = synapse_list[int(uniform(0,len(synapse_list)))] # Bad practice - This should be part of NetworkML - no random numbers in simulations, only in generators.
            #print "Connecting spikegen to ", synapse.path
            inhspikegen = spikegen_list[int(uniform(0,len(spikegen_list)))] # Bad practice - This should be part of NetworkML - no random numbers in simulations, only in generators.
            inhspikegen.threshold = THRESHOLD # from globalConstants.py
            inhspikegen.refractT = 0.25e-3 ## important to set this, else events are raised at every time step that Vm > Threshold
            inhspikegen.edgeTriggered = 1 # This ensures that spike is generated only on leading edge.
            inhspikegen.connect("event",synapse,"synapse")
            #### The delay and weight can be set only after connecting a spike event generator.
            #### delay and weight are arrays: multiple event messages can be connected to a single synapse
            # first argument below is the array index
            synapse.setWeight(0, 1)
            synapse.setDelay(0, 1.8e-3) # seconds

    def attachORNs(self):
        ##### Inhibitory Synapse NOT added as I am not sure if it is due to PG inhibition - Perhaps Cleland and Sethupathy also did not model it       
        ##### Excitatory + Inhibitory(not added) Synapse combo taken from the paper Djurisic etal 2008 JNeurosci.
        ##### Actually it's only an excitatory synapse, but they have used the inhibitory one to model the later time course.
        ##### Though this might be needed to account for PG cell inhibition?
        ##### Gbar-s have been set to ensure 16mV EPSP at the glom tuft as done in the paper.
        ##### Paper's defaults give only 8mV EPSP peak at glom tuft. So here multiplied by two. Cannot use supplemental table values as no cytoplasmic resistivity Ri in this model. Only axial resistance from glom to prim which doesn't matter much [maybe it does - loading rather than input resistance?]. Makes sense only if dendritic tree with many compartments.
        ##### Also no idea of how much to change the inhibitory part without Ri, so multiplied that also by 2

        synapse_list = []
        for compartment in self.cell.getChildren(self.cell.id): # compartments
            for child in moose.Neutral(compartment).getChildren(compartment): # channels and synapses
                if moose.Neutral(child).name in [self.ORNCellSynapse]:
                    synapse_list.append(moose.SynChan(child))
        for i in range(self.numORNSyns):
            synapse = synapse_list[int(uniform(0,len(synapse_list)))] # Bad practice - This should be part of NetworkML - no random numbers in simulations, only in generators.
            #print "Connecting ", synapse.path
            spiketable = moose.TimeTable(synapse.path+'/tt')
            #### SynChan's synapse MsgDest takes time as its argument. Thus spiketable should contain a list of spike times.
            spiketable.connect("event", synapse,"synapse")
            synapse.setWeight(0, 1) # 0th element in synaptic array set to weight 1
            self.spikeTableList.append(spiketable)

    def run(self):
        self.context.setClock(0, SIMDT, 0)
        self.context.setClock(1, SIMDT, 0) #### The hsolve and ee methods use clock 1
        self.context.setClock(2, SIMDT, 0) #### hsolve uses clock 2 for mg_block, nmdachan and others.
        self.context.setClock(PLOTCLOCK, PLOTDT, 0)
        
        firingrate = ORNfiringrates[mpirank-1]
        # Attach ORN firing files to PG and mitral
        for spike_table_list in [self.cellList[0][3],self.cellList[1][3]]: # spikeTableLists for PG and mitral.
            for i,spiketable in enumerate(spike_table_list):
                spiketable.filename = '../firefiles/firetimes_constrate_'+str(firingrate)+'_'+str(i%NUM_SPIKEFILES)+'.txt'
                #### To replicate ON shock
                #if i<100:
                #    spiketable.filename = '../firefiles/firetimes_ON_shock.txt' # just a single spike at 0.5s to all ORN-mitral/PG synapses
        self.context.reset() # A second reset will not work properly with tables!
        self.context.step(RUNTIME)
        print "sending for value ",firingrate,' from ',mpirank
        if self.outputType in ['mitral','both']:
            mpicomm.send( array(self.cellList[1][2]), dest=boss, tag=0 ) # send the mitral's soma Vm / spike time table
        else:
            mpicomm.send( array(self.cellList[0][2]), dest=boss, tag=0 ) # send the granule's soma Vm / spike time table
        ### testing - channels
        #mpicomm.send( array(self.var), dest=boss, tag=1 )

def collate():
    responsefreqlist = []
    for i,firingrate in enumerate(ORNfiringrates): # firingrate in Hz.
        soma_vm = mpicomm.recv(source=i+1, tag=0)
        print "received from ",i+1," on boss machine ",mpiname
        while soma_vm[-1] == 0.0: # MOOSE inserts crap zeros at the end in TAB_SPIKE table
            soma_vm = delete(soma_vm,-1) # numpy's delete 0th element
            if len(soma_vm) == 0: break
        responsefreqlist.append(len(soma_vm)/RUNTIME)
        #### testing - channels
        #var = mpicomm.recv(source=i+1, tag=1)
        #figure()
        #plot(var,',-b')
    figure()
    plot(ORNfiringrates,responsefreqlist,'+-r')
    show()
        
if __name__ == "__main__":
    seed([100.0])
    if mpirank==boss:
        collate()
    else:
        cell = CellTest(sys.argv[1])
        cell.run()
