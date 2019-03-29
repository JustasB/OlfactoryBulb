#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../channels","../synapses","../networks","../simulations","../generators"])

from networkConstants import *
from synapseConstants import *
from stimuliConstants import * # has SETTLETIME and NUM_ORN_FILES_PER_GLOM
from simset_activinhibition import * # has REALRUNTIME
from moose.utils import * # has printCellTree(cell)
from moose.neuroml import *

from load_channels import *
from load_synapses import *
from moose_utils import * # has attach_spikes()

from pylab import * # part of matplotlib that depends on numpy but not scipy

RUNTIME = REALRUNTIME + SETTLETIME

########## Need to run (from node000):
## mpiexec -machinefile ~/hostfile -n <length of ORNfiringrates+1> ~/Python-2.6.4/bin/python2.6 <full script path> <cellname>
## example: (from node000)
## mpiexec -machinefile ~/hostfile -n 41 ~/Python-2.6.4/bin/python2.6 CellTest_cleland_sethupathy.py PG
## OR
## ~/Python-2.6.4/bin/python2.6 CellTest_cleland_sethupathy.py PG
#### 0 rank process is for collating all jobs. (rank starts from 0)
## rank 0 process should run on the machine whose X window system has a Display connected and can show the graphs!!
## The rank 0 stdout is always directed to the terminal from which mpiexec was run. But X Display will be on node of boss process.
##### For long simulations save results in a text file for replotting later and avoid above ambiguity.
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size() # Total number of processes
mpirank = mpicomm.Get_rank() # Number of my process
mpiname = MPI.Get_processor_name() # Name of my node
# The 0th process is the boss who collates/receives all data from workers
boss = 0
print 'Process '+str(mpirank)+' on '+mpiname+'.'

ORNmax=20.0
if mpisize>1: ORNfiringrates = arange(0.0,ORNmax,ORNmax/float(mpisize-1))
else: ORNfiringrates = [5.0]
NUM_SPIKEFILES = 40

class CellTest:

    def __init__(self, cellname):
        self.context = moose.PyMooseBase.getContext()
        self.cellname = cellname
        load_channels()
        load_synapses(synchan_activation_correction)
        MML = MorphML({'temperature':CELSIUS})
        if cellname == 'PG':
            filename = 'PG_aditya2010_neuroML_L1_L2_L3.xml'
            self.ORNCellSynapse = 'ORN_PG'
            self.somaName = 'soma_0'
            self.numORNSyns = NUM_ORN_PG_SYNS
        elif cellname == 'mitral':
            filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3_mod_withspikeinit.xml'
            self.ORNCellSynapse = 'ORN_mitral'
            self.somaName = 'Seg0_soma_0'
            self.numORNSyns = NUM_ORN_MITRAL_SYNS
        self.cellSegmentDict = MML.readMorphMLFromFile(filename,{})
        self.cell = moose.Cell(self.context.deepCopy(self.context.pathToId('/library/'+cellname),\
            self.context.pathToId('/'),cellname))
        self.soma = moose.Compartment('/'+self.cellname+'/'+self.somaName)
        self.soma_vm = setupTable('soma_vm', self.soma,'Vm')
        if mpisize>1:
            self.soma_vm.stepMode = TAB_SPIKE # store only spikes
            self.soma_vm.stepSize = -0.030 # that cross -30 mV
        self.spikeTableList = []
        self.attachORNs()
        
        ### testing - channels
        # setup only one of the Tables below
        #self.soma_var = setupTable('soma_kca', moose.HHChannel2D(self.soma.path+'/Kca_mit_usb'), 'X')
        #self.soma_var = setupTable('soma_ca', moose.CaConc(self.soma.path+'/Ca_mit_conc'), 'Ca')

        ### For testing the mitral cell wrt the original mit.p, comment the line above for the ORNs
        ### and inject current below and compare with original mit.p
        ### comment out the test lines in self.run() and self.collate() also
        #inject = 0.5e-9 # 500pA injection
        #self.soma.inject = inject
        #printCellTree(self.cell)
        #print "______________________________________________________________"
        #print "The mit.p cell for comparison with the xml cell"
        #self.context.readCell('mit.p','/mitral2')
        #self.cell2 = moose.Cell("/mitral2")
        ## Since .p file has no way of knowing which channels are Ca-based or Ca-dependent, we need to do these two connections:
        #connect_CaConc( [moose.Neutral(comp) for comp in self.cell2.getChildren(self.cell2.id)] ) # assume all children of self.cell2 are compartments
        #printCellTree(self.cell2)
        #self.soma2 = moose.Compartment('/mitral2/soma')
        #self.soma2.inject = inject
        #self.soma2_vm = setupTable('soma2_vm', self.soma2,'Vm')
        
    def attachORNs(self):
        ##### Inhibitory Synapse NOT added as I am not sure if it is due to PG inhibition - Perhaps Cleland and Sethupathy also did not model it       
        ##### Excitatory + Inhibitory(not added) Synapse combo taken from the paper Djurisic etal 2008 JNeurosci.
        ##### Actually it's only an excitatory synapse, but they have used the inhibitory one to model the later time course.
        ##### Though this might be needed to account for PG cell inhibition?
        ##### Gbar-s have been set to ensure 16mV EPSP at the glom tuft as done in the paper.
        ##### Paper's defaults give only 8mV EPSP peak at glom tuft. So here multiplied by two.
        #####   Cannot use supplemental table values as no cytoplasmic resistivity Ri in this model.
        #####   Only axial resistance from glom to prim which doesn't matter much [maybe it does - loading rather than input resistance?].
        #####   Makes sense only if dendritic tree with many compartments.
        ##### Also no idea of how much to change the inhibitory part without Ri, so multiplied that also by 2

        synapse_list = []
        NML = NetworkML({'temperature':CELSIUS}) # has make_new_synapse()
        syn_name = self.ORNCellSynapse
        for segment in self.cellSegmentDict[self.cellname].values():
            # segment = [ segname,(proximalx,proximaly,proximalz),(distalx,distaly,distalz),\
            #    diameter,length,[potential_syn1, ... ] ]
            if syn_name in segment[5]: # is our ORN->cell synapse one of the potential synapses in this segment?
                syn_path = self.cell.name+'/'+segment[0]+'/'+syn_name
                if not self.context.exists(syn_path):
                    NML.make_new_synapse(syn_name, moose.Compartment(self.cell.name+'/'+segment[0]), syn_name)
                syn = moose.SynChan(syn_path)
                synapse_list.append(syn)

        for i in range(self.numORNSyns):
            ## Bad practice below - should be in NetworkML - no random numbers in simulations, only in generators.
            synapse = synapse_list[int(uniform(0,len(synapse_list)))]
            print "Connecting ", synapse.path
            spiketable = moose.TimeTable(synapse.path+'/tt'+str(i)) # unique spiketable (could be to same synapse)
            #### SynChan's synapse MsgDest takes time as its argument. Thus spiketable should contain a list of spike times.
            spiketable.connect("event", synapse,"synapse")
            synapse.setWeight(synapse.numSynapses-1, 2) # above added connection in synaptic array is set to weight 1
            self.spikeTableList.append(spiketable)
            ## constrate firefiles are NUM_SPIKEFILES in number, each of which has NUM_ORN_FILES_PER_GLOM lines.
            ## Each line can connect to a spiketable. Choose one randomly and connect. fileNumbers field is used by attach_spikes()
            spiketable.addField('fileNumbers')
            spiketable.setField( 'fileNumbers', str(int(uniform(0,NUM_ORN_FILES_PER_GLOM))) )

    def run(self):
        self.context.setClock(0, SIMDT, 0)
        self.context.setClock(1, SIMDT, 0) #### The hsolve and ee methods use clock 1
        self.context.setClock(2, SIMDT, 0) #### hsolve uses clock 2 for mg_block, nmdachan and others.
        self.context.setClock(PLOTCLOCK, PLOTDT, 0)
        
        firingrate = ORNfiringrates[mpirank-1]
        for i,spiketable in enumerate(self.spikeTableList):
            attach_spikes( '../firefiles/firefiles_constrate/firetimes_constrate'\
                +str(firingrate)+'_trial'+str(i%NUM_SPIKEFILES),\
                spiketable, 'cellORNtest'+str(mpirank) )
        self.context.reset() # A second reset will not work properly with tables!
        self.context.step(RUNTIME)
        if mpisize>1:
            print "sending for value ",firingrate,' from ',mpirank
            spiketimes = array(self.soma_vm)
            spiketimes_post = spiketimes[spiketimes>SETTLETIME]
            mpicomm.send( spiketimes_post, dest=boss, tag=0 )
            ### testing - channels
            #mpicomm.send( array(self.soma_var), dest=boss, tag=1 )
            ### testing - comparison with mit.p
            #print "sending 2 for value ",firingrate,' from ',mpirank
            #mpicomm.send( array(self.soma2_vm), dest=boss, tag=1 )
        else:
            figure()
            plot(array(self.soma_vm),'-k')
            show()

def collate():
    responsefreqlist = []
    for i,firingrate in enumerate(ORNfiringrates): # firingrate in Hz.
        soma_times = mpicomm.recv(source=i+1, tag=0)
        print "received from ",i+1," on boss machine ",mpiname
        if len(soma_times)>0:
            while soma_times[-1] == 0.0: # MOOSE inserts crap zeros at the end in TAB_SPIKE table
                soma_times = delete(soma_times,-1) # numpy's delete 0th element
                if len(soma_times) == 0: break
        responsefreqlist.append(len(soma_times)/REALRUNTIME)
        ## plot only if table's stepMode is not TAB_SPIKE
        #figure()
        #plot(soma_times,',r')
        ### testing - channels
        #soma_var = mpicomm.recv(source=i+1, tag=1)
        #figure()
        #plot(soma_var,'-b')
        ### testing - comparison with mit.p
        #soma2_vm = mpicomm.recv(source=i+1, tag=1)
        #print "received 2 from ",i+1," on boss machine ",mpiname
        #plot(soma2_vm,',b')
    figure()
    plot(ORNfiringrates,responsefreqlist,'+-r')
    show()
        
if __name__ == "__main__":
    seed([100.0])
    cellname = sys.argv[1]
    if mpisize>1:
        if mpirank==boss:
            collate()
        else:
            cell = CellTest(cellname)
            cell.run()
    else:
        cell = CellTest(cellname)
        cell.run()
        
