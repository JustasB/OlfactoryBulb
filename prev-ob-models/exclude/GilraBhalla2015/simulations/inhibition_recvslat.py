#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import pickle

sys.path.extend(["..","../networks","../generators","../simulations"])

from OBNetwork import *

from stimuliConstants import * # has SETTLETIME
from simset_inhibition import * # has REALRUNTIME
from sim_utils import * # has setup_tables(), plot_extras() and build_tweaks()
from data_utils import *

SETTLETIME = 500e-3 # increasing to 500ms from 250ms as PG cells spike I think right at the start.
RUNTIME = REALRUNTIME + SETTLETIME

# Whether to inject the current into tuft or soma, typically soma
TUFT_INJECT = False

from pylab import * # part of matplotlib that depends on numpy but not scipy

## from node000:
## mpiexec -machinefile ~/hostfile -n 3 ~/Python-2.6.4/bin/python2.6 inhibition_recvslat.py
## or from node059:
## mpiexec -machinefile ~/hostfile_end5 -n 3 ~/Python-2.6.4/bin/python2.6 inhibition_recvslat.py

#-----------------------------------------------------------

class lateral_granule_inhibition:
    
    def __init__(self):
        pass

    def setup_stim(self,network,rank):
        if SPIKEBLOCK:
            print "Blocking Na in mitral cells only"\
                " (don't have granules, else inconsistent) ..."
            for cell in network.mitralTable.values():
                #blockChannels(cell, ['Na','K2']) # only K2 i.e. Kfast, not all K
                blockChannels(cell, ['Na']) # only TTX, only Na blocked not K2
        iAinject = offInject
        if rank == 1:
            iBinject = 0.0
        else:
            iBinject = onInject
        ipulse_duration = 400e-3 # seconds ## for self and rec inh: Arevian et al
        #ipulse_duration = 125e-3 # seconds ## only for self inh: 40Hz for 125ms -- Friedman & Strowbridge 2000
        ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
        ## I need to use much larger currents (inhibition folded into mitral model)
        if TUFT_INJECT:
            tuftseg = 'Seg0_glom_81_102'
            print "injecting into tuft segment",tuftseg
            compA = moose.Compartment(network.mitralTable[mitralidx].path+'/'+tuftseg)
            compB = moose.Compartment(network.mitralTable[mitralsidekickidx].path+'/'+tuftseg)
        else:
            compA = network.mitralTable[mitralidx].soma
            compB = network.mitralTable[mitralsidekickidx].soma
        iA = setup_iclamp(compA, '_mitralA', SETTLETIME, ipulse_duration, iAinject)
        ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
        ## slightly stagger the start of current pulses in the two cells
        ## so that the mitrals do not continuously co-fire.
        iB = setup_iclamp(compB,'_mitralB', SETTLETIME-5e-3, ipulse_duration+5e-3, iBinject)
        ## To check backpropagating AP, check at various distances
        ## from the soma on the lateral dendrites.
        dend_seg = moose.Compartment(network.mitralTable[mitralidx].path+\
            #'/Seg0_sec_dendp4_0_254') # at 43.86 microns from soma
            #'/Seg0_sec_dendd3_0_204') # at 190.56 microns from soma
            '/Seg0_sec_dendd4_2_269') # at 1004.47 microns from soma
        self.mitAdendTable = setupTable('mitdendTable',dend_seg,'Vm')
        dend_seg = moose.Compartment(network.mitralTable[mitralsidekickidx].path+\
            #'/Seg0_sec_dendp4_0_254') # at 43.86 microns from soma
            #'/Seg0_sec_dendd3_0_204') # at 190.56 microns from soma
            '/Seg0_sec_dendd4_2_269') # at 1004.47 microns from soma
        self.mitBdendTable = setupTable('mitdendTable',dend_seg,'Vm')
        print "Glomerulus segment of mitral A = ", network.mitralTable[mitralidx].glom.path
        print "Glomerulus segment of mitral B = ", network.mitralTable[mitralsidekickidx].glom.path
        print 'Injecting mitral A with '+str(iAinject)+' and B with '+str(iBinject)+' at rank '+str(rank)

    def run_inhibition(self,network):
        resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
        network.context.step(RUNTIME)
        return (array(network.mitralTable[mitralidx]._vmTableSoma),\
            array(network.mitralTable[mitralsidekickidx]._vmTableSoma))

if __name__ == "__main__":
    includeProjections = ['granule_baseline']

    if mpirank == boss:
        mit_responses = []
        for procnum in [1,2]:
            mitral_Vm_both = mpicomm.recv(source=procnum, tag=0)
            mit_responses.append(mitral_Vm_both)

        ## write results to a file
        outfilename = 'data_recvslat.pickle'
        f = open(outfilename,'w')
        pickle.dump(mit_responses, f)
        f.close()
        print "Wrote", outfilename

        timevec = arange(0.0,RUNTIME+1e-10,SIMDT)*1e3
        fig = figure(facecolor='w')
        ax = fig.add_subplot(111)
        title('mitral A',fontsize=30)
        #plot(plotBins(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        #plot(plotSpikes(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        ax.plot(timevec,mit_responses[0][0]*1e3,color=(0,0,0),linewidth=2,label='Recurrent')
        ax.plot(timevec,mit_responses[1][0]*1e3,color=(0,0,1),linewidth=2,label='Lateral')
        biglegend()
        axes_labels(ax,'time (ms)','Vm (ms)')

        fig = figure(facecolor='w')
        ax = fig.add_subplot(111)
        title('mitral B',fontsize=30)
        #plot(plotBins(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        #plot(plotSpikes(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        ax.plot(timevec,mit_responses[0][1]*1e3,color=(0,0,0),linewidth=2,label='Recurrent')
        ax.plot(timevec,mit_responses[1][1]*1e3,color=(0,0,1),linewidth=2,label='Lateral')
        biglegend()
        axes_labels(ax,'time (ms)','Vm (ms)')
        
        show()
    else:
        ## includeProjections gets used only if ONLY_TWO_MITS is True:
        ## Keep below projections to 'second order cells'
        ## i.e. to cells (granules) connected to mits0&1.
        ## The connections between second order cell
        ## and mits0&1 are automatically retained of course.
        ## no need for 'PG' below as 'ORN_PG' and 'SA_PG' are not needed,
        ## and 'PG_mitral', 'mitral_PG' connections to/from mits0&1 are kept automatically.
        includeProjections = ['granule_baseline']
        tweaks = build_tweaks( CLUB_MITRALS, NO_SPINE_INH,\
            NO_SINGLES, NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
            includeProjections, (mitralidx,mitralsidekickidx) )
        network = OBNetwork(OBNet_file, synchan_activation_correction,\
            tweaks, mpirank, "recvslat", granfilebase, spiketable=False)
        #printNetTree() # from moose_utils.py

        sim =  lateral_granule_inhibition()
        sim.setup_stim(network, mpirank)
        mitral_responses_both = sim.run_inhibition(network)
        mpicomm.send( mitral_responses_both, dest=boss, tag=0 )
