#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle

## python2.6 tuft_nonlinearity.py
## Remove all interneurons in simset_odor

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has REALRUNTIME, NUMBINS
from generate_constrate_files import * # has the firingrates list

RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy
from plot_odor_morphs import *

#-----------------------------------------------------------

class constOdorResponse:
    
    def __init__(self):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self, network, avgnum):
        ### connect the ORNs
        for projname in network.projectionDict.keys():
            #### Calling attach_spikes() for each projection,
            #### would reconnect files to the same segment multiple times.
            #### But attach_files_uniquely() checks whether timetable.tableSize is zero or not
            #### i.e. files already attached or not.
            ### connect ORNs to mitrals
            if 'ORN_mitral' in projname:
                print "connecting ORN files to mitrals"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    #get the glomnum from the post path proj[2]
                    mitname = string.split(proj[2],'/')[1] # name of the mitral cell from '/mitrals_2/...'
                    mitnum = int(string.split(mitname,'_')[1]) # mitral number from 'mitrals_2'
                    firingrate = firingrates[mitnum]
                    filebase = '../firefiles/firetimes_constrate_'+str(firingrate)+'.txt'
                    self.attach_files_uniquely(filebase,proj[0],proj[2],avgnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,avgnum=None):
        ttpath = postsegpath+'/'+synname+'_tt'
        if self.context.exists(ttpath):
            # timetable already created by networkML reader - just wrap it below.
            tt = moose.TimeTable(ttpath) # post_segment_path+'/'+syn_name+'_tt'
        else:
            ## if timetable was not already created by networkML reader,
            ## it means that the synaptic weights must be zero!
            ## (no extra inhibition - only main inhibition)
            ## hence do not attach spikefiles
            return
        if tt.tableSize != 0: return # if files are already attached, do nothing!
        filebase += '_odor_'+str(odorA)+'_'+str(odorB)
        if avgnum is not None: filebase += '_avgnum'+str(avgnum)
        ## attach_spikes() accesses filenumbers to this segment
        ## from 'fileNumbers' field (of the timetable object in MOOSE)
        ## which is created while reading in networkML.
        attach_spikes(filebase, tt, self.mpirank)

    def connect_granule_baselines_to_files(self, network, avgnum):
        for projname in network.projectionDict.keys():
            if 'granule_baseline' in projname:
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    ## just wrap the existing timetable created when loading in neuroml
                    ## it has a field fileNumbers that has which lines to take from the firingfile
                    tt = moose.TimeTable(proj[2]+'/'+proj[0]+'_tt') # post_segment_path+'/'+syn_name+'_tt'
                    if IN_VIVO:
                        if gran_base_resp_tuned:
                            filebase = '../firefiles/firetimes_gran_baseline_'+str(avgnum)
                        else:
                            filebase = '../firefiles/firetimes_gran_baseline_noresp_'+str(avgnum)
                    else:
                        filebase = '../firefiles/firetimes_gran_baseline_invitro_'+str(avgnum)
                    attach_spikes(filebase, tt, self.mpirank)

    def run(self,network, binned):
        print "Resetting MOOSE."
        # from moose_utils.py sets clocks and resets
        resetSim(network.context, SIMDT, PLOTDT)
        print "Running at",self.mpirank
        network.context.step(RUNTIME)
        mitral_responses = []
        mitral_responses_binned = []
        if ONLY_TWO_MITS: num_mits = MIT_SISTERS
        else: num_mits = NUM_GLOMS*MIT_SISTERS
        for mitnum in range(num_mits):
            mitral = network.mitralTable[mitnum]
            ## only the last respiration cycle is taken
            if binned: mitral_responses_binned.append(
                plotBins(mitral._vmTableSoma, NUMBINS, RUNTIME, (NUM_RESPS-1)*RESPIRATION+SETTLETIME) )
            ## need to convert to numpy's array(),
            ## else MOOSE table cannot be pickled for mpi4py send()
            mitral_responses.append(array(mitral._vmTableSoma))
        return (mitral_responses,mitral_responses_binned)

#----------------------------------------------------------------

if __name__ == "__main__":

    sim =  constOdorResponse()
    includeProjections = []
    tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES,\
        NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections)
    BINNED = True
    ## if not BINNED, save the full mitral Vm-s
    ## and not just their spiketimes by setting spiketable = False below.
    network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
        mpirank=mpirank, invivo=IN_VIVO, spiketable=BINNED)
    #printNetTree() # from moose_utils.py
    ## monitor those interneurons that are connected to mitral indices 0 and 1
    ## save only spiketimes by setting extras_spikes_only=True
    sim.setupStim(network, avgnum=0)
    mitral_responses,mitral_responses_binned = sim.run(network,BINNED)


    figure()
    title('Glomerulus 0')
    if BINNED:
        deltabin = RESPIRATION/NUMBINS
        # Take only the last respiration cycle
        timevec = arange(SETTLETIME+(NUM_RESPS-1)*RESPIRATION+deltabin/2,RUNTIME,deltabin)
        mitral_responses = mitral_responses_binned
    else:
        timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
    plot(timevec,mitral_responses[0],color=(0.0,1.0,0.0))
    plot(timevec,mitral_responses[1],color=(0.0,1.0,0.5))
    show()
