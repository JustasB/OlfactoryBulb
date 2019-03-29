#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle
 
## from node000:
## mpiexec -machinefile ~/hostfile -n <numtrains*numtrials+1> ~/Python-2.6.4/bin/python2.6 odor_whitenoise.py <simtype>
## nohup mpiexec -machinefile ~/hostfile -n 251 ~/Python-2.6.4/bin/python2.6 odor_whitenoise.py -1 < /dev/null &
## typical value for numtrains = 250
## typical value for num of trials = 1
## (depends on number of available processing nodes and number of odorfiles generated)

## Need to set the type of simulation i.e. whether for mainmitral kernel (simtype = -1)
## or simtype = 0..6 index of the varied_mainrate, etc. variable in stimuliConstants.py

## See varied_mainrate, varied_distance, varied_gran_baseline in stimuliConstants.py after
## setting the option of 'varied' to one of ('mainrate', 'distance', 'gran_baseline') )
## Set various option like NO_PGs or ONLY_TWO_MITS in simset_odor.py

## single simulation USAGE:
## python2.6 odor_whitenoise.py -1

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, varied...
from simset_odor import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self, simtype):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()
        self.simtype = int(simtype)

    def setupStim(self,network,trainnum,trialnum):
        self.setupOdor(network, trainnum, trialnum)        
        print "Setup trainnum =",trainnum,"at",self.mpirank

    def setupOdor(self, network, trainnum, trailnum):
        ## first figure out which PG belongs to which glom
        ## PG_glom_map[pgname] returns the glom num of the PG:
        ## needed for ORN to PG connections.
        PG_glom_map = {}
        for projname in network.projectionDict.keys():
            if 'PG_mitral' in projname:
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    ## get the glomnum from the post path proj[2]
                    ## name of the mitral cell from '/mitrals_2/...'
                    mitname = string.split(proj[2],'/')[1]
                    ## glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    glomnum = int(string.split(mitname,'_')[1]) / 2
                    ## name of the PG cell from '/PGs_2/...'
                    pgname = string.split(proj[1],'/')[1]
                    PG_glom_map[pgname] = glomnum
        ## Now connect the ORNs
        for projname in network.projectionDict.keys():
            #### Calling attach_spikes() for each projection,
            #### would reconnect files to the same segment multiple times.
            #### But attach_files_uniquely() checks whether timetable.tableSize is zero or not
            #### i.e. files already attached or not.
            ## connect ORNs to mitrals
            if 'ORN_mitral' in projname:
                print "connecting ORN files to mitrals"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    ## get the glomnum from the post path proj[2]
                    ## name of the mitral cell from '/mitrals_2/...'
                    mitname = string.split(proj[2],'/')[1]
                    ## glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    glomnum = int(string.split(mitname,'_')[1]) / 2
                    self.attach_appropriate_files_to_glom(proj[0],proj[2],trainnum,trialnum,glomnum)
            ## connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    self.attach_appropriate_files_to_glom(proj[0],proj[2],trainnum,trialnum,glomnum)

    def attach_appropriate_files_to_glom(self,proj1,proj2,trainnum,trialnum,glomnum):
        ## for the excitatory kernel
        if self.simtype == -1:
            ## only connect glom0 to noise train to get exc kernel
            if glomnum == central_glom:
                filebase = ORNpathINHstr+'firetimes_whitenoise_glom'+str(glomnum)
                self.attach_files_uniquely(filebase,proj1,proj2,trainnum,trialnum)
        ## for inhibitory kernel (simtype != -1)
        else:
            ## connect glom0 to a constant firing rate
            if glomnum == central_glom:
                if varied == 'mainrate':
                    firingrate = varied_mainrate[self.simtype]
                    filebase = ORNpathINHstr+'firetimes_constrate'\
                        +str(firingrate)+'_avg'+str(trainnum)
                    self.attach_files_uniquely(filebase,proj1,proj2)                                
            ## connect other glom to noise train to get inh kernel
            else:
                filebase = ORNpathINHstr+'firetimes_whitenoise_glom'+str(glomnum)
                self.attach_files_uniquely(filebase,proj1,proj2,trainnum,trialnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,trainnum=None,trialnum=None):
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
        if trainnum is not None: filebase += '_train'+str(trainnum)
        if trialnum is not None: filebase += '_trial'+str(trialnum)
        ## attach_spikes() accesses filenumbers to this segment
        ## from 'fileNumbers' field (of the timetable object in MOOSE)
        ## which is created while reading in networkML.
        attach_spikes(filebase, tt, self.mpirank)

    def run(self,network, binned):
        print "Resetting MOOSE."
        # from moose_utils.py sets clocks and resets
        resetSim(network.context, SIMDT, PLOTDT)
        print "Running at",self.mpirank
        network.context.step(PULSE_RUNTIME)
        mitral_responses = []
        if ONLY_TWO_MITS: mits = [mitralidx, mitralsidekickidx]
        else: mits = range(NUM_GLOMS*MIT_SISTERS)
        ## network.mitralTable is a dictionary.
        for mitnum in mits:
            mitral = network.mitralTable[mitnum]
            ## need to convert to numpy's array(),
            ## else MOOSE table cannot be pickled for mpi4py send()
            mitral_responses.append(array(mitral._vmTableSoma))
        return mitral_responses

#----------------------------------------------------------------

if __name__ == "__main__":
    ## simtype == '-1': main mitral kernel
    ## other simtype gives the index of varied... variable
    ## to set the conditions for inh kernel due to nearby glomerulus.
    simtype = sys.argv[1]

    #### if only one process is called, plot one sim directly
    if mpisize == 1:
        #### run the slave processes
        sim = odorResponse(simtype)
        ## includeProjections gets used only if ONLY_TWO_MITS is True:
        ## Keep below projections to 'second order cells'
        ## i.e. to cells (granules) connected to mits0&1.
        ## The connections between second order cell
        ## and mits0&1 are automatically retained of course.
        # 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
        includeProjections = ['PG','granule_baseline']
        tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES, NO_JOINTS,\
            NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections, mitralsidekickidx)
        network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
            mpirank, granfilebase+'_noresp', spiketable=False)
        #printNetTree() # from moose_utils.py

        trialnum = 0
        trainnum = 0
        sim.setupStim(network, trainnum, trialnum)
        mitral_responses = sim.run(network, binned=False)
        timevec = arange(SETTLETIME,SETTLETIME+PULSE_RUNTIME+3*PLOTDT/2.0,PLOTDT)
        plot(timevec, mitral_responses[0])
        show()

    #### multiple processes
    else:
        if mpirank == boss:
            #### collate at boss process
            ## mitral_responses_list[avgnum][trainnum][mitnum][spikenum]
            mitral_responses_list = []
            numavgs = (mpisize-1)/NUMWHITETRAINS
            for avgnum in range(numavgs):
                response_set = []
                for trainnum in range(NUMWHITETRAINS):
                    procnum = avgnum*NUMWHITETRAINS + trainnum + 1
                    print 'waiting for process '+str(procnum)+'.'
                    ## below: you get a numpy array of 
                    ## rows=NUM_GLOMS*MIT_SISTERS and cols=spike times
                    ## mitral responses has spike times
                    ## we calculate STA to get kernel from spike times.
                    mitral_responses = mpicomm.recv(source=procnum, tag=0)
                    response_set.append( mitral_responses )
                mitral_responses_list.append(response_set)
            
            # write results to a file
            today = datetime.date.today()
            if NO_SINGLES: singles_str = '_NOSINGLES'
            else: singles_str = '_SINGLES'
            if NO_JOINTS: joints_str = '_NOJOINTS'
            else: joints_str = '_JOINTS'
            if NO_PGS: pgs_str = '_NOPGS'
            else: pgs_str = '_PGS'
            now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
            outfilename = '../results/odor_whitenoise/'+now+'_whitenoise'+singles_str+\
                joints_str+pgs_str+'_numgloms'+str(NUM_GLOMS)+'_simtype'+simtype+'.pickle'
            f = open(outfilename,'w')
            pickle.dump(mitral_responses_list, f)
            f.close()
            print "Wrote", outfilename

        else:
            #### run the slave processes
            sim = odorResponse(simtype)
            ## includeProjections gets used only if ONLY_TWO_MITS is True:
            ## Keep below projections to 'second order cells'
            ## i.e. to cells (granules) connected to mits0&1.
            ## The connections between second order cell
            ## and mits0&1 are automatically retained of course.
            # 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
            includeProjections = ['PG','granule_baseline']
            tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES, NO_JOINTS,\
                NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections, mitralsidekickidx)
            network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
                mpirank, granfilebase+'_noresp', spiketable=True)
            #printNetTree() # from moose_utils.py

            trialnum = (mpirank-1)/NUMWHITETRAINS
            trainnum = (mpirank-1)%NUMWHITETRAINS
            sim.setupStim(network, trainnum, trialnum)
            mitral_responses = sim.run(network, binned=True)
            mpicomm.send( mitral_responses, dest=boss, tag=0 )
            print 'sent from process',mpirank
