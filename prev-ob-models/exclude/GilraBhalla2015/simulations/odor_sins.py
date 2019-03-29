#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle
 
## from node000:
## mpiexec -machinefile ~/hostfile -n <num_sins*numtrials+1> ~/Python-2.6.4/bin/python2.6 odor_sins.py
## nohup mpiexec -machinefile ~/hostfile -n 201 ~/Python-2.6.4/bin/python2.6 odor_sins.py < /dev/null &
## typical value for num_sins = 5
## typical value for num of trials = 40
## (depends on number of available processing nodes and number of odorfiles generated)

## Set various option like NO_PGs or ONLY_TWO_MITS in simset_odor.py

## single simulation USAGE:
## python2.6 odor_sins.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, varied...
from simset_odor import *

from pylab import * # part of matplotlib that depends on numpy but not scipy
from plot_odor_sins import *

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self,network,sinnum,trialnum):
        self.setupOdor(network, sinnum, trialnum)        
        print "Setup sin frequency =",sine_frequencies[sinnum],"at",self.mpirank

    def setupOdor(self, network, sinnum, trailnum):
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
                    self.attach_appropriate_files_to_glom(proj[0],proj[2],sinnum,trialnum,glomnum)
            ## connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    self.attach_appropriate_files_to_glom(proj[0],proj[2],sinnum,trialnum,glomnum)

    def attach_appropriate_files_to_glom(self,proj1,proj2,sinnum,trialnum,glomnum):
        if LAT_SINS and glomnum == central_glom:
            filebase = '../firefiles/firefiles_constrate/firetimes_constrate'+str(SIN_MAIN_CONSTRATE)
            self.attach_files_uniquely(filebase,proj1,proj2,-1,trialnum) # sinnum=-1 to not use sinusoid
        else:
            filebase = '../firefiles/firefiles_sins/firetimes_sin_glom'+str(glomnum)
            self.attach_files_uniquely(filebase,proj1,proj2,sinnum,trialnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,sinnum=None,trialnum=None):
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
        ## If sinnum == -1, it means do not use sinusoids.
        if sinnum is not None and sinnum != -1: filebase += '_fnum'+str(sinnum)
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
        network.context.step(SIN_RUNTIME)
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

    #### if only one process is called, plot one sim directly
    if mpisize == 1:
        #### run the slave processes
        sim = odorResponse()
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
        sinnum = 0
        sim.setupStim(network, sinnum, trialnum)
        mitral_responses = sim.run(network, binned=False)
        timevec = arange(SETTLETIME,SETTLETIME+SIN_RUNTIME+3*PLOTDT/2.0,PLOTDT)
        plot(timevec, mitral_responses[0])
        show()

    #### multiple processes
    else:
        if mpirank == boss:
            #### collate at boss process
            ## mitral_responses_list[avgnum][sinnum][mitnum][spikenum]
            mitral_responses_list = []
            numavgs = (mpisize-1)/num_sins
            for avgnum in range(numavgs):
                response_set = []
                for sinnum in range(num_sins):
                    procnum = avgnum*num_sins + sinnum + 1
                    print 'waiting for process '+str(procnum)+'.'
                    ## below: you get a numpy array of 
                    ## rows=NUM_GLOMS*MIT_SISTERS and cols=spike times
                    ## mitral responses has spike times
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
            if NO_LATERAL: lat_str = '_NOLAT'
            else: lat_str = ''
            now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
            outfilename = '../results/odor_sins/'+now+'_sins'+singles_str+\
                joints_str+pgs_str+lat_str+'_numgloms'+str(NUM_GLOMS)+'.pickle'
            f = open(outfilename,'w')
            pickle.dump(mitral_responses_list, f)
            f.close()
            print "Wrote", outfilename
            
            plot_sin_responses(outfilename)
            show()

        else:
            #### run the slave processes
            sim = odorResponse()
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

            trialnum = (mpirank-1)/num_sins
            sinnum = (mpirank-1)%num_sins
            sim.setupStim(network, sinnum, trialnum)
            mitral_responses = sim.run(network, binned=True)
            mpicomm.send( mitral_responses, dest=boss, tag=0 )
            print 'sent from process',mpirank
