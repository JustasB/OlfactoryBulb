#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle

## from node000:
## mpiexec -machinefile ~/hostfile -n <numavgs*numvarinhs+1> ~/Python-2.6.4/bin/python2.6 odor_variedinhibition.py
## nohup mpiexec -machinefile ~/hostfile -n 61 ~/Python-2.6.4/bin/python2.6 odor_variedinhibition.py < /dev/null &
## typical value for numavgs = 10
## (depends on number of available processing nodes and number of odorfiles generated)
## typical value for numvariedinh = 6 (see generate_firerates_variedinhibition.py).
## Set various options like NO_PGs or ONLY_TWO_MITS in simset_odor
## For a single run:
## python2.6 odor_variedinhibition.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has REALRUNTIME, NUMBINS

RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy
from plot_odor_variedinhibition import *

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self,network,inhnum,avgnum):
        self.setupOdor(network, inhnum, avgnum)        
        print "Setup inhnum =",inhnum,"at",self.mpirank

    def setupOdor(self, network, inhnum, avgnum):
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
                    filebase = ORNpathINHstr+'firetimes_2sgm_glom_'+str(glomnum)
                    ## same input to glom0, but varied input to other gloms
                    if glomnum==0: thisinhnum = 0
                    else: thisinhnum = inhnum
                    self.attach_files_uniquely(filebase,proj[0],proj[2],thisinhnum,avgnum)
            ## connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathINHstr+'firetimes_2sgm_glom_'+str(glomnum)
                    ## same input to glom0, but varied input to other gloms
                    if glomnum==0: thisinhnum = 0
                    else: thisinhnum = inhnum
                    self.attach_files_uniquely(filebase,proj[0],proj[2],thisinhnum,avgnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,inhnum,avgnum=None):
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
        filebase += '_inhnum'+str(inhnum)
        if avgnum is not None: filebase += '_avgnum'+str(avgnum)
        ## attach_spikes() accesses filenumbers to this segment
        ## from 'fileNumbers' field (of the timetable object in MOOSE)
        ## which is created while reading in networkML.
        attach_spikes(filebase, tt, self.mpirank)

    def run(self,network, binned):
        print "Resetting MOOSE."
        # from moose_utils.py sets clocks and resets
        resetSim(network.context, SIMDT, PLOTDT)
        print "Running at",self.mpirank
        network.context.step(RUNTIME)
        mitral_responses = []
        mitral_responses_binned = []
        self.mitseg_responses = []
        if ONLY_TWO_MITS: mits = [mitralidx, mitralsidekickidx]
        else: mits = range(NUM_GLOMS*MIT_SISTERS)
        ## network.mitralTable is a dictionary.
        for mitnum in mits:
            ## mitralTable is a dict with mitnum as key
            mitral = network.mitralTable[mitnum]
            ## only the last respiration cycle is taken
            if binned: mitral_responses_binned.append(
                plotBins(mitral._vmTableSoma, NUMBINS, RUNTIME,\
                (NUM_RESPS-1)*RESPIRATION+SETTLETIME) )
            ## need to convert to numpy's array(),
            ## else MOOSE table cannot be pickled for mpi4py send()
            mitral_responses.append(array(mitral._vmTableSoma))
            #self.mitseg_responses.append(\
            #    (array(mitral._vmTableSoma),array(mitral._vmTableDend)))
        return (mitral_responses,mitral_responses_binned)

#----------------------------------------------------------------

if __name__ == "__main__":

    #### if only one process is called, plot one odor directly
    if mpisize == 1:
        sim =  odorResponse()
        ## includeProjections gets used only if ONLY_TWO_MITS is True:
        ## Keep below projections to 'second order cells'
        ## i.e. to cells (granules) connected to mits0&1.
        ## The connections between second order cell
        ## and mits0&1 are automatically retained of course.
        ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
        includeProjections = ['PG','granule_baseline']
        tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES, NO_JOINTS,\
            NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections, mitralsidekickidx)
        BINNED = False
        ## if not BINNED, save the full mitral Vm-s
        ## and not just their spiketimes by setting spiketable = False below.
        network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
            mpirank, granfilebase, spiketable=BINNED)
        #printNetTree() # from moose_utils.py

        ## monitor those interneurons that are connected to mitral indices 0 and 1
        ## save only spiketimes by setting extras_spikes_only=True
        #extras_spikes_only = True
        #tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS,\
        #    {'mitrals':[0,1]}, spikes=extras_spikes_only)

        avgnum = 0
        inhnum = 5
        sim.setupStim(network, inhnum, avgnum)
        mitral_responses,mitral_responses_binned = sim.run(network,binned=BINNED)
        #if not extras_spikes_only:
        #    timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
        #    plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS)
        #else:
        #    deltabin = RESPIRATION/NUMBINS
        #    ## Only the last respiration cycle
        #    timevec = arange(SETTLETIME+(NUM_RESPS-1)*RESPIRATION+deltabin/2,RUNTIME,deltabin)
        #    plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS,\
        #        NUMBINS, RUNTIME, SETTLETIME)
        figure()
        title('mitrals 0 and 2')
        if BINNED:
            deltabin = RESPIRATION/NUMBINS
            # Take only the last respiration cycle
            timevec = arange(SETTLETIME+(NUM_RESPS-1)*RESPIRATION+deltabin/2,RUNTIME,deltabin)
            mitral_responses = mitral_responses_binned
        else:
            timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
        plot(timevec,mitral_responses[0],color=(1.0,0.0,0.0))
        plot(timevec,mitral_responses[1],color=(0.0,1.0,0.0))
        #figure()
        #title('mitral 2')
        #plot(timevec,sim.mitseg_responses[1][0],color=(1.0,0.0,0.0))
        #plot(timevec,sim.mitseg_responses[1][1],color=(0.0,1.0,0.0))
        show()

    #### if multiple processes are called, average over odor morphs
    else:

        numodors = len(inputList)
        if mpirank == boss:
            #### collate at boss process
            mitral_responses_list = []
            mitral_responses_binned_list = []
            numavgs = (mpisize-1)/NUMINHS
            for avgnum in range(numavgs):
                response_odorset = []
                response_odorset_binned = []
                for inhnum in range(NUMINHS):
                    procnum = avgnum*NUMINHS + inhnum + 1
                    print 'waiting for process '+str(procnum)+'.'
                    ## below: you get a numpy array of 
                    ## rows=NUM_GLOMS*MIT_SISTERS and cols=NUMBINS
                    ## mitral responses has spike times,
                    ## mitral_responses_binned has binned firing rates
                    mitral_responses,mitral_responses_binned = \
                        mpicomm.recv(source=procnum, tag=0)
                    response_odorset.append( mitral_responses )
                    response_odorset_binned.append( mitral_responses_binned )
                mitral_responses_list.append(response_odorset)
                mitral_responses_binned_list.append(response_odorset_binned)
            
            # write results to a file
            today = datetime.date.today()
            if NO_SINGLES: singles_str = '_NOSINGLES'
            else: singles_str = '_SINGLES'
            if NO_JOINTS: joints_str = '_NOJOINTS'
            else: joints_str = '_JOINTS'
            if NO_PGS: pgs_str = '_NOPGS'
            else: pgs_str = '_PGS'
            now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
            outfilename = '../results/odor_varinh/'+now+'_odorvarinh'+singles_str+\
                joints_str+pgs_str+'_numgloms'+str(NUM_GLOMS)+'.pickle'
            f = open(outfilename,'w')
            pickle.dump((mitral_responses_list,mitral_responses_binned_list), f)
            f.close()
            print "Wrote", outfilename
            
            plot_varinh(outfilename)
            show()

        else:
            #### run the slave processes
            sim =  odorResponse()
            ## includeProjections gets used only if ONLY_TWO_MITS is True:
            ## Keep below projections to 'second order cells'
            ## i.e. to cells (granules) connected to mits0&1.
            ## The connections between second order cell
            ## and mits0&1 are automatically retained of course.
            ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
            includeProjections = ['PG','granule_baseline']
            tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES, NO_JOINTS,\
                NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections, mitralsidekickidx)
            network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
                mpirank, granfilebase, spiketable=True)
            #printNetTree() # from moose_utils.py

            avgnum = (mpirank-1)/NUMINHS
            inhnum = (mpirank-1)%NUMINHS
            sim.setupStim(network, inhnum, avgnum)
            mitral_responses_both = sim.run(network, binned=True)
            mpicomm.send( mitral_responses_both, dest=boss, tag=0 )
            print 'sent from process',mpirank
