#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle

## from node000:
## mpiexec -machinefile ~/hostfile -n <numavgs*numodors+1> ~/Python-2.6.4/bin/python2.6 odor_morphs.py
## nohup mpiexec -machinefile ~/hostfile -n 57 ~/Python-2.6.4/bin/python2.6 odor_morphs.py < /dev/null &
## typical value for numavgs = 8 
## (depends on number of available processing nodes and number of odorfiles generated)
## typical value for numodors = len(inputList) = 7.
## OR for a single odor run; from any node (or laptop -- can use system python (2.7)):
## python2.6 odor_morphs.py [SAVEFORMOVIE]
## Set various option like NO_PGs or ONLY_TWO_MITS in simset_odor

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has REALRUNTIME, NUMBINS

RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy
from plot_odor_morphs import *

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self,mpirank=mpirank): # mpirank is defined in data_utils.py
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self,network,args,avgnum):
        odorA, odorB = args[0]
        self.setupOdor(network, odorA, odorB, avgnum)        
        print "Setup odorA =",odorA,"odorB =",odorB,"at",self.mpirank

    def setupOdor(self, network, odorA, odorB, avgnum):
        ### first figure out which PG belongs to which glom
        ### PG_glom_map[pgname] returns the glom num of the PG: needed for ORN to PG connections.
        PG_glom_map = {}
        for projname in network.projectionDict.keys():
            if 'PG_mitral' in projname:
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    # get the glomnum from the post path proj[2]
                    # name of the mitral cell from '/mitrals_2/...'
                    mitname = string.split(proj[2],'/')[1]
                    # glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    glomnum = int(string.split(mitname,'_')[1]) / 2
                    # name of the PG cell from '/PGs_2/...'
                    pgname = string.split(proj[1],'/')[1]
                    PG_glom_map[pgname] = glomnum
        ### Now connect the ORNs
        for projname in network.projectionDict.keys():
            #### Calling attach_spikes() for each projection,
            #### would reconnect files to the same segment multiple times.
            #### But attach_files_uniquely() checks whether timetable.tableSize is zero or not
            #### i.e. files already attached or not.
            ### connect ORNs to mitrals
            if 'ORN_mitral' in projname:
                print "connecting ORN files to mitrals"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    # get the glomnum from the post path proj[2]
                    mitname = string.split(proj[2],'/')[1] # name of the mitral cell from '/mitrals_2/...'
                    glomnum = int(string.split(mitname,'_')[1]) / 2 # glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    filebase = ORNpathseedstr+'firetimes_2sgm_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],odorA,odorB,avgnum)
            ### connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathseedstr+'firetimes_2sgm_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],odorA,odorB,avgnum)
            ### connect SAs to PG
            if 'SA_PG' in projname:
                print "connecting SA files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathseedstr+'firetimes_SA'
                    self.attach_files_uniquely(filebase,proj[0],proj[2],odorA,odorB)
            ###### I am back to 'extra-connecting' modelled mitral as extra sister mitrals excitation to granules
            ###### Previously, as below, I was connecting ORNs of the glom to granules
            ###### which caused inhibition even when the sister mitrals were not even firing!
            ### connect unmodelled extra sister mitrals as files to granules
            #if 'mitral_granule_extra' in projname:
            #    print "Connecting unmodelled sister excitation files to granules"
            #    for i,proj in enumerate(network.projectionDict[projname][2]):
            #        granulename = string.split(proj[2],'/')[1] # name of the granule cell from '/granules_singles_2/...'
            #        # glomnum from pre_path = proj[1] = 'file[+<glomnum>]_<filenumber1>[_<filenumber2>...]'
            #        glomstr = proj[1].split('+')[1].split('_',1)[0]
            #        filebase = ORNpathseedstr+'firetimes_2sgm_glom_'+glomstr
            #        self.attach_files_uniquely(filebase,proj[0]+'_'+glomstr,proj[2],odorA,odorB,avgnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,odorA,odorB,avgnum=None):
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
        attach_spikes(filebase, tt, uniquestr+str(self.mpirank))

    def run(self,network, binned):
        print "Resetting MOOSE."
        # from moose_utils.py sets clocks and resets
        resetSim(network.context, SIMDT, PLOTDT)
        print "Running at",self.mpirank
        network.context.step(RUNTIME)
        mitral_responses = []
        mitral_responses_binned = []
        if ONLY_TWO_MITS or NO_LATERAL: num_mits = MIT_SISTERS
        else: num_mits = NUM_GLOMS*MIT_SISTERS
        for mitnum in range(num_mits):
            mitral = network.mitralTable[mitnum]
            ## only the last respiration cycle is taken
            if binned: mitral_responses_binned.append(
                plotBins(mitral._vmTableSoma, NUMBINS, RUNTIME,\
                (NUM_RESPS-1)*RESPIRATION+SETTLETIME) )
            ## need to convert to numpy's array(),
            ## else MOOSE table cannot be pickled for mpi4py send()
            mitral_responses.append(array(mitral._vmTableSoma))
        return (mitral_responses,mitral_responses_binned)

#----------------------------------------------------------------

if __name__ == "__main__":

    ## uniquestr to put in every temp filename to avoid clashing with other processes
    if len(sys.argv)>2: uniquestr = sys.argv[2]+'_' # _ necessary, else say 'morphs2'+mpirank is screwed
    else: uniquestr = 'morphs_'
    
    numodors = len(inputList)
    #### if only one process is called, plot one odor directly
    if mpisize == 1:
        ## for ref: inputList = [ (0.0,1.0), (0.2,0.8), (0.4,0.6), (0.6,0.4), (0.8,0.2), (1.0,0.0), (0.0,0.0) ]
        odorA = 0.0
        odorB = 1.0
        avgnum = 0
        ## OBNetwork.py uses the proxympirank passed here
        proxympirank = avgnum*numodors + 0 + 1 # avgnum*numodors + odornum + 1
        sim =  odorResponse(proxympirank)
        ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
        if ONLY_TWO_MITS and not NO_PGS: includeProjections = ['PG']
        else: includeProjections = []
        tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES,\
            NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
            includeProjections=includeProjections, nolateral=NO_LATERAL)
        BINNED = True#False # for mitrals only
        ## if not BINNED, save the full mitral Vm-s
        ## and not just their spiketimes by setting spiketable = False below.
        network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
            proxympirank, 'morphs', granfilebase, spiketable=BINNED)
        #printNetTree() # from moose_utils.py
        ## monitor those interneurons that are connected to mitral indices 0 and 1
        ## save only spiketimes by setting extras_spikes_only=True
        extras_spikes_only = True # for interneurons
        ## choose one of the below for only interneurons connected to mits0/1 vs all.
        interneurons_args = {}
        #extras_args = {'mitrals':[0,1]}
        tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,\
            interneurons_args, spikes=extras_spikes_only)
        ## To watch the pre compartment of mit2 that inhibits soma of mit 1
        #mit2 = moose.Cell('/mitrals_2')
        #mit2.precomp = moose.Compartment(get_matching_children(mit2, ['Seg0_sec_dendd4_4_278'])[0])
        #mit2._vmTablePrecomp = setupTable("vmTablePrecomp",mit2.precomp,'Vm')
        sim.setupStim(network, ((odorA,odorB),), avgnum=avgnum)
        ## widely different resting potentials of mit0 and mit1
        if VARY_MITS_RMP:
            tweak_field('/mitrals_0/##[TYPE=Compartment]', 'Em', '-58e-3')
            tweak_field('/mitrals_1/##[TYPE=Compartment]', 'Em', '-70e-3')
        if "SAVEFORMOVIE" in sys.argv:
            SAVEFORMOVIE = True
            mitTables = setupMitralTables(network,BINNED)
        else: SAVEFORMOVIE = False
        mitral_responses,mitral_responses_binned = sim.run(network,BINNED)

        ## Save data for movie
        if SAVEFORMOVIE:
            moviedatafile = 'movie_data_netseed'+netseedstr+'_stimseed'+rateseedstr+'_directed'+str(frac_directed)
            timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
            ## number of colours in mitcolours must match number of mitral cells
            ## Each colour-entry below is a tuple of (baseline/initial colour, spiking/peak colour, colourmap)
            ## Each colour is a tuple of (r,g,b,a)
            ## Set in Moogli's config file, whether to change color, and/or alpha, or use colourmap.
            mitcolours = [ 
                ((0.3,0,0,0.3),(1,0,0,1),'jet'), ((0.3,0,0.3,0.3),(1,0,1,1),'jet'),
                ((0,0,0.3,0.3),(0,0,1,1),'jet'), ((0,0,0.3,0.3),(0,0,0.3,0.3),'jet'),
                ((0,0.3,0,0.3),(0,1,0,1),'jet'), ((0,0.3,0,0.3),(0,0.3,0,0.3),'jet') ]                
            dataTables = exportTables(network, \
                NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS, BINNED, \
                interneurons_args,mitcolours) ## args tells: take only interneurons connected to mits 0,1
            mitDataTables = exportMitralTables(mitTables,mitcolours,BINNED)
            ### Save colour info -- not implemented yet
            #f = open(moviedatafile+'_colours.pickle','w')
            #pickle.dump( colourTable, f )
            #f.close()
            ## Save movie data
            if BINNED: moviedatafile += '_mitspiketimes.pickle'
            else: moviedatafile += '_mitVms.pickle'
            f = open(moviedatafile,'w')
            pickle.dump( {'projections':network.projectionDict,\
                'sim_data':(timevec,dataTables,mitDataTables)}, f )
            f.close()
            print "Saved output file",moviedatafile
            sys.exit(0)
            
        ## plot and display
        if not extras_spikes_only:
            timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
            plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS)
        else:
            deltabin = RESPIRATION/NUMBINS
            ## Only the last respiration cycle
            timevec = arange(SETTLETIME+(NUM_RESPS-1)*RESPIRATION+deltabin/2,RUNTIME,deltabin)
            plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS,\
                NO_MULTIS, NUMBINS, RUNTIME, SETTLETIME)
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
        ## plot soma; and precompartment of mit2 that inhibits mit0.
        #figure()
        #title('mitral 2')
        #plot(timevec,mitral_responses[2],color=(1,0,0))
        #plot(timevec,mit2._vmTablePrecomp,color=(0,0,0))
        show()

    #### if multiple processes are called, average over odor morphs
    else:
        ## construct the results filename
        today = datetime.date.today()
        if NO_SINGLES: singles_str = '_NOSINGLES'
        else: singles_str = '_SINGLES'
        if NO_JOINTS: joints_str = '_NOJOINTS'
        else: joints_str = '_JOINTS'
        if NO_PGS: pgs_str = '_NOPGS'
        else: pgs_str = '_PGS'
        if NO_LATERAL: lat_str = '_NOLAT'
        else: lat_str = '_LAT'
        if VARY_MITS_RMP: varmitstr = '_VARMIT'
        else: varmitstr = '_NOVARMIT'
        ## stable enough that time tags are not needed
        now =  ''#datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")+'_'
        outfilename = '../results/odor_morphs/'+now+'odormorph'+\
            '_netseed'+netseedstr+'_stimseed'+rateseedstr
        if NONLINEAR_ORNS: outfilename += '_NL'+NONLINEAR_TYPE
        outfilename += singles_str+joints_str+pgs_str+lat_str+varmitstr+\
            '_numgloms'+str(NUM_GLOMS)
        if DIRECTED: outfilename += '_directed'+str(FRAC_DIRECTED)
        outfilename += '.pickle'
        
        ## if NOSHOW, then check if resultfile exists, proceed only if non-existent.
        if 'NOSHOW' in sys.argv:
            NOSHOW = True
            ## If NOSHOW, then automatic mode, hence don't overwrite resultfile, if exists beforehand.
            if os.path.exists(outfilename):
                ## activdep_inhibition_repeats.py searches for Wrote in first word,
                ## and filename in second word. so output that even if not simulating.
                if mpirank==boss:
                    for procnum in range(1,mpisize):
                        mpicomm.recv(source=procnum,tag=10)
                    print "ExistsSoNotWrote",outfilename
                else:
                    mpicomm.send('done',dest=boss,tag=10)
                sys.exit()
        else: NOSHOW = False

        if mpirank == boss:
            #### collate at boss process
            mitral_responses_list = []
            mitral_responses_binned_list = []
            numavgs = (mpisize-1)/numodors
            for avgnum in range(numavgs):
                response_odorset = []
                response_odorset_binned = []
                for odornum in range(numodors):
                    procnum = avgnum*numodors + odornum + 1
                    print 'waiting for process '+str(procnum)+'.'
                    #### you get a numpy array of rows=NUM_GLOMS*MIT_SISTERS and cols=NUMBINS
                    #### mitral responses has spike times, mitral_responses_binned has binned firing rates
                    mitral_responses,mitral_responses_binned = mpicomm.recv(source=procnum, tag=0)
                    response_odorset.append( mitral_responses )
                    response_odorset_binned.append( mitral_responses_binned )
                mitral_responses_list.append(response_odorset)
                mitral_responses_binned_list.append(response_odorset_binned)
            
            ## write results to a file
            f = open(outfilename,'w')
            pickle.dump((mitral_responses_list,mitral_responses_binned_list), f)
            f.close()
            print "Wrote", outfilename
            
            if not NOSHOW:
                plot_morphs(outfilename)
                show()

        else:
            #### run the slave processes
            sim =  odorResponse()
            avgnum = (mpirank-1)/numodors
            odorset = inputList[(mpirank-1)%numodors]
            odorA, odorB = odorset
            ## If CLUB_MITRAL=False, then extra exc from mitral sisters
            ## (to certain connected granules as proxy for unmodelled sisters) does NOT get used.
            ## Instead, here I connect extra baseline excitation to ALL granules if odor is non-zero.
            if not CLUB_MITRALS and not (odorA==0.0 and odorB==0.0):
                granfilebase += '_extra'
            ## includeProjections gets used only if ONLY_TWO_MITS is True:
            ## Keep below projections to 'second order cells'
            ## i.e. to cells (granules) connected to mits0&1.
            ## The connections between second order cell
            ## and mits0&1 are automatically retained of course.
            ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
            includeProjections = ['PG','granule_baseline']
            tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES,\
                NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
                includeProjections=includeProjections, nolateral=NO_LATERAL)
            ## unique str = 'morphs_', etc so that temp files of morphs and pulses etc do not overlap
            network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
                mpirank, uniquestr, granfilebase, spiketable=True)
            ## widely different resting potentials of mit0 and mit1
            if VARY_MITS_RMP:
                tweak_field('/mitrals_0/##[TYPE=Compartment]', 'Em', '-58e-3')
                tweak_field('/mitrals_1/##[TYPE=Compartment]', 'Em', '-70e-3')
            #printNetTree() # from moose_utils.py

            sim.setupStim(network, (odorset,), avgnum)
            mitral_responses_both = sim.run(network, binned=True)
            mpicomm.send( mitral_responses_both, dest=boss, tag=0 )
            print 'sent from process',mpirank
