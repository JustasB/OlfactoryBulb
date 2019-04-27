#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle

######################### ONLY odorA is currently simulated, see sim.setupStim() at the very end.
######################### Firefiles are available for both odors, so can set odorB in sim.setupStim().
## from node000:
## mpiexec -machinefile ~/hostfile -n <numavgs*numscalings+1> ~/Python-2.6.4/bin/python2.6 odor_morphs.py
## nohup mpiexec -machinefile ~/hostfile -n 55 ~/Python-2.6.4/bin/python2.6 odor_scaledpulses.py < /dev/null &
## typical value for numavgs = 9
## (depends on number of available processing nodes and number of odorfiles generated)
## typical value for numscalings = 6. number of items in scaledList (air + 5 scalings).
## OR for a single odor run; from any node:
## python2.6 odor_scaledpulses.py
## Set various option like NO_PGs, etc in simset_odor_minimal.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, SCALED_RUNTIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has REALRUNTIME, NUMBINS
## NUMBINS=10 for respiration, wrong numbins to use her,
##  rebinned later in analysis, so not an issue while saving here though.

RUNTIME = SCALED_RUNTIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self,network,args,avgnum):
        scalestr = args[0]
        self.setupOdor(network, scalestr, avgnum)        
        print "Setup odor scaling",scalestr,"at",self.mpirank

    def setupOdor(self, network, scalestr, avgnum):
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
                    filebase = ORNpathseedstr+'firetimes_scaledpulses_width'+str(scaledWidth)+'_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],scalestr,avgnum)
            ### connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathseedstr+'firetimes_scaledpulses_width'+str(scaledWidth)+'_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],scalestr,avgnum)
            ### connect SAs to PG
            if 'SA_PG' in projname:
                print "SA not implemented for scaled pulses."
                #print "connecting SA files to PGs"
                #for i,proj in enumerate(network.projectionDict[projname][2]):
                #    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                #    glomnum = PG_glom_map[pgname]
                #    filebase = ORNpathseedstr+'firetimes_SA'
                #    self.attach_files_uniquely(filebase,proj[0],proj[2],odorA,odorB)
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

    def attach_files_uniquely(self,filebase,synname,postsegpath,scalestr,avgnum=None):
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
        filebase += '_odor'+scalestr
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
            ## NUMBINS=10 for respiration, wrong numbins to use her,
            ##  rebinned later in analysis, so not an issue while saving here though.
            ## BAD! this sim is not for respiration; but rebinned later, hence saved.
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
    else: uniquestr = 'scaledpulses_'

    #### if only one process is called, plot one odor directly
    if mpisize == 1:
        sim =  odorResponse()
        ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
        if ONLY_TWO_MITS and not NO_PGS: includeProjections = ['PG']
        else: includeProjections = []
        tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES,\
            NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
            includeProjections=includeProjections, nolateral=NO_LATERAL)
        BINNED = False # for mitrals
        ## if not BINNED, save the full mitral Vm-s
        ## and not just their spiketimes by setting spiketable = False below.
        network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
            mpirank, 'scaledpulses', granfilebase, spiketable=BINNED)
        #printNetTree() # from moose_utils.py
        ## monitor those interneurons that are connected to mitral indices 0 and 1
        ## save only spiketimes by setting extras_spikes_only=True
        extras_spikes_only = False#True # for interneurons
        tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,\
            {'mitrals':[0,1]}, spikes=extras_spikes_only)
        ### To watch the pre compartment of mit2 that inhibits soma of mit 1
        #mit2 = moose.Cell('/mitrals_2')
        #mit2.precomp = moose.Compartment(get_matching_children(mit2, ['Seg0_sec_dendd4_4_278'])[0])
        #mit2._vmTablePrecomp = setupTable("vmTablePrecomp",mit2.precomp,'Vm')
        ## To watch the inactivation yGate of Na channel of mit0
        mit0 = moose.Cell('/mitrals_0')
        #printCellTree(mit0)
        mit0.soma_Na = moose.HHChannel('/mitrals_0/Seg0_soma_0/Na_mit_usb')
        mit0._ygatetable = setupTable("ygatetable",mit0.soma_Na,'Y')
        chanIs = []
        channames = ['Na_mit_usb','K_mit_usb','K2_mit_usb',\
            'KA_bsg_yka','Kca_mit_usb','LCa3_mit_usb']
        for channame in channames:
            chan = moose.HHChannel('/mitrals_0/Seg0_soma_0/'+channame)
            chanIs.append(setupTable(channame+'I',chan,'Ik'))
        mit1 = moose.Cell('/mitrals_1')
        mit1.soma_Na = moose.HHChannel('/mitrals_1/Seg0_soma_0/Na_mit_usb')
        mit1._ygatetable = setupTable("ygatetable",mit1.soma_Na,'Y')
        
        sim.setupStim(network, ('A_scale4',), avgnum=0)
        ## widely different resting potentials of mit0 and mit1
        if VARY_MITS_RMP:
            tweak_field('/mitrals_0/##[TYPE=Compartment]', 'Em', '-58e-3')
            tweak_field('/mitrals_1/##[TYPE=Compartment]', 'Em', '-70e-3')
        mitral_responses,mitral_responses_binned = sim.run(network,BINNED)
        if not extras_spikes_only:
            timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
            plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS)
        else:
            deltabin = RUNTIME/50e-3
            timevec = arange(SETTLETIME+deltabin/2,RUNTIME,deltabin)
            numberofbins = len(timevec)
            plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS,\
                NO_MULTIS, numberofbins, RUNTIME, SETTLETIME)
        figure()
        title('Glomerulus 0')
        if BINNED:
            deltabin = RUNTIME/50e-3
            timevec = arange(SETTLETIME+deltabin/2,RUNTIME,deltabin)
            mitral_responses = mitral_responses_binned
        else:
            timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
        plot(timevec,mitral_responses[0],color=(0.0,1.0,0.0))
        plot(timevec,mitral_responses[1],color=(0.0,1.0,0.5))
        ### plot soma; and precompartment of mit2 that inhibits mit0.
        #figure()
        #title('mitral 2')
        #plot(timevec,mitral_responses[2],color=(1,0,0))
        #plot(timevec,mit2._vmTablePrecomp,color=(0,0,0))
        ## plot yGate of Na channel in soma of mit0.
        figure()
        title('mitrals 0 & 1 soma Na inactivation gate')
        plot(timevec,mit0._ygatetable,color=(1,0,0))
        plot(timevec,mit1._ygatetable,color=(0,0,1))
        figure()
        for i,chanI in enumerate(chanIs):
            plot(timevec,chanI,label=channames[i])
        title("mit 0 channel Is")
        legend()
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
        outfilename = '../results/odor_pulses/'+now+'scaledpulses_width'+str(scaledWidth)+\
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

        numodors = len(scaledList)
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
                figure()
                show()

        else:
            #### run the slave processes       
            sim =  odorResponse()
            avgnum = (mpirank-1)/numodors
            scalenum = (mpirank-1)%numodors
            ## If CLUB_MITRAL=False (in simset_odor.py), then extra exc from mitral sisters
            ## (to certain connected granules as proxy for unmodelled sisters) does NOT get used.
            ## Instead, here I connect extra baseline excitation to ALL granules
            ## Don't set this True ever, as the baseline should scale with odor which it does not
            if not CLUB_MITRALS:
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

            sim.setupStim(network, ('A_scale'+str(scalenum),), avgnum)
            mitral_responses_both = sim.run(network, binned=True)
            mpicomm.send( mitral_responses_both, dest=boss, tag=0 )
            print 'sent from process',mpirank
