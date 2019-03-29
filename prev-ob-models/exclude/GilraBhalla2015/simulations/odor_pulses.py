#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import datetime
import pickle

## from node000:
## mpiexec -machinefile ~/hostfile -n <numavgs*numpulses+1> ~/Python-2.6.4/bin/python2.6 odor_pulses.py
## typical value for numavgs = 9
## (depends on number of available processing nodes and number of odorfiles generated)
## typical value for numpulses = 2*RANDOM_NUM_SEQ = 6.
## nohup mpiexec -machinefile ~/hostfile -n 55 ~/Python-2.6.4/bin/python2.6 odor_pulses.py < /dev/null &
## OR for a single odor run; from any node:
## python2.6 odor_pulses.py
## Set various option like NO_PGs or ONLY_TWO_MITS in simset_odor

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from data_utils import * # has mpi import and variables also
from OBNetwork import *
from sim_utils import *

from stimuliConstants import * # has SETTLETIME, PULSE_RUNTIME, pulsebins
from simset_odor import *

from pylab import * # part of matplotlib that depends on numpy but not scipy
from plot_odor_pulses import *

#-----------------------------------------------------------

class odorResponse:
    
    def __init__(self):
        self.mpirank = mpirank
        self.context = moose.PyMooseBase.getContext()

    def setupStim(self,network,pulse_i,avgnum):
        self.setupOdor(network, pulse_i, avgnum)
        print "Setup pulse number =",pulse_i,"at",self.mpirank

    def setupOdor(self, network, pulse_i, avgnum):
        ### first figure out which PG belongs to which glom
        ### PG_glom_map[pgname] returns the glom num of the PG: needed for ORN to PG connections.
        PG_glom_map = {}
        for projname in network.projectionDict.keys():
            if 'PG_mitral' in projname:
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    #get the glomnum from the post path proj[2]
                    # name of the mitral cell from '/mitrals_2/...'
                    mitname = string.split(proj[2],'/')[1]
                    # glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    glomnum = int(string.split(mitname,'_')[1]) / 2
                    # name of the PG cell from '/PGs_2/...'
                    pgname = string.split(proj[1],'/')[1]
                    PG_glom_map[pgname] = glomnum
        ### Now connect the ORNs
        for projname in network.projectionDict.keys():
            #### Calling attach_spikes() for each projection, would reconnect files to the same segment multiple times.
            #### But attach_files_uniquely() checks whether timetable.tableSize is zero or not
            #### i.e. files already attached or not.
            ### connect ORNs to mitrals
            if 'ORN_mitral' in projname:
                print "connecting ORN files to mitrals"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    #get the glomnum from the post path proj[2]
                    mitname = string.split(proj[2],'/')[1] # name of the mitral cell from '/mitrals_2/...'
                    glomnum = int(string.split(mitname,'_')[1]) / 2 # glomerulus number from 'mitrals_2' by integer division i.e. 2/2 = 1
                    filebase = ORNpathseedstr+'firetimes_rndpulse_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],pulse_i,avgnum)
            ### connect ORNs to PG
            if 'ORN_PG' in projname:
                print "connecting ORN files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathseedstr+'firetimes_rndpulse_glom_'+str(glomnum)
                    self.attach_files_uniquely(filebase,proj[0],proj[2],pulse_i,avgnum)
            ### connect SAs to PG
            if 'SA_PG' in projname:
                print "connecting SA files to PGs"
                for i,proj in enumerate(network.projectionDict[projname][2]):
                    pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                    glomnum = PG_glom_map[pgname]
                    filebase = ORNpathseedstr+'firetimes_SA_rndpulse'
                    self.attach_files_uniquely(filebase,proj[0],proj[2],pulse_i)
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
            #        filebase = '../firefiles/firetimes_rndpulse_glom_'+glomstr
            #        self.attach_files_uniquely(filebase,proj[0]+'_'+glomstr,proj[2],pulse_i,avgnum)

    def attach_files_uniquely(self,filebase,synname,postsegpath,pulse_i,avgnum=None):
        ttpath = postsegpath+'/'+synname+'_tt'
        if self.context.exists(ttpath):
            # timetable already created by networkML reader - just wrap it below.
            tt = moose.TimeTable(ttpath) # post_segment_path+'/'+syn_name+'_tt'
        else:
            # if timetable was not already created by networkML reader,
            # it means that the synaptic weights must be zero! (no extra inhibition - only main inhibition)
            # hence do not attach spikefiles
            return
        if tt.tableSize != 0: return # if files are already attached, do nothing!
        filebase += '_pulse_'+str(pulse_i)
        #### I have separate firefiles for each odor, each glom and each avgnum!
        #### avgnum is what decides variability.
        if avgnum is not None: filebase += '_avgnum'+str(avgnum)
        ## attach_spikes() accesses filenumbers to this segment
        ## from 'fileNumbers' field (of the timetable object in MOOSE)
        ## which is created while reading in networkML.
        attach_spikes(filebase, tt, uniquestr+str(self.mpirank))

    def run(self,network, binned):
        print "Resetting MOOSE."
        ## from moose_utils.py sets clocks and resets
        resetSim(network.context, SIMDT, PLOTDT)
        print "Running at",self.mpirank
        network.context.step(PULSE_RUNTIME)
        mitral_responses = []
        mitral_responses_binned = []
        if ONLY_TWO_MITS or NO_LATERAL: num_mits = MIT_SISTERS
        else: num_mits = NUM_GLOMS*MIT_SISTERS
        for mitnum in range(num_mits):
            mitral = network.mitralTable[mitnum]
            if binned: mitral_responses_binned.append(
                plotBins(mitral._vmTableSoma, pulsebins, PULSE_RUNTIME, 0.0) )
            ## need to convert to numpy's array(),
            ## else MOOSE table cannot be pickled for mpi4py send()
            mitral_responses.append(array(mitral._vmTableSoma))
        return (mitral_responses,mitral_responses_binned)

#----------------------------------------------------------------

if __name__ == "__main__":

    if mpisize == 1:
        ## if only one process is called, plot one odor directly
        sim =  odorResponse()
        ## 'PG' includes 'ORN_PG', 'PG_mitral', 'mitral_PG' and 'SA_PG'
        if ONLY_TWO_MITS and not NO_PGS: includeProjections = ['PG']
        else: includeProjections = []
        tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES,\
            NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
            includeProjections=includeProjections, nolateral=NO_LATERAL)
        BINNED = True#False
        ## if not BINNED, save the full mitral Vm-s and
        ## not just their spiketimes by setting spiketable = False below.
        uniquestr = 'pulses_'
        network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,\
            mpirank, uniquestr, granfilebase+'_noresp_extra', spiketable=BINNED)
        #printNetTree() # from moose_utils.py
        ## monitor those interneurons that are connected to mitral indices 0 and 1
        ## save only spiketimes by setting spikes=True
        extras_spikes_only = False#True
        tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,\
            {'mitrals':[0,1]}, spikes=extras_spikes_only)
        ## pulse 0 is random air pulse, pulse 1 is odor A,
        ## pulse 2 is odor B, so on alternately, last pulse is A&B.
        pulse_i = 1
        sim.setupStim(network, pulse_i, avgnum=0)
        ## widely different resting potentials of mit0 and mit1
        if VARY_MITS_RMP:
            tweak_field('/mitrals_0/##[TYPE=Compartment]', 'Em', '-58e-3')
            tweak_field('/mitrals_1/##[TYPE=Compartment]', 'Em', '-70e-3')
        mit0soma_Im = setupTable('mit0somaIm',network.mitralTable[0].soma,'Im')
        ## setup a separate Vm table, since the one setup by OBNetwork() stores only spikes.
        mit0soma_Vm = setupTable('mit0somaVm',network.mitralTable[0].soma,'Vm')
        mitral_responses,mitral_responses_binned = sim.run(network,BINNED)
        if not extras_spikes_only:
            timevec = arange(0.0,PULSE_RUNTIME+PLOTDT+1e-12,PLOTDT)
            plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS)
        else:
            pulsebindt = PULSE_RUNTIME/pulsebins
            timevec = arange(pulsebindt/2.0,PULSE_RUNTIME,pulsebindt)
            plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES,
                NO_JOINTS, NO_MULTIS, pulsebins, PULSE_RUNTIME, 0.0)
        ## Plot the Vm, total current, membrane current and axial current in mit0's soma
        ## Total current in mitral soma = C dV/dt
        mit0soma_C = network.mitralTable[0].soma.Cm
        mit0soma_Vm = array(mit0soma_Vm)
        mit0soma_Itotal = mit0soma_C * (mit0soma_Vm[1:]-mit0soma_Vm[:-1])/PLOTDT
        mit0soma_Im = array(mit0soma_Im)
        mit0soma_Iaxial = mit0soma_Itotal - mit0soma_Im[1:]
        timevec_unbinned = arange(0.0,PULSE_RUNTIME+1e-12,PLOTDT)
        ## plot mitral 0 soma Vm on left-side y axis
        figure(facecolor='w')
        title('Vm in mitral 0 soma')
        plot(timevec_unbinned,mit0soma_Vm[1:],label='Vm',color='k')
        ## plot mitral 0 soma I-s on right-side y axis
        figure(facecolor='w')
        title('Currents in mitral 0 soma')
        plot(timevec_unbinned,mit0soma_Itotal,label='Itotal',color='r')
        plot(timevec_unbinned,mit0soma_Im[1:],label='Im',color='g')
        plot(timevec_unbinned,mit0soma_Iaxial,label='Iaxial',color='b')
        legend()
        ## plot the responses of a glomerulus
        figure(facecolor='w')
        title('Glomerulus 0')
        if BINNED:
            mitral_responses = mitral_responses_binned
        #plot(timevec,mitral_responses[0],color=(0.0,1.0,0.0))
        #plot(timevec,mitral_responses[1],color=(0.0,1.0,0.5))
        show()
    else:
        ## if multiple processes are called, average over odor morphs
        ## pulse 0 is random air pulse, pulses 1 & 2, 3 & 4
        ## are odors A & B respectively, last pulse 5 is A&B.
        numpulses = RANDOM_PULSE_NUMS*2

        ## construct the results filename
        now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
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
        outfilename = '../results/odor_pulses/'+now+'odorpulses'+\
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
            ## mitral_responses_list[avgnum][pulsenum][mitnum]
            mitral_responses_list = []
            mitral_responses_binned_list = []
            numavgs = (mpisize-1)/numpulses
            for avgnum in range(numavgs):
                response_odorset = []
                response_odorset_binned = []
                for pulsenum in range(numpulses):
                    procnum = avgnum*numpulses + pulsenum + 1
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
                plot_pulse_responses(outfilename)
                show()

        else:
            #### run the slave processes
            ## uniquestr to put in every temp filename to avoid clashing with other processes
            if len(sys.argv)>2: uniquestr = sys.argv[2]+'_' # _ necessary, else say 'pulses2'+mpirank is screwed
            else: uniquestr = 'pulses_'
            sim =  odorResponse()

            ## no respiration tuning
            granfilebase += '_noresp'
            ## If CLUB_MITRAL=False, then extra exc from mitral sisters
            ## (to certain connected granules as proxy for unmodelled sisters) does NOT get used.
            ## Instead, here I connect extra baseline excitation to ALL granules if odor is non-zero.
            ## Ideally, this extra exc should be only during odor period, but
            ## one doesn't know if suction is comparable to respiration, and also
            ## I'm only fitting the odor period for kernels,
            ## so extra exc was generated for the full duration in generate_firefiles_gran_baseline.py.
            if not CLUB_MITRALS: granfilebase += '_extra'

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
            ## unique str = 'pulses_' so that temp files of morphs and pulses etc do not overlap
            network = OBNetwork(OBNet_file, synchan_activation_correction, tweaks,
                mpirank, uniquestr, granfilebase, spiketable=True)
            ## widely different resting potentials of mit0 and mit1
            if VARY_MITS_RMP:
                tweak_field('/mitrals_0/##[TYPE=Compartment]', 'Em', '-58e-3')
                tweak_field('/mitrals_1/##[TYPE=Compartment]', 'Em', '-70e-3')
            #printNetTree() # from moose_utils.py

            avgnum = (mpirank-1)/numpulses
            pulsenum = (mpirank-1)%numpulses
            sim.setupStim(network, pulsenum, avgnum)
            mitral_responses_both = sim.run(network, binned=True)
            mpicomm.send( mitral_responses_both, dest=boss, tag=0 )
            print 'sent from process',mpirank
