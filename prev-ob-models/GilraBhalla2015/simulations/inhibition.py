#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose
from OBNetwork import *

from stimuliConstants import * # has SETTLETIME
from simset_inhibition import * # has REALRUNTIME
from sim_utils import * # has setup_tables(), plot_extras() and build_tweaks()
from moose_utils import *

RUNTIME = REALRUNTIME + SETTLETIME

# Whether to inject the current into tuft or soma, typically soma
TUFT_INJECT = False

from pylab import * # part of matplotlib that depends on numpy but not scipy

## Run: python2.6 inhibition.py <PGvsGranule|PGvsGranuleIPSC|lateral_granule>

#-----------------------------------------------------------

class lateral_granule_inhibition:
    
    def __init__(self):
        pass

    def setup_stim(self,network,args):
        if SPIKEBLOCK:
            print "Blocking Na in mitral cells only"\
                " (don't have granules, else inconsistent) ..."
            for cell in network.mitralTable.values():
                #blockChannels(cell, ['Na','K2']) # only K2 i.e. Kfast, not all K
                blockChannels(cell, ['Na']) # only TTX, only Na blocked not K2
        iAinject = offInject
        iBinject = onInject
        ipulse_duration = 400e-3#150e-3#400e-3 # seconds
        ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
        ## I seem to be using much larger currents - increase the inhibition
        if TUFT_INJECT:
            tuftsegname = 'Seg0_prim_dend_5_20'
            #tuftsegname = 'Seg0_glom_81_102'
            print "injecting into tuft/tuft-base segment",tuftsegname
            self.compA = moose.Compartment(network.mitralTable[mitralidx].path+'/'+tuftsegname)
            self.compB = moose.Compartment(network.mitralTable[mitralsidekickidx].path+'/'+tuftsegname)
            self.mitAtuftbaseTable = setupTable('mittuftbaseTable',self.compA,'Vm')
        else:
            self.compA = network.mitralTable[mitralidx].soma
            self.compB = network.mitralTable[mitralsidekickidx].soma
        iA = setup_iclamp(self.compA, '_mitralA', SETTLETIME, ipulse_duration, iAinject)
        ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
        ## slightly stagger the start of current pulses in the two cells
        ## so that the mitrals do not continuously co-fire.
        iB = setup_iclamp(self.compB,'_mitralB', SETTLETIME-50e-3, ipulse_duration+50e-3, iBinject)
        ## To check backpropagating AP, check at various distances from the soma on the lateral dendrites.
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
        print 'Injecting mitral A with '+str(iAinject)+' and B with '+str(iBinject)

    def run_inhibition(self,network,args):
        resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
        network.context.step(RUNTIME)
        
        timevec = arange(0.0,RUNTIME+1e-10,SIMDT)*1000
        fig = figure(facecolor='w')
        ax = fig.add_subplot(111)
        #ax.set_title('mitral A',size='large')
        #plot(plotBins(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        #plot(plotSpikes(network.mitralTable[str(mitralidx)]._vmTableSoma),'r-,')
        plot(timevec,array(network.mitralTable[mitralidx]._vmTableSoma)*1000,\
            'r',linestyle='-',linewidth=2.0,label='soma')
        if TUFT_INJECT:
            plot(timevec,array(self.mitAtuftbaseTable)*1000,\
                'b',linestyle='--',linewidth=2.0,label='tuft base')
        else:
            plot(timevec,array(network.mitralTable[mitralidx]._vmTableGlom)*1000,\
                'b',linestyle='--',linewidth=2.0,label='tuft')
        plot(timevec,array(self.mitAdendTable)*1000,\
            'g',linestyle='-.',linewidth=2.0,label='lat dend')
        biglegend()
        axes_labels(ax,"time (ms)","Vm (mV)")
        
        figure(facecolor='w')
        title('mitral B')
        #plot(plotBins(network.mitralTable[str(mitralsidekickidx)]._vmTableSoma),'g-,')
        #plot(plotSpikes(network.mitralTable[str(mitralsidekickidx)]._vmTableSoma),'g-,')
        plot(timevec,network.mitralTable[mitralsidekickidx]._vmTableSoma,'r-,',label='Soma')
        plot(timevec,network.mitralTable[mitralsidekickidx]._vmTableGlom,'b-,',label='Glom')
        plot(timevec,self.mitBdendTable,'g-,',label='Dend')
        legend()

#----------------------------------------------------------------

class PGvsGranule:
    
    def __init__(self):
        pass

    def setup_stim(self,network,args):
        mitstr = "mitrals_"+str(mitralidx) # mitralidx from simset_inhibition.py
        #### Connect a short current injector into a granule_single to generate an AP.
        #### This granule_single must connect to our mitral number mitralidx
        found_granmitsyn = False
        # the post_segment on the mitral must be ~= required_distance from the soma.
        self.required_distance = args[0] # meters
        dist_list = []

        for proj in network.projectionDict.values(): # see NetworkML_reader.py for projectionDict format
            # is 'granules_singles' part of the source population name and 'mitral' part of target:
            if 'granules_singles' in proj[0] and 'mitrals' in proj[1]:
                #### The same granule_mitral_GABA SynChan is used for multiple connections
                #### (each SynChan has an array of inputs)
                #### I don't want to change the Gbar for a SynChan multiple times, so maintain a unique list.
                conn_changed = {}
                for conn in proj[2]: # loop through the connections that are part of this projection
                    if mitstr in conn[2]: # if mitral_<mitralidx> is in the post_segment_path
                        post_segment = moose.Compartment(conn[2])
                        ### I correct ALL the granule-|mitral Gbar for different clubbing of PGs vs granules
                        ### I only activate the synapses I am interested in so it doesn't matter what I do to others!       
                        granmitsyn = moose.SynChan(get_matching_children(post_segment, ['granule_mitral_GABA'])[0])
                        if granmitsyn.path in conn_changed: continue
                        conn_changed[granmitsyn.path] = True
                        print 'original', granmitsyn.path, granmitsyn.Gbar
                        mitsoma = network.mitralTable[mitralidx].soma
                        distance = sqrt( (post_segment.x-mitsoma.x)**2 \
                            + (post_segment.y-mitsoma.y)**2 + (post_segment.z-mitsoma.z)**2 )
                        dist_list.append( (abs(distance-self.required_distance), conn, distance) )
                        # linear scaling by factor 20 - suspect!!!
                        granmitsyn.Gbar = granmitsyn.Gbar*PG_CLUB/GRANS_CLUB_SINGLES
                dist_list.sort()
                distdiff,conn,self.dist = dist_list[0]
                print "Lateral compartment is at",self.dist,"meters."
                self.dendsegpath = conn[2]
                post_seg = moose.Compartment(self.dendsegpath) # wrap post_segment_path conn[2]
                #printRecursiveTree(post_seg.id,2)
                ## take the cellname part [-2] of pre_segment_path conn[1],
                ## then take the last [-1] id part of cellname
                granidx = string.split( string.split(conn[1],'/')[-2], '_' )[-1]
                gran = network.populationDict['granules_singles'][1][int(granidx)]
                 ## assumes at least one soma and takes the first!
                gran.soma = moose.Compartment(get_matching_children(gran, ['Soma','soma'])[0])
                gran_inject = 0.25e-9 # Amperes
                ipulse_duration = 5e-3 # seconds
                igran = setup_iclamp(gran.soma, '_gran', SETTLETIME+100e-3, ipulse_duration, gran_inject)
                self.granTable = setupTable('granTable',gran.soma,'Vm')
                self.granDendTable = setupTable('granDendTable',moose.Compartment(conn[1]),'Vm') # pre_segment
                ### measure the Vm in the post_segment of the mitral
                self.mitdendTable = setupTable('mitdendTable',post_seg,'Vm')
                break # only one segment needed.

        if args[1] == 'IPSC': # voltage clamp mitral soma
            #### Since Dong et al use a Ca chelator and TEA (K blocker) in the Vclamp pipette,
            #### I block Ca and K channels in the mitral cell
            blockChannels(network.mitralTable[mitralidx], ['K','Ca']) # in moose_utils.py
            ## PID gain: I think for Davison 4-comp-mitral/granule: 0.5e-5 # optimal gain
            ## too high 0.5e-4 drives it to oscillate at high frequency,
            ## too low 0.5e-6 makes it have an initial overshoot (due to Na channels?)
            ## But for BBmit1993, gain of 1e-6 is optimal
            self.PIDITable = setup_vclamp(mitsoma, '_somavclamp', 0, RUNTIME+150e-3, 0.0e-3, gain=1e-6)
            self.mitITable = setupTable('mitITable',mitsoma,'Im')

        #### Find any one PG that connects to our mitral number mitralidx
        #### Connect a short current injector into this PG to generate an AP.
        #### Find the number of synapses between this PG and this mitral mitralidx
        ####  and normalize by it to compare with granule.
        found_PGmitsyn = False
        mitSegsList = []
        for proj in network.projectionDict.values(): # see NetworkML_reader.py for projectionDict format
            # is 'PG' part of the source population name and 'mitral' part of target:
            if 'PG' in proj[0] and 'mitrals' in proj[1]:
                for conn in proj[2]: # loop through the connections that are part of this projection
                    if mitstr in conn[2]: # if mitral_<mitralidx> is in the post_segment_path
                        ### select the first PG that connects to this mitral mitralidx
                        ###  and connect a current clamp to it. also set up tables for the PG dend and mitral tuft.
                        if not found_PGmitsyn:
                            PG_path = conn[1]
                            # take the cellname part [-2] of PG_path, then take the last [-1] id part of cellname
                            PGidx = string.split( string.split(PG_path,'/')[-2], '_' )[-1]
                            PG = network.populationDict['PGs'][1][int(PGidx)]
                            PG.soma = moose.Compartment(get_matching_children(PG, ['Soma','soma'])[0]) # assumes at least one soma and takes the first!
                            self.PGTable = setupTable('PGTable',PG.soma,'Vm')
                            self.PGDendTable = setupTable('PGDendTable',moose.Compartment(conn[1]),'Vm') # pre_segment
                            PG_inject = 0.25e-9 # Amperes
                            ipulse_duration = 5e-3 # seconds
                            iPG = setup_iclamp(PG.soma, '_PG', SETTLETIME+RUNTIME/2.0, ipulse_duration, PG_inject)
                            self.tuftsegpath = conn[2]
                            post_seg = moose.Compartment(self.tuftsegpath)
                            ### measure the Vm in the post_segment of the mitral by wrapping post_segment_path conn[2]
                            self.mittuftTable = setupTable('mittuftTable',post_seg,'Vm')
                            found_PGmitsyn = True
                            numsyns_thisPG_to_thismitral = 1
                            mitSegsList.append(self.tuftsegpath)
                        ### find out how many times this selected PG connects to this mitral mitralidx
                        else:
                            comparedPG_path = conn[1]
                            # take the cellname part [-2] of PG_path, then take the last [-1] id part of cellname
                            comparedPGidx = string.split( string.split(comparedPG_path,'/')[-2], '_' )[-1]
                            if PGidx == comparedPGidx:
                                numsyns_thisPG_to_thismitral += 1
                                mitSegsList.append(conn[2])
            if found_PGmitsyn: break

        ### I correct ALL the PG-|mitral Gbar for multiple PG-|mit synapses of selected PG to desired mitral
        ### I only activate the synapses I am interested in so it doesn't matter what I do to others!       
        ### make a set (keep unique elements only) of the post-segments, so as not to correct any segment twice.
        for segpath in set(mitSegsList):
            post_seg = moose.Compartment(segpath)
            PGmitsyn = moose.SynChan(get_matching_children(post_seg, ['PG_mitral'])[0])
            PGmitsyn.Gbar = PGmitsyn.Gbar/numsyns_thisPG_to_thismitral
        print "Corrected PG Gbar for selected PG",PGidx,"making",\
            numsyns_thisPG_to_thismitral,"synapses to mitral",mitralidx
        print "Also corrected granule Gbar for",GRANS_CLUB_SINGLES,\
            "single grans clubbed and",PG_CLUB,"PGs clubbed."
        print "Finally PG-|mit and gran-|mit synapses are each PG_CLUB =",\
            PG_CLUB,"times the single synapse."

        print 'Injecting granule_single',granidx,'with',gran_inject,\
            'at',SETTLETIME,'s and PG',PGidx,'with',PG_inject,'at',SETTLETIME+RUNTIME/2.0,'s.'        

    def run_inhibition(self,network,args):
        resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
        network.context.step(RUNTIME)
        timevec = arange(0.0,RUNTIME+1e-10,SIMDT)
        figure(facecolor='w')
        plot(timevec,network.mitralTable[mitralidx]._vmTableSoma,'r-,',label='mitral soma Vm')
        plot(timevec,self.mitdendTable,',-g',label=self.dendsegpath+' Vm')
        plot(timevec,self.mittuftTable,',-b',label=self.tuftsegpath+' Vm')
        legend(loc='lower right')
        ylabel('SecDend is at %0.2f microns'%(self.dist*1e6,),fontsize='large')
        xlabel('time (s)', fontsize='large')
        figure()
        plot(timevec,self.granTable,',-g',label='granule soma Vm')
        plot(timevec,self.granDendTable,',-y',label='granule dend Vm')
        plot(timevec,self.PGTable,',-b',label='PG soma Vm')
        plot(timevec,self.PGDendTable,',-c',label='PG dend Vm')
        legend(loc='upper right')
        if args[1] == 'IPSC':
            figure(facecolor='w')
            #plot(timevec,self.mitITable,',-r',label='mitral soma Im')
            plot(timevec,self.PIDITable,',-g',label='PID I to mitral')
            legend(loc='upper right')
            ylabel('I (A); SecDend is at %0.2f microns'%(self.dist*1e6,),fontsize='large')
            xlabel('time (s)', fontsize='large')

#----------------------------------------------------

if __name__ == "__main__":
    err_string = "You must tell me the type of simulation\
 to run as a command line parameter:\n\
 'lateral_granule' or 'PGvsGranule' or 'PGvsGranuleIPSC'"
    if len(sys.argv) < 2:
        print err_string
        sys.exit(1)
    sim_option = sys.argv[1]
    if sim_option == 'lateral_granule':
        sim =  lateral_granule_inhibition()
        includeProjections = ['granule_baseline']
        args = []
    elif sim_option == 'PGvsGranule':
        sim =  PGvsGranule()
        includeProjections = []
        # the post_segment on the mitral must be ~ this distance from the soma
        args = [100e-6,'IPSP'] # meters
    elif sim_option == 'PGvsGranuleIPSC':
        sim =  PGvsGranule()
        includeProjections = []
        # the post_segment on the mitral must be ~ this distance from the soma
        args = [100e-6,'IPSC'] # meters
    else:
        print "Sorry wrong option.", err_string
        sys.exit(1)

    tweaks = build_tweaks( CLUB_MITRALS, NO_SPINE_INH, NO_SINGLES, NO_JOINTS,\
        NO_MULTIS, NO_PGS, ONLY_TWO_MITS, includeProjections, (mitralidx,mitralsidekickidx) )
    network = OBNetwork(OBNet_file, synchan_activation_correction,\
        tweaks, mpirank, granfilebase, spiketable=False)
    #printNetTree() # from moose_utils.py

    sim.setup_stim(network, args)
    spikes = True
    tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, spikes=spikes)
    sim.run_inhibition(network, args) # tests simset_inhibition i.e. inhibition between two mitrals
    if spikes:
        bins = 20
        timevec = arange(0.0,REALRUNTIME,REALRUNTIME/float(bins))
        plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, bins, RUNTIME, SETTLETIME)
    else:
        timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
        plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, '')
    show()
