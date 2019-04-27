#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import pickle
import datetime

sys.path.extend(["..","../networks","../generators","../simulations"])

from OBNetwork import *

from stimuliConstants import * # has SETTLETIME
from simset_activinhibition import * # has REALRUNTIME
from sim_utils import * # has build_tweaks(), and print_extras_activity()
from data_utils import *

########## You need to run:
## From any node on the gj cluster:
## python2.6 activdep_inhibition_mitspikeinit.py [SAVEFIG]

SPIKERUNTIME = 50e-3
RUNTIME = SPIKERUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

## 33 for first soma, then tuft; 35 onwards for first tuft, then soma
## Weak shock -- 33 mit, 4 PG; Transition shock -- 35 mit, 4 PG; Strong shock -- 66 mit, 8 PG.
#num_tuft_comps_shock_mit = 33
##num_tuft_comps_shock_mit = 40
num_tuft_comps_shock_mit = 66
#num_tuft_comps_shock_pg = 4 # keep this a multiple of 2, distributed on 2 PG dendrites
num_tuft_comps_shock_pg = 8 # keep this a multiple of 2, distributed on 2 PG dendrites
PLOT_EXTRAS = True

## set lateral_mitnum to 1 for 2MITS / 2 for 2GLOMS option set in generate_neuroML.py
## Note that for directed connectivity, mit 3 is used for directed conn to mit 0 in generate_neuroml.py,
## thus mit 2 represents a non-directed conn cell.
## If you want to show asymm inhibition between directed cells, you should use mit 3 below.
lateral_mitnum = 2#3
if REVERSED_ADI:
    mitralmainidx = lateral_mitnum
    mitralsidekickidx = 0
else:
    mitralmainidx = 0
    mitralsidekickidx = lateral_mitnum

from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size() # Total number of processes
mpirank = mpicomm.Get_rank() # Number of my process
mpiname = MPI.Get_processor_name() # Name of my node
# The 0th process is the boss who collates/receives all data from workers
boss = 0
print 'Process '+str(mpirank)+' on '+mpiname+'.'

if PLOT_EXTRAS and mpisize>11:
    print "You want to plot Vm-s of mitrals, singles and granules for",mpisize,"processes."
    print "To avoid you the embarrassment of lots of figures, I'm aborting."
    sys.exit(1)


def onespike2synapse(network,compartmentobj,synName,synNum):
    ## connectSynapse is in moose.utils in MOOSE beta 1.4 -- deepcopies a new synapse
    ## so, if synapse exists, don't try to deepcopy over it again.
    if moose.context.exists(compartmentobj.path+'/'+synName):
        synapse = moose.SynChan(compartmentobj.path+'/'+synName)
    else:
        synapse = connectSynapse(network.context, compartmentobj, synName, 1.0)
    ## put in i also, because randomly a synapse might get connected twice.
    spiketable = moose.TimeTable(synapse.path+'/tt_'+str(synNum))
    #### SynChan's synapse MsgDest takes time as its argument.
    ## Thus spiketable should contain a list of spike times.
    spiketable.connect("event", synapse,"synapse")
    ## Presently only method 4 i.e. loading from file is supported for TimeTable.
    ## Hence the need to write a single value in a file.
    fn = 'temp_spikefile_mitspikeinit.txt'
    f = open(fn,'w')
    f.write(str(SETTLETIME+uniform(0,4e-3))) ## 4ms jitter in nerve shock
    f.close()
    spiketable.filename = fn
    os.remove(fn)


def setup_stim(network):
    ## network.populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
    ## network.projectionDict = { 'projectionname1':(source,target,[(syn_name1,pre_seg_path,post_seg_path),...]) , ... }

    ## We want a single sharp current injection at SETTLETIME, so cannot use soma.inject
    
    ## Only connect to num_tuft_comps_shock_mit number of tuft compartments
    mitcell = network.mitralTable[mitralmainidx] ## mitralTable is same as populationDict['mitrals'][1]
    for segi,seginfo in enumerate((network.cellSegmentDict['mitral'].values())[0:num_tuft_comps_shock_mit]):
        ## segment info - see MorphML_reader.py
        ## cellSegmentDict = { segid1 : [ segname,(proximalx,proximaly,proximalz),
        ##    (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] , ... }
        ## seginfo[5] is a list of potential synapses at this segment
        if 'ORN_mitral' in seginfo[5]:
            ## wrap this segment, seginfo[0] is segment name
            tuftcomp = moose.Compartment(mitcell.path+'/'+seginfo[0])
            onespike2synapse(network,tuftcomp,'ORN_mitral',segi)
            ### compartment, name_extn, start_time, duration, current (all SI)
            #setup_iclamp(tuftcomp, '_nerveshock'+str(segi),\
            #    SETTLETIME, 1e-3, tuft_shock_inject)
    
    ## Only connect to num_tuft_comps_shock_pg number of PG compartments
    for pgcell in network.populationDict['PGs'][1].values():
        for segi,seginfo in enumerate((network.cellSegmentDict['PG'].values())): # 3 compartments
            ## cellSegmentDict = { segid1 : [??, [ segname,(proximalx,proximaly,proximalz),
            ##    (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] ] , ... }
            ## seginfo[5] is a list of potential synapses at this segment
            if 'ORN_PG' in seginfo[5]: # 2 compartments
                for i in range(num_tuft_comps_shock_pg/2):
                    ## wrap this segment, seginfo[0] is segment name
                    tuftcomp = moose.Compartment(pgcell.path+'/'+seginfo[0])
                    onespike2synapse(network,tuftcomp,'ORN_PG',i*2+segi)
                    ### compartment, name_extn, start_time, duration, current (all SI)
                    #setup_iclamp(tuftcomp, '_nerveshock'+str(i*2+segi),\
                    #    SETTLETIME, 1e-3, tuft_shock_inject)
                    
    print "Set up simultaneous synapses on mainidx mitral and each PG:",\
        num_tuft_comps_shock_mit,num_tuft_comps_shock_pg


def run_inhibition(network, tables):
    resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
    network.context.step(RUNTIME)
    ## get mitral A's firing rate


#----------------------------------------------------

if __name__ == "__main__":
    uniquestr = 'spikeinit_'
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
        includeProjections, (mitralmainidx,mitralsidekickidx) )
    ## send mpirank to put in ORN filenames / gran baseline temp files
    ## so they do not clash between mpi processes
    ## also, unique str, so that temp files of morphs, pulses, etc do not overlap
    network = OBNetwork(OBNet_file, synchan_activation_correction,\
        tweaks, mpirank, uniquestr, granfilebase, spiketable=SPIKETABLE)
    #printNetTree() # from moose_utils.py

    ## setup a separate Vm table, since the one setup by OBNetwork() above stores only spikes.
    mitAsoma_Vm = setupTable('mitAsomaVm',network.mitralTable[mitralmainidx].soma,'Vm')

    ## monitor Vm at the base of the tuft
    tuftbase_seg = moose.Compartment(network.mitralTable[mitralmainidx].path+\
        #'/Seg0_tuftden_19_23') # for migliore and shepherd 2007 cell
        '/Seg0_prim_dend_5_20') # tuft base
    tuftBaseTable_Vm = setupTable('mitdendTable',tuftbase_seg,'Vm')

    ## monitor Vm in the tuft compartment
    tuftcomp = moose.Compartment(network.mitralTable[mitralmainidx].path+'/Seg0_glom_81_102')
    tuftTable_Vm = setupTable('mitTuftTable',tuftcomp,'Vm')

    setup_stim(network)
    SPIKETABLE = False
    ## if not SPIKETABLE: record the Vm-s of a few interneurons
    ## else: record spiketimes of all interneurons
    tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,
        args={'mitrals':(mitralmainidx,)}, spikes=SPIKETABLE)
    run_inhibition(network, tables)

    timevec_unbinned = arange(0.0,RUNTIME+1e-12,PLOTDT)
    ## Paper figure: mitA's soma, tuft base and tuft
    fig2 = figure(figsize=(columnwidth/2.0,linfig_height/2.0),dpi=1200,facecolor='w') # 'none' is transparent
    ax = fig2.add_subplot(111)
    plot(timevec_unbinned*1000,array(mitAsoma_Vm)*1000,'-r', label='soma',\
        linewidth=plot_linewidth, linestyle='-') # ms and mV
    plot(timevec_unbinned*1000,array(tuftBaseTable_Vm)*1000,'-b', label='base',\
        linewidth=plot_linewidth, dashes=(0.5,0.5)) # ms and mV
    plot(timevec_unbinned*1000,array(tuftTable_Vm)*1000,'-k', label='tuft',\
        linewidth=plot_linewidth, dashes=(1.5,0.5)) # ms and mV
    beautify_plot(ax,x0min=False,drawxaxis=True,drawyaxis=True,\
        xticks=[SETTLETIME*1000,SETTLETIME*1000+20],yticks=[-75,0,40])
    #add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
    #    sizex=50,labelx='50 ms',sizey=20,labely='20 mV',label_fontsize=5,\
    #    bbox_to_anchor=[0.7,0.6],bbox_transform=ax.transAxes)
    ax.set_xlim(SETTLETIME*1000,SETTLETIME*1000+20)
    ax.set_xticklabels(['0','20'])
    ax.set_ylim(-75,40)
    axes_labels(ax,"time (ms)","Vm (mV)",fontsize=label_fontsize,xpad=-3,ypad=-4)
    if 'SAVEFIG' in sys.argv:
        fig2.tight_layout()
        fig2.savefig('../figures/connectivity/cells/mitral_spikeinit_'+\
            str(num_tuft_comps_shock_mit)+'ORNmits.png',dpi=fig2.dpi)
        fig2.savefig('../figures/connectivity/cells/mitral_spikeinit_'+\
            str(num_tuft_comps_shock_mit)+'ORNmits.svg',dpi=fig2.dpi)
        ## Paper figure done

    if PLOT_EXTRAS:
        timevec = arange(0.0,RUNTIME+2*PLOTDT+1e-12,PLOTDT)
        titlestr = ''
        #figure()
        #title('red:mitA, green:mitB, '+titlestr)
        #plot(timevec, mitA, 'r,')
        #plot(timevec, mitB, 'g,')
        plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS, titlestr)
    else:
        spikestables = \
            print_extras_activity(tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS, '')
            
    show()
