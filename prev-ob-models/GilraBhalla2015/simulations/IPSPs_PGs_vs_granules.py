#!/usr/bin/env python
# -*- coding: utf-8 -*-

########## You need to run:
## python2.6 IPSPs_PGs_vs_granules.py

import os
import sys
import math
import pickle
import datetime

sys.path.extend(["..","../networks","../generators","../simulations"])

from OBNetwork import *

from stimuliConstants import *
## Set OBNet_file in simset_activinhibition to have 10:1 singles, so that IPSCs are same height
from simset_activinhibition import *
from sim_utils import * # has build_tweaks(), and print_extras_activity()
from data_utils import *

RUNTIME = 3.0#6.0 # s # this must be a float, else takes as time steps

## Ensure that you use 10:1 aggregation of singles, if multiplicity = 'singles'
## so that with SYNS_PER_CLUBBED_SINGLE=10, there are no synchronous large IPSCs/IPSPs.
## But there will still be 10 IPSCs of same height for every single spike!
## For joints, 1:1 but there are 4x directed joints which'll show up say 1 in 4.
multiplicity = 'joints' # 'joints' / 'singles'

VOLTAGE_CLAMP = True

mitralmainidx = 0
## set lateral_mitnum to 1 for 2MITS / 2 for 2GLOMS option set in generate_neuroML.py
## Note that for directed connectivity, mit 3 is used for directed conn to mit 0 in generate_neuroml.py,
## thus mit 2 represents a non-directed conn cell.
## If you want to show asymm inhibition between directed cells, you should use mit 3 below.
lateral_mitnum = 2
mitralsidekickidx = lateral_mitnum

from pylab import * # part of matplotlib that depends on numpy but not scipy

## Output file name
today = datetime.date.today()
if NO_SINGLES: singles_str = '_NOSINGLES'
else: singles_str = '_SINGLES'
if NO_JOINTS: joints_str = '_NOJOINTS'
else: joints_str = '_JOINTS'
if NO_PGS: pgs_str = '_NOPGS'
else: pgs_str = '_PGS'
if IN_VIVO: invivo_str = '_invivo'
else: invivo_str = ''
if DIRECTED: dirt_str = '_directed'+str(FRAC_DIRECTED)
else: dirt_str = ''
if REVERSED_ADI: rev_str = '_reversed'
else: rev_str = ''
if ASYM_TEST: asym_str = '_asym'
else: asym_str = ''
if VOLTAGE_CLAMP: vstr = '_vclamp'
else: vstr = ''
#now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")+'_'
now = '' # stable enough to not bother about the date-time of simulation
outfilename = '../results/ADI/'+now+'IPSPs_PGvsgran_seed'+netseedstr+mitdistancestr+\
    singles_str+joints_str+pgs_str+invivo_str+dirt_str+rev_str+asym_str+vstr+'.pickle'
 
#-----------------------------------------------------------

def setup_stim(network):
    num_spikes = int(10*RUNTIME)
    if VOLTAGE_CLAMP:
        mitcell = network.mitralTable[mitralmainidx]
        mitsoma = moose.Compartment(get_matching_children(mitcell, ['Soma','soma'])[0])
        #### I block Ca and K channels in the cell
        #### We don't want these channels to be active when cell is held at 0mV; Na inactivates.
        blockChannels(mitcell, ['K','Ca']) # in moose_utils.py
        ## PID gain: I think for Davison 4-comp-mitral/granule: 0.5e-5 # optimal gain
        ## too high 0.5e-4 drives it to oscillate at high frequency,
        ## too low 0.5e-6 makes it have an initial overshoot (due to Na channels?)
        ## But for BBmit1993, gain of 1e-6 is optimal
        ## from moose.utils, setup_vclamp(compartment, name, delay1, width1, level1, gain=0.5e-5)
        network.PIDITable = setup_vclamp(mitsoma, '_somavclamp', 0, RUNTIME, 0e-3, gain=1e-6)
    ## network.populationDict = 
    ##      { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
    ## network.projectionDict = 
    ##      { 'projectionname1':(source,target,[(syn_name1,pre_seg_path,post_seg_path),...]) , ... }
    ## PGs are aggregated 5:1 with async synapses. So for every PG spike, there will be 5 IPSCs/IPSPs.
    if 'PGs' in network.populationDict:
        PGs = network.populationDict['PGs'][1]
        PG_conns = network.projectionDict['PG_mitral'][2]
        num_PG_conns = len(PG_conns)
        i = 0
        while i<num_spikes/PG_CLUB:
            ## find a PG-|mitral synapse on the main mitral, and make only that PG spike
            conn_id = int(uniform(num_PG_conns))
            conn = PG_conns[conn_id]
            ## get the mitnum from the post_seg_path
            mitname = string.split(conn[2],'/')[1] # name of the mitral cell from '/mitrals_2/...'
            mitnum = int(string.split(mitname,'_')[1]) # mit number from 'mitrals_2'
            if mitnum != mitralmainidx: continue
            
            ## get the PGnum from the post_seg_path
            PGname = string.split(conn[1],'/')[1] # name of the PG cell from '/PGs_2/...'
            PGnum = int(string.split(PGname,'_')[1]) # PG number from 'PGs_2'
            PG = PGs[PGnum]
            PGsoma = moose.Compartment(get_matching_children(PG, ['Soma','soma'])[0])
            ## We want a single spike at given time, so cannot use soma.inject
            setup_iclamp(PGsoma, '_PG'+str(PGnum),\
                uniform(0,RUNTIME), 1e-3, 200e-12) # start_time, duration, current (all SI)
            i += 1
    ## Ensure that you use 10:1 aggregation of singles, if multiplicity = 'singles'
    ## so that with SYNS_PER_CLUBBED_SINGLE=10, there are no synchronous large IPSCs/IPSPs.
    ## But there will still be 10 IPSCs of same height for every single spike!
    ## For joints, 1:1 but there are 4x directed joints which'll show up say 1 in 4.
    if 'granules_'+multiplicity in network.populationDict:
        grans = network.populationDict['granules_'+multiplicity][1]
        gran_conns = network.projectionDict['granule_mitral_inh_'+multiplicity][2]
        num_gran_conns = len(gran_conns)
        i = 0
        if multiplicity=='singles': club_factor = SYNS_PER_CLUBBED_SINGLE
        else: club_factor = GRANS_CLUB_JOINTS_2GLOMS
        while i<num_spikes/club_factor:
            conn_id = int(uniform(num_gran_conns))
            conn = gran_conns[conn_id]
            ## get the mitnum from the post_seg_path
            mitname = string.split(conn[2],'/')[1] # name of the mitral cell from '/mitrals_2/...'
            mitnum = int(string.split(mitname,'_')[1]) # mit number from 'mitrals_2'
            if mitnum != mitralmainidx: continue
            
            ## get the grannum from the post_seg_path
            granname = string.split(conn[1],'/')[1] # name of the granule cell from '/granules_joints_2/...'
            grannum = int(string.split(granname,'_')[2]) # granule number from 'granules_joints_2'
            gran = grans[grannum]
            gransoma = moose.Compartment(get_matching_children(gran, ['Soma','soma'])[0])
            ## We want a single spike at given time, so cannot use soma.inject
            setup_iclamp(gransoma, '_gran_'+multiplicity+str(grannum),\
                uniform(0,RUNTIME), 1e-3, 100e-9) # start_time, duration, current (all SI)
            i += 1

def run_inhibition(network, tables):
    resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
    network.context.step(RUNTIME)

def save_and_plot(network):
    if VOLTAGE_CLAMP:
        Vclamp_I = array(network.PIDITable)
    mit_Vm = array(network.mitralTable[mitralmainidx]._vmTableSoma)
    ## sometimes mitVm has one or two extra value
    ## so generate extra time values and discard if mitVm is not that long
    tlist = arange(0.0,RUNTIME+3*PLOTDT,PLOTDT)
    tlist = tlist[:len(Vclamp_I)]

    fvsifile = open(outfilename,'w')
    pickle.dump((tlist,Vclamp_I), fvsifile)
    fvsifile.close()
    print "Wrote",outfilename

    mainfig = figure(facecolor='w')
    mainaxes = mainfig.add_subplot(111)
    if VOLTAGE_CLAMP:
        starti = int(50e-3/PLOTDT)
        mainaxes.plot(tlist[starti:],Vclamp_I[starti:]*1e12,color='k')
        axes_labels(mainaxes,'time (s)','I (pA)')

    mainfig = figure(facecolor='w')
    mainaxes = mainfig.add_subplot(111)
    mainaxes.plot(tlist,mit_Vm*1e3,color='k')
    axes_labels(mainaxes,'time (s)','V (mV)')

    show()

#----------------------------------------------------

if __name__ == "__main__":
    seed([100.0])
    if os.path.exists(outfilename):
        print "Have to implement purely plotting if output file exists:", outfilename
        #sys.exit()
    uniquestr = 'pgvsgran_'
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
    ## spiketable is False that mitral Vm-s are recorded not spiketimes.
    network = OBNetwork(OBNet_file, synchan_activation_correction,\
        tweaks, mpirank, uniquestr, granfilebase, spiketable=False)
    #printNetTree() # from moose_utils.py

    setup_stim(network)
    ## if SPIKETABLE: record the Vm-s of a few interneurons
    ## else: record spiketimes of all interneurons
    tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,
        args={'mitrals':(mitralmainidx,mitralsidekickidx)}, spikes=SPIKETABLE)
    run_inhibition(network, tables)
    save_and_plot(network)
