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
from sim_utils import * # has build_tweaks()
from data_utils import *

## constrate firefiles start firing immediately (no SETTLETIME)
## Hence just have 20ms as SETTLETIME. overrides that from stimuliConstants.py.
## But PGs have an initial spike on a high initial settling transient, so 100ms:
SETTLETIME = 100e-3#20e-3 # s
RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

## set lateral_mitnum to 1 for 2MITS / 2 for 2GLOMS option set in generate_neuroML.py
## Note that for directed connectivity, mit 3 is used for directed conn to mit 0 in generate_neuroml.py,
## thus mit 2 represents a non-directed conn cell.
lateral_mitnum = 3#2#3
if REVERSED_ADI:
    mitralmainidx = lateral_mitnum
    mitralsidekickidx = 0
else:
    mitralmainidx = 0
    mitralsidekickidx = lateral_mitnum

########## You need to run:
## From gj:
## ./restart_mpd_static
## The 0th (boss process) will always be node000 as it is the first node in ~/hostfile.
## HENCE FROM node000: cd to the working directory simulations/ (so that sys.path has accurate relative paths)
## mpiexec -machinefile ~/hostfile -n <numplotpts*2+1> ~/Python-2.6.4/bin/python2.6 inhibition_tuftinput.py
## I typically take numplots = 20 (generate_firefiles_constrate.py should also be called accordingly first)
## nohup mpiexec -machinefile ~/hostfile -n 41 ~/Python-2.6.4/bin/python2.6 inhibition_tuftinput.py &> nohup_tuftinh.out < /dev/null &
## For not showing the plots, append NOSHOW as a commandline argument.

## for lateral propagation of AP :
## You need to set seed,mitdist,mitdiststr (simset_activinhibition_minimal.py)
## and segnames (in code below), and run :
## python2.6 inhibition_tuftinput.py

##### 0 rank process is for collating all jobs. (rank starts from 0)
##### I assume rank 0 process always runs on the machine whose X window system has a Display connected
##### and can show the graphs!!!!!!
##### The rank 0 stdout is always directed to the terminal from which mpiexec was run.
##### I hope X output also works the same way.
##### For long simulations save results in a text file for replotting later and avoid above ambiguity.
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
    print "Check if PLOT_EXTRAS in simset_activinhibition.py is True."
    sys.exit(1)

## half jobs run with mitral B off, half with on, hence twice the step
if mpisize>1: Ainjectarray = arange(0.0,ORNmax,ORNmax/float(mpisize-1)*2)
else: Ainjectarray = [15.0]
numpts = len(Ainjectarray)
onInject = oninject_ext # Hz # from simset_activinhibition_minimal.py
ORN_constrate_str = '../firefiles/firefiles_constrate/firetimes_constrate'

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
if ASYM_TEST: asym_str = '_asym' ## asym test acts as same ORN frate into A and B
else: asym_str = ''
if ODORINH: odorinh_str = '_odorinh'
else: odorinh_str = '_airinh'
#now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")+'_'
now = '' # stable enough to not bother about the date-time of simulation
outfilename = '../results/tuftADI/'+now+'tuftADI_'+str(lateral_mitnum)+'_seed'+netseedstr+mitdistancestr+\
    singles_str+joints_str+pgs_str+invivo_str+dirt_str+rev_str+asym_str+odorinh_str+'.pickle'

#-----------------------------------------------------------

def setup_stim(network,mpirank_given):
    iAinject = Ainjectarray[(mpirank_given-1)%numpts]
    ## half jobs run with mitral B off, half with on
    if ASYM_TEST: # for asymmetry, inject same current in mitA and mitB
        iBinject = ((mpirank_given-1)/numpts) * iAinject # integer division
    else: # for ADI, inject current to get 80Hz in mitB
        iBinject = ((mpirank_given-1)/numpts) * onInject # integer division
    ipulse_duration = REALRUNTIME # seconds

    setup_tuftinject(network, mitralmainidx, iAinject, mpirank_given)        
    setup_tuftinject(network, mitralsidekickidx, iBinject, mpirank_given)        
    print 'Tuft ORN injecting mitral A with '+str(iAinject)+' Hz and B with '+\
        str(iBinject)+' Hz at process = '+str(mpirank_given)+'.'

def setup_tuftinject(network, mitnum, frate, mpirank_given=mpirank):
    ## first figure out those PGs connected to this mitral mitnum
    PGs_to_thismit = []
    for projname in network.projectionDict.keys():
        if 'PG_mitral' in projname:
            for i,proj in enumerate(network.projectionDict[projname][2]):
                ## get the glomnum from the post path proj[2]
                ## name of the mitral cell from '/mitrals_2/...'
                mitname = string.split(proj[2],'/')[1]
                ## mit number from 'mitrals_2'
                conn_mitnum = int(string.split(mitname,'_')[1])
                if conn_mitnum == mitnum:
                    ## name of the PG cell from '/PGs_2/...'
                    pgname = string.split(proj[1],'/')[1]
                    PGs_to_thismit.append(pgname)
    ## Now connect the ORNs
    ORNfilebase = ORN_constrate_str + str(frate)+'_trial'+str(mpirank_given-1)
    for projname in network.projectionDict.keys():
        #### Calling attach_spikes() for each projection,
        #### would reconnect files to the same segment multiple times.
        #### But attach_files_uniquely() checks whether timetable.tableSize is zero or not
        #### i.e. files already attached or not.
        ## connect ORNs to mitrals
        if 'ORN_mitral' in projname:
            print "connecting ORN files to mitrals for mitnum",mitnum
            for i,proj in enumerate(network.projectionDict[projname][2]):
                ## get the mitnum from the post path proj[2]
                ## name of the mitral cell from '/mitrals_2/...'
                mitname = string.split(proj[2],'/')[1]
                ## mit number from 'mitrals_2'
                conn_mitnum = int(string.split(mitname,'_')[1])
                if conn_mitnum == mitnum:
                    attach_files_uniquely(ORNfilebase,proj[0],proj[2])
        ## connect ORNs to PG
        if 'ORN_PG' in projname:
            print "connecting ORN files to PGs for mitnum",mitnum
            for i,proj in enumerate(network.projectionDict[projname][2]):
                pgname = string.split(proj[2],'/')[1] # name of the PG cell from '/PGs_2/...'
                if pgname in PGs_to_thismit:
                    attach_files_uniquely(ORNfilebase,proj[0],proj[2])

def attach_files_uniquely(filebase,synname,postsegpath):
    ttpath = postsegpath+'/'+synname+'_tt'
    if moose.context.exists(ttpath):
        # timetable already created by networkML reader - just wrap it below.
        tt = moose.TimeTable(ttpath) # post_segment_path+'/'+syn_name+'_tt'
    else:
        ## if timetable was not already created by networkML reader,
        ## it means that the synaptic weights must be zero!
        ## (no extra inhibition - only main inhibition)
        ## hence do not attach spikefiles
        return
    if tt.tableSize != 0: return # if files are already attached, do nothing!
    ## attach_spikes() accesses filenumbers to this segment
    ## from 'fileNumbers' field (of the timetable object in MOOSE)
    ## which is created while reading in networkML.
    attach_spikes(filebase, tt, uniquestr+str(mpirank))

def run_inhibition(network, tables):
    resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
    network.context.step(RUNTIME)
    # get mitral A's firing rate
    oneoverISI, meanfreq, events = \
        calcFreq(network.mitralTable[mitralmainidx]._vmTableSoma,\
            RUNTIME, SETTLETIME, PLOTDT, THRESHOLD, SPIKETABLE)
    mpicomm.send( meanfreq, dest=boss, tag=0 ) # frequency tag
    oneoverISI, meanfreq, events = \
        calcFreq(network.mitralTable[mitralsidekickidx]._vmTableSoma,\
            RUNTIME, SETTLETIME, PLOTDT, THRESHOLD, SPIKETABLE)
    mpicomm.send( meanfreq, dest=boss, tag=1 ) # frequency of mit B tag
    #print 'Firing rate = '+str(result[3])+'Hz on injecting mitral A with '+str(iAinject)+\
    #    ' and B with '+str(iBinject)+' at process = '+str(mpirank)+'.'
    ##mpicomm.send( array(network.mitralTable[mitralmainidx]._vmTableSoma), dest=boss, tag=1 )
    ##mpicomm.send( array(network.mitralTable[mitralsidekickidx]._vmTableSoma), dest=boss, tag=2 )
    mpicomm.send( numpy_convert_tables(tables), dest=boss, tag=3 ) # extra tables tag   
    print 'Sent output from process '+str(mpirank)+'.'


def collate():
    if not NOSHOW:
        mainfig = figure(facecolor='w')
        mainaxes = mainfig.add_subplot(111)
    dual_firingratearray = []
    for mitralBinject in [0,1]:
        firingratearray = []
        for A,mitralAinject in enumerate(Ainjectarray):
            procnum = mitralBinject*numpts + A + 1
            print 'waiting for process '+str(procnum)+'.'
            Afiringrate = mpicomm.recv(source=procnum, tag=0)
            firingratearray.append(Afiringrate)
            Bfiringrate = mpicomm.recv(source=procnum, tag=1)
            print "mitral B firing at",Bfiringrate
            ##mitA = mpicomm.recv(source=procnum, tag=1)
            ##mitB = mpicomm.recv(source=procnum, tag=2)
            tables = mpicomm.recv(source=procnum, tag=3)
            if PLOT_EXTRAS:
                timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
                titlestr = 'Ainject='+str(mitralAinject)+'Hz Binject='+str(mitralBinject*onInject)+'Hz'
                if not NOSHOW:
                    #figure()
                    #title('red:mitA, green:mitB, '+titlestr)
                    #plot(timevec, mitA, 'r,')
                    #plot(timevec, mitB, 'g,')
                    plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, titlestr)
            else:
                iAinject = mitralAinject
                iBinject = mitralBinject * onInject
                activity_table = print_extras_activity(tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,\
                    'I_A='+str(iAinject)+'Hz & I_B='+str(iBinject)+'Hz.')
                ## annotate each point with #joints;#singles firing
                firingcells = ''
                if not NO_JOINTS:
                    firingcells += str(activity_table['joints'][0])
                    if not NO_SINGLES: firingcells += ';'+str(activity_table['singles'][0])
                elif not NO_SINGLES: firingcells += str(activity_table['singles'][0])
                if not NOSHOW:
                    mainaxes.annotate(firingcells,xy=(mitralAinject,Afiringrate))
        if not NOSHOW:
            mainaxes.plot(Ainjectarray, firingratearray, color=(mitralBinject,1-mitralBinject,0),\
                marker=['+','x'][mitralBinject], label="mitral B ORNs = "+str(mitralBinject*onInject)+" Hz.")
        dual_firingratearray.append(firingratearray)

    fvsifile = open(outfilename,'w')
    pickle.dump((Ainjectarray, dual_firingratearray), fvsifile)
    fvsifile.close()
    print "Wrote",outfilename

    if not NOSHOW:
        mainaxes.legend(loc="lower right")
        mainaxes.set_xlabel("mit A tuft ORNs mean firing rate (Hz)",fontsize=16)
        mainaxes.set_ylabel("mitral A firing rate (Hz)",fontsize=16)
        show()

#----------------------------------------------------

def varyVrest_inject(network,popname,raiseRMP,inject):
    for cell in network.populationDict[popname][1].values():
        tweak_field(cell.path+'/##[TYPE=Compartment]', 'Em', 'Em+'+str(raiseRMP))
        tweak_field(cell.path+'/##[TYPE=Compartment]', 'inject', str(inject))

if __name__ == "__main__":
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
    if mpisize == 1:
        if len(sys.argv)>2: uniquestr = sys.argv[2]+'_' # _ necessary, else say 'adi2'+mpirank is screwed
        else: uniquestr = 'tuft_'
        ## includeProjections gets used only if ONLY_TWO_MITS is True:
        ## Keep below projections to 'second order cells'
        ## i.e. to cells (granules) connected to mits0&1.
        ## The connections between second order cell
        ## and mits0&1 are automatically retained of course.
        ## important to include 'PG' below as 'ORN_PG' is needed,
        ## though 'PG_mitral', 'mitral_PG' connections to/from mits0&1 are kept automatically.
        includeProjections = ['granule_baseline','PG']
        tweaks = build_tweaks( CLUB_MITRALS, NO_SPINE_INH,\
            NO_SINGLES, NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
            includeProjections, (mitralmainidx,mitralsidekickidx) )
        ## send mpirank to put in ORN filenames
        ## so they do not clash between processes - not needed here
        network = OBNetwork(OBNet_file, synchan_activation_correction,\
            tweaks, mpirank, uniquestr, granfilebase, spiketable=SPIKETABLE)
        #printNetTree() # from moose_utils.py

        ## Vary the resting potential of PGs
        if not NO_PGS: varyVrest_inject(network,'PGs',PG_raiseRMP,PG_inject)
        ## Vary the resting potential of granules
        if not NO_SINGLES: varyVrest_inject(network,'granules_singles',granule_raiseRMP,granule_inject)
        if not NO_JOINTS: varyVrest_inject(network,'granules_joints',granule_raiseRMP,granule_inject)
        if not NO_MULTIS: varyVrest_inject(network,'granules_multis',granule_raiseRMP,granule_inject)

        setup_stim(network,2) # 1 to give no i/p, 2 to give oninject_ext to mitB.
        
        ## setup a separate Vm table, since the one setup by OBNetwork() above stores only spikes.
        mitAsoma_Vm = setupTable('mitAsomaVm',network.mitralTable[mitralmainidx].soma,'Vm')
        #mitBsoma_Vm = setupTable('mitBsomaVm',network.mitralTable[mitralsidekickidx].soma,'Vm')
        mitBsoma_Vm = setupTable('mitBsomaVm',moose.Compartment('/mitrals_3/Seg0_soma_0'),'Vm')
        ## get the nearA and afterA ids for each netfile instance -- see my notes of 17Apr2013 and 5Feb2013.
        mitB_nearA_Vm = setupTable('mitB_nearA_Vm',moose.Compartment('/mitrals_3/Seg0_sec_dendd4_0_261'),'Vm')
        mitB_afterA_Vm = setupTable('mitB_afterA_Vm',moose.Compartment('/mitrals_3/Seg0_sec_dendd4_1_265'),'Vm')

        ## if SPIKETABLE (in simset_activinhibition.py): record the Vm-s of a few interneurons
        ## else: record spiketimes of all interneurons
        tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,
            args={'mitrals':(mitralmainidx,)}, spikes=SPIKETABLE)

        resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
        network.context.step(RUNTIME)

        mitAsoma_Vm = array(mitAsoma_Vm)
        mitBsoma_Vm = array(mitBsoma_Vm)
        timevec_unbinned = arange(0.0,RUNTIME+1e-12,PLOTDT)
        ## plot mitral A/B soma Vm
        figure(facecolor='w')
        title('Vm in mitral A soma')
        plot(timevec_unbinned,mitAsoma_Vm,label='Vm',color='k')
        axes_labels(gca(),"time (s)","Vm (V)",fontsize=label_fontsize+2)
        figure(facecolor='w')
        title('Vm in mitral B soma')
        plot(timevec_unbinned,mitBsoma_Vm,label='Vm',color='k')
        plot(timevec_unbinned,mitB_nearA_Vm,label='Vm',color='b')
        plot(timevec_unbinned,mitB_afterA_Vm,label='Vm',color='r')
        axes_labels(gca(),"time (s)","Vm (V)",fontsize=label_fontsize+2)
        
        ## Paper figure: mitB i.e. directed mitral 3's soma vs inhibitory sites
        fig2 = figure(figsize=(columnwidth/2.0,linfig_height/2.0),dpi=1200,facecolor='w') # 'none' is transparent
        ax = fig2.add_subplot(111)
        plot(timevec_unbinned*1000,array(mitBsoma_Vm)*1000,'-r', label='soma',\
            linewidth=0.5, linestyle='-') # ms and mV
        plot(timevec_unbinned*1000,array(mitB_nearA_Vm)*1000,'-b', label='at inh site',\
            linewidth=0.5, dashes=(1.5,0.5)) # ms and mV
        plot(timevec_unbinned*1000,array(mitB_afterA_Vm)*1000,'-k', label='after inh site',\
            linewidth=0.5, dashes=(0.5,0.5)) # ms and mV
        beautify_plot(ax,x0min=False,drawxaxis=True,drawyaxis=True,\
            xticks=[150,200],yticks=[-75,0,40])
        #add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        #    sizex=50,labelx='50 ms',sizey=20,labely='20 mV',label_fontsize=5,\
        #    bbox_to_anchor=[0.7,0.6],bbox_transform=ax.transAxes)
        ax.set_xlim(150,200)
        ax.set_xticklabels(['150','200'])
        ax.set_ylim(-75,40)
        axes_labels(ax,"time (ms)","Vm (mV)",fontsize=label_fontsize,xpad=-4,ypad=-3)
        fig2.tight_layout()
        fig2.savefig('../figures/connectivity/cells/mitral_spikeprop.png',dpi=fig2.dpi)
        fig2.savefig('../figures/connectivity/cells/mitral_spikeprop.svg',dpi=fig2.dpi)
        ## Paper figure done

        timevec = arange(0.0,RUNTIME+PLOTDT+1e-12,PLOTDT)
        if SPIKETABLE:
            ## 50ms bins for firingrate
            timevec = arange(50e-3/2.0,REALRUNTIME,50e-3)
            plot_extras_spikes(timevec, tables, NO_PGS, NO_SINGLES, \
                NO_JOINTS, NO_MULTIS, len(timevec), RUNTIME, SETTLETIME)
        else:
            plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS)
        
        show()

    else:
        if mpirank==boss:
            collate()
        else:
            if len(sys.argv)>2: uniquestr = sys.argv[2]+'_' # _ necessary, else say 'adi2'+mpirank is screwed
            else: uniquestr = 'tuft_'
 
            ## If CLUB_MITRAL=False, then extra exc from mitral sisters
            ## (to certain connected granules as proxy for unmodelled sisters) does NOT get used.
            ## Instead, here I connect extra baseline excitation to ALL granules if odor is non-zero.
            ## This extra basline is only if ODORINH i.e. we are simulating odor inh and not air inh
            if ODORINH and not CLUB_MITRALS:
                granfilebase += '_extra'

            ## includeProjections gets used only if ONLY_TWO_MITS is True:
            ## Keep below projections to 'second order cells'
            ## i.e. to cells (granules) connected to mits0&1.
            ## The connections between second order cell
            ## and mits0&1 are automatically retained of course.
            ## important to include 'PG' below as 'ORN_PG' is needed,
            ## though 'PG_mitral', 'mitral_PG' connections to/from mits0&1 are kept automatically.
            includeProjections = ['granule_baseline','PG']
            tweaks = build_tweaks( CLUB_MITRALS, NO_SPINE_INH,\
                NO_SINGLES, NO_JOINTS, NO_MULTIS, NO_PGS, ONLY_TWO_MITS,\
                includeProjections, (mitralmainidx,mitralsidekickidx) )
            ## send mpirank to put in ORN filenames
            ## so they do not clash between processes - not needed here
            network = OBNetwork(OBNet_file, synchan_activation_correction,\
                tweaks, mpirank, uniquestr, granfilebase, spiketable=SPIKETABLE)
            #printNetTree() # from moose_utils.py

            ### Vary the resting potential of PGs
            #if not NO_PGS: varyVrest_inject(network,'PGs',PG_raiseRMP,PG_inject)
            ### Vary the resting potential of granules
            #if not NO_SINGLES: varyVrest_inject(network,'granules_singles',granule_raiseRMP,granule_inject)
            #if not NO_JOINTS: varyVrest_inject(network,'granules_joints',granule_raiseRMP,granule_inject)
            #if not NO_MULTIS: varyVrest_inject(network,'granules_multis',granule_raiseRMP,granule_inject)

            setup_stim(network, mpirank)
            ## if SPIKETABLE: record the Vm-s of a few interneurons
            ## else: record spiketimes of all interneurons
            tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,
                args={'mitrals':(mitralmainidx,)}, spikes=SPIKETABLE)
            run_inhibition(network, tables)
