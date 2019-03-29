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

RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

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

########## You need to run:
## From gj:
## ./restart_mpd_static
## The 0th (boss process) will always be node000 as it is the first node in ~/hostfile.
## HENCE FROM node000: cd to the working directory simulations/ (so that sys.path has accurate relative paths)
## mpiexec -machinefile ~/hostfile -n <numplotpts*2+1> ~/Python-2.6.4/bin/python2.6 activdep_inhibition.py
## I typically take numplotpts = 30
## nohup mpiexec -machinefile ~/hostfile -n 61 ~/Python-2.6.4/bin/python2.6 activdep_inhibition.py < /dev/null &
## For not showing the plots, append NOSHOW as a commandline argument.
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
    sys.exit(1)

## half jobs run with mitral B off, half with on, hence twice the step
Ainjectarray = arange(0.0,Imax,Imax/(mpisize-1)*2)
numpts = len(Ainjectarray)

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
#now =  datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")+'_'
now = '' # stable enough to not bother about the date-time of simulation
outfilename = '../results/ADI/'+now+'ADI_'+str(lateral_mitnum)+'_seed'+netseedstr+mitdistancestr+\
    singles_str+joints_str+pgs_str+invivo_str+dirt_str+rev_str+asym_str+'.pickle'

#-----------------------------------------------------------

def setup_stim(network):
    iAinject = Ainjectarray[(mpirank-1)%numpts]
    ## half jobs run with mitral B off, half with on
    if ASYM_TEST: # for asymmetry, inject same current in mitA and mitB
        iBinject = ((mpirank-1)/numpts) * iAinject # integer division
    else: # for ADI, inject current to get 80Hz in mitB
        iBinject = ((mpirank-1)/numpts) * onInject # integer division
    ipulse_duration = REALRUNTIME # seconds
    ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
    ## I seem to be using much larger currents - increase the inhibition
    iA = setup_iclamp(network.mitralTable[mitralmainidx].soma, '_mitralA',\
        SETTLETIME, ipulse_duration, iAinject)
    ## 1-1200pA for 400ms was used by Arevian et al to generate FvsI curves.
    ## slightly stagger the start of current pulses in the two cells
    ## so that the mitrals do not continuously co-fire.
    iB = setup_iclamp(network.mitralTable[mitralsidekickidx].soma, '_mitralB',\
        SETTLETIME-5e-3, ipulse_duration+5e-3, iBinject)
    print 'Injecting mitral A with '+str(iAinject)+' and B with '+\
        str(iBinject)+' at process = '+str(mpirank)+'.'
        

def run_inhibition(network, tables):
    resetSim(network.context, SIMDT, PLOTDT) # from moose_utils.py sets clocks and resets
    network.context.step(RUNTIME)
    ## get mitral A's firing rate
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
            if not NOSHOW:
                if PLOT_EXTRAS:
                    timevec = arange(0.0,RUNTIME+1e-12,PLOTDT)
                    titlestr = 'Ainject='+str(mitralAinject)+' Binject='+str(mitralBinject*onInject)
                    #figure()
                    #title('red:mitA, green:mitB, '+titlestr)
                    #plot(timevec, mitA, 'r,')
                    #plot(timevec, mitB, 'g,')
                    plot_extras(timevec, tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS, titlestr)
                else:
                    iAinject = mitralAinject
                    iBinject = mitralBinject * onInject
                    spikestables = \
                        print_extras_activity(tables, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS, \
                        'I_A='+str(iAinject)+'A & I_B='+str(iBinject)+'A.')
                    ## annotate each point with # spikes for singles;joints;multis
                    ## use index 2 for #spikes , index 0 for #cells spiking; index 1 for total #cells.
                    firingcells = ''
                    if not NO_SINGLES: firingcells = str(spikestables['singles'][2])
                    if not NO_JOINTS: firingcells += ';'+str(spikestables['joints'][2])
                    if not NO_MULTIS: firingcells += ';'+str(multistables['multis'][2])
                    mainaxes.annotate(firingcells,xy=(mitralAinject,Afiringrate))
        if not NOSHOW:
            mainaxes.plot(Ainjectarray, firingratearray, color=(mitralBinject,1-mitralBinject,0),\
                marker=['+','x'][mitralBinject], label="mitral B inject = "+str(mitralBinject*onInject)+" A.")
        dual_firingratearray.append(firingratearray)

    fvsifile = open(outfilename,'w')
    pickle.dump((Ainjectarray, dual_firingratearray), fvsifile)
    fvsifile.close()
    print "Wrote",outfilename

    if not NOSHOW:
        mainaxes.legend(loc="lower right")
        mainaxes.set_xlabel("mitral A inject (A)",fontsize=14)
        mainaxes.set_ylabel("mitral A firing rate (Hz)",fontsize=14)
        show()

#----------------------------------------------------

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
    if mpirank==boss:
        collate()
    else:
        if len(sys.argv)>2: uniquestr = sys.argv[2]+'_' # _ necessary, else say 'adi2'+mpirank is screwed
        else: uniquestr = 'adi_'
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

        setup_stim(network)
        ## if SPIKETABLE: record the Vm-s of a few interneurons
        ## else: record spiketimes of all interneurons
        tables = setupTables(network, NO_PGS, NO_SINGLES, NO_JOINTS, NO_MULTIS,
            args={'mitrals':(mitralmainidx,)}, spikes=SPIKETABLE)
        run_inhibition(network, tables)
