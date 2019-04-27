#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import pickle

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose.utils import * # imports moose

from stimuliConstants import * # has SETTLETIME
from simset_activinhibition import * # has REALRUNTIME
from data_utils import *
from plot_mitral_connectivity_NetworkML import * # get connectivity from netfile

RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

ASYM_TEST = False
DIRECTED = True
FRAC_DIRECTED = 0.01#0.003#0.03
IN_VIVO = False#True
if IN_VIVO: invivo_str = '_invivo'
else: invivo_str = ''
if DIRECTED: dirt_str = '_directed'+str(FRAC_DIRECTED)
else: dirt_str = ''

USE_BEST_HALF = True

## inh =  (no_singles,no_joints,no_PGs)
inh = (False,False,False)
netseeds = arange(100.0,5000.0,100.0)

if IN_VIVO:
    distance_list = [50.0, 200.0, 400.0, 600.0]#, 800.0] # microns
    reverse_list = [True,False]
else:
    distance_list = [50.0] # microns
    reverse_list = [False]
## whether to calculate and print the number of joint granules connected
SHOW_CONN = False

########## You need to run: python2.6 average_activdep_inhibition.py

## Read in the .csv file having the average ADI of Arevian et al.
## numpy's genfromtxt() to read csv files
## Arevian_mean_ADI['x'], Arevian_mean_ADI['Curve1'] contain x and y lists
Arevian_mean_ADI = genfromtxt('arevianetal_average_ADI.csv', delimiter=',', dtype=None, names=True)


def get_filename(netseed,inh,mitdistance,reverse,latmitnum=2):
    ### read filename from the output file of automated run
    ## construct the filename
    if inh[0]: singles_str = '_NOSINGLES'
    else: singles_str = '_SINGLES'
    if inh[1]: joints_str = '_NOJOINTS'
    else: joints_str = '_JOINTS'
    if inh[2]: pgs_str = '_NOPGS'
    else: pgs_str = '_PGS'
    if reverse: rev_str = '_reversed'
    else: rev_str = ''
    if ASYM_TEST: asym_str = '_asym'
    else: asym_str = ''
    ## stable enough that time tags are not needed
    outfilename = '../results/ADI/ADI_'+str(latmitnum)+\
        '_seed'+str(netseed)+'_mitdist'+str(mitdistance)+\
        singles_str+joints_str+pgs_str+invivo_str+dirt_str+\
        rev_str+asym_str+'.pickle'
    return outfilename, (singles_str, joints_str, pgs_str)

def get_rates_inh(filename):
    f = open(filename,'r')
    Ainjectarray, both_firingratearrays = pickle.load(f)
    f.close()
    mit_alone = array(both_firingratearrays[0])
    mit_inhibited = array(both_firingratearrays[1])
    only_reasonable_inh = where(Ainjectarray<1.5e-9)
    diff_frates = mit_alone[only_reasonable_inh] - mit_inhibited[only_reasonable_inh]
    max_redux = max(diff_frates)
    avg_redux = mean(diff_frates)
    ## take max f rate redux from 6 to 55Hz
    only_reasonable_55Hz = where(Ainjectarray<1.2e-9)
    mit_alone_55Hz = mit_alone[only_reasonable_55Hz]
    only_to_55Hz = where(mit_alone_55Hz<=55.0)
    only_6Hz_to_55Hz = where(mit_alone_55Hz[only_to_55Hz]>=6.0)
    diff_frates = mit_alone_55Hz[only_6Hz_to_55Hz] - mit_inhibited[only_6Hz_to_55Hz]
    max_asym_redux = max(diff_frates)
    return Ainjectarray, both_firingratearrays, max_redux, avg_redux, max_asym_redux

def plot_avg_inh(reverse, mit_distance, title_str):
    ## load the data
    full_firingratearray = []
    redux_frate_array = []
    ## 0 to 140 Hz into 10 Hz bins
    maxfrate = 140
    fratebinsize = 10
    redux_vs_frate_array = [ [] for i in range(maxfrate/fratebinsize) ]
    max_redux_frate_array = []
    avg_redux_frate_array = []
    max_asym_redux_frate_array = []
    mit_somadend_joints_avg = array([0,0])
    num_available = 0
    for netseed in netseeds:
        filename,_ = get_filename(netseed,inh,mit_distance,reverse)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        print 'checking for',filename
        if not os.path.exists(filename): continue
        num_available += 1
        print 'found',filename

        Ainjectarray, both_firingratearrays, max_redux, avg_redux, max_asym_redux = \
            get_rates_inh(filename)
        full_firingratearray.append(both_firingratearrays)
        redux_frate_array_this = array(both_firingratearrays[1])-array(both_firingratearrays[0])
        redux_frate_array.append(redux_frate_array_this)
        max_redux_frate_array.append(max_redux)
        avg_redux_frate_array.append(avg_redux)
        max_asym_redux_frate_array.append(max_asym_redux)
        
        if SHOW_CONN:
            ## get connectivity details about this netfile
            bits = filename.split('_')
            netseedstr = bits[1]
            mitdistancestr = bits[2]
            OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed20_jointsclubbed1'\
                '_numgloms2_'+netseedstr+'_'+mitdistancestr
            if directed: OBNet_file += '_directed0.0_proximal'
            OBNet_file += '_2GLOMS'
            if not IN_VIVO: OBNet_file += '_INVITRO.xml'
            else: OBNet_file += '.xml'
            netml = NetworkML(twogloms=True)
            tweaks = {}
            ## mit_posn_grannos = 
            ## (mit_prim_joints,mit_somadend_joints,mit_sec_joints,mit_prim_singles,mit_sec_singles)
            (mitralsDict,connectivityDict,mit_posn_grannos) = \
                netml.readNetworkMLFromFile(OBNet_file,params=tweaks)
            mit_somadend_joints = mit_posn_grannos[1] # dict {0:0,mitB:0}, mitB=2 if TWOGLOMS
            mit_somadend_joints_avg[0] += mit_somadend_joints[0]
            mit_somadend_joints_avg[1] += mit_somadend_joints[2]

    ## select 50% of the strongest ADI : may or may not use for plots: see below
    ## Currently not using only strongest.
    ## get indices of the best ADI, sorted according to max redux in frate, to keep them
    low_idx = [i[0] for i in sorted(enumerate(max_redux_frate_array), key=lambda x:x[1], reverse=True)]
    low_idx = low_idx[0:num_available/2]
    keep_arrays = []
    keep_redux_frate_array = []
    for idx in low_idx:
        keep_arrays.append(full_firingratearray[idx])
        keep_redux_frate_array.append(redux_frate_array[idx])
    keep_arrays = array(keep_arrays)
    full_firingratearray = array(full_firingratearray)

    ## calculate means and std errs
    ## comment/uncomment one of the two blocks below to take mean of top 50% OR all
    ## mean inh trough of 50% is approx 25% lower than that of all
    ## mean of 50%
    if USE_BEST_HALF:
        print "Using means of top 50% of the mitral pairs having most inh redux."
        mean_firingratearray = mean(keep_arrays, axis=0)
        std_firingratearray = std(keep_arrays, axis=0)
        redux_frate_array = keep_redux_frate_array
        nummitpairs = len(keep_arrays)
    else:
    ## mean of all
        print "Using means of all mitral pairs, rather than top half."
        mean_firingratearray = mean(full_firingratearray, axis=0)
        std_firingratearray = std(full_firingratearray, axis=0)
        nummitpairs = len(full_firingratearray)

    print "Total number of mitral pairs used for means =",nummitpairs
    ## bin the reduction so that you get reduction vs firing rate of A (rather than I_A)
    ## some bins will have more elements than others
    for redux_frate_array_this in redux_frate_array:
        for i,frateA in enumerate(both_firingratearrays[0]):
            redux_vs_frate_array[int(frateA/fratebinsize)].append(redux_frate_array_this[i])

    if SHOW_CONN:
        mit_somadend_joints_avg /= float(len(filelist))
        print "Avg # of granules to soma and nearest segments to mits A and B are",\
            mit_somadend_joints_avg

    ## plot the frate redux and percentage redux in firing rates
    frateA_array = []
    frateA_running = fratebinsize/2.0
    frate_redux_mean_array = []
    frate_redux_std_array = []
    for arr in redux_vs_frate_array:
        if arr:
            frate_redux_mean_array.append( mean(arr) )
            frate_redux_std_array.append( std(arr) )
            frateA_array.append(frateA_running)
        frateA_running += fratebinsize
    redux_perc_mean_array = array(frate_redux_mean_array)/array(frateA_array)*100
    redux_perc_std_array = array(frate_redux_std_array)/array(frateA_array)*100
    
    ## simulation of fig 2e of Arevian et al
    fig = figure(figsize=(columnwidth*3/4, linfig_height/2.0), dpi=fig_dpi, facecolor='w')
    ### plot reduction in percentage
    #ax2 = fig.add_subplot(111,frameon=False)
    #ax2.fill_between(frateA_array, \
    #    redux_perc_mean_array + redux_perc_std_array/sqrt(nummitpairs), \
    #    redux_perc_mean_array - redux_perc_std_array/sqrt(nummitpairs), \
    #    color='b',alpha=0.4,linewidth=0)
    #ax2.plot(frateA_array, redux_perc_mean_array, \
    #    color='b',marker='s',ms=marker_size,linewidth=linewidth)
    #ax2.set_ylim(-30,0)
    #ax2.set_yticks([-30,-20,-10,0])
    #ax2.add_artist(Line2D((0, 0), (-30, 0), \
    #    color='black', linewidth=axes_linewidth))
    #ax2.get_xaxis().set_ticks_position('bottom')
    #axes_labels(ax2,"mitral A firing rate (Hz)","change in rate of A (%)")
    
    ## plot delta rate in Hz
    #ax = ax2.twinx()
    ax = fig.add_subplot(111)
    ax.plot(Arevian_mean_ADI['x'],Arevian_mean_ADI['Curve1'],\
        color='b',dashes=(1.5,0.5),linewidth=linewidth)
    ax.fill_between(frateA_array, \
        array(frate_redux_mean_array) + array(frate_redux_std_array)/sqrt(nummitpairs), \
        array(frate_redux_mean_array) - array(frate_redux_std_array)/sqrt(nummitpairs), \
        color='r',alpha=0.4,linewidth=0)
    ax.plot(frateA_array, frate_redux_mean_array, \
        color='r',marker='o',ms=marker_size,linewidth=linewidth)
    ax.set_xlim(0,maxfrate)
    ax.set_xticks([0,maxfrate/4,maxfrate/2,maxfrate*3/4,maxfrate])
    ax.set_ylim(-15,0)
    ax.set_yticks([-15,-10,-5,0])
    axes_labels(ax,"rate of A (Hz)","$\Delta$ rate of A (Hz)")
    ax.xaxis.set_label_position('top')
    ax.get_yaxis().set_ticks_position('left')
    ax.get_xaxis().set_ticks_position('top')
    ## set xticks text sizes
    for label in ax.get_xticklabels():
        label.set_fontsize(label_fontsize)
    ## hide the top and left axes
    for loc, spine in ax.spines.items(): # items() returns [(key,value),...]
        spine.set_linewidth(axes_linewidth)
        if loc in ['right','bottom']:
            spine.set_color('none') # don't draw spine
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig('../figures/ADI/ADI'+invivo_str+dirt_str+'.svg',dpi=fig.dpi)
    
    ## Plot the best ADI-s separately also:
    for ploti,idx in enumerate(low_idx):
        mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
        mainaxes = mainfig.add_subplot(111,frameon=False)
        best_firingratearray = full_firingratearray[idx]
        for mitralBinject in [0,1]:
            if ASYM_TEST:
                if mitralBinject == 0: mitB_str = '0.0 A'
                else: mitB_str = 'A inject'
            else: mitB_str = str(mitralBinject*onInject)+' A'
            mainaxes.plot(Ainjectarray*1e9, best_firingratearray[mitralBinject],\
            color=(mitralBinject,0,0), marker=['s','o'][mitralBinject],\
            label="mitral B inject = "+mitB_str, linewidth=linewidth, markersize=marker_size)
        mainaxes.get_yaxis().set_ticks_position('left')
        mainaxes.get_xaxis().set_ticks_position('bottom')
        mainaxes.set_xticks([0,1,2.5])
        mainaxes.set_ylim(0,120)
        mainaxes.set_yticks([0,20,40,60,80,100,120])
        axes_labels(mainaxes,"mitral A inject (nA)","mitral A firing rate (Hz)") # enlarges x and y ticks and labels
        xmin, xmax = mainaxes.get_xaxis().get_view_interval()
        ymin, ymax = mainaxes.get_yaxis().get_view_interval()
        mainaxes.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
        mainaxes.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
        fig_clip_off(mainfig)
        mainfig.tight_layout()
        mainfig.savefig('../figures/ADI/ADI' \
            +invivo_str+dirt_str+rev_str+'_'+str(idx)+'.png', \
            dpi=mainfig.dpi)
        mainfig.savefig('../figures/ADI/ADI' \
            +invivo_str+dirt_str+rev_str+'_'+str(idx)+'.svg', \
            dpi=mainfig.dpi)

    ## plot the mean ADI firing rate vs I curves with B on and off
    mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
    mainaxes = mainfig.add_subplot(111,frameon=False)
    for mitralBinject in [0,1]:
        if ASYM_TEST:
            if mitralBinject == 0: mitB_str = '0.0 A'
            else: mitB_str = 'A inject'
        else: mitB_str = str(mitralBinject*onInject)+' A'
        mainaxes.errorbar(x=Ainjectarray*1e9,\
        y=mean_firingratearray[mitralBinject],yerr=std_firingratearray[mitralBinject],\
        color=(mitralBinject,0,0), marker=['s','o'][mitralBinject], capsize=cap_size,\
        label="mitralB I = "+mitB_str, linewidth=linewidth, markersize=marker_size)
    #if IN_VIVO: biglegend("lower right")#"upper left")
    #else: biglegend("lower right")
    axes_labels(mainaxes,"mitral A inject (nA)","mitral A firing rate (Hz)")
    mainaxes.set_title('mean '+title_str,fontsize=label_fontsize)
    mainaxes.get_yaxis().set_ticks_position('left')
    mainaxes.get_xaxis().set_ticks_position('bottom')
    mainaxes.set_xticks([0,1,2.5])
    mainaxes.set_ylim(0,120)
    mainaxes.set_yticks([0,20,40,60,80,100,120])
    xmin, xmax = mainaxes.get_xaxis().get_view_interval()
    ymin, ymax = mainaxes.get_yaxis().get_view_interval()
    mainaxes.set_xlim(0,xmax)
    mainaxes.set_ylim(0,ymax)
    mainaxes.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
    mainaxes.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
    fig_clip_off(mainfig)
    mainfig.tight_layout()
    mainfig.savefig('../figures/ADI/ADI'+invivo_str+'_'+title_str+'.png',\
        bbox_inches='tight',dpi=mainfig.dpi)

    ## print the max reduction in firing rates
    print 'Max redux '+title_str+' till IA=1.5nA =',max_redux_frate_array
    print 'Avg redux '+title_str+' till IA=1.5nA =',avg_redux_frate_array
    print 'Max redux '+title_str+' 6Hz to 55Hz =',max_asym_redux_frate_array
    return max_redux_frate_array, avg_redux_frate_array, max_asym_redux_frate_array

if __name__ == "__main__":
    redux = []
    redux_std = []
    for i,reverse in enumerate(reverse_list):
        print "Reversed =", reverse
        if reverse: rev_str = "Reversed"
        else: rev_str = "Forward"
        redux_directed = []
        redux_directed_std = []
        for j,mit_distance in enumerate(distance_list):
            print 'Distance =',mit_distance,'microns:'
            max_redux_frate_array, avg_redux_frate_array, max_asym_redux_frate_array = \
                plot_avg_inh(reverse, mit_distance, rev_str+' @ '+str(mit_distance)+' microns')
            redux_directed.append(mean(avg_redux_frate_array))
            redux_directed_std.append(std(avg_redux_frate_array))
        redux.append((rev_str,redux_directed))
        redux_std.append(redux_directed_std)
    if IN_VIVO:
        fig = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w')
        ax = fig.add_subplot(111,frameon=False)
        errorbar(distance_list,redux[0][1],yerr=redux_std[0],\
            color=(1,0,0),linewidth=linewidth,marker='o',label=redux[0][0])
        errorbar(distance_list,redux[1][1],yerr=redux_std[1],\
            color=(0,0,1),linewidth=linewidth,marker='x',label=redux[1][0])
        biglegend('upper right')#'lower left')
        axes_labels(ax,'separation (microns)','mean frate redux (Hz)') # (0 to 1.5nA) in 400ms
        ax.set_xticks(distance_list)
        ax.get_yaxis().set_ticks_position('left')
        ax.get_xaxis().set_ticks_position('bottom')
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.set_ylim(0,ymax)
        ax.set_yticks([0,ymax])
        ax.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
        ax.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
        fig.tight_layout()
        fig.savefig('../figures/ADI/ADI_asym_invivo.png',bbox_inches='tight',dpi=fig.dpi)
    if ASYM_TEST:
        max_redux_frate_array_rev, max_asym_redux_frate_array_rev =  \
            plot_avg_inh(filelist_invitro_reversed, 'reversed ')
    
        mainfig = figure(facecolor='w')
        mainaxes = mainfig.add_subplot(111)
        for i,redux in enumerate(max_asym_redux_frate_array):
            mainaxes.plot([1,2],[-max_asym_redux_frate_array_rev[i],-redux],
                marker='+',linewidth='2.0')
        mainaxes.set_xticks([0.8,1,2,2.2])
        mainaxes.set_xticklabels(['','A','B',''])
        axes_labels(mainaxes,"mitral cell","reduction in mitral firing (Hz)")

    ## plot the fig 2c of Arevian et al.
    expdata_Boff = [(0,0),(79.429,0),(159.286,5.067),(239.143,9.966),(317.464,17.748),\
    (400.393,19.765),(481.786,30.141),(560.107,44.839),(639.964,59.826),(719.821,82.306),\
    (801.214,97.292),(879.536,109.973),(959.393,116.890),(1039.250,129.571),(1122.179,136.776)]
    expdata_Bon = [(0,0),(79.429,4.779),(160.821,7.661),(239.143,12.272),(322.071,17.460),\
    (398.857,25.241),(480.250,29.853),(560.107,39.940),(639.964,52.332),(719.821,65.013),\
    (801.214,84.899),(879.536,104.497),(959.393,119.196),(1040.786,127.265),(1120.643,139.658)]
    ## remove the two highest points
    expdata_Boff = expdata_Boff[:-2]
    expdata_Bon = expdata_Bon[:-2]
    fig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='none') # transparent
    ax = fig.add_subplot(111,frameon=False)
    expdata_BoffX,expdata_BoffY = zip(*expdata_Boff)
    ax.plot(array(expdata_BoffX)*1e-3,expdata_BoffY,\
        color='k', linewidth=linewidth, marker='s',markersize=marker_size)
    expdata_BonX,expdata_BonY = zip(*expdata_Bon)
    ax.plot(array(expdata_BonX)*1e-3,expdata_BonY,\
        color='r', linewidth=linewidth, marker='o',markersize=marker_size)
    axes_labels(ax,'mitral A inject (nA)','mitral A firing rate (Hz)')
    ax.get_yaxis().set_ticks_position('left')
    ax.get_xaxis().set_ticks_position('bottom')
    ax.set_ylim(0,120)
    ax.set_yticks([0,20,40,60,80,100,120])
    ax.set_xticks([0,0.5,1.2])
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
    ax.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
    fig_clip_off(fig)
    fig.tight_layout()
    fig.savefig('../figures/ADI/ADI_exp.png',dpi=fig.dpi,transparent=True)
    fig.savefig('../figures/ADI/ADI_exp.svg',dpi=fig.dpi,transparent=True)

    show()
