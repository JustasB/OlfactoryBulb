#!/usr/bin/env python
# -*- coding: utf-8 -*-

########## You need to run: python2.6 average_tuft_inhibition.py <tuftADIresults_foldername>

import os
import sys
import math
import pickle

sys.path.extend(["..","../networks","../generators","../simulations"])

from moose.utils import * # imports moose
from data_utils import *

from stimuliConstants import * # has SETTLETIME
from simset_activinhibition import * # has oninject_ext
from plot_mitral_connectivity_NetworkML import * # get connectivity from netfile

from pylab import * # part of matplotlib that depends on numpy but not scipy

## override sim time settings of stimuliConstants.py (has stimuliConstantsMinimal.py)
## Here in tuftADI, unlike activdep inh,
## ASYM_TEST puts same ORN frate into B as A, rather than const.
ASYM_TEST = False#True
DIRECTED = True
FRAC_DIRECTED = 0.01#0.03
IN_VIVO = True
onInject = oninject_ext # in Hz, from simset_activinhibition_minimal
REVERSED_ADI = False
NONLINEAR_ORNS = False
## if ODORINH, sims used 1x granule bgnd (ensure extraexc_factor=1 in stimuliConstants.py)
## and 10Hz ORN to mitB, else 1x and 5Hz.
## I'm using this to just deliver 0.5x and 1x lat i/p, to show linearity of lat inh.
ODORINH = True#False

if IN_VIVO: invivo_str = '_invivo'
else: invivo_str = ''
if DIRECTED: dirt_str = '_directed'+str(FRAC_DIRECTED)
else: dirt_str = ''
if REVERSED_ADI: rev_str = '_reversed'
else: rev_str = ''
if ASYM_TEST: asym_str = '_asym'
else: asym_str = ''
if ODORINH: odorinh_str = '_odorinh'
else: odorinh_str = '_airinh'

## inh =  (no_singles,no_joints,no_PGs)
inh = (False,False,False)
print "Using NOSINGLES =",inh[0],"NOJOINTS =",inh[1],"NOPGS =",inh[2]
netseeds = [100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0]

if IN_VIVO:
    mit_distance_list = [50.0,200.0,400.0,600.0,800.0,1000.0,1200.0,1400.0,1600.0,1800.0,1850.0,1900.0] # microns
else:
    mit_distance_list = [50.0] # microns

## whether to calculate and print the number of joint granules connected
SHOW_CONN = False

if len(sys.argv)<2:
    print 'Specify folder name where tuftADI results reside.'
    sys.exit(1)
else:
    filebase = sys.argv[1]
    ## eg: choose '_aug17_2' from '../results/tuftADI_aug17_2/'
    dirextns = filebase.split('/')
    if dirextns[-1]=='': dirextn = dirextns[-2]
    else: dirextn = dirextns[-1]
    dirextn = dirextn.split('tuftADI')[-1]

def get_filename(netseed,inh,mitdistance,lat_mitnum=3,
        _directed=DIRECTED,_frac_directed=FRAC_DIRECTED,
        _asym_str=asym_str,_rev_str=rev_str,_odorinh_str=odorinh_str,
        directed_str=dirt_str,filebase='../results/tuftADI'):
    ### read filename from the output file of automated run
    ## construct the filename
    if inh[0]: singles_str = '_NOSINGLES'
    else: singles_str = '_SINGLES'
    if inh[1]: joints_str = '_NOJOINTS'
    else: joints_str = '_JOINTS'
    if inh[2]: pgs_str = '_NOPGS'
    else: pgs_str = '_PGS'
    ## stable enough that time tags are not needed
    outfilename = filebase+'/tuftADI_'+str(lat_mitnum)+\
        '_seed'+str(netseed)+'_mitdist'+str(mitdistance)+\
        singles_str+joints_str+pgs_str+invivo_str+directed_str+\
        _rev_str+_asym_str+_odorinh_str+'.pickle'
    return outfilename, (singles_str, joints_str, pgs_str)

def get_rates_inh(filename):
    f = open(filename,'r')
    Ainjectarray, both_firingratearrays = pickle.load(f)
    f.close()
    mit_alone = array(both_firingratearrays[0])
    mit_inhibited = array(both_firingratearrays[1])
    diff_frates = mit_alone - mit_inhibited
    max_redux = max(diff_frates)
    avg_redux = mean(diff_frates)
    return Ainjectarray, both_firingratearrays, max_redux, avg_redux

def plot_avg_inh_lat(mit_distance,_rev_str,manyplots=False,directed_str=dirt_str):
    ## load the data
    full_firingratearray = []
    half_firingratearray = []
    redux_frate_array = []
    max_redux_frate_array = []
    avg_redux_frate_array = []
    for netseed in netseeds:
        filename,_ = get_filename(netseed,inh,mit_distance,
            _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
            directed_str=directed_str,filebase=filebase)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        if not os.path.exists(filename):
            print "File not found",filename
            continue
            
        Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
            get_rates_inh(filename)
        full_firingratearray.append(both_firingratearrays)
        max_redux_frate_array.append(max_redux)
        avg_redux_frate_array.append(avg_redux)

        ## load for 0.5x (_airinh) firing rate in B too.
        if not ASYM_TEST:
            filename,_ = get_filename(netseed,inh,mit_distance,
                _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_airinh',
                directed_str=directed_str,filebase=filebase)
            ## if the result file for these seeds & tweaks doesn't exist,
            ## then carry on to the next.
            if not os.path.exists(filename): continue
            print filename
            Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
                get_rates_inh(filename)
            half_firingratearray.append(both_firingratearrays)

    ## calculate means and std errs
    full_firingratearray = array(full_firingratearray)
    mean_firingratearray = mean(full_firingratearray, axis=0)
    std_firingratearray = std(full_firingratearray, axis=0)
    half_firingratearray = array(half_firingratearray)
    mean_half_firingratearray = mean(half_firingratearray, axis=0)
    std_half_firingratearray = std(half_firingratearray, axis=0)

    if manyplots:
        ## plot the lone vs inhibited firing rate vs input (Hz) curves
        mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
        mainaxes = mainfig.add_subplot(111,frameon=False)
        for mitralBinject in [0,1]:
            if ASYM_TEST:
                if mitralBinject == 0: mitB_str = '0.0 Hz'
                else: mitB_str = 'as for A'
            else: mitB_str = str(mitralBinject*onInject)+' Hz'
            mainaxes.errorbar(x=Ainjectarray,\
            y=mean_firingratearray[mitralBinject],\
            yerr=std_firingratearray[mitralBinject],\
            color=['k','b'][mitralBinject], marker=['s','o'][mitralBinject],\
            markersize=marker_size, label="mitB i/p = "+mitB_str, linewidth=linewidth)
        if not ASYM_TEST and 'reverse' not in _rev_str:
            mitB_str = str(0.5*onInject)+' Hz'
            mainaxes.errorbar(x=Ainjectarray,\
            y=mean_half_firingratearray[1],\
            yerr=std_half_firingratearray[1],\
            color='r', marker='d', markersize=marker_size,\
            label="mitB i/p = "+mitB_str, linewidth=linewidth)    
        mainaxes.get_yaxis().set_ticks_position('left')
        mainaxes.get_xaxis().set_ticks_position('bottom')
        xmin, xmax = mainaxes.get_xaxis().get_view_interval()
        ymin, ymax = mainaxes.get_yaxis().get_view_interval()
        mainaxes.set_xlim([0,xmax])
        mainaxes.set_ylim([0,ymax])
        mainaxes.set_xticks([0,xmax])
        mainaxes.set_yticks([0,ymax])
        mainaxes.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
        mainaxes.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
        ##biglegend("lower right")
        #biglegend("upper left")
        axes_labels(mainaxes,"mitral A ORN rate (Hz)","mitral A firing rate (Hz)")
        mainfig.tight_layout()
        figfilename = '../figures/tuftADI/tuftADI'+'_mitdist'+str(mit_distance) \
            +invivo_str+dirt_str+_rev_str+asym_str+dirextn
        mainfig.savefig(figfilename+'.png', dpi=mainfig.dpi)
        mainfig.savefig(figfilename+'.svg', dpi=mainfig.dpi)

    return max_redux_frate_array, avg_redux_frate_array

def plot_avg_inh_lat_all(mit_distance,_rev_str,manyplots=False):
    mainfig = figure(figsize=(columnwidth*2.0,linfig_height),dpi=300,facecolor='w')
    ax1 = mainfig.add_subplot(1,5,1)
    ax2 = mainfig.add_subplot(1,5,2)
    ax3 = mainfig.add_subplot(1,5,3)
    ax4 = mainfig.add_subplot(1,5,4)
    ax5 = mainfig.add_subplot(1,5,5)
    ## inh =  (no_singles,no_joints,no_PGs)
    inh_options = [
        ('../results/tuftADI/',(False,False,False),ax1),
        ('../results/tuftADI/',(False,False,True),ax2),
        ('../results/tuftADI/',(True,True,False),ax3),
        ('../results/tuftADI_aug17_2_PGmod/',(False,False,False),ax4),
        ('../results/tuftADI_aug17_4_PGmod/',(False,False,False),ax5)]
    maxy = 0.0
    for axnum,(filebase,inh,mainaxes) in enumerate(inh_options):
        ## load the data
        full_firingratearray = []
        half_firingratearray = []
        redux_frate_array = []
        max_redux_frate_array = []
        avg_redux_frate_array = []
        for netseed in netseeds:
            filename,_ = get_filename(netseed,inh,mit_distance,
                _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
                filebase=filebase)
            ## if the result file for these seeds & tweaks doesn't exist,
            ## then carry on to the next.
            if not os.path.exists(filename): continue
            print filename
                
            Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
                get_rates_inh(filename)
            full_firingratearray.append(both_firingratearrays)
            max_redux_frate_array.append(max_redux)
            avg_redux_frate_array.append(avg_redux)

            ## load for 0.5x (_airinh) firing rate in B too.
            if not ASYM_TEST:
                filename,_ = get_filename(netseed,inh,mit_distance,
                    _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_airinh',
                    filebase=filebase)
                ## if the result file for these seeds & tweaks doesn't exist,
                ## then carry on to the next.
                if not os.path.exists(filename): continue
                print filename
                Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
                    get_rates_inh(filename)
                half_firingratearray.append(both_firingratearrays)

        ## calculate means and std errs
        full_firingratearray = array(full_firingratearray)
        mean_firingratearray = mean(full_firingratearray, axis=0)
        std_firingratearray = std(full_firingratearray, axis=0)
        half_firingratearray = array(half_firingratearray)
        mean_half_firingratearray = mean(half_firingratearray, axis=0)
        std_half_firingratearray = std(half_firingratearray, axis=0)

        if manyplots:
            ## plot the lone vs inhibited firing rate vs input (Hz) curves
            for mitralBinject in [0,1]:
                if ASYM_TEST:
                    if mitralBinject == 0: mitB_str = '0.0 Hz'
                    else: mitB_str = 'as for A'
                else: mitB_str = str(mitralBinject*onInject)+' Hz'
                mainaxes.errorbar(x=Ainjectarray,\
                y=mean_firingratearray[mitralBinject],\
                yerr=std_firingratearray[mitralBinject],\
                color=['k','b'][mitralBinject], marker=['s','o'][mitralBinject],\
                markersize=marker_size, label="mitB i/p = "+mitB_str, linewidth=linewidth)
            if not ASYM_TEST and 'reverse' not in _rev_str:
                mitB_str = str(0.5*onInject)+' Hz'
                mainaxes.errorbar(x=Ainjectarray,\
                y=mean_half_firingratearray[1],\
                yerr=std_half_firingratearray[1],\
                color='r', marker='d', markersize=marker_size,\
                label="mitB i/p = "+mitB_str, linewidth=linewidth)    

            xmin,xmax,ymin,ymax = \
                beautify_plot(mainaxes,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
            maxy = max(maxy,ymax)
            if axnum==0:
                axes_labels(mainaxes,"","mitral A firing rate (Hz)")
            elif axnum==2:
                axes_labels(mainaxes,"mitral A ORN rate (Hz)","")
    for ax in [ax1,ax2,ax3,ax4,ax5]:
        ax.set_ylim([0,maxy])
        ax.set_yticks([0,maxy])
    #mainfig.tight_layout()
    mainfig.subplots_adjust(top=0.95,left=0.07,right=0.95,bottom=0.15,wspace=0.25,hspace=0.4)
    figfilename = '../figures/tuftADI/tuftADI_alllat_'+'_mitdist'+str(mit_distance) \
        +invivo_str+dirt_str+_rev_str+asym_str+dirextn
    mainfig.savefig(figfilename+'.png', dpi=mainfig.dpi)
    mainfig.savefig(figfilename+'.svg', dpi=mainfig.dpi)

    return

def plot_avg_inh_self(mit_distance,_rev_str,manyplots=False):
    ## plot the lone vs self-inhibited firing rate vs input (Hz) curves
    mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
    mainaxes = mainfig.add_subplot(111)
    ## inh =  (no_singles,no_joints,no_PGs)
    inh_list = [(True,True,True),(False,False,True),(True,True,False),(False,False,False)]
    labels_list = ['no PGs+granules','only granules','only PGs','both']
    for inh_i,inh in enumerate(inh_list):
        full_firingratearray = []
        dashes = [(1.5,0.5),(0.5,0.5),(1.5,0.5,0.5,0.5),(1.0,0.0)]
        for netseed in netseeds:
            filename,_ = get_filename(netseed,inh,mit_distance,
                _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
                filebase=filebase)
            ## if the result file for these seeds & tweaks doesn't exist,
            ## then carry on to the next.
            if not os.path.exists(filename): continue
                
            ## load the data
            Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
                get_rates_inh(filename)
            ## not interested in lateral inhibition here, so leave that out
            full_firingratearray.append(both_firingratearrays[0])

        ## calculate means and std errs
        full_firingratearray = array(full_firingratearray)
        mean_firingratearray = mean(full_firingratearray, axis=0)
        std_firingratearray = std(full_firingratearray, axis=0)

        mainaxes.errorbar(x=Ainjectarray,\
            y=mean_firingratearray,\
            yerr=std_firingratearray,\
            color=['k','b','c','r'][inh_i], marker=['s','^','o','d'][inh_i],\
            markersize=marker_size, linewidth=linewidth,\
            label=labels_list[inh_i]) # removed dashes=dashes[inh_i]

    beautify_plot(mainaxes,x0min=True,y0min=True,\
        xticksposn='bottom',yticksposn='left',drawxaxis=True,drawyaxis=True)
    ### Legend too big
    #biglegend('lower right',mainaxes,fontsize=label_fontsize-2,\
    #    labelspacing=0.1,borderpad=0.1,markerscale=0.1,columnspacing=0.1,frameon=False)
    ## Legend without the error bars
    ## get handles
    handles, labels = mainaxes.get_legend_handles_labels()
    ## remove the errorbars
    handles = [h[0] for h in handles]
    ## use them in the legend
    leg = mainaxes.legend(handles, labels, loc='upper left',numpoints=1,\
        labelspacing=0.0,borderpad=0.01,markerscale=1,columnspacing=0.0,\
        handletextpad=0.3,prop={'size':label_fontsize-2},frameon=False)
    ### modify legend text sizes -- use prop={'size':...} above,
    ### else below method causes alignment issues of handles with labels
    #for t in leg.get_texts():
    #    t.set_fontsize(label_fontsize-2)
    axes_labels(mainaxes,"mitral A ORN rate (Hz)","mitral A firing rate (Hz)",ypad=-6)
    mainfig.tight_layout()
    figfilename = '../figures/tuftADI/tuftADI'+'_mitdist'+str(mit_distance) \
        +invivo_str+dirt_str+_rev_str+asym_str
    mainfig.savefig(figfilename+'_selfinh.png', dpi=mainfig.dpi)
    mainfig.savefig(figfilename+'_selfinh.svg', dpi=mainfig.dpi)

def plot_avg_diff_inh(mit_distance,_rev_str):
    ## load the data
    full_firingratearray = []
    noinh_firingratearray = []
    noself_firingratearray = []
    redux_frate_array = []
    max_redux_frate_array = []
    avg_redux_frate_array = []
    for netseed in netseeds:
        filename,_ = get_filename(netseed,inh,mit_distance,
            _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
            filebase=filebase)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        if not os.path.exists(filename): continue
        print filename
            
        Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
            get_rates_inh(filename)
        full_firingratearray.append(both_firingratearrays)
        max_redux_frate_array.append(max_redux)
        avg_redux_frate_array.append(avg_redux)

        filename,_ = get_filename(netseed,(True,False,True),mit_distance,
            _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
            filebase=filebase)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        if not os.path.exists(filename): continue
        print filename
            
        Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
            get_rates_inh(filename)
        noself_firingratearray.append(both_firingratearrays)

        filename,_ = get_filename(netseed,(True,True,True),mit_distance,
            _asym_str=asym_str,_rev_str=_rev_str,_odorinh_str='_odorinh',
            filebase=filebase)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        if not os.path.exists(filename): continue
        print filename
            
        Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
            get_rates_inh(filename)
        noinh_firingratearray.append(both_firingratearrays)

    ## calculate means and std errs
    full_firingratearray = array(full_firingratearray)
    mean_firingratearray = mean(full_firingratearray, axis=0)
    std_firingratearray = std(full_firingratearray, axis=0)
    noinh_firingratearray = array(noinh_firingratearray)
    mean_noinh_firingratearray = mean(noinh_firingratearray, axis=0)
    std_noinh_firingratearray = std(noinh_firingratearray, axis=0)
    noself_firingratearray = array(noself_firingratearray)
    mean_noself_firingratearray = mean(noself_firingratearray, axis=0)
    std_noself_firingratearray = std(noself_firingratearray, axis=0)

    ############## plotting
    ## plot the lone vs inh firing rate vs I curves
    mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
    mainaxes = mainfig.add_subplot(111,frameon=False)
    ## plot mit cell alone
    mainaxes.errorbar(x=Ainjectarray,\
        y=mean_noinh_firingratearray[0], yerr=std_noinh_firingratearray[0],\
        color='k', marker='s',\
        markersize=marker_size, linewidth=linewidth)
    ## plot mit cell without self inhibition
    mainaxes.errorbar(x=Ainjectarray,\
        y=mean_noself_firingratearray[0], yerr=std_noself_firingratearray[0],\
        color='r', marker='o',\
        markersize=marker_size, linewidth=linewidth)
    ## plot mit cell without self inhibition + lateral odorinh
    mainaxes.errorbar(x=Ainjectarray,\
        y=mean_noself_firingratearray[1], yerr=std_noself_firingratearray[1],\
        color='r', marker='d',\
        markersize=marker_size, linewidth=linewidth)
    for mitralBinject in [0,1]:
        if ASYM_TEST:
            if mitralBinject == 0: mitB_str = '0.0 Hz'
            else: mitB_str = 'as for A'
        else: mitB_str = str(mitralBinject*onInject)+' Hz'
        mainaxes.errorbar(x=Ainjectarray,\
            y=mean_firingratearray[mitralBinject],\
            yerr=std_firingratearray[mitralBinject],\
            color=['b','b'][mitralBinject], marker=['x','+'][mitralBinject],\
            markersize=marker_size, label="mitB i/p = "+mitB_str, linewidth=linewidth)
    mainaxes.get_yaxis().set_ticks_position('left')
    mainaxes.get_xaxis().set_ticks_position('bottom')
    xmin, xmax = mainaxes.get_xaxis().get_view_interval()
    ymin, ymax = mainaxes.get_yaxis().get_view_interval()
    mainaxes.set_xlim([0,xmax])
    mainaxes.set_ylim([0,ymax])
    mainaxes.set_xticks([0,xmax])
    mainaxes.set_yticks([0,ymax])
    mainaxes.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
    mainaxes.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
    ##biglegend("lower right")
    #biglegend("upper left")
    axes_labels(mainaxes,"mitral A ORN frate (Hz)","mitral A firing rate (Hz)")
    mainfig.tight_layout()
    figfilename = '../figures/tuftADI/tuftADI'+'_mitdist'+str(mit_distance) \
        +invivo_str+dirt_str+_rev_str+asym_str
    #mainfig.savefig(figfilename+'.png', dpi=mainfig.dpi)
    #mainfig.savefig(figfilename+'.svg', dpi=mainfig.dpi)

    return max_redux_frate_array, avg_redux_frate_array

def plot_avg_inh_Bconst_BasA(mit_distance,_rev_str=rev_str,asym_list=['_asym','']):
    ## plot the firing rate vs I curves for lone, BasA, Bconst
    mainfig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
    mainaxes = mainfig.add_subplot(111,frameon=False)

    curvei = 0
    max_redux_frate_array = []
    avg_redux_frate_array = []
    ## asym (True/False) specifies whether mitB gets sameasA / const tuft input.
    for _asym_str in asym_list:
        ## load the data
        full_firingratearray = []
        for netseed in netseeds:
            filename,_ = get_filename(netseed,inh,mit_distance,\
                _asym_str=_asym_str,_rev_str=_rev_str,filebase=filebase)
            ## if the result file for these seeds & tweaks doesn't exist,
            ## then carry on to the next.
            if not os.path.exists(filename): continue
            print 'Using',filename

            Ainjectarray, both_firingratearrays, max_redux, avg_redux = \
                get_rates_inh(filename)
            full_firingratearray.append(both_firingratearrays)
            ## use only those frate reductions that are with same input as A in mitB.
            if not _asym_str:
                max_redux_frate_array.append(max_redux)
                avg_redux_frate_array.append(avg_redux)

        ## calculate means and std errs
        full_firingratearray = array(full_firingratearray)
        mean_firingratearray = mean(full_firingratearray, axis=0)
        std_firingratearray = std(full_firingratearray, axis=0)

        for mitralBinject in [0,1]:
            if _asym_str: ## if '', then const i/p to mitB
                mitB_str = str(mitralBinject*onInject)+' Hz'
            else:
                if mitralBinject == 0: mitB_str = '0.0 Hz'
                else: mitB_str = 'as for A'
            ## condition to not plot the lone mitA f-vs-I twice
            #if mitralBinject!=0 or _asym_str:
            if True:
                mainaxes.errorbar(x=Ainjectarray,\
                y=mean_firingratearray[mitralBinject],\
                yerr=std_firingratearray[mitralBinject],\
                color=['k','b','r'][curvei], marker=['s','o','v'][curvei],\
                markersize=marker_size, label="mitB i/p = "+mitB_str, linewidth=linewidth)
                curvei += 1

    mainaxes.get_yaxis().set_ticks_position('left')
    mainaxes.get_xaxis().set_ticks_position('bottom')
    xmin, xmax = mainaxes.get_xaxis().get_view_interval()
    ymin, ymax = mainaxes.get_yaxis().get_view_interval()
    mainaxes.set_xlim([0,xmax])
    mainaxes.set_ylim([0,ymax])
    mainaxes.set_xticks([0,xmax])
    mainaxes.set_yticks([0,ymax])
    mainaxes.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
    mainaxes.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
    ##biglegend("lower right")
    #biglegend("upper left")
    axes_labels(mainaxes,"mitral A ORN frate (Hz)","mitral A firing rate (Hz)")
    mainfig.tight_layout()
    mainfig.savefig('../figures/tuftADI/tuftADI_both' \
        +invivo_str+dirt_str+rev_str+odorinh_str+'.png', \
        dpi=mainfig.dpi)
    mainfig.savefig('../figures/tuftADI/tuftADI_both' \
        +invivo_str+dirt_str+rev_str+odorinh_str+'.svg', \
        dpi=mainfig.dpi)
    return max_redux_frate_array, avg_redux_frate_array

def plot_avg_inh_distancedep(mit_distance_list):
    redux = []
    redux_std = []
    for i,reverse in enumerate([False,True]):
        print "Reversed =", reverse
        if reverse:
            rev_label = "Reversed"
            _rev_str = '_reversed'
        else:
            rev_label = "Forward"
            _rev_str = ''
        redux_directed = []
        redux_directed_std = []
        for mit_distance in mit_distance_list:
            print 'Distance =',mit_distance,'microns:'
            ## the former below, plots for only one of asym=True/False
            ## (i.e. sameasA/const i/p to mitB) set at the top.
            ## whereas the latter below, plots for both in the same graph.
            max_redux_frate_array, avg_redux_frate_array = \
                plot_avg_inh_lat(mit_distance,_rev_str,manyplots=False)
            #max_redux_frate_array, avg_redux_frate_array = \
            #    plot_avg_inh_Bconst_BasA(mit_distance,_rev_str)
            redux_directed.append(mean(avg_redux_frate_array))
            redux_directed_std.append(std(avg_redux_frate_array))
        redux.append((rev_label,redux_directed))
        redux_std.append(redux_directed_std)

    ## plot the reduction-in-frate vs separation-of-mitrals for forward & reverse tuft inh.
    fig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
    ax = fig.add_subplot(111)
    mit_distance_list = array(mit_distance_list)/1000.0
    errorbar(mit_distance_list,-array(redux[0][1]),yerr=redux_std[0],\
        color=(1,0,0),linewidth=linewidth,marker='o',label=redux[0][0])
    errorbar(mit_distance_list,-array(redux[1][1]),yerr=redux_std[1],\
        color=(0,0,1),linewidth=linewidth,marker='x',label=redux[1][0])
    #biglegend('upper right')#'lower left')
    axes_labels(ax,'separation (mm)','mean $\Delta$ rate (Hz)',ypad=-3) # (0 to 1.5nA) in 400ms
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
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.set_xlim(0,xmax)
    ax.set_ylim(ymin,0)
    ax.set_yticks([ymin,0])
    ax.set_xticks([0,1,2])
    fig_clip_off(fig)
    fig.tight_layout()
    fig.savefig('../figures/tuftADI/tuftADI_redux'+invivo_str+asym_str+'.svg',\
        bbox_inches='tight',dpi=fig.dpi)
    fig.savefig('../figures/tuftADI/tuftADI_redux'+invivo_str+asym_str+'.png',\
        bbox_inches='tight',dpi=fig.dpi)

def plot_avg_inh_distancedep_v2(mit_distance_list):
    redux = []
    redux_std = []
    mit_distance_list = array(mit_distance_list)
    for i,reverse in enumerate([False,True]):
        ## plot the reduction-in-frate vs separation-of-mitrals for forward & reverse tuft inh.
        fig = figure(figsize=(columnwidth/2.0,linfig_height),dpi=300,facecolor='w')
        ax = fig.add_subplot(111)
        print "Reversed =", reverse
        if reverse:
            rev_label = "Reversed"
            _rev_str = '_reversed'
        else:
            rev_label = "Forward"
            _rev_str = ''
        for dirti,directed_str in enumerate(['','_directed0.0','_directed0.01']):
            redux_directed = []
            redux_directed_std = []
            print "Directed str set as",directed_str,"."
            for mit_distance in mit_distance_list:
                print 'Distance =',mit_distance,'microns:'
                ## the former below, plots for only one of asym=True/False
                ## (i.e. sameasA/const i/p to mitB) set at the top.
                ## whereas the latter below, plots for both in the same graph.
                max_redux_frate_array, avg_redux_frate_array = \
                    plot_avg_inh_lat(mit_distance,_rev_str,manyplots=False,\
                    directed_str=directed_str)
                #max_redux_frate_array, avg_redux_frate_array = \
                #    plot_avg_inh_Bconst_BasA(mit_distance,_rev_str)
                redux_directed.append(mean(avg_redux_frate_array))
                redux_directed_std.append(std(avg_redux_frate_array))
            redux.append((rev_label,redux_directed))
            redux_std.append(redux_directed_std)

            errorbar(mit_distance_list/1000.,-array(redux_directed),yerr=redux_directed_std,\
                color=['g','b','r'][dirti],linewidth=linewidth,\
                marker=['o','x','s'][dirti],label=redux[0][0])
        #biglegend('upper right')#'lower left')
        axes_labels(ax,'separation (mm)','mean $\Delta$ rate (Hz)',ypad=-3) # (0 to 1.5nA) in 400ms
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
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.set_xlim(0,xmax)
        ymin = -35 ## hard-coded for my simulations
        ax.set_ylim(ymin,1)
        ax.set_yticks([ymin,0,2])
        ax.set_xticks([0,1,2])
        fig_clip_off(fig)
        fig.tight_layout()
        fig.savefig('../figures/tuftADI/tuftADI_redux'+invivo_str+asym_str+rev_label+'.svg',\
            bbox_inches='tight',dpi=fig.dpi)
        fig.savefig('../figures/tuftADI/tuftADI_redux'+invivo_str+asym_str+rev_label+'.png',\
            bbox_inches='tight',dpi=fig.dpi)

if __name__ == "__main__":
    ### for single / few test runs...
    #plot_avg_inh_Bconst_BasA(600.0,_rev_str='',asym_list=[''])
    #plot_avg_diff_inh(400.0,_rev_str='')
    ## PAPER figure 4: example tuft inh -- lateral
    #plot_avg_inh_lat(400.0,_rev_str='',manyplots=True)
    ## PAPER figure 5 (diff pg inputs) : tuft inh -- lateral
    #plot_avg_inh_lat_all(400.0,_rev_str='',manyplots=True)
    ## For a single plot of lat inh -- air and odor
    ## set the directory on commandline and the inh_options in the function:
    #plot_avg_inh_lat(400.0,_rev_str='',manyplots=True)
    ## PAPER figure 2: tuft inh -- self
    plot_avg_inh_self(400.0,_rev_str='',manyplots=True)
    ## OBSOLETE PAPER figure 4: tuft inh -- distance dependence
    #plot_avg_inh_distancedep(mit_distance_list)
    ## PAPER figure 4: tuft inh -- distance dependence
    #plot_avg_inh_distancedep_v2(mit_distance_list)

    show()
