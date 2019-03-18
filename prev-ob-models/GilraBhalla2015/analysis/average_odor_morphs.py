# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: python2.6 average_odor_morphs.py [morphs results directory]

import os,sys,math,string
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has RESPIRATION
from pylab import * # part of matplotlib that depends on numpy but not scipy
from scipy import stats
from sim_utils import * # has rebin() and imports data_utils.py for axes_off()
from calc_corrs import * # has calc_corrs()
## Choose whether to have linear weights and rectifier: FULLlin = True
## or just linear weights and sigmoidal non-linearity: FULLlin = False
FULLlin = True
if FULLlin:
    import fit_odor_morphs_withFULLlin as fit_om # has fit_morphs()
    linextn = 'FULLlin'
else:
    import fit_odor_morphs as fit_om # has fit_morphs()
    linextn = 'lin'
## I only need residual2noise, but average_odor_pulses imports average_odor_morphs,
## hence import <...> as <...> to avoid conflict.
import average_odor_pulses as forR2N # has function to calculate residual to noise

IN_VIVO = True
## these below are often set in the respective function arguments in __main__
DIRECTED = True ## overrides the one in networkConstants.py (came in via sim_utils.py)
FRAC_DIRECTED = 0.01 ## overrides the one in networkConstants.py (came in via sim_utils.py)
## Below two overide the variables that came in from stimuliConstantsMinimal.py via stimuliConstants.py
NONLINEAR_ORNS = False
NONLINEAR_TYPE = 'P' # P for primary glom non-linear, L for lateral gloms non-linear

fullfilename = '../results/odor_morphs/morphs_random'
if NONLINEAR_ORNS: fullfilename += 'NL'+NONLINEAR_TYPE
fullfilename += '.pickle'
#fullfile = open(fullfilename,'r')
#morphs = pickle.load(fullfile)
#fullfile.close()

PLOTRESP_NUM = 1 # whether to plot 2 respiration cycles or 1
NUMBINS = 5 ## overriding that in calc_corrs.py
BIN_WIDTH_TIME = RESPIRATION/NUMBINS
bindt = RESPIRATION/float(NUMBINS)
## I take the last PLOTRESP_NUM of respiration cycles out of NUM_RESPS simulated
responsetlist = arange( SETTLETIME+(NUM_RESPS-PLOTRESP_NUM)*RESPIRATION+bindt/2.0, \
    SETTLETIME+NUM_RESPS*RESPIRATION, RESPIRATION/NUMBINS )

## IMPORTANT!!!!!! num_gloms defined below is only for some functions.
## Typically used plot_directed() has num_gloms defined in __main__ below.
salient = False#True
if salient:
    stim_seeds = [-1,-10,-19,-28]#range(-1,-37,-1)#[-1,-2,-3,-4,-5,-6,-7,-8,-9]
    num_gloms_list = [3] # used only by some functions, see __main__ for others
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order, below options are:
    ## all cells; no lat; no joints, varyRMP; no PGs; no singles + no joints, only mitrals
    #inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
    #    (2,(False,True,False,False,False)), (3,(False,False,False,False,True)), \
    #    (4,(False,False,False,True,False)), (5,(True,True,False,False,False)), (6,(True,True,False,True,False))]
    inh_options = [ (0,(False,False,False,False,False)), (0,(False,False,True,False,False)) ]
else:
    #stim_seeds = append(stim_seeds,[-1,-10,-19,-28])
    
    ##################################### IMPORTANT ###################################
    ###### USE ONE OR THE OTHER FOR FIGURES
    ## I used all simulations for decorrelation FIGURES
    stim_seeds = arange(750.0,1100.0,1.0)
    ## I used only 50 seeds for fittng Adil morphs.
    #stim_seeds = arange(750.0,800.0,1.0)
    
    num_gloms_list = [3] # used only by some functions, see __main__ for others
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order,below options are: all cells; no lat; no joints, varyRMP; no PGs; only mitrals
    #inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
    #    (2,(False,True,False,False,False)), (3,(False,False,False,False,True)), \
    #    (6,(False,False,False,True,False)), (8,(True,True,False,True,False))]
    ###### THESE inh_options ARE USED BY ONLY SOME FUNCTIONS, SEE __main__ FOR OTHERS.
    inh_options = [ (0,(False,False,False,False,False))]#, (1,(False,False,True,False,False))]
        #(2,(False,False,False,False,True)) ]
net_seeds = [100.0,200.0]

def get_filename(netseed,stimseed,inh,numgloms,stimi,neti,inhi,
        _directed=DIRECTED,_frac_directed=FRAC_DIRECTED,\
        resultsdir='../results/odor_morphs'):
    ### read filename from the output file of automated run
    #filename = morphs[ngi,stimi,neti,inhi]
    ## construct the filename
    if inh[0]: singles_str = '_NOSINGLES'
    else: singles_str = '_SINGLES'
    if inh[1]: joints_str = '_NOJOINTS'
    else: joints_str = '_JOINTS'
    if inh[3]: pgs_str = '_NOPGS'
    else: pgs_str = '_PGS'
    if inh[2]: lat_str = '_NOLAT'
    else: lat_str = '_LAT'
    if inh[4]: varmitstr = '_VARMIT'
    else: varmitstr = '_NOVARMIT'
    ## stable enough that time tags are not needed
    filename = resultsdir+'/odormorph'+\
        '_netseed'+str(netseed)+'_stimseed'+str(stimseed)
    if NONLINEAR_ORNS: filename += '_NL'+NONLINEAR_TYPE
    filename += singles_str+joints_str+pgs_str+lat_str+varmitstr+\
        '_numgloms'+str(numgloms)
    if _directed: filename += '_directed'+str(_frac_directed)
    filename += '.pickle'
    return filename, (singles_str, joints_str, pgs_str, lat_str, varmitstr)

def read_morphfile_allmits(filename):
    f = open(filename,'r')
    (mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
    f.close()
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, numbins=PLOTRESP_NUM*NUMBINS,\
            bin_width_time=BIN_WIDTH_TIME, numresps=PLOTRESP_NUM)
    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    return numavgs, mitral_responses_avg, mitral_responses_std

def plot_responses():
    for stimi,stimseed in enumerate(stim_seeds):
        if not salient: net_seeds_here = [stimseed]
        else: net_seeds_here = net_seeds
        numsubfigs_rows = len(net_seeds_here)*len(num_gloms_list)
        numsubfigs_cols = len(inh_options)
        fig = figure(facecolor='w')
        ax = fig.add_subplot(111)
        figtext(0.1,0.94,'A,B,air: mit0 vs mit1 : stim ='+str(stimseed),fontsize=20)
        delaxes()

        for neti,netseed in enumerate(net_seeds_here):
            for ngi,num_gloms in enumerate(num_gloms_list):
                ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
                for plotyi,(inhi,inh) in enumerate(inh_options):

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,num_gloms,stimi,neti,inhi)
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    switches_str = string.join(switch_strs,'')
                    print filename

                    ## read in the spiketimes/binned responses for this run
                    numavgs, mitral_responses_avg, mitral_responses_std = \
                        read_morphfile_allmits(filename)

                    ## calc phase corr-s for printing on the title
                    (air_corr,odorA_corr,odorB_corr),\
                        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
                        Dfrates = \
                        calc_corrs(filename, norm_str="overall", \
                            numbins=NUMBINS, bin_width_time=BIN_WIDTH_TIME, printinfo=False)

                    ## plot the binned firing rates highlighting VARMIT differences
                    ## add_subplot(rows,cols,fignum)
                    ## fignum fills first row to end, then next row
                    ax = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                        (neti+ngi)*numsubfigs_cols+plotyi+1)
                    ## air
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[6,0],\
                        yerr=mitral_responses_std[6,0],color=(0,0,0),linewidth=2)
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[6,1],\
                        yerr=mitral_responses_std[6,1],color=(0,0,0.5),linewidth=2)                    
                    ## odorA
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,0],\
                        yerr=mitral_responses_std[5,0],color=(1,0,0),linewidth=2)
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,1],\
                        yerr=mitral_responses_std[5,1],color=(1,0,0.5),linewidth=2)
                    ## odorB
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[0,0],\
                        yerr=mitral_responses_std[0,0],color=(0,1,0),linewidth=2)
                    ax.errorbar(x=responsetlist,y=mitral_responses_avg[0,1],\
                        yerr=mitral_responses_std[0,1],color=(0,1,0.5),linewidth=2)
                    ax.set_title(
                        str(netseed)+str(num_gloms)+'G'+'_%1.1f_%1.1f_%1.1f'\
                        %(air_corr,odorA_corr,odorB_corr)+\
                        switches_str,fontsize=8)

def plot_decorr_single(resultsdir,stimseed=stim_seeds[0],\
        numgloms=num_gloms_list[0],inh=inh_options[0][1],odornum=None,splax=None):
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    figtext(0.12,0.90,'sisters: mit0 (r) vs mit1 (g)',fontsize=34)
    netseed = stimseed
    if stimseed<0:
        stimseed = int(stimseed)
        netseed = net_seeds[0]

    filename, switch_strs \
        = get_filename(netseed,stimseed,inh,numgloms,\
            None,None,None,resultsdir=resultsdir)
    switches_str = string.join(switch_strs,'')
    print filename

    ## read in the spiketimes/binned responses for this run
    numavgs, mitral_responses_avg, mitral_responses_std = \
        read_morphfile_allmits(filename)
    ## since I plot the mean response, I must plot standard error of the mean
    ## = standard deviation of a repeat / sqrt(num of repeats).
    ## NOTE: our plot depends on number of trials.
    mitral_responses_se = mitral_responses_std/sqrt(numavgs)

    ## calc phase corr-s for printing on the title
    (air_corr,odorA_corr,odorB_corr),\
        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
        Dfrates = \
        calc_corrs(filename, norm_str="overall", \
            numbins=NUMBINS, bin_width_time=BIN_WIDTH_TIME, printinfo=False)

    odor_corr = odorA_corr
    ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,0],\
        yerr=mitral_responses_std[5,0],color=(1,0,0),linewidth=4)
    ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,1],\
        yerr=mitral_responses_std[5,1],color=(0,1,0),linewidth=4)
    ax.set_title('corr = %1.2f'%(odor_corr,),fontsize=30)
    #ax.set_xlim(0.75,1.25)
    #ax.set_xticks([0.75,1.25])
    ax.set_xticks([0.8,0.9,1.0,1.1,1.2])
    ax.set_xticklabels(['1','2','3','4','5'])
    ## just to scale up the ticks fontsize.
    axes_labels(ax,'phase bin','firing rate (Hz)',adjustpos=False,fontsize=34)
    
    if odornum is not None and splax is not None:
        ## mit0 (red) vs mit1 (blue)
        simresponse = mitral_responses_avg[odornum,0]
        simerr = mitral_responses_se[odornum,0]
        plottimes = arange(BIN_WIDTH_TIME/2.0,RESPIRATION,BIN_WIDTH_TIME)
        #splax.fill_between(plottimes, simresponse+simerr, simresponse-simerr,
        #    color='r',alpha=0.4)
        #splax.plot(plottimes,simresponse,color='r',linewidth=linewidth)
        splax.errorbar(x=plottimes,y=simresponse,color='r',\
            marker='s',ms=marker_size,linewidth=linewidth)
        simresponse = mitral_responses_avg[odornum,1]
        simerr = mitral_responses_se[odornum,1]
        #splax.fill_between(plottimes, simresponse+simerr, simresponse-simerr,
        #    color='b',alpha=0.4)
        #splax.plot(plottimes,simresponse,color='b',linewidth=linewidth)
        splax.errorbar(x=plottimes,y=simresponse,color='b',\
            marker='o',ms=marker_size,linewidth=linewidth,linestyle='dashed')
        splax.set_xlim(0,RESPIRATION)
        beautify_plot(splax,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')
        axes_labels(splax,'time (s)','firing rate (Hz)',adjustpos=False,fontsize=label_fontsize)
        #splax.text(-0.3,0.9,'c',fontweight='bold',fontsize=label_fontsize,transform=splax.transAxes)

def plot_responses_mits_paperfigure(resultsdir,odornum,stimseed=stim_seeds[0],\
        numgloms=num_gloms_list[0],inh=inh_options[0][1]):
    netseed = stimseed
    filename, switch_strs \
        = get_filename(netseed,stimseed,inh,numgloms,\
            None,None,None,resultsdir=resultsdir)
    ## read in the spiketimes/binned responses for this run
    numavgs, mitral_responses_avg, mitral_responses_std = \
        read_morphfile_allmits(filename)
    ## plot separate figures for responses of mitral cells
    plottimes = arange(BIN_WIDTH_TIME/2.0,RESPIRATION,BIN_WIDTH_TIME)
    for miti,mitnum in enumerate([0,1,2,4]):
        fig = figure(figsize=(columnwidth/5.0,linfig_height/4.0),\
            dpi=300,facecolor='none') # none is transparent
        ax = fig.add_subplot(111)
        simresponse = mitral_responses_avg[odornum,mitnum]
        ax.errorbar(x=plottimes,y=simresponse,color=['r','m','g','b'][miti],\
            marker='o',ms=marker_size,linewidth=linewidth,\
            linestyle=['solid','dashed','solid','dashed'][miti])
        if mitnum==4:
            add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=True,\
                sizex=0.2,labelx='0.2 s',sizey=8,labely='8 Hz',\
                bbox_to_anchor=[0.7,-0.3],bbox_transform=ax.transAxes)
        beautify_plot(ax,x0min=False,y0min=True,\
            drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
        ax.set_ylim(0,40) # same ymax to set common scale bar
        fig.tight_layout()
        fig_clip_off(fig)
        fig.savefig('../figures/decorr/response_mit'+str(mitnum)+'.svg',\
            dpi=fig.dpi,transparent=True)
        
def plot_decorrs_special_paperfigure(resultsdir,numglomslist,_inh_options=inh_options,
        _directed=DIRECTED,_frac_directed=FRAC_DIRECTED,graph=True):
    if graph:
        fig = figure(figsize=(8, 3), dpi=150, facecolor='w')
        ax = fig.add_subplot(111)
        #figtext(0.1,0.94,'sisters: mit0 (r) vs mit1 (g): frate(Hz) vs phase bin',fontsize=34)
        delaxes()
        figair = figure(facecolor='w')
        axair = figair.add_subplot(111)
        #figtext(0.1,0.94,'Air: mit0 (r) vs mit1 (g): frate(Hz) vs phase bin',fontsize=34)
        delaxes()
    all_odor_corrs = [] # all odor corrs including nan-s
    odor_corrs = [] # not nan odor corrs
    air_corrs = []
    good_odor_corrs = []
    good_air_corrs = []
    total_frate = array([0]*NUMBINS)
    num_responses = 0
    total_airfrate = array([0]*NUMBINS)
    deltafrate_sis0 = []
    deltafrate_sis1 = []
    for stimi,stimseed in enumerate(stim_seeds):
        netseed = stimseed
        if stimseed<0:
            stimseed = int(stimseed)
            netseed = net_seeds[0]
        for ngi,num_gloms in enumerate(numglomslist):
            ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
            for plotyi,(inhi,inh) in enumerate(_inh_options):

                filename, switch_strs \
                    = get_filename(netseed,stimseed,inh,num_gloms,None,None,None,\
                        _directed,_frac_directed,resultsdir)
                ## if the result file for these seeds & tweaks doesn't exist,
                ## then carry on to the next.
                if not os.path.exists(filename): continue
                switches_str = string.join(switch_strs,'')
                #if graph: print filename

                ## read in the spiketimes/binned responses for this run
                f = open(filename,'r')
                (mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
                f.close()
                ## rebin the spikes
                mitral_responses_binned_list = \
                    rebin(mitral_responses_list, numbins=PLOTRESP_NUM*NUMBINS,\
                        bin_width_time=BIN_WIDTH_TIME, numresps=PLOTRESP_NUM)
                numavgs = len(mitral_responses_list)
                mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
                mitral_responses_std = std(mitral_responses_binned_list, axis=0)
                ## mean 1% odor firing rate
                total_frate += mitral_responses_avg[0,0]+\
                    mitral_responses_avg[0,1]+\
                    mitral_responses_avg[5,0]+\
                    mitral_responses_avg[5,1]
                num_responses += 4
                total_airfrate += mitral_responses_avg[6,0] + mitral_responses_avg[6,1]
                #if graph:
                #    print "Odors:"
                #    print mitral_responses_avg[0,0]
                #    print mitral_responses_avg[0,1]
                #    print mitral_responses_avg[5,0]
                #    print mitral_responses_avg[5,1]
                #    print "Airs:"
                #    print mitral_responses_avg[6,0]
                #    print mitral_responses_avg[6,1]
                
                ## Calc change in firing rate for each odor, for each sister
                numphasebins = float(len(mitral_responses_avg[0,0]))
                for odornum in [0,5]:
                    deltafrate0 = sum(mitral_responses_avg[odornum,0])/numphasebins \
                        - sum(mitral_responses_avg[6,0])/numphasebins
                    deltafrate1 = sum(mitral_responses_avg[odornum,1])/numphasebins \
                        - sum(mitral_responses_avg[6,1])/numphasebins
                    deltafrate_sis0.append(deltafrate0)
                    deltafrate_sis1.append(deltafrate1)

                ## calc phase corr-s for printing on the title
                (air_corr,odorA_corr,odorB_corr),\
                    (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
                    Dfrates = \
                    calc_corrs(filename, norm_str="overall", \
                        numbins=NUMBINS, bin_width_time=BIN_WIDTH_TIME, printinfo=False)
                all_odor_corrs.append(odorA_corr)
                all_odor_corrs.append(odorB_corr)
                ## If one of the phasic responses is all zeroes, stats.pearsonr gives a one-time warning
                ## and returns nan-s, but nan-s are removed by below, for obtaining the distribution
                if not isnan(odorA_corr): odor_corrs.append(odorA_corr)
                if not isnan(odorB_corr): odor_corrs.append(odorB_corr)
                if not isnan(air_corr): air_corrs.append(air_corr)
                print stimseed, odorA_corr, odorB_corr, air_corr

                corr_cutoff = -0.5
                numrows = 4
                numcols = 6
                numplots = numrows*numcols
                ## nan-s always return false on comparison with non-nan, hence ignored
                ### fignum fills first row to end, then next row
                ## odorA
                if odorA_corr<corr_cutoff:
                    good_odor_corrs.append(odorA_corr)
                    if graph and len(good_odor_corrs)<=numplots:
                        ax = fig.add_subplot(numrows,numcols,len(good_odor_corrs))
                        ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,0],\
                            yerr=mitral_responses_std[5,0],color=(1,0,0),linewidth=4)
                        ax.errorbar(x=responsetlist,y=mitral_responses_avg[5,1],\
                            yerr=mitral_responses_std[5,1],color=(0,1,0),linewidth=4)
                        ax.set_title('%1.2f'%(good_odor_corrs[-1]),fontsize=16)
                        ymax = ax.get_ylim()[1]
                        ax.set_yticks([0,ymax])
                        ax.set_yticklabels(['0',str(ymax)])
                        #ax.set_xlim(0.75,1.25)
                        #ax.set_xticks([0.75,1.25])
                        ax.set_xticks([])
                        #ax.set_xticklabels(['0.75','1.25'])
                        ## just to scale up the ticks fontsize.
                        axes_labels(ax,'','',adjustpos=False,fontsize=18)
                        #print filename, "corr =", good_odor_corrs[-1]
                ## odorB
                if odorB_corr<corr_cutoff:
                    good_odor_corrs.append(odorB_corr)
                    if graph and len(good_odor_corrs)<=numplots:
                        ax = fig.add_subplot(numrows,numcols,len(good_odor_corrs))
                        ax.errorbar(x=responsetlist,y=mitral_responses_avg[0,0],\
                            yerr=mitral_responses_std[0,0],color=(1,0,0),linewidth=4)
                        ax.errorbar(x=responsetlist,y=mitral_responses_avg[0,1],\
                            yerr=mitral_responses_std[0,1],color=(0,1,0),linewidth=4)
                        ax.set_title('%1.2f'%(good_odor_corrs[-1]),fontsize=16)
                        ymax = ax.get_ylim()[1]
                        ax.set_yticks([0,ymax])
                        ax.set_yticklabels(['0',str(ymax)])
                        #ax.set_xlim(0.75,1.25)
                        #ax.set_xticks([0.75,1.25])
                        ax.set_xticks([])
                        #ax.set_xticklabels(['0.75','1.25'])
                        ## just to scale up the ticks fontsize.
                        axes_labels(ax,'','',adjustpos=False,fontsize=18)
                        #print filename, "corr =", good_odor_corrs[-1]
                ## air
                if air_corr<corr_cutoff:
                    good_air_corrs.append(air_corr)
                    if graph and len(good_air_corrs)<=numplots:
                        axair = figair.add_subplot(numrows,numcols,len(good_air_corrs))
                        axair.errorbar(x=responsetlist,y=mitral_responses_avg[6,0],\
                            yerr=mitral_responses_std[6,0],color=(1,0,0),linewidth=4)
                        axair.errorbar(x=responsetlist,y=mitral_responses_avg[6,1],\
                            yerr=mitral_responses_std[6,1],color=(0,1,0),linewidth=4)
                        axair.set_title('%1.2f'%(good_air_corrs[-1]),fontsize=16)
                        ymax = ax.get_ylim()[1]
                        ax.set_yticks([0,ymax])
                        ax.set_yticklabels(['0',str(ymax)])
                        #axair.set_xlim(0.75,1.25)
                        #axair.set_xticks([0.75,1.25])
                        axair.set_xticks([])
                        #axair.set_xticklabels(['0.75','1.25'])
                        ## just to scale up the ticks fontsize.
                        axes_labels(axair,'','',adjustpos=False,fontsize=18)
                        #print filename, "corr =", good_air_corrs[-1]

    ## correlation in change in firing rate
    corr_deltafrate = \
        dot(deltafrate_sis0,deltafrate_sis1)/ \
            sqrt(dot(deltafrate_sis0,deltafrate_sis0)*dot(deltafrate_sis1,deltafrate_sis1))
    ## * is element wise multiplication.
    ## contribution of each mitral-pair--odor combo to the deltafrate corr.
    ## normalized by the normalization of the dot-product.
    corr_deltafrate_eachmit = array(deltafrate_sis0)*array(deltafrate_sis1) \
            / sqrt(dot(deltafrate_sis0,deltafrate_sis0)*dot(deltafrate_sis1,deltafrate_sis1))
                    
    if num_responses==0:
        print "No result files found ..."
        sys.exit(1)
    overall_odor_mean = sum(total_frate/float(num_responses))/float(len(total_frate))
    overall_air_mean = sum(total_airfrate/float(num_responses/2.0))/float(len(total_airfrate))
    print "Phasic odor mean firing rate =",total_frate/float(num_responses)
    print "Overall odor mean = ", overall_odor_mean
    print "Phasic air mean firing rate =",total_airfrate/float(num_responses/2.0)
    print "Overall air mean = ", overall_air_mean
    print "delta frates of sisters :",zip(deltafrate_sis0,deltafrate_sis1)
    print "Correlation between sisters' change in firing rate with odor",corr_deltafrate

    if graph:
        fig.tight_layout() # for previous figure

        ## paper figure of decorr example and correlation distribution
        fig = figure(figsize=(columnwidth,linfig_height/2.0),dpi=300,facecolor='w')
        splax = fig.add_subplot(1,2,1)
        ###### paper figure for example decorrelation.
        ## odornum should be 5 (odor A) or 0 (odor B).
        plot_decorr_single(resultsdir=resultsdir,stimseed=844.0,numgloms=3,\
            inh=(False,False,False,False,False),odornum=0,splax=splax)

        ax = fig.add_subplot(1,2,2)
        #title("PDF of correlations between non-zero sister responses",fontsize=30)
        ax.hist(air_corrs,10,range=(-1.0,1.0),normed=True,histtype='step',\
            linewidth=linewidth,color='b',ls='dotted')
        ax.hist(odor_corrs,10,range=(-1.0,1.0),normed=True,histtype='step',\
            linewidth=linewidth,color='r',ls='solid')
        #ax.text(-0.3, 0.9, 'd', fontweight='bold', fontsize=label_fontsize, transform=ax.transAxes)
        beautify_plot(ax,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')
        axes_labels(ax,'correlation','prob. density',adjustpos=False)
        ax.set_xticks([-1,0,1])
        fig.tight_layout()
        fig_clip_off(fig)
        fig.savefig('../figures/decorr/sim_decorr.svg', dpi=fig.dpi)
        fig.savefig('../figures/decorr/sim_decorr.png', dpi=fig.dpi)
    
        fig = figure(figsize=(columnwidth/3.0, linfig_height/2.0), dpi=300, facecolor='w')
        ax = fig.add_subplot(111)
        #title("PDF of correlations between non-zero sister responses",fontsize=30)
        scatter(all_odor_corrs,corr_deltafrate_eachmit,marker='o',s=marker_size)
        axes_labels(ax,'correlation','deltafrate',adjustpos=False)
        #biglegend('upper left')
        fig.tight_layout()

    return corr_deltafrate, odor_corrs, air_corrs, overall_odor_mean, overall_air_mean

def plot_xcorrgrams():
    ## plot the average xcorrgrams over all stimuli
    ## for each netseed (figure) and inhibition type (panel)
    for neti,netseed in enumerate(net_seeds):
        
        numsubfigs_rows = len(num_gloms_list)
        numsubfigs_cols = len(inh_options)
        ## fig with subpanels for spike corr
        figxcorr = figure(facecolor='w')
        figtext(0.1,0.94,'odor A,B,air xcorrs between sisters: net'+str(netseed),fontsize=20)
        delaxes() # delete the main axes to accomodate subpanels
        ## printed output header
        print 'odors A, B, air phase similarity / corr between sisters: netseed = '+str(netseed)
        ## fig with subpanels for phase corr / similarity
        figcorr = figure(facecolor='w')
        figtext(0.1,0.94,'odorA,B,air ph sim/corr between sisters: net'+str(netseed),fontsize=20)
        delaxes() # delete the main axes to accomodate subpanels
        
        for ngi,num_gloms in enumerate(num_gloms_list):

            ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
            for plotyi,(inhi,inh) in enumerate(inh_options):

                ## spike corrs for all odors and airs to this network instance
                air_xcorrgrams = []
                odorA_xcorrgrams = []
                odorB_xcorrgrams = []
                ## phase tuning corrs for all odors and airs to this network instance
                air_corrs = []
                odorA_corrs = []
                odorB_corrs = []
                ## average change in frate (compared to air) for all odors
                mit0_Dfrates = []
                mit1_Dfrates = []
                mit2_Dfrates = []
                mit3_Dfrates = []

                for stimi,stimseed in enumerate(stim_seeds):

                    filename, switch_strs = \
                        get_filename(netseed,stimseed,inh,num_gloms,stimi,neti,inhi)
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    switches_str = string.join(switch_strs,'')

                    ## calculate spike time based xcorrs:
                    (air_corr,odorA_corr,odorB_corr),\
                        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
                        Dfrates = \
                        calc_corrs(filename, norm_str="overall", \
                            numbins=NUMBINS, bin_width_time=BIN_WIDTH_TIME, printinfo=False)
                    if air_corr<0 or odorA_corr<0 or odorB_corr<0:
                        print filename
                        print "air corr =",air_corr,"odor A corr =",odorA_corr,"odor B corr =",odorB_corr
                    ## if first value in list is nan, so are all the others; exclude them
                    if not math.isnan(airxcorrgram[0]): air_xcorrgrams.append(airxcorrgram)
                    if not math.isnan(odorAxcorrgram[0]): odorA_xcorrgrams.append(odorAxcorrgram)
                    if not math.isnan(odorBxcorrgram[0]): odorB_xcorrgrams.append(odorBxcorrgram)
                    ## if the phase similarity / correlation is nan, exclude
                    if not math.isnan(air_corr): air_corrs.append(air_corr)
                    if not math.isnan(odorA_corr): odorA_corrs.append(odorA_corr)
                    if not math.isnan(odorB_corr): odorB_corrs.append(odorB_corr)
                    ## average Delta frates for the sisters 0 and 1 -- odor A and B
                    mit0_Dfrates.extend((Dfrates[0],Dfrates[2]))
                    mit1_Dfrates.extend((Dfrates[1],Dfrates[3]))
                    #mit2_Dfrates.extend((Dfrates[4],Dfrates[6]))
                    #mit3_Dfrates.extend((Dfrates[5],Dfrates[7]))
                
                #### plot various analyses for this network instance + conditions

                ## plot spike time based xcorrs:
                ax = figxcorr.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                    (ngi)*numsubfigs_cols+plotyi+1)
                if len(air_xcorrgrams)>0:
                    ax.plot(tcorrlist, mean(air_xcorrgrams,axis=0), color=(0,0,0), label='air')
                if len(odorA_xcorrgrams)>0:
                    ax.plot(tcorrlist, mean(odorA_xcorrgrams,axis=0), color=(1,0,0), label='odor A')
                if len(odorB_xcorrgrams)>0:
                    ax.plot(tcorrlist, mean(odorB_xcorrgrams,axis=0), color=(0,1,0), label='odor B')
                ax.set_title(
                    str(num_gloms)+'G'+switches_str,fontsize=8)

                ## print phase similarity / corrs
                print '\t'+str(num_gloms)+'Glomeruli'+switches_str
                print "\t\tAir phase similarity =",mean(air_corrs),'SD',std(air_corrs)
                print "\t\tOdorA phase similarity =",mean(odorA_corrs),'SD',std(odorA_corrs)
                print "\t\tOdorB phase similarity =",mean(odorB_corrs),'SD',std(odorB_corrs)

                ## plot histogram of phase corr / similarity between sisters for air & odors:
                ## similar to plots in fig 7e of Ashesh Dhawale et al
                ax = figcorr.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                    (ngi)*numsubfigs_cols+plotyi+1)
                ## typically air_corrs will have less data points as many are nan-s
                ## how to normalize? normed=True is not what I want.
                if len(air_corrs)>0:
                    ax.hist(air_corrs, bins=5, histtype='step', normed=False, color=(0,0,0), label='air')
                if len(odorA_corrs)>0:
                    ax.hist(odorA_corrs, bins=5, histtype='step', normed=False, color=(1,0,0), label='odor A')
                if len(odorB_corrs)>0:
                    ax.hist(odorB_corrs, bins=5, histtype='step', normed=False, color=(0,1,0), label='odor B')
                ax.set_title(str(num_gloms)+'G'+switches_str,fontsize=8)
                ax.set_xlim(-1.0,1.0)

                ## print corr between mit0 and mit1 mean Dfrates:
                #print "\t\tmit0 mean Delta firing rate change (Hz) for various odors =",mit0_Dfrates
                #print "\t\tmit1 mean Delta firing rate change (Hz) for various odors =",mit1_Dfrates
                print "\t\tCorrelation of 'mean change in firing rate' Odor Response Spectrum (Ashesh et al)"
                print "\t\t\tbetween sisters 0 and 1 =",\
                    stats.pearsonr(mit0_Dfrates,mit1_Dfrates)[0]
                #print "\t\t\tbetween non-sisters 0 and 2 =",\
                #    stats.pearsonr(mit0_Dfrates,mit2_Dfrates)[0]

def plot_directed(glomnums):
    """ Be sure to set frac_directed=0.0/0.05 above at the very beginning. """
    odor_corrs_means = []
    odor_corrs_SDs = []
    air_corrs_means = []
    air_corrs_SDs = []
    corrs_deltafrate = []
    fig = figure()
    for gni,glomnum in enumerate(glomnums):
        print "Computing phasic and deltafrate correlations for # of gloms =",glomnum
        ## Set graph=True below to plot neg corr-ed responses too.
        corr_deltafrate, odor_corrs, air_corrs, overall_odor_mean, overall_air_mean = \
            plot_decorrs_special([glomnum],graph=True)
        ax = fig.add_subplot(len(glomnums),1,gni+1)
        #hist(air_corrs,20,range=(-1.0,1.0),normed=True,histtype='step',\
        #    color='b',linewidth=2,label='air %2.1f'%overall_air_mean+'Hz')
        hist(odor_corrs,20,range=(-1.0,1.0),normed=True,histtype='step',\
            color='r',linewidth=2,label='odor %2.1f'%overall_odor_mean+'Hz')
        ax.set_xticks([])
        #ax.set_xticklabels(['0.75','1.25'])
        ## just to scale up the ticks fontsize.
        axes_labels(ax,'','',adjustpos=False,fontsize=34)

        corrs_deltafrate.append(corr_deltafrate)
        ## mean and SD of phasic correlations of odor and air
        odor_corrs_means.append(mean(odor_corrs))
        odor_corrs_SDs.append(std(odor_corrs))
        air_corrs_means.append(mean(air_corrs))
        air_corrs_SDs.append(std(air_corrs))

        ax.set_yticks([])
        #biglegend(legendlocation='upper left')
        if gni == len(glomnums)-1:
            ax.set_xticks([-1.0,0.0,1.0])
            ax.set_xticklabels(['-1','0','1'])
            axes_labels(ax,'phase correlation','',adjustpos=False,fontsize=30)
    plt.tight_layout()

    ## mean phase corr vs number of connected gloms
    fig=figure()
    ax=fig.add_subplot(111)
    #plot(glomnums,air_corrs_means,color='b',linewidth=2,label='air')
    plot(glomnums,odor_corrs_means,color='r',linewidth=2,label='odor')
    ax.set_xticks(glomnums)
    ax.set_xticklabels([str(glomnum) for glomnum in glomnums])
    axes_labels(ax,'# of connected glomeruli','phase correlation mean',\
        adjustpos=False,fontsize=30)
    #biglegend(legendlocation='lower left')
    plt.tight_layout()
    ## spread of phase corr vs number of connected gloms
    fig=figure()
    ax=fig.add_subplot(111)
    #errorbar(glomnums,air_corrs_SDs,color='b',linewidth=2,label='air')
    errorbar(glomnums,odor_corrs_SDs,color='r',linewidth=2,label='odor')
    ax.set_xticks(glomnums)
    ax.set_xticklabels([str(glomnum) for glomnum in glomnums])
    axes_labels(ax,'# of connected glomeruli','phase correlation spread',\
        adjustpos=False,fontsize=30)
    #biglegend(legendlocation='upper left')
    plt.tight_layout()
    ## delta frate corr vs number of connected gloms
    fig=figure()
    ax=fig.add_subplot(111)
    plot(glomnums,corrs_deltafrate,color='b',linewidth=2)
    ax.set_xticks(glomnums)
    ax.set_xticklabels([str(glomnum) for glomnum in glomnums])
    axes_labels(ax,'# of connected glomeruli','$\Delta$frate correlation',\
        adjustpos=False,fontsize=30)
    tight_layout()
    
def plot_across_sims_paperfigure():
    odor_corrs_means = []
    odor_corrs_SDs = []
    air_corrs_means = []
    air_corrs_SDs = []
    corrs_deltafrate = []
    fig = figure(figsize=(columnwidth/2.0,linfig_height*2.0),dpi=300,facecolor='w')
    ## [(dirextn, directed, frac_directed, inh_options, numgloms, color),...]
    ## [ ('', directed, 0.0, novarmit, [...]), ('', directed, 0.0, varmit, [...]), ('', directed, 0.05, novarmit, [...]) ]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    #decorr_net_options = [ (True,0.0,(0,(False,False,False,False,False)),[2,3,6],'b'),
    #    (True,0.0,(0,(False,False,False,False,True)),[2,3],'g'),
    #    (True,0.05,(0,(False,False,False,False,False)),[2,3,4,5,6,7],'r') ]
    ####### CAUTION: EVEN THOUGH folder suffix is _0.33x_may14; the scaling is 0.5x odor, 1x air, 1x background
    decorr_net_options = [ \
        ('_0.33x_may14',False,0.0,(0,(False,False,False,False,False)),[3],'k'),\
        ('_0.33x_may14',True,0.0,(1,(False,False,False,False,False)),[3],'k'),\
        ('_0.33x_may14',True,0.0,(2,(False,False,False,False,True)),[3],'k'),\
        ('',True,0.01,(3,(False,False,False,False,False)),[3],'k'),\
        ('_2x4x_may16',True,0.01,(4,(False,False,False,False,False)),[7],'k') ]
    #decorr_net_options = [ (True,0.0,(0,(False,False,False,False,False)),[2],'b'),
    #    (True,0.05,(0,(False,False,False,False,False)),[3],'r') ]
    neti=0
    indexstrs = []
    colorslist = []
    numsubplots = sum([1 for net_option in decorr_net_options for glomnum in net_option[4]])
    for (dirextn,directed,frac_directed,inh_option,numglomslist,color) in decorr_net_options:
        for numgloms in numglomslist:
            print "Computing phasic and deltafrate correlations for # of gloms =",numgloms
            print "The inh tweak (no_singles,no_joints,no_lat,no_PGs,varyRMP) is :",inh_option
            ## Set graph=True below to plot neg corr-ed responses too.
            corr_deltafrate, odor_corrs, air_corrs, overall_odor_mean, overall_air_mean = \
                plot_decorrs_special_paperfigure('../results/odor_morphs'+dirextn, \
                    [numgloms],[inh_option],directed,frac_directed,graph=False)
            ax = fig.add_subplot(numsubplots,1,neti+1)
            #hist(air_corrs,20,range=(-1.0,1.0),normed=True,histtype='step',\
            #    color='b',linewidth=2,label='air %2.1f'%overall_air_mean+'Hz')
            counts,bins,patches = hist(odor_corrs,10,range=(-1.0,1.0),normed=True,histtype='step',\
                color=color,linewidth=linewidth,label='odor %2.1f'%overall_odor_mean+'Hz')
            bar(bins[:5],counts[:5],width=bins[1]-bins[0],facecolor='r',edgecolor='r')
            ax.set_xticks([])
            ax.set_yticks([])

            corrs_deltafrate.append(corr_deltafrate)
            ## mean and SD of phasic correlations of odor and air
            odor_corrs_mean = mean(odor_corrs)
            odor_corrs_SD = std(odor_corrs)
            odor_corrs_means.append(odor_corrs_mean)
            odor_corrs_SDs.append(odor_corrs_SD)
            air_corrs_means.append(mean(air_corrs))
            air_corrs_SDs.append(std(air_corrs))
            if neti==0:
                ymax = ax.get_ylim()[1]*1.2 # hard coded to 1.2x initial-ymax, to fit all subplots!
                ylim = (0.0,ymax)
                ax.set_ylim(ylim)
                ymid = ymax/2.0
            else:
                ax.set_ylim(ylim)
            ### The mean and SD are obvious from the figure.
            ### Drawing extra lines/bars like below are distracting -- hence commented out.
            ### draw an arrow of default 'Curve' style (only line) for the mean on the correlation distribution
            #arrmean = matplotlib.patches.FancyArrowPatch(
            #    (odor_corrs_mean, ymax*0.75), (odor_corrs_mean, ymax*0.25), linewidth=linewidth)
            #ax.add_patch(arrmean)
            ### draw an arrow of style BarAB for the SD on the correlation distribution
            #arrSD = matplotlib.patches.FancyArrowPatch(
            #    (odor_corrs_mean-odor_corrs_SD, ymid), (odor_corrs_mean+odor_corrs_SD, ymid), linewidth=linewidth )
            #arrSD.set_arrowstyle('|-|',widthA=10,angleA=90,widthB=10,angleB=90)
            #ax.add_patch(arrSD)
            beautify_plot(ax,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')
            if neti==2:
                axes_labels(ax,'','probability density',adjustpos=False,fontsize=label_fontsize)
            if neti>=numsubplots-1:
                ax.set_xticks([-1,0,1])
                ax.set_xticklabels(['-1','0','1'])
                ## scale up the ticks and axes labels fontsize.
                axes_labels(ax,'correlation','',adjustpos=False,fontsize=label_fontsize)
            else: ax.set_xticks([])
            #ax.text(0.15, 1.0, ['d','e','f','g','h','i'][neti], \
            #    fontweight='bold', fontsize=label_fontsize, transform=ax.transAxes)
            ## Print the delta-frate correlation on the plot/histogram
            ax.text(0.15, 0.75, '%1.3f'%(corr_deltafrate), fontsize=label_fontsize, transform=ax.transAxes)

            indexstr = ''
            if directed: indexstr += 'directed, '
            if frac_directed>0.0: indexstr += 'differential, '
            if inh_option[1][4]: indexstr += 'varymits, '
            indexstr += 'latgloms '+str(numgloms-1)
            indexstrs.append(indexstr)
            colorslist.append(color)
            neti+=1

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig_clip_off(fig)
    fig.savefig('../figures/decorr/decorr_conns.svg',dpi=fig.dpi)
    fig.savefig('../figures/decorr/decorr_conns.png',dpi=fig.dpi)

    ## figure for plotting mean and SD of the distribution
    mainfig = figure(figsize=(8, 6), dpi=150, facecolor='w')
    ax=mainfig.add_subplot(111)
    indices = range(neti)
    ax.bar(indices, odor_corrs_means, color=colorslist,
       align='center', yerr=odor_corrs_SDs, ecolor='black', width=0.3)
    ax.set_ylim([-0.25,1.0])
    ax.set_yticks([-0.25,0,0.5,1.0])
    ax.set_yticklabels(['-0.25','0','0.5','1'])
    ax.set_xticks(indices)
    ax.set_xticklabels(indexstrs)
    mainfig.autofmt_xdate() # rotates xlabels
    axes_labels(ax,'','',\
        adjustpos=False,fontsize=16)
    mainfig.tight_layout()

def plot_peaks_tufted_vs_mitrals_paper_figure():
    fig = figure(figsize=(columnwidth/2.0,linfig_height*2.0),dpi=300,facecolor='w')
    ## [(dirextn, directed, frac_directed, inh_options, numglomslist, color),...]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ####### CAUTION: EVEN THOUGH folder suffix is _0.33x_may14; the scaling is 0.5x odor, 1x air, 1x background
    decorr_net_options = [ \
        ('',True,0.01,(1,(False,False,False,False,False)),[2],'k'),\
        ('',True,0.01,(2,(False,False,False,True,False)),[2],'k') ]
    neti=0
    indexstrs = []
    colorslist = []
    numsubplots = sum([1 for net_option in decorr_net_options for glomnum in net_option[4]])
    for (dirextn,directed,frac_directed,inh_option,numglomslist,color) in decorr_net_options:
        for numgloms in numglomslist:
            print "Computing phasic and deltafrate correlations for # of gloms =",numgloms
            print "The inh tweak (no_singles,no_joints,no_lat,no_PGs,varyRMP) is :",inh_option
            ## Set graph=True below to plot neg corr-ed responses too.
            corr_deltafrate, odor_corrs, air_corrs, overall_odor_mean, overall_air_mean = \
                plot_decorrs_special_paperfigure('../results/odor_morphs'+dirextn, \
                    [numgloms],[inh_option],directed,frac_directed,graph=False)
            ax = fig.add_subplot(numsubplots,1,neti+1)
            #hist(air_corrs,20,range=(-1.0,1.0),normed=True,histtype='step',\
            #    color='b',linewidth=2,label='air %2.1f'%overall_air_mean+'Hz')
            counts,bins,patches = hist(odor_corrs,10,range=(-1.0,1.0),normed=True,histtype='step',\
                color=color,linewidth=linewidth,label='odor %2.1f'%overall_odor_mean+'Hz')
            bar(bins[:5],counts[:5],width=bins[1]-bins[0],facecolor='r',edgecolor='r')
            ax.set_xticks([])
            ax.set_yticks([])

            corrs_deltafrate.append(corr_deltafrate)
            ## mean and SD of phasic correlations of odor and air
            odor_corrs_mean = mean(odor_corrs)
            odor_corrs_SD = std(odor_corrs)
            odor_corrs_means.append(odor_corrs_mean)
            odor_corrs_SDs.append(odor_corrs_SD)
            air_corrs_means.append(mean(air_corrs))
            air_corrs_SDs.append(std(air_corrs))
            if neti==0:
                ymax = ax.get_ylim()[1]*1.2 # hard coded to 1.2x initial-ymax, to fit all subplots!
                ylim = (0.0,ymax)
                ax.set_ylim(ylim)
                ymid = ymax/2.0
            else:
                ax.set_ylim(ylim)
            ### The mean and SD are obvious from the figure.
            ### Drawing extra lines/bars like below are distracting -- hence commented out.
            ### draw an arrow of default 'Curve' style (only line) for the mean on the correlation distribution
            #arrmean = matplotlib.patches.FancyArrowPatch(
            #    (odor_corrs_mean, ymax*0.75), (odor_corrs_mean, ymax*0.25), linewidth=linewidth)
            #ax.add_patch(arrmean)
            ### draw an arrow of style BarAB for the SD on the correlation distribution
            #arrSD = matplotlib.patches.FancyArrowPatch(
            #    (odor_corrs_mean-odor_corrs_SD, ymid), (odor_corrs_mean+odor_corrs_SD, ymid), linewidth=linewidth )
            #arrSD.set_arrowstyle('|-|',widthA=10,angleA=90,widthB=10,angleB=90)
            #ax.add_patch(arrSD)
            beautify_plot(ax,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')
            if neti==2:
                axes_labels(ax,'','probability density',adjustpos=False,fontsize=label_fontsize)
            if neti>=numsubplots-1:
                ax.set_xticks([-1,0,1])
                ax.set_xticklabels(['-1','0','1'])
                ## scale up the ticks and axes labels fontsize.
                axes_labels(ax,'correlation','',adjustpos=False,fontsize=label_fontsize)
            else: ax.set_xticks([])
            #ax.text(0.15, 1.0, ['d','e','f','g','h','i'][neti], \
            #    fontweight='bold', fontsize=label_fontsize, transform=ax.transAxes)
            ## Print the delta-frate correlation on the plot/histogram
            ax.text(0.15, 0.75, '%1.3f'%(corr_deltafrate), fontsize=label_fontsize, transform=ax.transAxes)

            indexstr = ''
            if directed: indexstr += 'directed, '
            if frac_directed>0.0: indexstr += 'differential, '
            if inh_option[1][4]: indexstr += 'varymits, '
            indexstr += 'latgloms '+str(numgloms-1)
            indexstrs.append(indexstr)
            colorslist.append(color)
            neti+=1

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig_clip_off(fig)
    fig.savefig('../figures/decorr/decorr_conns.svg',dpi=fig.dpi)
    fig.savefig('../figures/decorr/decorr_conns.png',dpi=fig.dpi)

    ## figure for plotting mean and SD of the distribution
    mainfig = figure(figsize=(8, 6), dpi=150, facecolor='w')
    ax=mainfig.add_subplot(111)
    indices = range(neti)
    ax.bar(indices, odor_corrs_means, color=colorslist,
       align='center', yerr=odor_corrs_SDs, ecolor='black', width=0.3)
    ax.set_ylim([-0.25,1.0])
    ax.set_yticks([-0.25,0,0.5,1.0])
    ax.set_yticklabels(['-0.25','0','0.5','1'])
    ax.set_xticks(indices)
    ax.set_xticklabels(indexstrs)
    mainfig.autofmt_xdate() # rotates xlabels
    axes_labels(ax,'','',\
        adjustpos=False,fontsize=16)
    mainfig.tight_layout()

## Obsolete: This is called by fit_odor_morphs.py as a panel in Figure.
def plot_chisq_hist_paperfigure(ax1,ax2,resultsdir='../results/odor_morphs'):
    """ Plot chi-sq histogram of the morph fits. """
    #fig = figure(figsize=(columnwidth/3.0,columnwidth/2.0),dpi=300,facecolor='w') # 'none' is transparent
    #ax = fig.add_subplot(111,frameon=False)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [ (0,(False,False,False,False,False),'lat inh') ]
    for ploti,(inhi,inh,inhstr) in enumerate(inh_options):
        chisqs = []
        lin_chisqs = []
        n_accept = 0
        for stimi,stimseed in enumerate(stim_seeds):
            if not salient: net_seeds = [stimseed]
            for neti,netseed in enumerate(net_seeds):
                for ngi,num_gloms in enumerate([3]):

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,num_gloms,stimi,neti,inhi,resultsdir=resultsdir)
                    switches_str = string.join(switch_strs,'')
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    #print filename
                    for fitted_mitral in [0,1]:
                        ## First the weighted-linear sigmoid:
                        ## If the fitted params file does not exist, create it (them).
                        if not os.path.exists(filename+'_params'+str(fitted_mitral)):
                            print "fitting file",filename
                            ## fit the responses for this result file
                            params,chisq,inputsA,inputsB,fitted_responses,\
                                numavgs,firingbinsmeanList,firingbinserrList\
                                = fit_om.fit_morphs(filename, fitted_mitral, 'arb')
                        else:
                            f = open(filename+'_params'+str(fitted_mitral),'r')
                            params,chisq = pickle.load(f)
                            f.close()
                        chisqs.append(chisq)
                        ## linear-sigmoid [perhaps get rid of sigmoid also a la Priyanka?]
                        ## If the fitted params file does not exist, create it (them).
                        if not os.path.exists(filename+'_paramslin'+str(fitted_mitral)):
                            ## fit the responses for this result file
                            params,chisq,inputsA,inputsB,fitted_responses,\
                                numavgs,firingbinsmeanList,firingbinserrList\
                                = fit_om.fit_morphs(filename, fitted_mitral, 'lin')
                        else:
                            f = open(filename+'_paramslin'+str(fitted_mitral),'r')
                            params,chisq = pickle.load(f)
                            f.close()
                        lin_chisqs.append(chisq)
                        n_accept += 1

        ax1.hist(chisqs,20,histtype='step',linewidth=linewidth,label=inhstr,color='k')
        ax2.hist(lin_chisqs,20,histtype='step',linewidth=linewidth,label=inhstr,color='k')
        print "Number of mitral cells accepted =",n_accept
        
        ## beautify plots
        for axnum,ax in enumerate([ax1,ax2]):
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            #text(0.15, 1.0, ['E','F'][axnum], fontsize=label_fontsize, transform = ax.transAxes)
            ax.get_yaxis().set_ticks_position('left')
            ax.get_xaxis().set_ticks_position('bottom')
            xmin, xmax = ax.get_xaxis().get_view_interval()
            ymin, ymax = ax.get_yaxis().get_view_interval()
            ax.set_xlim(0,xmax)
            ax.set_ylim(0,ymax)
            ax.set_xticks([0,xmax])
            ax.set_yticks([0,ymax])
            ax.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
            ax.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
            axes_labels(ax,'','',adjustpos=False) # sets font-size for tick labels also
        ax2.text(-0.4,1.3,'count',fontsize=label_fontsize, rotation='vertical', transform=ax.transAxes)
        ax2.text(0.3,-0.38,'chi-sq',fontsize=label_fontsize, transform=ax.transAxes)
    #fig.tight_layout()
    #subplots_adjust(top=0.90,wspace=1.0)
    #fig.savefig('../figures/morph_chisqs.svg', bbox_inches='tight',dpi=fig.dpi)
    #fig.savefig('../figures/morph_chisqs.png', bbox_inches='tight',dpi=fig.dpi)

def plot_onemitexample_R2N_hist_paperfigure(eg_netseed,eg_mitnum,resultsdir='../results/odor_morphs'):
    """ Plot residual to noise histogram of the morph fits. """
    fig = figure(figsize=(columnwidth,columnwidth/2.0),dpi=300,facecolor='w') # 'none' is transparent
    ax3 = fig.add_subplot(2,3,1)
    ax4 = fig.add_subplot(2,3,2)
    ax5 = fig.add_subplot(2,3,4)
    ax6 = fig.add_subplot(2,3,5)
    ax1 = fig.add_subplot(2,3,3)
    ax2 = fig.add_subplot(2,3,6)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [ (0,(False,False,False,False,False),'lat inh') ]
    for ploti,(inhi,inh,inhstr) in enumerate(inh_options):
        R2Ns = []
        lin_R2Ns = []
        chilist = []
        n_accept = 0
        for stimi,stimseed in enumerate(stim_seeds):
            if not salient: net_seeds = [stimseed]
            for neti,netseed in enumerate(net_seeds):
                for ngi,num_gloms in enumerate([3]):

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,num_gloms,stimi,neti,inhi,resultsdir=resultsdir)
                    switches_str = string.join(switch_strs,'')
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    print filename
                    for fitted_mitral in [0,1]:
                        ## First the weighted-linear sigmoid:
                        ## If the fitted params file does not exist, create it (them).
                        if not os.path.exists(filename+'_params'+str(fitted_mitral)):
                            print "fitting file",filename
                            refit = True
                        else: refit = False
                        ## read in params & responses for this result file
                        mit_fit_params = \
                            fit_om.fit_morphs(filename, fitted_mitral, 'arb', refit=refit)
                        params,chisq,inputsA,inputsB,fitted_responses,\
                                numavgs,firingbinsmeanList,firingbinserrList = mit_fit_params
                        S2N,S2R = forR2N.residual2noise(fitted_responses[-2],firingbinsmeanList[-2],\
                            firingbinserrList[-2]*sqrt(numavgs),starti=0) # odor A
                        R2N_A = S2N/S2R
                        if isnan(R2N_A): continue
                        S2N,S2R = forR2N.residual2noise(fitted_responses[0],firingbinsmeanList[0],\
                            firingbinserrList[0]*sqrt(numavgs),starti=0) # odor B
                        R2N_B = S2N/S2R
                        if isnan(R2N_B): continue
                        R2Ns.append(R2N_A)
                        R2Ns.append(R2N_B)
                        if netseed == eg_netseed and fitted_mitral == eg_mitnum:
                            fit_om.plot_example_onemit(ax3,ax4,eg_mitnum,mit_fit_params)
                        
                        ## Linear-rectifier or Linear-sigmoid depending on FULLlin variable above.
                        ## If the fitted params file does not exist, create it (them).
                        if not os.path.exists(filename+'_params'+linextn+str(fitted_mitral)):
                            print "fitting FULLlin file",filename
                            refit = True
                        else: refit = False
                        ## fit/get the params and responses for this result file
                        mit_fit_params = \
                            fit_om.fit_morphs(filename, fitted_mitral, 'lin', refit=refit)
                        params,chisq,inputsA,inputsB,fitted_responses,\
                            numavgs,firingbinsmeanList,firingbinserrList = mit_fit_params
                        S2N,S2R = forR2N.residual2noise(fitted_responses[-2],firingbinsmeanList[-2],\
                            firingbinserrList[-2]*sqrt(numavgs),starti=0) # odor A
                        R2N_A = S2N/S2R
                        if isnan(R2N_A): continue
                        S2N,S2R = forR2N.residual2noise(fitted_responses[0],firingbinsmeanList[0],\
                            firingbinserrList[0]*sqrt(numavgs),starti=0) # odor B
                        R2N_B = S2N/S2R
                        if isnan(R2N_B): continue
                        lin_R2Ns.append(R2N_A)
                        lin_R2Ns.append(R2N_B)
                        chilist.append(sqrt(chisq))
                        if netseed == eg_netseed and fitted_mitral == eg_mitnum:
                            fit_om.plot_example_onemit(ax5,ax6,eg_mitnum,mit_fit_params)

                        n_accept += 1

        R2N_max = 1.0
        ax1.hist(clip(R2Ns,0,R2N_max),20,normed=True,edgecolor='b',facecolor='b')
        _,y1 = ax1.get_ylim()
        ax2.hist(clip(lin_R2Ns,0,R2N_max),20,normed=True,edgecolor='b',facecolor='b')
        #ax2.hist(clip(chilist,0,R2N_max),20,normed=True,edgecolor='b',facecolor='b')
        _,y2 = ax2.get_ylim()
        yR2Nmax = max(y1,y2)
        print "Number of mitral cells accepted =",n_accept
        
        ## beautify plots
        for axnum,ax in enumerate([ax1,ax2]):
            xmin,xmax,ymin,ymax = \
                beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
            ax.set_xlim([0,R2N_max])
            ax.set_xticks([0,R2N_max])
            ax.set_ylim([0,yR2Nmax])
            ax.set_yticks([0,yR2Nmax])
        for ax in [ax1,ax3,ax4]:
            ax.set_xticklabels(['',''])
        ## axes_labels() sets sizes of tick labels too.
        axes_labels(ax1,'','prob. density',adjustpos=False,xpad=0,ypad=0)
        ax1.yaxis.set_label_coords(-0.29,-0.3)
        axes_labels(ax2,'$\sqrt{residual/noise}$','',adjustpos=False,xpad=1,ypad=0)

        axes_labels(ax3,'','firing rate (Hz)',adjustpos=False,xpad=0,ypad=0)
        ax3.yaxis.set_label_coords(-0.29,-0.3)
        axes_labels(ax5,'time (s)','',adjustpos=False,xpad=3,ypad=0)

        axes_labels(ax4,'','fitted weight',adjustpos=False,xpad=0,ypad=0)
        ax4.yaxis.set_label_coords(-0.24,-0.3)
        axes_labels(ax6,'conc (% SV)','',adjustpos=False,xpad=3,ypad=0)

        fig_clip_off(fig)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.3,wspace=0.5) # has to be after tight_layout()
        fig.savefig('../figures/morphs_R2Ns.svg',dpi=fig.dpi)
        fig.savefig('../figures/morphs_R2Ns.png',dpi=fig.dpi)

if __name__ == "__main__":
    if len(stim_seeds)<5:
        plot_responses()
        #plot_xcorrgrams()
        #plot_decorr_single()
    else:
        ## default network: set directed=True, frac_directed=0.05, varmit False (in inh_options)
        #glomnums = range(1,6)
        ## non-decorrelating networks:
        ## set directed=True, frac_directed=0.0, varmit True/False (in inh_options)
        ## OR directed=False; varmit False (in inh_options)
        #glomnums = [2,3,6]
        #plot_directed(glomnums)
        
        ## PAPER figure 7: plot the correlation distributions for various network connectivities
        ## CAUTION: SET stim_seeds AT THE TOP to (750.0,1100.0)
        #plot_across_sims_paperfigure()
        
        ###### for the default network, plot the best decorr-ed responses:
        ###### PAPER figure 6: for distribution of air and odor distributions & example.
        ## CAUTION: SET stim_seeds AT THE TOP to (750.0,1100.0)
        if len(sys.argv)>1: resultsdir = sys.argv[1]
        else: resultsdir = '../results/odor_morphs/'
        print 'Using results directory:',resultsdir
        #plot_decorrs_special_paperfigure(resultsdir,[3],[(0,(False,False,False,False,False))],\
        #    _directed=True,_frac_directed=0.01,graph=True)
        ## PAPER figure 6: example decorr: choose odornum as 5 for odorA and 0 for odor B
        ## CAUTION: SET stim_seeds AT THE TOP to (750.0,1100.0)
        #plot_responses_mits_paperfigure( resultsdir,odornum=0,stimseed=844.0,\
        #    numgloms=3,inh=(False,False,False,False,False) )
        
        ############ ODOR MORPHS -- Adil style fits
        ## OBSOLETE -- odor morphs chisq histogram
        #plot_chisq_hist()
        
        ## PAPER FIGURE supplementary figure 5 : Adil's morph fits
        ## CAUTION: SET stim_seeds AT THE TOP to (750.0,800.0)
        ## plots one mitral morphs fits example
        ## AND histogram of overall fits
        ## pass example netseed, example mitnum and resultsdir
        #plot_onemitexample_R2N_hist_paperfigure(754.0,0,resultsdir) # 752.0,0 is an alternative eg

        ## PAPER FIGURE supplementary figure 6: tufted (no PG) vs mitral (with PG) phase difference
        ## JUST STARTED, INCOMPLETE
        plot_peaks_tufted_vs_mitrals_paper_figure()

    show()
