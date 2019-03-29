# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: python2.6 average_odor_scaledpulses.py

import os,sys,math,string
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has RESPIRATION
from pylab import * # part of matplotlib that depends on numpy but not scipy
from scipy import stats
from scipy import signal # for gaussian smoothing
from sim_utils import * # has rebin() and imports data_utils.py for axes_off()
from calc_corrs import * # has calc_corrs()
## BE CAREFUL: use _constair_separateodors below, others are obsolete.
import fit_scaledpulses as fit_sp # has fit_plot_scaledpulses and ref pulse constants
import average_odor_morphs as corr_utils # has get_filename() for morphs

IN_VIVO = True
directed = True
FRAC_DIRECTED = 0.01 ## overrides the one in networkConstants.py (came in via sim_utils.py)
## Below two overide the variables that came in from stimuliConstantsMinimal.py via stimuliConstants.py
NONLINEAR_ORNS = False
NONLINEAR_TYPE = 'P' # P for primary glom non-linear, L for lateral gloms non-linear
_scaledWidth = 0.2 # s # overrides that in stimuliConstantsMinimal.py
num_scalings = len(pulseList)-1

#fullfilename = '../results/odor_morphs/morphs_random'
#if NONLINEAR_ORNS: fullfilename += 'NL'+NONLINEAR_TYPE
#fullfilename += '.pickle'
#fullfile = open(fullfilename,'r')
#morphs = pickle.load(fullfile)
#fullfile.close()

PLOTRESP_NUM = 1 # whether to plot 2 respiration cycles or 1
NUMBINS = 5
BIN_WIDTH_TIME = RESPIRATION/NUMBINS
bindt = RESPIRATION/float(NUMBINS)
## I take the last PLOTRESP_NUM of respiration cycles out of NUM_RESPS simulated
responsetlist = arange( SETTLETIME+(NUM_RESPS-PLOTRESP_NUM)*RESPIRATION+bindt/2.0, \
    SETTLETIME+NUM_RESPS*RESPIRATION, RESPIRATION/NUMBINS )

min_frate_cutoff = 0.25 # Hz

salient = False#True
if salient:
    stim_seeds = [-1,-19]#range(-1,-37,-1)#[-1,-2,-3,-4,-5,-6,-7,-8,-9]
    num_gloms_list = [3]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order,below options are:
    ## all cells; no lat; no joints, varyRMP; no PGs; no singles + no joints, only mitrals
    #inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
    #    (2,(False,True,False,False,False)), (3,(False,False,False,False,True)), \
    #    (4,(False,False,False,True,False)), (5,(True,True,False,False,False)), \
    #    (6,(True,True,False,True,False))]
    inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)) ]
else:
    #stim_seeds = [157.0,160.0,190.0,191.0,212.0,441.0]
    #num_gloms_list = [5,2]
    stim_seeds = arange(750.0,800.0,1.0)#[157.0,160.0,190.0,191.0]
    num_gloms_list = [3]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order,below options are: all cells; no lat; no joints, varyRMP; no PGs; only mitrals
    inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
        (2,(True,False,False,False,False)), (3,(True,True,False,False,True)), \
        (6,(False,False,False,True,False)), (8,(True,True,False,True,False))]
net_seeds = [200.0]

def get_filename(netseed,stimseed,inh,ngi,\
        nl_orns=NONLINEAR_ORNS,nl_type=NONLINEAR_TYPE,\
        resultsdir='../results/odor_pulses'):
    ### read filename from the output file of automated run
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
    filename = resultsdir+'/scaledpulses_width'+str(_scaledWidth)+\
        '_netseed'+str(netseed)+'_stimseed'+str(stimseed)
    if nl_orns: filename += '_NL'+nl_type
    filename += singles_str+joints_str+pgs_str+lat_str+varmitstr+\
        '_numgloms'+str(num_gloms_list[ngi])
    if directed: filename += '_directed'+str(FRAC_DIRECTED)
    filename += '.pickle'
    return filename, (singles_str, joints_str, pgs_str, lat_str, varmitstr)

def plot_lin_contribs_paperfigure():
    """ plot linearity scores for all options
    """
    fig1 = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w')
    ax1a = plt.subplot2grid((2,5),(0,0)) # full
    ax1b = plt.subplot2grid((2,5),(1,0)) # full
    ax2a = plt.subplot2grid((2,5),(0,1)) # no lat
    ax2b = plt.subplot2grid((2,5),(1,1)) # no lat
    ax3a = plt.subplot2grid((2,5),(0,2)) # no self
    ax3b = plt.subplot2grid((2,5),(1,2)) # no self
    ax4a = plt.subplot2grid((2,5),(0,3)) # no inh
    ax4b = plt.subplot2grid((2,5),(1,3)) # no inh
    ax5a = plt.subplot2grid((2,5),(0,4)) # non-lin
    ax5b = plt.subplot2grid((2,5),(1,4))
    ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [
        ('',0,(False,False,False,False,False),'all',False,None,(ax1a,ax1b)), \
        ('',1,(False,False,True,False,False),'no lat-inh',False,None,(ax2a,ax2b)), \
        ('',2,(True,False,False,True,False),'no self-inh',False,None,(ax3a,ax3b)), \
        ('',3,(True,True,False,True,False),'no inh',False,None,(ax4a,ax4b)), \
        ('',4,(False,False,False,False,False),'non-lin ORNs',True,'P',(ax5a,ax5b)) ]
    Rsqs_all = [ [] for i in range(len(inh_options)) ]
    R_all = [ [] for i in range(len(inh_options)) ]
    maxy_hist = 0.0
    for ploti,(dirextn,inhi,inh,inhstr,nl_orns,nl_type,(axa,axb)) in enumerate(inh_options):
        slopes_all = []
        peaks_all = []
        avg_frates_all = []
        for stimi,stimseed in enumerate(stim_seeds):
            filename, switch_strs \
                = get_filename(stimseed,stimseed,inh,0,nl_orns,nl_type,\
                    resultsdir='../results/odor_pulses'+dirextn)
            ## if the result file for these seeds & tweaks doesn't exist,
            ## then carry on to the next.
            if not os.path.exists(filename): continue
            worker = fit_sp.fit_plot_scaledpulses(filename,stimseed,_scaledWidth)
            fits_2mits,peaks_2mits = worker.fit_pulses(dirextn,True,False,False) # (dirextn,noshow,savefig,test)
            for mitrali in [0,1]:
                _,slopes,intercepts,r_values,_,_,se_slope,avg_frates = zip(*fits_2mits[mitrali])
                ## leave out the fit for reference with itself
                r_values = delete(r_values,fit_sp.ref_response_scalenum-1)
                R_all[ploti].extend(r_values)
                Rsqs_all[ploti].extend(r_values**2)
                nan_present = False
                for slope in slopes:
                    if isnan(slope): nan_present = True
                intercepts_bad = False
                #for intercept in intercepts:
                #    if abs(intercept)>2.0: intercepts_bad = True
                if not nan_present and not intercepts_bad:
                    slopes_all.append(slopes)
                    peaks_all.append(peaks_2mits[mitrali])
                    print filename
                else:
                    print 'Left out :',filename
                avg_frates_all.append(avg_frates)
        slopes_mean = append([0],mean(slopes_all,axis=0))
        slopes_std = append([0],std(slopes_all,axis=0))
        peaks_mean = append([0],mean(peaks_all,axis=0))
        peaks_std = append([0],std(peaks_all,axis=0))
        avg_frates_mean = mean(avg_frates_all,axis=0)
        print "The average firing rates for the different scaling is",avg_frates_mean
        print "Peaks mean is",peaks_mean
        
        ## plots
        axa.hist(R_all[ploti],10,range=(0,1.0),normed=True,histtype='step',\
            linewidth=linewidth,color='k',ls='solid')
        xmin,xmax,ymin,ymax = \
            beautify_plot(axa,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left',yticks=[])
        axa.set_ylim(0,ymax) # later this is set to the max over all histograms maxy_hist
        maxy_hist = max(maxy_hist,ymax)
        ## x-axis label for full first row
        if inhi==2:
            axes_labels(axa,"R (correlation)","",xpad=2)
            #axa.xaxis.set_label_coords(1.3,-0.25)
        ### inset plots within subplots
        ### passing transform=ax.transAxes to add_axes() doesn't work, hence jugglery from
        ### http://matplotlib.1069221.n5.nabble.com/Adding-custom-axes-within-a-subplot-td20316.html
        #Bbox = matplotlib.transforms.Bbox.from_bounds(1.6, 0.4, 0.6, 0.6) 
        #trans = ax.transAxes + fig1.transFigure.inverted() 
        #l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox,trans).bounds
        #axinset = fig1.add_axes([l, b, w, h])
        xconcs = array(scaledList[1:])
        ref_scale = scaledList[fit_sp.ref_response_scalenum]
        print "Slopes are =",slopes_mean*ref_scale
        axb.plot(range(6),range(6),color=(0,0,0.7,0.5),dashes=(2.0,1.0)) # linear reference
        axb.errorbar(x=append([0],xconcs),y=slopes_mean*ref_scale,yerr=slopes_std*ref_scale,\
            color='b',linewidth=linewidth,capsize=cap_size,dashes=(2.0,1.0))
        axb_twin = axb.twinx()
        axb_twin.errorbar(x=append([0],xconcs),y=peaks_mean,yerr=peaks_std,\
            color='k',linewidth=linewidth,capsize=cap_size,dashes=(0.5,1.0))
        #for slopes in slopes_all:
        #    axb.plot(xconcs,array(slopes)*ref_scale,linewidth=linewidth)
        _,_,_,ymax = beautify_plot(axb,x0min=True,y0min=True,\
            xticksposn='bottom',yticksposn='left',yticks=[])
        axb.set_xlim(0,5)
        axb.set_ylim(0,5)
        axb.set_xticks([0,1,5])
        axb_twin.set_ylim(0,180)
        axb_twin.set_yticks([])
        ## Draw the twin y axis (turned off always by beautify_plot)
        for loc, spine in axb.spines.items(): # items() returns [(key,value),...]
            spine.set_linewidth(axes_linewidth)
            if loc in ['right']:
                spine.set_color('k') # draw spine in black

        ## x-axis label for full second row
        if inhi==2:
            axes_labels(axb,"ORN scaling","",xpad=2)
            #axb.xaxis.set_label_coords(1.3,-0.25)
        if inhi==4: # laebl twin y axes of rightmost plot
            axb_twin.set_yticks([0,180])
            axes_labels(axb_twin,'','peak (Hz)',ypad=-6) # to set font size for twin y ticklabels

    for i,(_,_,_,_,_,_,(axa,axb)) in enumerate(inh_options):
        axa.set_ylim(0,maxy_hist)
        if i==0: # label y axis for left-most plots
            axa.set_yticks([0,maxy_hist])
            axb.set_yticks([0,1,5])
            axes_labels(axa,'','density',ypad=2)
            axes_labels(axb,'','mitral scaling',ypad=2)
    fig1.tight_layout()
    fig_clip_off(fig1)
    fig1.subplots_adjust(top=0.95,left=0.1,bottom=0.15,right=0.91,wspace=0.4,hspace=0.4)
    #fig1.text(0.31,0.65,'density',fontsize=label_fontsize,\
    #        rotation='vertical', transform=fig1.transFigure)
    #fig1.text(0.6,0.025,'$R^2$',fontsize=label_fontsize,transform=fig1.transFigure)
    fig1.savefig('../figures/lin_contribs_scaledpulses.svg',dpi=fig1.dpi)
    fig1.savefig('../figures/lin_contribs_scaledpulses.png',dpi=fig1.dpi)

if __name__ == "__main__":
    ## plot linearity with parts of the network removed / tweaked
    plot_lin_contribs_paperfigure()
    
    show()
