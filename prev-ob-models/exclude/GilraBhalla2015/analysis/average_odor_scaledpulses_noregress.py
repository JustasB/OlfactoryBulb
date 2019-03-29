# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: python2.6 average_odor_scaledpulses_noregress.py

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
from scipy import interpolate # for interp1d
from scipy import optimize # for fsolve
from sim_utils import * # has rebin() and imports data_utils.py for axes_off()
from analysis_utils import * # has rebin_mean() and read_pulsefile()
from calc_corrs import * # has calc_corrs()
## BE CAREFUL: use _constair_separateodors below, others are obsolete.
import fit_scaledpulses_noregress as fit_sp # has fit_plot_scaledpulses and ref pulse constants
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

def plot_lin_contribs_oldpaperfigure():
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
            axes_labels(axa,"correlation with 1x","",xpad=2)
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
        #axb.plot(range(6),range(6),color=(0,0,0.7,0.5),dashes=(2.0,1.0)) # linear reference
        axb.errorbar(x=append([0],xconcs),y=peaks_mean,yerr=peaks_std,\
            color='k',linewidth=linewidth,capsize=cap_size)
        #for slopes in slopes_all:
        #    axb.plot(xconcs,array(slopes)*ref_scale,linewidth=linewidth)
        _,_,_,ymax = beautify_plot(axb,x0min=True,y0min=True,\
            xticksposn='bottom',yticksposn='left',yticks=[])
        axb.set_xlim(0,5)
        axb.set_ylim(0,180)
        axb.set_xticks([0,1,5])
        axb.set_yticks([])
        ### Draw the twin y axis (turned off always by beautify_plot)
        #for loc, spine in axb.spines.items(): # items() returns [(key,value),...]
        #    spine.set_linewidth(axes_linewidth)
        #    if loc in ['right']:
        #        spine.set_color('k') # draw spine in black

        ## x-axis label for full second row
        if inhi==2:
            axes_labels(axb,"ORN scaling","",xpad=2)
            #axb.xaxis.set_label_coords(1.3,-0.25)
        if inhi==0: # laebl twin y axes of rightmost plot
            axb.set_yticks([0,180])
            axes_labels(axb,'','peak (Hz)',ypad=-6) # to set font size for twin y ticklabels

    for i,(_,_,_,_,_,_,(axa,axb)) in enumerate(inh_options):
        axa.set_ylim(0,maxy_hist)
        if i==0: # label y axis for left-most plots
            axa.set_yticks([0,maxy_hist])
            axes_labels(axa,'','density',ypad=2)
    fig1.tight_layout()
    fig_clip_off(fig1)
    fig1.subplots_adjust(top=0.95,left=0.1,bottom=0.15,right=0.98,wspace=0.4,hspace=0.4)
    #fig1.text(0.31,0.65,'density',fontsize=label_fontsize,\
    #        rotation='vertical', transform=fig1.transFigure)
    #fig1.text(0.6,0.025,'$R^2$',fontsize=label_fontsize,transform=fig1.transFigure)
    fig1.savefig('../figures/lin_contribs_scaledpulses.svg',dpi=fig1.dpi)
    fig1.savefig('../figures/lin_contribs_scaledpulses.png',dpi=fig1.dpi)

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
    maxy_R = 0.0
    for ploti,(dirextn,inhi,inh,inhstr,nl_orns,nl_type,(axa,axb)) in enumerate(inh_options):
        slopes_all = []
        peaks_all = []
        avg_frates_all = []
        R_thisinh = []
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
                R_thisinh.append(r_values)
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
        R_mean = mean(R_thisinh,axis=0)
        R_std = std(R_thisinh,axis=0)
        xconcs = scaledList[1:]
        axa.errorbar(xconcs,R_mean,R_std,linewidth=linewidth,color='b',\
            ls='solid',marker='o',ms=marker_size,capsize=cap_size)
        #axa.plot(xconcs,R_std,linewidth=linewidth,color='r',ls='solid',marker='x',ms=marker_size)
        xmin,xmax,ymin,ymax = \
            beautify_plot(axa,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left',yticks=[])
        axa.set_ylim(0,ymax) # later this is set to the max over all histograms maxy_hist
        axa.set_xticks([0,1,2,5])
        maxy_R = max(maxy_R,ymax)
        ## x-axis label for full first row
        if inhi==2:
            axes_labels(axa,"","",xpad=2) # ORN scaling is on the next row x-axis
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
        #axb.plot(range(6),range(6),color=(0,0,0.7,0.5),dashes=(2.0,1.0)) # linear reference
        axb.errorbar(x=append([0],xconcs),y=peaks_mean,yerr=peaks_std,\
            color='b',linewidth=linewidth,capsize=cap_size)
        #for slopes in slopes_all:
        #    axb.plot(xconcs,array(slopes)*ref_scale,linewidth=linewidth)
        _,_,_,ymax = beautify_plot(axb,x0min=True,y0min=True,\
            xticksposn='bottom',yticksposn='left',yticks=[])
        axb.set_xlim(0,5)
        axb.set_ylim(0,180)
        axb.set_xticks([0,1,2,5])
        axb.set_yticks([])
        ### Draw the twin y axis (turned off always by beautify_plot)
        #for loc, spine in axb.spines.items(): # items() returns [(key,value),...]
        #    spine.set_linewidth(axes_linewidth)
        #    if loc in ['right']:
        #        spine.set_color('k') # draw spine in black

        ## x-axis label for full second row
        if inhi==2:
            axes_labels(axb,"conc (% SV)","",xpad=2)
            #axb.xaxis.set_label_coords(1.3,-0.25)
        if inhi==0: # laebl twin y axes of rightmost plot
            axb.set_yticks([0,180])
            axes_labels(axb,'','peak (Hz)',ypad=-6) # to set font size for twin y ticklabels

    for i,(_,_,_,_,_,_,(axa,axb)) in enumerate(inh_options):
        axa.set_ylim(0,1)
        if i==0: # label y axis for left-most plots
            axa.set_yticks([0,1])
            axes_labels(axa,'','corr to 1%',ypad=2)
    fig1.tight_layout()
    fig_clip_off(fig1)
    fig1.subplots_adjust(top=0.95,left=0.1,bottom=0.15,right=0.98,wspace=0.4,hspace=0.4)
    #fig1.text(0.31,0.65,'density',fontsize=label_fontsize,\
    #        rotation='vertical', transform=fig1.transFigure)
    #fig1.text(0.6,0.025,'$R^2$',fontsize=label_fontsize,transform=fig1.transFigure)
    fig1.savefig('../figures/lin_contribs_scaledpulses.svg',dpi=fig1.dpi)
    fig1.savefig('../figures/lin_contribs_scaledpulses.png',dpi=fig1.dpi)

def plot_responseshapes_vs_conc():
    fig = figure(figsize=(columnwidth/3.*2,linfig_height/2.),dpi=300,facecolor='w')
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    xconcs = array(scaledList[1:])
    ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_option = ('',0,(False,False,False,False,False),'all',False,None)
    tpeaks_all = []
    tpeaks_normed_all = []
    durations_normed_all = []
    risetimes_all = []
    durations_all = []
    deltapeaks_all = []
    latencies_all = []
    latencies_normed_all = []
    (dirextn,inhi,inh,inhstr,nl_orns,nl_type) = inh_option
    for stimi,stimseed in enumerate(stim_seeds):
        filename, switch_strs \
            = get_filename(stimseed,stimseed,inh,0,nl_orns,nl_type,\
                resultsdir='../results/odor_pulses'+dirextn)
        ## if the result file for these seeds & tweaks doesn't exist,
        ## then carry on to the next.
        if not os.path.exists(filename): continue
        print filename
        ######### load in the mitral pulse responses
        #### mitral_responses_list[avgnum][scalenum][mitralnum][spikenum]
        #### mitral_responses_binned_list[avgnum][scalenum][mitralnum][binnum]
        mitral_responses_list, mitral_responses_binned_list = read_pulsefile(filename)
        ##-------------------------- rebin the responses and pulses ------------------------------
        ## rebin sim responses to pulserebindt=50ms, then take mean
        numtrials,mitral_responses_mean,mitral_responses_std = \
                rebin_mean(mitral_responses_list,fit_sp.pulserebindt,SCALED_RUNTIME)
        for mitrali in [0,1]:
            ## take the odor responses of the mitral to be fitted
            firingbinsmeanList = mitral_responses_mean[:,mitrali]
            
            starti = int(PULSE_START/fit_sp.pulserebindt)
            endi = int((PULSE_START+scaledWidth+kernel_time)/fit_sp.pulserebindt)
            air_bgnd = firingbinsmeanList[0]
            air_bgnd_relevant = firingbinsmeanList[0][starti:endi]
            tpeaks = []
            risetimes = []
            durations = []
            deltapeaks = []
            latencies = []
            for scalenum in [1,2,3,4,5]: ## conc scaled pulses
                ## find the mean first spike latency after odor onset
                latency = 0
                for avgi in range(len(mitral_responses_list)):
                    for spiketime in mitral_responses_list[avgi][scalenum][mitrali]:
                        if spiketime > PULSE_START:
                            latency += spiketime - PULSE_START
                            break
                latency_mean = latency / float(len(mitral_responses_list))
                latencies.append(latency_mean)
                ## find peak times and durations of binned responses
                scale = scaledList[scalenum]
                response = firingbinsmeanList[scalenum][starti:endi]-air_bgnd_relevant
                tpeak = argmax(response)
                tpeaks.append(tpeak)
                ##--------- Calculate half-max duration of the response
                ## f is an interpolated 1D function that gives
                ## down-shifted response that is zero at halfmax for fsolve/bisect below
                f = interpolate.interp1d(range(len(response)),response-max(response)/2.,\
                        bounds_error=False,fill_value=0.0)
                ### fsolve to find half-max points on either side of tpeak
                ### but often gives negative values!! so ditch this
                #halfresp_low = optimize.fsolve(f,tpeak-1)[0]
                #halfresp_high = optimize.fsolve(f,tpeak+1)[0]
                ## use bisect to find halfmax pt within an interval
                ## but bisect need opposite signs at the ends of the interval
                ## else ValueError, so increase interval, till you get it
                tmin = tpeak-1
                while True:
                    try:
                        halfresp_low = optimize.bisect(f,tpeak,tmin-1)
                        break
                    except ValueError:
                        tmin -= 1
                tmax = tpeak+1
                while True:
                    try:
                        halfresp_high = optimize.bisect(f,tpeak,tmax)
                        break
                    except ValueError:
                        tmax += 1
                risetimes.append(halfresp_low)
                duration = (halfresp_high-halfresp_low)
                durations.append(duration)
                
            ## take only the reduction in the time to peak vs concentration
            tpeaks = array(tpeaks)*fit_sp.pulserebindt*1000
            risetimes = array(risetimes)*fit_sp.pulserebindt*1000
            durations = array(durations) *fit_sp.pulserebindt*1000
            tpeaks_all.append(tpeaks)
            tpeaks_normed_all.append(tpeaks/mean(tpeaks))
            risetimes_all.append(risetimes)
            durations_all.append(durations)
            durations_normed_all.append(durations/mean(durations))
            deltapeaks = [tpeaks[-1]-tpeaks[-2],tpeaks[-2]-tpeaks[-3]]
            deltapeaks_all.extend(deltapeaks)
            latencies_all.append(latencies)
            latencies_normed_all.append(latencies/mean(latencies))
            #ax1.plot(xconcs,tpeaks)
            #ax2.plot(xconcs,durations)
    
    tpeaks_mean = mean(tpeaks_all,axis=0)
    tpeaks_std = std(tpeaks_all,axis=0)
    tpeaks_normed_mean = mean(tpeaks_normed_all,axis=0)
    tpeaks_normed_se = std(tpeaks_normed_all,axis=0) / sqrt(len(tpeaks_normed_all))
    risetimes_mean = mean(risetimes_all,axis=0)
    risetimes_std = std(risetimes_all,axis=0)
    durations_mean = mean(durations_all,axis=0)
    durations_std = std(durations_all,axis=0)
    durations_normed_mean = mean(durations_normed_all,axis=0)
    durations_normed_se = std(durations_normed_all,axis=0) / sqrt(len(durations_normed_all))
    latencies_mean = mean(latencies_all,axis=0)
    latencies_std = std(latencies_all,axis=0)
    latencies_normed_mean = mean(latencies_normed_all,axis=0)
    latencies_normed_se = std(latencies_normed_all,axis=0) / sqrt(len(latencies_normed_all))
    ## plots
    ax1.errorbar(xconcs,latencies_normed_mean,latencies_normed_se,color='b',linewidth=linewidth)
    #ax1.errorbar(xconcs,risetimes_mean,risetimes_std,color='b',linewidth=linewidth)
    ax2.errorbar(xconcs,tpeaks_normed_mean,tpeaks_normed_se,color='b',linewidth=linewidth)
    #ax2.errorbar(xconcs,durations_normed_mean,durations_normed_se,color='b',linewidth=linewidth)
    for i,ax in enumerate([ax1,ax2]):
        _,_,ymin,ymax = beautify_plot(ax,x0min=True,y0min=False,\
            xticksposn='bottom',yticksposn='left')
        ax.set_xlim(0,5)
        ax.set_xticks([0,1,2,5])
        ax.set_yticks([ymin,1,ymax])
        axes_labels(ax,'conc (% SV)',['norm-t-spike','norm-t-peak'][i])
    ### Draw the twin y axis (turned off always by beautify_plot)
    #for loc, spine in axb.spines.items(): # items() returns [(key,value),...]
    #    spine.set_linewidth(axes_linewidth)
    #    if loc in ['right']:
    #        spine.set_color('k') # draw spine in black

    fig_clip_off(fig)
    fig.tight_layout()
    fig.savefig('../figures/scaledpulses_latency.svg',dpi=fig.dpi)
    fig.savefig('../figures/scaledpulses_latency.png',dpi=fig.dpi)
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.hist(deltapeaks_all)
    axes_labels(ax,'delta peak times','#')

if __name__ == "__main__":
    ## plot linearity with parts of the network removed / tweaked
    #plot_lin_contribs_oldpaperfigure()
    ## PAPER figure. Similar to above, but plotting mean and sd of corr with 1x
    ## versus concentration instead of a combined distrib.
    plot_lin_contribs_paperfigure()
    ## Compare pulse shapes versus concentration
    plot_responseshapes_vs_conc()
    
    show()
