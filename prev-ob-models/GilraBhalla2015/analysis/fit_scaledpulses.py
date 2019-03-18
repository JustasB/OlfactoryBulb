# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
########## This variant does not use an air kernel, rather a constant air rate.
## USAGE1: python2.6 fit_scaledpulses.py <../results/odor_pulses/scaledpulses_ ... .pickle> <stimseed>

from scipy import optimize
from scipy import special
import scipy.interpolate
from scipy import signal
import scipy.io
from scipy.stats import linregress
from pylab import *
import pickle
import sys
import math
import copy as cp

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, SCALED_RUNTIME
from networkConstants import * # has central_glom
from data_utils import * # has axes_labels()
from analysis_utils import * # has load_morph...(), predict_respresp() and set_min_err()

## index of the reference response
ref_response_scalenum = 3
## if peak detect is True, peak scaling is plotted,
## rather than the scaling wrt reference response above
peak_detect = False#True

## rebin the spikes as below, irrespective of previous binning
pulserebindt = 50e-3#fitting_dt # 50ms as Priyanka uses
pulserebins = int(SCALED_RUNTIME/pulserebindt)
#bin_width_time = 2*pulserebindt ## overlapping bins with this bin width
pulsetlist = arange(pulserebindt/2.0,SCALED_RUNTIME,pulserebindt)

NOISE_ANALYSIS = False

## number of the mitral to be fitted.
fitted_mitrals = [2*central_glom,2*central_glom+1]

class fit_plot_scaledpulses():
    def __init__(self,filename,stimseed,scaledWidth=scaledWidth):
        self.filename = filename
        self.stimseed = stimseed
        self.scaledWidth = scaledWidth

    def fit_pulses(self,dirextn='',_noshow=True,_savefig=False, _test=False):
        if _test:
            ########## load in the stimuli
            ## scaledPulseList[glomnum][scalenum][odornum][binnum]
            scaledPulsesList,genORNkernels \
                = read_scaledpulses_stimuli(self.stimseed,self.scaledWidth)
            ## decimate pulseList by taking every convolutiondt/FIRINGFILLDT bins
            ## and decimate ORNkernels by taking every kerneldt/FIRINGFILLDT=50/1=50 bins
            ## pulserebinList[scalenum,binnum] numpy array
            pulserebinList,linORNkernelA,linORNkernelB = \
                decimate(scaledPulsesList[0,:,0],pulserebindt,genORNkernels,kerneldt,kernel_size)       
            numtrials = 1
        else:
            ######### load in the mitral pulse responses
            #### mitral_responses_list[avgnum][scalenum][mitralnum][spikenum]
            #### mitral_responses_binned_list[avgnum][scalenum][mitralnum][binnum]
            mitral_responses_list, mitral_responses_binned_list = read_pulsefile(self.filename)

            if NOISE_ANALYSIS:
                ##---------- Print/plot for best bin size, but this bin size is not used presently
                ## From Neural Computation 19, 1503–1527 (2007) Shimazaki & Shinomoto, pg 6
                DeltaList = arange(25.0e-3,500.0e-3,25.0e-3)
                SSvals = []
                SSvals_extended = []
                for Delta in DeltaList:
                    numtrials,mitral_responses_mean1,mitral_responses_std1 = \
                            rebin_mean(mitral_responses_list,Delta,SCALED_RUNTIME)
                    fitted_mitral = 0
                    firingbinsmeanList1 = mitral_responses_mean1[:,fitted_mitral]
                    firingbinserrList1 = mitral_responses_std1[:,fitted_mitral]
                    ## each val in firingbinsmeanList[scalenum][binnum] is k_i/(n*Delta)
                    kbar = mean(firingbinsmeanList1)*numtrials*Delta
                    kvariance = (std(firingbinsmeanList1)*numtrials*Delta)**2
                    shimazaki_shinomoto_val = (2*kbar-kvariance)/((numtrials*Delta)**2)
                    SSvals.append(shimazaki_shinomoto_val)
                    print 'Delta t =',Delta,"shimazaki_shinomoto_val =",shimazaki_shinomoto_val
                    ## Cost function i.e. shimazaki_shinomoto_val for different # of trials
                    Priyanka_numtrials = 12
                    ssval_extended = (1.0/Priyanka_numtrials - 1.0/numtrials)*kbar/numtrials/Delta**2 + \
                                        shimazaki_shinomoto_val
                    SSvals_extended.append(ssval_extended)
                if not _noshow:
                    fig = figure()
                    ax = fig.add_subplot(111)
                    ax.plot(DeltaList,SSvals,label='#='+str(numtrials))
                    ax.plot(DeltaList,SSvals_extended,label='#='+str(Priyanka_numtrials))
                    ax.legend()
                    axes_labels(ax,'Delta t (s)','SS Cost function (Hz^2)',fontsize=20)

            ##-------------------------- rebin the responses and pulses ------------------------------
            ## rebin sim responses to pulserebindt=50ms, then take mean
            numtrials,mitral_responses_mean,mitral_responses_std = \
                    rebin_mean(mitral_responses_list,pulserebindt,SCALED_RUNTIME)

        ## full fitting data for both mitrals
        fits_2mits = []
        peak_scales_2mits = []
        for mit_i,fitted_mitral in enumerate(fitted_mitrals):
            if _test:
                firingbinsmeanList = pulserebinList
                firingbinsmeanList += uniform(-0,0,shape(firingbinsmeanList))
                firingbinserrList = zeros(shape(firingbinsmeanList))
            else:
                ## take the odor responses of the mitral to be fitted
                firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
                ## The model predicts the individual response not the mean.
                ## Hence below fitting uses standard deviation, not standard error of the mean.
                firingbinserrList = mitral_responses_std[:,fitted_mitral]
            
            starti = int(PULSE_START/pulserebindt)
            endi = int((PULSE_START+scaledWidth+kernel_time)/pulserebindt)
            air_bgnd = firingbinsmeanList[0]
            air_bgnd_relevant = firingbinsmeanList[0][starti:endi]

            ##---------------------------- fit scaled pulse responses ------------------------------------------
                            
            ## define the reference response / scaling
            ref_scale = scaledList[ref_response_scalenum]
            ref_response = firingbinsmeanList[ref_response_scalenum][starti:endi]-air_bgnd_relevant
            fits = []
            peak_scales = []
            for scalenum in [1,2,3,4,5]: ## conc scaled pulses
                scale = scaledList[scalenum]
                response = firingbinsmeanList[scalenum][starti:endi]-air_bgnd_relevant
                ## http://www.jerrydallal.com/LHSP/slrout.htm for defn of std error of the estimate
                ## SEE is std error of data about the regression line
                slope, intercept, r_value, p_value, see = linregress(ref_response,response)
                ## SEE is called \hat{\sigma}_\epsilon i.e. sqrt(MSE) here:
                ## http://en.wikipedia.org/wiki/Regression_analysis
                ## formula for std error of slope is also from above Wikipedia article
                se_slope = see/std(ref_response)
                avgfrate = sum(firingbinsmeanList[scalenum][starti:endi])/float(endi-starti)
                fits.append((scale,slope,intercept,r_value,p_value,see,se_slope,avgfrate))
                ## peak scaling
                peak_scales.append(max(response))
            peak_scales =  array(peak_scales)#/peak_scales[ref_response_scalenum-1]*scaledList[ref_response_scalenum]
            peak_scales_2mits.append(peak_scales)
            fits_2mits.append(fits)

            ##---------------------------- plot scaled pulse responses -----------------------------------------

            if not _noshow:
                if peak_detect: print "BEWARE. Using peak scaling"
                else: print "BEWARE. Using fitted scaling"
                ############################### plot the responses and the fits
                fig = figure(figsize=(columnwidth,linfig_height/2),dpi=300,facecolor='w') # 'none' is transparent
                ## conc scaled pulses, leave the 0th pulse which is air_bgnd
                ax = plt.subplot2grid((1,3),(0,0),rowspan=1)
                for scaleiter,scale in enumerate(scaledList[1:]): ## conc scaled pulses
                    sister_ratio = (fitted_mitral%MIT_SISTERS)/float(MIT_SISTERS)
                    scaledpulsetime = array(pulsetlist[starti:endi]) - pulsetlist[starti] # start from t=0
                    ################### Plot the simulated responses
                    ## smooth the simulated response
                    ## windowsize=5 and SD=0.65 are defaults from matlab's smoothts() for gaussian smoothing
                    Gwindow = signal.gaussian(5,0.65)
                    ## help from http://www.scipy.org/Cookbook/SignalSmooth
                    simresponse = convolve(Gwindow/Gwindow.sum(),\
                        firingbinsmeanList[scaleiter+1]-air_bgnd,mode='same')
                    ## ditch the smoothing above for scaled pulses
                    simresponse = firingbinsmeanList[scaleiter+1][starti:endi]-air_bgnd_relevant
                    ## numpy array, hence adds element by element
                    scale_color = ['r','b','g','m','c'][scaleiter]
                    fill_between(scaledpulsetime,
                        simresponse+firingbinserrList[scaleiter+1][starti:endi]/sqrt(numtrials),
                        simresponse-firingbinserrList[scaleiter+1][starti:endi]/sqrt(numtrials),
                        color=scale_color,alpha=0.3)
                    plot(scaledpulsetime,simresponse,linewidth=plot_linewidth,color=scale_color)
                xmin,xmax,ymin,ymax = beautify_plot(ax,x0min=True,y0min=False,xticksposn='bottom',yticksposn='left')
                ax.set_xticks([0,xmax])
                if ymin<10: ax.set_yticks([ymin,0,ymax])
                else: ax.set_yticks([ymin,ymax])
                plot([0,scaledWidth],[ymin+2,ymin+2],linewidth=plot_linewidth*3,color='r')
                axes_labels(ax,'s','Hz',adjustpos=False,xpad=0,ypad=-3)
                #ax.yaxis.set_label_coords(-0.4,1.2)

                ## plot the scaling
                ax = plt.subplot2grid((1,3),(0,1),rowspan=1)
                minx = min(ref_response)
                maxx = max(ref_response)
                xlist = arange(minx,maxx,(maxx-minx)/50.0) # this is in Hz
                for scalenum in [1,2,3,4,5]: ## conc scaled pulses
                    scale = scaledList[scalenum]
                    response = firingbinsmeanList[scalenum][starti:endi]-air_bgnd_relevant
                    scale_normed,slope,intercept,r_value,_,_,_,_ = fits[scalenum-1]
                    print "stimseed =",self.stimseed,", scale =",scale_normed,"r^2 = ",r_value**2
                    if scalenum != ref_response_scalenum:
                        color4scale = ['r','b','g','m','c'][scalenum-1]
                        marker4scale = ['s','+','d','x','.'][scalenum-1]
                        scatter(ref_response,response,s=marker_size,marker=marker4scale,\
                            color=color4scale,edgecolor=color4scale)
                        ylist = slope*xlist+intercept
                        plot(xlist,ylist,color=color4scale,linewidth=linewidth)
                xmin,xmax,ymin,ymax = beautify_plot(ax,x0min=False,y0min=False,\
                    xticksposn='bottom',yticksposn='left')
                ax.set_xticks([xmin,0,xmax])
                ax.set_yticks([ymin,0,ymax])
                axes_labels(ax,'Hz','',adjustpos=False,xpad=1)
                #ax.yaxis.set_label_coords(-0.5,1.2)
                    
                ## plot response scaling vs conc scaling
                ax = plt.subplot2grid((1,3),(0,2),rowspan=1)
                plot(range(6),range(6),color=(0,0,0.7,0.5),dashes=(2.0,1.0)) ## linear reference
                concratios,slopevsconc,_,_,_,_,se_slope,avg_frates = zip(*fits)
                print "Average firing rates for mitral",fitted_mitral,\
                    "for different scales is",avg_frates
                errorbar(append([0],concratios),y=append([0],array(slopevsconc))*ref_scale,\
                    yerr=append([0],array(se_slope))*ref_scale,\
                    color='b',linewidth=linewidth,marker='o',ms=marker_size,dashes=(2.0,1.0))
                ax_twin = ax.twinx()
                ax_twin.plot(scaledList,append([0],peak_scales),color='k',linewidth=linewidth,\
                    marker='o',ms=marker_size,dashes=(0.5,1.0))
                beautify_plot(ax,x0min=False,y0min=False,xticksposn='bottom',yticksposn='left')
                ## Draw the twin y axis (turned off always by beautify_plot)
                for loc, spine in ax.spines.items(): # items() returns [(key,value),...]
                    spine.set_linewidth(axes_linewidth)
                    if loc in ['right']:
                        spine.set_color('k') # draw spine in black

                ax.set_ylim(0,5)
                ax.set_xticks([0,1,5])
                ax.set_yticks([0,1,5])
                ax_twin.set_ylim(0,80)
                ax_twin.set_yticks([0,80])
                axes_labels(ax,'ORN scaling','mitral scaling',adjustpos=False,xpad=1,ypad=0)
                axes_labels(ax_twin,'','mitral peak',adjustpos=False,ypad=-1)
                #ax.yaxis.set_label_coords(-0.3,1.2)

                fig.tight_layout()
                fig_clip_off(fig)
                fig.subplots_adjust(top=0.94,left=0.1,right=0.91,hspace=0.4,wspace=0.5)

                if _savefig:
                    fig.savefig('../figures/scalelinearity_example_'+str(self.stimseed)+\
                        '_mit'+str(mit_i)+'.svg',dpi=fig.dpi)
                    fig.savefig('../figures/scalelinearity_example_'+str(self.stimseed)+\
                        '_mit'+str(mit_i)+'.png',dpi=fig.dpi)

                if NOISE_ANALYSIS:
                    ## plot the variance vs firing rate mean for each mitral
                    ## variance = mean/bintime of firng rate for Poisson process
                    fig2 = figure()
                    ax2 = fig2.add_subplot(111)
                    for scaleiter,scale in enumerate(scaledList): ## conc scaled pulses including air
                        ax2.scatter( firingbinsmeanList[scaleiter], \
                            firingbinserrList[scaleiter]**2, color='r' )
                    beautify_plot(ax2)
                    axes_labels(ax2,'mean rate (Hz)','variance (Hz^2)',fontsize=14)
                    
                    ## plot individual trials for a given response
                    fig3 = figure()
                    ax3 = fig3.add_subplot(111)
                    for trialspikelist in mitral_responses_list: 
                        plot( plotBins( trialspikelist[0][fitted_mitral],\
                            pulserebins, SCALED_RUNTIME, 0.0) )

        return fits_2mits, peak_scales_2mits
        

if __name__ == "__main__":
    NOSHOW = False
    SAVEFIG = True#False
    if len(sys.argv) > 2:
        filename = sys.argv[1]
        stimseed = sys.argv[2]
        worker = fit_plot_scaledpulses(filename,stimseed)
        post_pulses = filename.split('odor_pulses')[1]
        dirextn = post_pulses.split('/')[0]
        print 'directory extension is',dirextn
        if 'TEST' in sys.argv: TEST=True
        else: TEST=False
        worker.fit_pulses(dirextn,NOSHOW,SAVEFIG,TEST)
        show()
    else:
        print "At least specify data file containing pickled mitral responses, and ORN frate seed."
        sys.exit(1)
