#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pickle
sys.path.extend(["..","../networks","../generators","../simulations","../analysis"])

from OBNetwork import *
from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has REALRUNTIME, NUMBINS
from sim_utils import * # has rebin() and imports data_utils.py for axes_off()
from data_utils import * # has mpi import and variables also

plot_images=False#True ## plot Adil style images?
## below requires scipy which requires lapack (present only on gj not on nodes)
## hence import only if not running in parallel
if plot_images and mpisize == 1:
    from fit_odor_morphs import * # has fit_morphs()

from pylab import * # part of matplotlib that depends on numpy but not scipy

PLOTRESP_NUM = 2 # whether to plot 2 respiration cycles or 1
BINS_PER_RESP = 17
NUMBINS = BINS_PER_RESP*PLOTRESP_NUM
BIN_WIDTH_TIME = RESPIRATION/float(BINS_PER_RESP)
bindt = RESPIRATION/float(NUMBINS)
fitted_mitral = 2*central_glom+0
## I take the last PLOTRESP_NUM of respiration cycles out of NUM_RESPS simulated
startplottime = SETTLETIME+(NUM_RESPS-PLOTRESP_NUM)*RESPIRATION
endplottime = SETTLETIME+NUM_RESPS*RESPIRATION
responsetlist = arange( startplottime+bindt/2.0, endplottime, RESPIRATION/NUMBINS)

def plot_morph_responses(mitral_responses_avg,mitral_responses_se):

    numodors = len(inputList)
    if ONLY_TWO_MITS: mitlist = range(MIT_SISTERS)
    else: mitlist = range(1)#range(6)#range(NUM_GLOMS*MIT_SISTERS)
    for mitnum in mitlist:
        if mitnum%MIT_SISTERS == 0:
            figure()
            title('Glomerulus '+str(mitnum/MIT_SISTERS))
        for odornum in range(numodors):
            odorA,odorB = inputList[odornum]
            sister_ratio = (mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            errorbar(x=responsetlist,y=mitral_responses_avg[odornum,mitnum],\
                yerr=mitral_responses_se[odornum,mitnum],color=(odorA,odorB,sister_ratio))

    #### plot to compare air decorr
    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    title('Air: mit0 (red) vs mit1 (blue)',fontsize=36)
    errorbar(x=responsetlist,y=mitral_responses_avg[6,0],\
        yerr=mitral_responses_se[6,0],color='r',linewidth=2)
    errorbar(x=responsetlist,y=mitral_responses_avg[6,1],\
        yerr=mitral_responses_se[6,1],color='b',linewidth=2)
    axes_labels(ax,'time (s)','firing rate (Hz)',adjustpos=True)

    #### plot to compare odor decorr
    glomnum = central_glom
    ## 1/3rd column 85mm/3/25.4 inches wide
    fig = figure(figsize=(columnwidth/3,linfig_height),dpi=300,facecolor='w')
    ax = fig.add_subplot(2,1,1,frameon=False)
    #title('OdorA: mit0 (red) vs mit1 (blue)',fontsize=36)
    simresponse = mitral_responses_avg[5,2*glomnum+0]
    simerr = mitral_responses_se[5,2*glomnum+0]
    fill_between(responsetlist, simresponse+simerr, simresponse-simerr,
        color='r',alpha=0.4)
    plot(responsetlist,simresponse,color='r',linewidth=linewidth)
    simresponse = mitral_responses_avg[5,2*glomnum+1]
    simerr = mitral_responses_se[5,2*glomnum+1]
    fill_between(responsetlist, simresponse+simerr, simresponse-simerr,
        color='b',alpha=0.4)
    plot(responsetlist,simresponse,color='b',linewidth=linewidth)
    ax.set_xlim(startplottime,endplottime)
    ax.get_xaxis().set_ticks_position('none')
    ax.get_yaxis().set_ticks_position('left')
    ## need to add the axes lines that I want, after deleting full frame.
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.set_xticks([])
    ax.set_yticks([ymin,ymax])
    ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=linewidth))
    ## separate respirations
    ax.add_artist(Line2D(((xmax+xmin)/2.0, (xmax+xmin)/2.0), (ymin, ymax), color='black', linewidth=linewidth))
    axes_labels(ax,'','') # sets default fontsize too

    ax = fig.add_subplot(2,1,2,frameon=False)
    #title('OdorB: mit0 (red) vs mit1 (blue)',fontsize=36)
    simresponse = mitral_responses_avg[0,2*glomnum+0]
    simerr = mitral_responses_se[0,2*glomnum+0]
    fill_between(responsetlist, simresponse+simerr, simresponse-simerr,
        color='r',alpha=0.4)
    plot(responsetlist,simresponse,color='r',linewidth=linewidth)
    simresponse = mitral_responses_avg[0,2*glomnum+1]
    simerr = mitral_responses_se[0,2*glomnum+1]
    fill_between(responsetlist, simresponse+simerr, simresponse-simerr,
        color='b',alpha=0.4)
    plot(responsetlist,simresponse,color='b',linewidth=linewidth)
    ax.set_xlim(startplottime,endplottime)
    ax.get_xaxis().set_ticks_position('bottom')
    ax.get_yaxis().set_ticks_position('left')
    ## need to add the axes lines that I want, after deleting full frame.
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.set_xticks([xmin,xmax])
    ax.set_yticks([ymin,ymax])
    ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=linewidth))
    ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=linewidth))
    ## separate respirations
    ax.add_artist(Line2D(((xmax+xmin)/2.0, (xmax+xmin)/2.0), (ymin, ymax), color='black', linewidth=linewidth))
    axes_labels(ax,'time (s)','firing rate (Hz)')

    ## Below text position is wrt to last axes drawn.
    ## transform = ax.transAxes sets the test position as axes units and not data units.
    #text(-0.42,1.4,'firing rate(Hz)', fontsize=label_fontsize, rotation='vertical', transform = ax.transAxes)
    #text(0.20,-0.27,'time (s)', fontsize=label_fontsize, transform = ax.transAxes)
    fig.tight_layout()
    fig.savefig('../figures/sistermorphs.png',dpi=fig.dpi)
    fig.savefig('../figures/sistermorphs.svg',dpi=fig.dpi)
    
def plot_morph_images(filename,mitral_responses_avg,mitral_responses_binned_list):
    params,chisq,inputsA,inputsB,fitted_responses,numavgs,firingbinsmeanList,firingbinserrList\
        = fit_morphs(filename, fitted_mitral)

    fig = figure(facecolor='none') # the string 'none' sets transparent facecolor
    mitnum = fitted_mitral
    numodors = len(inputList)
    mitral_responses_binned_list = array(mitral_responses_binned_list)
    ## air image will be padded on both sides of odor response below
    ## array[::-1] is for reversing elements in an array
    air_image = [ trialresponse[::-1] \
            for trialresponse in mitral_responses_binned_list[:,numodors-1,mitnum,:] ]
    numtrials = len(air_image)
    for odornum in range(numodors-1):
        # odors are mapped to odornum: B<->0, A<->5; so reverse the index = realodornum to show A to B
        realodornum = numodors-2-odornum
        ############## One image showing individual trials
        ax = fig.add_axes([0.15*(odornum)+0.025,0.471,0.15,0.2]) # [left,bottom,width,height]
        ## block of air responses are padded to odor on both sides
        ## array[::-1] is for reversing elements in an array
        response_image = []
        ## one has to do response_image = [] and then extend(),
        ## instead of response_image=air_image.
        ## the latter doesn't create a new array but is like a pointer to air_image
        ## so that later statements keep extending air_image
        response_image.extend(air_image)
        response_image.extend( [ trialresponse[::-1] \
            for trialresponse in mitral_responses_binned_list[:,realodornum,mitnum,:] ] )
        response_image.extend(air_image)
        im = ax.imshow(transpose(response_image), cmap=cm.jet, norm=None)
        title(['A','0.8A+\n0.6B','0.6A+\n0.4B','0.4A+\n0.6B','0.6A+\n0.8B','B','air'][odornum], fontsize=28)
        axes_off(ax)
        #cb = colorbar(im)
        #cb.set_clim(vmin=0,vmax=50)
        im.set_clim(vmin=0,vmax=50)
        ############### One image showing average trials
        ax = fig.add_axes([0.15*(odornum)+0.025,0.33,0.15,0.2]) # [left,bottom,width,height]
        ## block of air responses are padded to odor on both sides
        ## array[::-1] is for reversing elements in an array
        response_image = []
        response_image.extend( [fitted_responses[-1][::-1]]*numtrials )
        response_image.extend( [fitted_responses[realodornum][::-1]]*numtrials )
        response_image.extend( [fitted_responses[-1][::-1]]*numtrials )
        im = ax.imshow(transpose(response_image), cmap=cm.jet, norm=None)
        #title(['A','0.8A+\n0.6B','0.6A+\n0.4B','0.4A+\n0.6B','0.6A+\n0.8B','B','air'][odornum])
        axes_off(ax)
        #cb = colorbar(im)
        #cb.set_clim(vmin=0,vmax=50)
        im.set_clim(vmin=0,vmax=50)
    ax = fig.add_axes([0.94,0.38,0.01,0.24]) # [left,bottom,width,height]
    #ax.set_aspect(100.0)
    #ax.set_position([0.85,0.38,0.01,0.24]) # [left,bottom,width,height]
    cb = fig.colorbar(im, cax=ax)
    cb.set_clim(vmin=0,vmax=50)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(18)

def plot_morphs(filename):
    f = open(filename,'r')
    (mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
    f.close()
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, numbins=PLOTRESP_NUM*NUMBINS,\
            bin_width_time=BIN_WIDTH_TIME, numresps=PLOTRESP_NUM)

    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## since I plot the mean response, I must plot standard error of the mean
    ## = standard deviation of a repeat / sqrt(num of repeats).
    ## NOTE: our plot depends on number of trials.
    mitral_responses_se = mitral_responses_std/sqrt(numavgs)

    if plot_images:
        #### Fit the morphs, then plot image plots of responses and fits a la Khan et al 2008.
        plot_morph_images(filename,mitral_responses_avg,mitral_responses_binned_list)
    else:
        #### Usual firing rate vs time plots of responses
        plot_morph_responses(mitral_responses_avg,mitral_responses_se)

    show()

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "You need to specify the morph responses pickle filename."
        sys.exit(1)

    plot_morphs(sys.argv[1])
