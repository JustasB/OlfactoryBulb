# -*- coding: utf-8 -*-

import math, os
import pickle

from pylab import *

## USAGE: python2.6 plot_corrs.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from networkConstants import * # has central_glom
from data_utils import *
from calc_corrs import *

from results_catalogue import *

fig1 = figure(facecolor='none') # 'none' is transparent
## A super axes to set common x and y axes labels
bigAxes1 = fig1.add_axes([0.1,0.1,0.8,0.8],frameon=False) # hide frame
#bigAxes1.set_xticks([0,1,2,3])
#bigAxes1.set_xticklabels(['none', 'singles',\
#    's+joints', 's+j+PGs'],fontsize=20)
bigAxes1.set_xticks([])
bigAxes1.set_yticks([])
bigAxes1.text(-0.1,0.3,'spiking probability', fontsize=24, rotation='vertical')
bigAxes1.text(-0.1,-0.11,'inhibition: none, singles, s+joints, s+j+PGs',\
    fontsize=24, rotation='horizontal')
bigAxes1.set_title('Cross-Correlogram peak',fontsize=36)

fig2 = figure(facecolor='none') # 'none' is transparent
## A super axes to set common x and y axes labels
bigAxes2 = fig2.add_axes([0.1,0.1,0.8,0.8],frameon=False) # hide frame
#bigAxes2.set_xticks([0,1,2,3])
#bigAxes2.set_xticklabels(['none', 'singles',\
#    's+joints', 's+j+PGs'],fontsize=20)
bigAxes2.set_xticks([])
bigAxes2.set_yticks([])
bigAxes2.text(-0.1,-0.11,'inhibition: none, singles, s+joints, s+j+PGs',\
    fontsize=24, rotation='horizontal')
bigAxes2.set_title('Binned Cross-Correlation peak',fontsize=36)

fig3 = figure(facecolor='none') # 'none' is transparent
## A super axes to set common x and y axes labels
bigAxes3 = fig3.add_axes([0.1,0.1,0.8,0.8],frameon=False) # hide frame
#bigAxes3.set_xticks([0,1,2,3])
#bigAxes3.set_xticklabels(['none', 'singles',\
#    's+joints', 's+j+PGs'],fontsize=20)
bigAxes3.set_xticks([])
bigAxes3.set_yticks([])
bigAxes3.text(-0.1,0.3,'time(s)', fontsize=24, rotation='vertical')
bigAxes3.text(-0.1,-0.11,'inhibition: none, singles, s+joints, s+j+PGs',\
    fontsize=24, rotation='horizontal')
bigAxes3.set_title('Cross-Correlogram shift',fontsize=36)

plotnum = 1
for runnum in range(len(filelist[0])):
    corrgramlist = []
    corrgrammaxlist = []
    corrgrampeaklist = []
    corrlist = []
    print seeds[runnum]
    for filename in array(filelist)[:,runnum]:
        filenamefull = '../results/odor_morphs/'+filename
        (air_corr,odorA_corr,odorB_corr),\
            (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram)\
            = calc_corrs(filenamefull, norm_str='overall')
        corrgramlist.append([airxcorrgram,odorAxcorrgram,odorBxcorrgram])
        maxes = [max(airxcorrgram),max(odorAxcorrgram),max(odorBxcorrgram)]
        peaks = [ where(airxcorrgram==maxes[0])[0],\
            where(odorAxcorrgram==maxes[1])[0],
            where(odorBxcorrgram==maxes[2])[0] ]
        if len(peaks[0])>1 or len(peaks[1])>1 or len(peaks[2])>1:
            print "Multiple peaks for",filenamefull,peaks
        corrgrammaxlist.append(maxes)
        ## binwidth = 1e-3 for xcorrgram, convert bin index to time shift
        peaks = ( array([mean(peaks[0]),mean(peaks[1]),mean(peaks[2])]) - 1 ) * 1e-3 - 0.5
        corrgrampeaklist.append(peaks)
        corrlist.append([air_corr,odorA_corr,odorB_corr])
    corrlist = array(corrlist)
    corrgramlist = array(corrgramlist)
    corrgrammaxlist = array(corrgrammaxlist)
    corrgrampeaklist = array(corrgrampeaklist)
    
    ## hardcoded - I know there are 4 files
    ax1 = fig1.add_subplot(2,3,plotnum)
    ax1.set_xticks([])
    ax1.plot(corrgrammaxlist[:,0],color=(0,0,0),linewidth=2)
    ax1.plot(corrgrammaxlist[:,1],color=(1,0,0),linewidth=2)
    ax1.plot(corrgrammaxlist[:,2],color=(0,1,0),linewidth=2)
    #ymax = ax1.get_ylim()[1]+0.01
    ## very important to give it after the plot functions else autoscales
    #ax1.set_ylim(-0.01,ymax)
    #ax1.set_yticks([0,ymax])
    #ax1.set_yticklabels(['0','%0.2f'%ymax],size='large')
    ## or set autoscaling off - gca().set_autoscale_on(False)
    ## - taken from http://old.nabble.com/ylim-does-not-work-td19000814.html
    ax1.set_title( str(seeds[runnum]), size='large')

    ## hardcoded - I know there are 4 files
    ax3 = fig3.add_subplot(2,3,plotnum)
    ax3.set_xticks([])
    #ax3.plot(corrgrampeaklist[:,0],color=(0,0,0),linewidth=2)
    #ax3.plot(corrgrampeaklist[:,1],color=(1,0,0),linewidth=2)
    #ax3.plot(corrgrampeaklist[:,2],color=(0,1,0),linewidth=2)
    ax3.plot(corrgramlist[3,0],color=(0,0,0),linewidth=2)
    ax3.plot(corrgramlist[3,2],color=(0,1,0),linewidth=2)
    ax3.plot(corrgramlist[3,1],color=(1,0,0),linewidth=2)
    ax3.set_title( str(seeds[runnum]), size='large')

    ## hardcoded - I know there are 4 files
    ax2 = fig2.add_subplot(2,3,plotnum)
    ax2.set_xticks([])
    ax2.plot(corrlist[:,0],color=(0,0,0),linewidth=2)
    ax2.plot(corrlist[:,1],color=(1,0,0),linewidth=2)
    ax2.plot(corrlist[:,2],color=(0,1,0),linewidth=2)
    #ymax = ax1.get_ylim()[1]+0.01
    ## very important to give it after the plot functions else autoscales
    #ax2.set_ylim(-0.01,ymax)
    #ax2.set_yticks([0,ymax])
    #ax2.set_yticklabels(['0','%0.2f'%ymax],size='large')
    ## or set autoscaling off - gca().set_autoscale_on(False)
    ## - taken from http://old.nabble.com/ylim-does-not-work-td19000814.html
    ax2.set_title( str(seeds[runnum]), size='large')

    plotnum += 1

show()
