# -*- coding: utf-8 -*-

import math, os
import pickle

from pylab import *

## USAGE: python2.6 plot_chisquares.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from networkConstants import * # has central_glom
from data_utils import *
from fit_odor_morphs import *

from results_catalogue import *

fig1 = figure(facecolor='none') # 'none' is transparent
## A super axes to set common x and y axes labels
bigAxes1 = fig1.add_axes([0.1,0.1,0.8,0.8],frameon=False) # hide frame
#bigAxes1.set_xticks([0,1,2,3])
#bigAxes1.set_xticklabels(['none', 'singles',\
#    's+joints', 's+j+PGs'],fontsize=20)
bigAxes1.set_xticks([])
bigAxes1.set_yticks([])
bigAxes1.text(-0.1,0.3,'chisquare', fontsize=24, rotation='vertical')
bigAxes1.text(-0.1,-0.11,'inhibition: none, singles, s+joints, s+j+PGs',\
    fontsize=24, rotation='horizontal')
bigAxes1.set_title('Morph chisquare vs inhibition',fontsize=36)

plotnum = 1
for runnum in range(len(filelist[0])):
    chisqlist = []
    print seeds[runnum]
    ## loop over files with the same seeds but different inhibitory components
    for filename in array(filelist)[:,runnum]:
        filenamefull = '../results/odor_morphs/'+filename
        chisq_mitlist = []
        for fitted_mitral in fitted_mitral_list:
            if fit_type=='arb':
                paramsfilename = filenamefull+'_params'+str(fitted_mitral)
            elif fit_type=='lin':
                paramsfilename = filenamefull+'_paramslin'+str(fitted_mitral)
            print paramsfilename
            if not os.path.exists(paramsfilename):
                fit_morphs(filenamefull,fitted_mitral)
            f = open(paramsfilename,'r')
            params,chisq = pickle.load(f)
            f.close()
            chisq_mitlist.append(chisq)
        chisqlist.append(chisq_mitlist)

    chisqlist=array(chisqlist)
    ## hardcoded - I know there are 4 files for different inhibitions
    ax1 = fig1.add_subplot(2,3,plotnum)
    ax1.set_xticks([])
    ## plot for mit0 and mit1
    ax1.plot(chisqlist[:,0],color=(1,0,0),linewidth=2)
    ax1.plot(chisqlist[:,1],color=(0,0,1),linewidth=2)
    #ymax = ax1.get_ylim()[1]+0.01
    ## very important to give it after the plot functions else autoscales
    #ax1.set_ylim(-0.01,ymax)
    #ax1.set_yticks([0,ymax])
    #ax1.set_yticklabels(['0','%0.2f'%ymax],size='large')
    ## or set autoscaling off - gca().set_autoscale_on(False)
    ## - taken from http://old.nabble.com/ylim-does-not-work-td19000814.html
    ax1.set_title( str(seeds[runnum]), size='large')

    plotnum += 1

show()
