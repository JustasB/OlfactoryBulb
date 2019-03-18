# -*- coding: utf-8 -*-

import math, os
import pickle

from pylab import *

## USAGE: python2.6 plot_morph_weights.py
## You must run fit_odor_morphs.py on the below filelist
## before calling this script, so that _param0 and _param1 files are present.

sys.path.extend(["..","../networks","../generators","../simulations"])

from networkConstants import * # has central_glom

seeds = [
("100.0","160.0"),
("100.0","190.0"),
("500.0","160.0"),
("500.0","set"),
("300.0","157.0")
]
filelist = [
"2011_05_15_01_07_odormorph_SINGLES_JOINTS_PGS_numgloms10.pickle",
"2011_05_15_22_30_odormorph_SINGLES_JOINTS_PGS_numgloms10.pickle",
"2011_05_17_22_51_odormorph_SINGLES_JOINTS_PGS_numgloms10.pickle",
"2011_05_18_14_51_odormorph_SINGLES_JOINTS_PGS_numgloms10.pickle",
"2011_05_20_18_46_odormorph_SINGLES_JOINTS_PGS_numgloms10.pickle"
]

coeffs = [0,0.2,0.4,0.6,0.8,1]
NUMWTS = len(coeffs[1:-1])

def constrain0to1(x):
    return math.exp(x)/(1+math.exp(x))

def get_weights(params, NUMBINS):
    #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
    #### constrain weights to be between 0 and 1
    #### sort the weights to ensure monotonicity
    inputsA = [ constrain0to1(x) for x in params[3*NUMBINS:(3*NUMBINS+NUMWTS)] ]
    inputsA.extend([0.0,1.0])
    inputsA.sort() # in place sort
    inputsB = [ constrain0to1(x) for x in params[(3*NUMBINS+NUMWTS):(3*NUMBINS+2*NUMWTS)] ]
    inputsB.extend([0.0,1.0])
    inputsB.sort(reverse=True) # weights of odor B need to be used in reverse
    return (inputsA,inputsB)

fig = figure(facecolor='none') # 'none' is transparent
## A super axes to set common x and y axes labels
bigAxes = axes(frameon=False) # hide frame
xticks([]) # don't want to see any ticks on this axis
yticks([])
text(-0.05,0.3,'fitted weights', fontsize=24, rotation='vertical')
text(0.4,-0.11,'real weights', fontsize=24, rotation='horizontal')

plotnum = 1
for filenum,filename in enumerate(filelist):
    for mitnum in [central_glom,central_glom+1]:
        filenamefull = '../results/odor_morphs/'+filename+'_params'+str(mitnum)
        f = open(filenamefull,'r')
        params = pickle.load(f)
        f.close()
        NUMBINS = len(params[:-9])/3
        inputsA,inputsB = get_weights(params,NUMBINS)
        maxerror = sqrt(sum(array([0.8,0.6,0.4,0.2])**2)/4.0) # max rms error
        ## normalized score = 1 - norm-ed rms error
        scoreA = 1 - sqrt( sum( (inputsA[1:-1]-arange(0.2,0.81,0.2))**2 )/4.0 )/maxerror
        scoreB = 1 - sqrt( sum( (inputsB[1:-1]-arange(0.8,0.19,-0.2))**2 )/4.0 )/maxerror    
        
        fig.add_subplot(3,4,plotnum) ### hardcoded - I know there are 5 files each with 2 mits
        plotnum += 1
        xticks([0,1],size='large')
        yticks([0,1],size='large')
        plot(coeffs,inputsA,color=(1,0,0),linewidth=2)
        plot(coeffs,inputsB,color=(0,1,0),linewidth=2)
        xlim(0,1)
        ylim(-0.01,1.01) #### very important to give it after the plot functions else autoscales
        ## or set autoscaling off - gca().set_autoscale_on(False)
        ## - taken from http://old.nabble.com/ylim-does-not-work-td19000814.html
        title( str(seeds[filenum])+' '+'mit'+str(mitnum)+'  %.2f'%scoreA + ' '+'%.2f'%scoreB, size='large')

show()
