# -*- coding: utf-8 -*-

## USAGE: python2.6 calc_entropy_morphs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle

from scipy import optimize
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from networkConstants import * # has central_glom
from data_utils import * # has info th functions
from analysis_utils import * # has read_morphfile() and NUM_REBINS, etc.
from infoth_test import * # has plot_table()

info_dt = 10e-3
## last time bin if smaller than info_dt will get thrown away.
num_infobins = int((ODORRUNTIME-SETTLETIME)/info_dt)
timerange = num_infobins*info_dt

def plot_table(rasters,rowlabels,collabels,data,cellcolours,titlestr,figfilename):
    ## 'plot' a table
    fig = figure(figsize=(8, 6), dpi=100)
    ax = fig.add_axes([0.14, 0.85, 0.95, 0.1])
    axes_off(ax)
    ## loop over rasters in reverse order, as they are plotted from bottom upwards
    numrasters = len(rasters)
    for rasteri,raster in enumerate(rasters[::-1]):
        raster = array(raster)
        ## find out indices of 1-s and plot them:
        rasterindices = where(raster==1)[0]
        ax.plot(rasterindices,zeros(len(rasterindices))+rasteri,'|k',\
            markersize=40/2**(numrasters-1), markeredgewidth='2') # | is the marker
    ax.set_ylim(-0.5,rasteri+0.5)
    dirtItable = ax.table(cellText=data, cellColours=cellcolours, rowLoc='right',\
        rowLabels=rowlabels, colLabels=collabels, colLoc='center', loc='bottom')
    table_props = dirtItable.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells:
        cell.set_height(1.5)
        cell.set_fontsize(18)
    ax.set_title(titlestr,fontsize=14)
    ## tight_layout() doesn't seem to work with table
    #fig.tight_layout()
    #fig.savefig(figfilename,dpi=300)

def calc_morph_entropyrates(filename):
    f = open(filename,'r')
    #### each mitral_responses_list[avgnum][odornum][mitralnum][spikenum] stores spiketime
    #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()
    spiketrains_mits = []
    nummitrals = len(mitral_responses_list[0][0])
    for mitnum in range(nummitrals):
        spiketrains = []
        for belowtrials in mitral_responses_list:
            for belowodornums in belowtrials:
                spiketrain = \
                    get_spiketrain_from_spiketimes(belowodornums[mitnum],SETTLETIME,timerange,num_infobins)
                spiketrains.append(spiketrain)
        spiketrains_mits.append(spiketrains)
        print "Entropy rate of mitral num",mitnum,'=',calc_entropyrate(spiketrains,markovorder=5)

    collabels = ['Order 1','2','4']
    rowlabels = ['Delay 0','1','2','3']
    dirtIs = []
    cellcolours = []
    for delay in [0,1,2,3]:
        print "delay =",delay
        dirtIorders = []
        cellcoloursorders = []
        for order in [1,2,4]:
            print "order =",order
            I2to0 = calc_dirtinforate(spiketrains_mits[2],spiketrains_mits[0],order,order,delay,delay)
            I0to2 = calc_dirtinforate(spiketrains_mits[0],spiketrains_mits[2],order,order,delay,delay)
            I4to1 = calc_dirtinforate(spiketrains_mits[4],spiketrains_mits[1],order,order,delay,delay)
            I1to4 = calc_dirtinforate(spiketrains_mits[1],spiketrains_mits[4],order,order,delay,delay)
            I1to2 = calc_dirtinforate(spiketrains_mits[0],spiketrains_mits[1],order,order,delay,delay)
            dirtIstr = '2to0 = {:1.3f}, 0to2 = {:1.3f}\n4to1 = {:1.3f}, 1to4 = {:1.3f}\n0to1 = {:1.3f}'\
                .format(I2to0,I0to2,I4to1,I1to4,I1to2)
            print dirtIstr
            dirtIorders.append(dirtIstr)
            if I2to0>0.9: cellcoloursorders.append('r')
            else: cellcoloursorders.append('w')
        dirtIs.append(dirtIorders)
        cellcolours.append(cellcoloursorders)
    titlestr = ""
    plot_table([spiketrains_mits[0][19],spiketrains_mits[1][19],\
        spiketrains_mits[2][19],spiketrains_mits[4][19]],\
        rowlabels,collabels,dirtIs,cellcolours,titlestr,'copycat_mydefn.svg')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print "Specify data file containing pickled list."
        sys.exit(1)

    calc_morph_entropyrates(filename)
    show()
