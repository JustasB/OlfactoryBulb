# -*- coding: utf-8 -*-

## USAGE: python2.6 test_corrs.py

from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME
from data_utils import *
from simset_odor import* # has REALRUNTIME

f = open('../generators/firerates_2sgm_'+str(rate_seednum)+'.pickle','r')
frateOdorList,fratePulseList,randomPulseList,randomPulseStepsList,randomResponseList,kernels\
    = pickle.load(f)
f.close()

RUNTIME = REALRUNTIME + SETTLETIME
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
numt = len(firingtsteps)

def calc_corrs():

    ## frateOdorList[glomnum][odornum]
    ## ORN firing rates are very low - increase.
    
    ## compare increases in firing rate
    #frate1 = 0.5*50*frateOdorList[0][1] # glom0 0.8*odorB+0.2*odorA
    #frate2 = 50*frateOdorList[1][0] # glom1 odorB
    #frate3 = 0.5*3*frate1
    
    ## compare wide (v1) vs narrow (v2) vs phase shift (v3)
    frate1 = 50*frateOdorList[0][0] # glom0 odorB
    frate2 = 200*frateOdorList[1][5] # glom1 odorA
    frate3 = 50*frateOdorList[8][0] # glom8 odorB

    peakfrate = max(max(frate1),max(frate2),max(frate3))*2.0
    v1 = []
    v2 = []
    v3 = []
    for i in range(10):
        v1.append(poissonTrainVaryingRate(RUNTIME,peakfrate,REFRACTORY,firingtsteps,frate1))
        v2.append(poissonTrainVaryingRate(RUNTIME,peakfrate,REFRACTORY,firingtsteps,frate2))
        v3.append(poissonTrainVaryingRate(RUNTIME,peakfrate,REFRACTORY,firingtsteps,frate3))

    ## pairs of vectors varnames between which to calculate xcorrgrams.
    ## CHANGE THIS FOR DIFFERENT PAIRS:
    #corrgram_pairs = (('v1','v1'),('v1','v2'),('v1','v3'))
    corrgram_pairs = (('v3','v3'),('v3','v2'),('v3','v1'))
    #corrgram_pairs = (('v1','v3'),('v3','v1'))
    #corrgram_pairs = (('v1','v1'),('v3','v3'))

    starttime = REALRUNTIME+SETTLETIME-2*RESPIRATION
    endtime = REALRUNTIME+SETTLETIME
    T = endtime-starttime
    ## Dhawale et al 2010: 5 ms time bin, T=0.5s.
    ## I have T=1.0s as rat resp is 1.0s whereas mouse is 0.5s.
    ## refractory period in my poisson generator is 1ms, so have that as the bin size:
    ## must ensure that there are never more than one spike per bin per moving window.
    dt = 1e-3
    tcorrlist = arange(-T/4.0,T/4.0+1e-6,dt)

    colorlist = ['r','g','b','c','m','y','k']
    labellist = [va+' '+vb for (va,vb) in corrgram_pairs]
    
    ## calc xcorrgrams for each type of normalization.
    for (norm_str,title_str) in (('none','no norm'),('overall','integral norm'),\
            ('ref','ref norm'),('analogous','standard norm')):
        xcorrgrams = []
        for (va_str,vb_str) in corrgram_pairs:
            ## lookup the variable value given the variable name
            ## don't use eval() as it can evaluate any python expression
            ## hackable if user enters a string
            va = locals()[va_str]
            vb = locals()[vb_str]
            xcorrgrams.append(crosscorrgram( va, vb, dt, T/4.0, starttime, endtime, norm_str ))
        plot_corrs(tcorrlist,xcorrgrams,title_str,colorlist,labellist)

    return (v1,v2,v3)

def plot_corrs(tcorrlist,xcorrgrams,titlestr,colorlist,labellist):
    fig = figure(facecolor='none')
    ax = fig.add_subplot(111)
    for xcgnum,xcorrgram in enumerate(xcorrgrams):
        plot(tcorrlist, xcorrgram, color=colorlist[xcgnum], label=labellist[xcgnum])
    biglegend()
    axes_labels(ax,'time shift (s)','')
    title(titlestr+' sisters xcorrelogram', fontsize=24)

if __name__ == "__main__":
    rasters = calc_corrs()
    plot_rasters(rasters, RUNTIME)
    show()

