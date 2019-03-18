# -*- coding: utf-8 -*-

## USAGE: python2.6 calc_kernel_xcorrs.py

from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, PULSE_RUNTIME, pulsebins
from data_utils import * # has axes_labels()

pulsebindt = PULSE_RUNTIME/pulsebins
kernel_size = int(2.0/pulsebindt) # 2.0 seconds kernel size

#listtype = 'firingrate'
listtype = 'non-linear ORNs'

if listtype == 'firingrate':
    ## firingrate, filename
    filelist = [
    (2.0,'../results/odor_pulses/2011_08_19_14_46_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle'),
    (3.5,'../results/odor_pulses/2011_08_19_17_15_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle'),
    (10,'../results/odor_pulses/2011_08_19_17_19_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle')
    ]
elif listtype == 'inhibition':
    ## inhibition, filename
    filelist = [
    (8.0,'../results/odor_pulses/2011_08_20_19_38_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle'),
    (10.0,'../results/odor_pulses/2011_08_19_17_19_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle')
    ]
elif listtype == 'non-linear ORNs':
    ## ORN frate followed by static non-linearity
    filelist = [
    (2.0,'../results/odor_pulses/2011_08_19_14_46_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle'),
    (3.0,'../results/odor_pulses/2011_08_19_21_00_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle'),
    (8.0,'../results/odor_pulses/2011_08_19_20_51_odorpulses_SINGLES_JOINTS_PGS_numgloms3.pickle')
    ]


def calc_corrs_chisq(filenamelist):
    corrsR = []
    corrsA = []
    corrsB = []
    chisqs = []
    for i,(frate,filename) in enumerate(filenamelist):
        f = open(filename+'_params','r')
        sigmoid,chisq,params = pickle.load(f)
        f.close()
        kernelR = params[0:kernel_size]
        kernelA = params[kernel_size:2*kernel_size]
        kernelB = params[2*kernel_size:3*kernel_size]
        if i==0:
            kernelR1 = kernelR
            kernelA1 = kernelA
            kernelB1 = kernelB
        else:
            ## find cross-correlation of current kernel with first kernel
            corrsR.append( corrcoef([kernelR1,kernelR],rowvar=1)[0,1] )
            corrsA.append( corrcoef([kernelA1,kernelA],rowvar=1)[0,1] )
            corrsB.append( corrcoef([kernelB1,kernelB],rowvar=1)[0,1] )
        chisqs.append(chisq)
    return corrsR,corrsA,corrsB,chisqs

if __name__ == "__main__":
    
    corrsR,corrsA,corrsB,chisqs = calc_corrs_chisq(filelist)
    frates = [frate for (frate,filename) in filelist]
    #between_frates = [(frate+frates[i+1])/2.0 for i,frate in enumerate(frates[:-1])]

    fig = figure(facecolor='none')
    ax = fig.add_subplot(111)
    plot(frates, chisqs, color=(0,0,1), marker='o', label='chisqs')
    axes_labels(ax,'ORN mean firing rate (Hz)','Chisq')
    ax2 = ax.twinx()
    plot(frates[1:], corrsR, color=(0,0,0), marker='o', label='corrR')
    plot(frates[1:], corrsA, color=(1,0,0), marker='o', label='corrA')
    plot(frates[1:], corrsB, color=(0,1,0), marker='o', label='corrB')
    biglegend('lower left')
    axes_labels(ax2,'ORN mean firing rate (Hz)','Corr')
    title('Chisq/Corrs vs ORN firing rate', fontsize=24)
    xlim(0.0,12.0)

    show()
