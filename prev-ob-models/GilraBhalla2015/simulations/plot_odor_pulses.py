#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pickle
sys.path.extend(["..","../networks","../generators","../simulations"])

from OBNetwork import *
from stimuliConstants import * # has PULSE_SETTLETIME, pulsebins
from simset_odor import * # has ONLY_TWO_MITS

from pylab import * # part of matplotlib that depends on numpy but not scipy

def plot_pulse_responses(picklefile):
    f = open(picklefile,'r')
    (mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
    f.close()

    numpulses = RANDOM_PULSE_NUMS*2
    pulsebindt = PULSE_RUNTIME/pulsebins
    pulsetlist = arange(pulsebindt/2.0,PULSE_RUNTIME,pulsebindt)

    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    # since I plot the mean response, I must plot standard error of the mean
    # = standard deviation of a repeat / sqrt(num of repeats).
    mitral_responses_se = mitral_responses_std/sqrt(numavgs)

    if ONLY_TWO_MITS: mitlist = range(MIT_SISTERS)
    else: mitlist = range(NUM_GLOMS*MIT_SISTERS)
    for mitnum in mitlist:
        if mitnum%MIT_SISTERS == 0:
            figure()
            title('Glomerulus '+str(mitnum/MIT_SISTERS))
        for pulsenum in range(numpulses):
            sister_ratio = (mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            errorbar(pulsetlist,y=mitral_responses_avg[pulsenum,mitnum],\
                yerr=mitral_responses_se[pulsenum,mitnum],\
                color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio))

    ## compare responses of central sisters to odors A & B
    figure()
    mit0 = central_glom*2+0
    mit1 = central_glom*2+1
    title('Sisters of central glomerulus - odorA')
    ## first pulse
    errorbar(pulsetlist,y=mitral_responses_avg[1,mit0],\
        yerr=mitral_responses_se[1,mit0],\
        color=(1,0,0),label='odorA_mit0')
    errorbar(pulsetlist,y=mitral_responses_avg[1,mit1],\
        yerr=mitral_responses_se[1,mit1],\
        color=(1,0,1),label='odorA_mit1')
    ## second pulse
    figure()
    title('Sisters of central glomerulus - odorA')
    errorbar(pulsetlist,y=mitral_responses_avg[3,mit0],\
        yerr=mitral_responses_se[3,mit0],\
        color=(1,0.5,0),label='odorA_mit0')
    errorbar(pulsetlist,y=mitral_responses_avg[3,mit1],\
        yerr=mitral_responses_se[3,mit1],\
        color=(1,0.5,1),label='odorA_mit1')
    figure()
    title('Sisters of central glomerulus - odorB')
    ## first pulse
    errorbar(pulsetlist,y=mitral_responses_avg[2,mit0],\
        yerr=mitral_responses_se[2,mit0],\
        color=(0,1,0),label='odorB_mit0')
    errorbar(pulsetlist,y=mitral_responses_avg[2,mit1],\
        yerr=mitral_responses_se[2,mit1],\
        color=(0,1,1),label='odorB_mit1')
    ## second pulse
    figure()
    title('Sisters of central glomerulus - odorB')
    errorbar(pulsetlist,y=mitral_responses_avg[4,mit0],\
        yerr=mitral_responses_se[4,mit0],\
        color=(0.5,1,0),label='odorB_mit0')
    errorbar(pulsetlist,y=mitral_responses_avg[4,mit1],\
        yerr=mitral_responses_se[4,mit1],\
        color=(0.5,1,1),label='odorB_mit1')
    

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "You need to specify the morph responses pickle filename."
        sys.exit(1)
    plot_pulse_responses(sys.argv[1])
    show()
