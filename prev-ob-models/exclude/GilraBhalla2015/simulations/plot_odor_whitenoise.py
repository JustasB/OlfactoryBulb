#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pickle
sys.path.extend(["..","../networks","../generators","../simulations"])

from OBNetwork import *
from stimuliConstants import *
from simset_odor import * # has ONLY_TWO_MITS
from sim_utils import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

def plot_whitenoise_responses(picklefile):
    f = open(picklefile,'r')
    mitral_responses_list = pickle.load(f)
    f.close()

    bindt = 10e-3
    numbins = int(PULSE_RUNTIME/bindt)
    mitral_responses_binned_list = \
        rebin_pulses(mitral_responses_list, numbins, PULSE_RUNTIME, 0.0)
    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mean(mitral_responses_binned_list, axis=0),axis=0)

    pulsetlist = arange(0.0, PULSE_RUNTIME, bindt)
    if ONLY_TWO_MITS: mitlist = range(MIT_SISTERS)
    else: mitlist = range(NUM_GLOMS*MIT_SISTERS)
    for mitnum in mitlist:
        if mitnum%MIT_SISTERS == 0:
            figure()
            if ONLY_TWO_MITS: title('mit0 and mit2')
            else: title('Glomerulus '+str(mitnum/MIT_SISTERS))
        sister_ratio = (mitnum%MIT_SISTERS)/float(MIT_SISTERS)
        plot(pulsetlist,mitral_responses_avg[mitnum],color=(0,0,sister_ratio))

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "You need to specify the whitenoise responses pickle filename."
        sys.exit(1)
    plot_whitenoise_responses(sys.argv[1])
    show()
