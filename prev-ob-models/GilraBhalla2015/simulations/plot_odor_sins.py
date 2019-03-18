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

bindt = 2e-3

def plot_sin_responses(picklefile):
    f = open(picklefile,'r')
    mitral_responses_list = pickle.load(f)
    f.close()
    ## mitral_responses_list[avgnum][sinnum][mitnum][spikenum]

    numbins = int(SIN_RUNTIME/bindt)
    mitral_responses_binned_list = \
        rebin_pulses(mitral_responses_list, numbins, SIN_RUNTIME, 0.0)
    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)

    sintlist = arange(0.0, SIN_RUNTIME, bindt)
    mitnum = 0
    for sinnum in range(num_sins):
        figure(facecolor='w')
        sincolor = (sinnum+1) / float(num_sins)
        plot(sintlist,mitral_responses_avg[sinnum][mitnum],
            color=(0,1-sincolor,sincolor))

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "You need to specify the whitenoise responses pickle filename."
        sys.exit(1)
    plot_sin_responses(sys.argv[1])
    show()
