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

NUMBINS = 17
BIN_WIDTH_TIME = RESPIRATION/8.0
fitted_mitral = 2*central_glom+0

def plot_varinh_responses(mitral_responses_avg,mitral_responses_se):

    if ONLY_TWO_MITS: mitlist = range(MIT_SISTERS)
    else: mitlist = range(NUM_GLOMS*MIT_SISTERS)
    for mitnum in mitlist:
        figure()
        title('Mitral '+str(mitnum))
        for inhnum in range(NUMINHS):
            sister_ratio = (mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            errorbar(x=range(NUMBINS),y=mitral_responses_avg[inhnum,mitnum],\
                yerr=mitral_responses_se[inhnum,mitnum],\
                color=(inhnum/float(NUMINHS),1-inhnum/float(NUMINHS),0))

def read_odorresults_file(filename):
    f = open(filename,'r')
    (mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
    f.close()
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, numbins=NUMBINS, bin_width_time=BIN_WIDTH_TIME)

    numavgs = len(mitral_responses_list)
    mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## since I plot the mean response, I must plot standard error of the mean
    ## = standard deviation of a repeat / sqrt(num of repeats).
    mitral_responses_se = mitral_responses_std/sqrt(numavgs)
    
    return mitral_responses_avg, mitral_responses_se

def plot_varinh(filename):
    
    mitral_responses_avg, mitral_responses_se =\
        read_odorresults_file(filename)

    #### Usual firing rate vs time plots of responses
    plot_varinh_responses(mitral_responses_avg,mitral_responses_se)

    show()

def plot_varinh_diff(fn1,fn2):
    mitral_responses_avg1, mitral_responses_se1 =\
        read_odorresults_file(fn1)
    mitral_responses_avg2, mitral_responses_se2 =\
        read_odorresults_file(fn2)

    if ONLY_TWO_MITS:
        mitral_responses_avg = mitral_responses_avg1[:,:2] - mitral_responses_avg2[:,:2]
        mitral_responses_se = sqrt( mitral_responses_se1[:,:2]**2 + mitral_responses_se2[:,:2]**2 )
    else:
        mitral_responses_avg = mitral_responses_avg1 - mitral_responses_avg2
        mitral_responses_se = sqrt( mitral_responses_se2**2 + mitral_responses_se2**2 )

    #### plot (difference of firing rates) vs time
    plot_varinh_responses(mitral_responses_avg,mitral_responses_se)

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "You need to specify the morph responses pickle filename."
        sys.exit(1)
    elif len(sys.argv)==2:
        plot_varinh(sys.argv[1])
    else:
        plot_varinh_diff(sys.argv[1],sys.argv[2])
    show()
