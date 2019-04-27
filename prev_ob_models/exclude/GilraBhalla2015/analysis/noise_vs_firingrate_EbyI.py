# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO BE A CLONE OF MUKUND'S AND ADIL'S MATLAB ONE
## USAGE: python2.6 fit_odor_morphs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle

metainfo = [
"E=400*1nS*2, no inhibition",
"E=400*1nS*2, I=10000*0.02nS, E/I~=4.0",
"E=400*1nS*2, I=10000*0.10nS, E/I~=0.8",
"E=400*1nS*2, I=10000*0.20nS, E/I~=0.4",
"E=400*1nS*2, I=10000*0.50nS, E/I~=0.16",
"E=400*1nS*2, I=10000*1.00nS, E/I~=0.08",
]
filelist = [
"2011_05_25_17_56_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_05_25_22_29_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_05_25_18_22_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_05_25_18_27_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_05_25_18_36_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_05_25_18_46_odormorph_NOSINGLES_NOJOINTS_NOPGS_numgloms1.pickle"
]
colorlist = ['r','g','b','c','m','y','k']
markerlist = ['s','o','d','+','x','^','<','>','v']

from pylab import *
import pickle
import sys

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has NUMBINS
from networkConstants import * # has central_glom
from sim_utils import * # has rebin() to alter binsize
from data_utils import *

## I override the NUMBINS in simset_odor above, and I rebin() below
NUMBINS = 17
## smooths / overlapping bins. non-overlapping would be RESPIRATION/NUMBINS
bin_width_time = 2.0*RESPIRATION/NUMBINS 

NUMMIX = len(inputList)

if __name__ == "__main__":
    
    for filenum,filename in enumerate(filelist):
        fig = figure(facecolor='w') # 'none' is transparent
        filenamefull = "../results/odor_morphs/"+filename
        f = open(filenamefull,'r')
        #### mitral_responses_list[avgnum][odornum][mitralnum][tnum]
        #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
        mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
        f.close()

        ###################### Input conditioning
        mitral_responses_binned_list = \
            rebin(mitral_responses_list, NUMBINS, bin_width_time)
        #### very important to convert to numpy array,
        #### else where() below returns empty list.
        mitral_responses_binned_list = array(mitral_responses_binned_list)
        mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
        mitral_responses_std = std(mitral_responses_binned_list, axis=0)
        ## since I fit the mean response,
        ## I must use standard error/deviation of the _mean_
        ## = standard deviation of a repeat / sqrt(num of repeats).
        NUMAVGs = len(mitral_responses_binned_list)
        mitral_responses_se = mitral_responses_std/sqrt(NUMAVGs)
        
        ######## scatter plot of std error vs firing rate
        scatter(mitral_responses_mean,mitral_responses_se,\
            s=20, color=colorlist[filenum], marker=markerlist[filenum], label=metainfo[filenum])
        
        ## Poisson noise
        plot(arange(0,50,1),[sqrt(val) for val in arange(0,50,1)],color='k',label='Poisson')
        legend()

    show()
        
        
