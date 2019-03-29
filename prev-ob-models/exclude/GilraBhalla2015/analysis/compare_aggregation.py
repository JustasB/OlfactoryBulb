# -*- coding: utf-8 -*-

########## Compare the variation in each bin for two different aggregation ratios
## First edit the reference and compared filenums in the script below
## USAGE: python2.6 compare_aggregation

metainfo = [
"1:1 singles (1:1 syns), no joints",
"2:1 singles (1:1 syns), no joints",
"2:1 singles (2:1 syns), no joints",
"100:1 singles (10:1 syns), no joints",
"100:1 singles (10:1 syns), 1:1 joints",
"100:1 singles (10:1 syns), 2:1 joints",
"100:1 singles (10:1 syns), 2:1 joints, firefileseed=300.0"
]
filelist = [
"2011_06_03_05_08_odormorph_SINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_11_59_odormorph_SINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_15_34_odormorph_SINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_12_24_odormorph_SINGLES_NOJOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_14_56_odormorph_SINGLES_JOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_15_06_odormorph_SINGLES_JOINTS_NOPGS_numgloms1.pickle",
"2011_06_03_15_21_odormorph_SINGLES_JOINTS_NOPGS_numgloms1.pickle"
]
colorlist = ['r','g','b','c','m','y','k']
markerlist = ['s','o','d','+','x','^','<','>','v']

reference_filenum = 4
compared_filenum = 6

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

def get_mitral_responses(filename):
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
    return mitral_responses_mean, mitral_responses_se

if __name__ == "__main__":
    ref_response_mean, ref_response_se = get_mitral_responses(filelist[reference_filenum])
    comp_response_mean, comp_response_se = get_mitral_responses(filelist[compared_filenum])
    mitnum = 2*central_glom+0
    ## add 1e-10 to ensure no divide by zero issues.
    ref_response_mean = array(ref_response_mean[:,mitnum])+1e-10
    comp_response_mean = array(comp_response_mean[:,mitnum])+1e-10
    fig1 = figure(facecolor='w') # 'none' is transparent
    fig2 = figure(facecolor='w') # 'none' is transparent
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ## 0 is odor B, 5 is odorA, 6 is air
    numodors = len(inputList)
    for odornum in range(numodors):
        odorA,odorB = inputList[odornum]
        ## normalize bin-wise
        error_mean = (comp_response_mean[odornum] - ref_response_mean[odornum]) \
            / maximum(ref_response_mean[odornum],comp_response_mean[odornum]) # element-wise maximum (numpy)
        ax1.plot(range(NUMBINS),error_mean,color=(odorA,odorB,0),marker=markerlist[odornum])
        ## normalize with maximum of all responses
        error_mean = (comp_response_mean[odornum] - ref_response_mean[odornum]) \
            / max(max(flatten(ref_response_mean)),max(flatten(comp_response_mean)))
        ax2.plot(range(NUMBINS),error_mean,color=(odorA,odorB,0),marker=markerlist[odornum])
    metastr = '\n' + metainfo[compared_filenum] + ' ' + metainfo[reference_filenum]
    ax1.set_title('response error normed bin-wise (Hz) vs bin'+metastr)
    ax2.set_title('response error normed with max (Hz) vs bin'+metastr)
    show()
