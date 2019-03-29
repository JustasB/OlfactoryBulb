# -*- coding: utf-8 -*-

## USAGE: python2.6 calc_corrs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle

from scipy import stats
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has NUMBINS, overridden below / in function call
from networkConstants import * # has central_glom
from sim_utils import * # has rebin() to alter binsize
from data_utils import * # has crosscorrgram

## I override the NUMBINS (and bin_width_time) defined in simset_odor above, and I rebin() below
NUMBINS = 17
## smooths / overlapping bins. non-overlapping would be RESPIRATION/NUMBINS
bin_width_time = 2.0*RESPIRATION/NUMBINS

NUMMIX = len(inputList)

def get_tvecs(flist, morphnum, mitnum):
    """ returns a NUMAVGS array of firing times for this morphnum and mitnum. """
    ## cannot convert to numpy array(),
    ## as number of firing times is not fixed
    ## thus cannot use flist[:,morphnum,mitnum]
    return [flist[i][morphnum][mitnum] for i in range(len(flist))]

def calc_corrs(filename, norm_str, numbins=NUMBINS, bin_width_time=bin_width_time, printinfo = True):
    """
    Plots correlations between air and odor responses of sister cells.
    """
    f = open(filename,'r')
    #### mitral_responses_list[avgnum][odornum][mitralnum][tnum]
    #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()

    ###################### Input conditioning
    ## By default, rebin takes only the last respiratory cycle i.e. last numbins
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, numbins, bin_width_time)
    #### very important to convert to numpy array,
    #### else where() below returns empty list.
    mitral_responses_binned_list = array(mitral_responses_binned_list)
    mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## since I fit the mean response, I must use standard error/deviation of the _mean_
    ## = standard deviation of a repeat / sqrt(num of repeats).
    NUMAVGs = len(mitral_responses_binned_list)
    mitral_responses_se = mitral_responses_std/sqrt(NUMAVGs)

    ## mean frate change sisters in glom0 and non-sisters (glom1) for each odor
    odorA_Dfrate_mit0 = mean(mitral_responses_mean[5,central_glom] - mitral_responses_mean[6,central_glom])
    odorA_Dfrate_mit1 = mean(mitral_responses_mean[5,central_glom+1] - mitral_responses_mean[6,central_glom+1])
    odorB_Dfrate_mit0 = mean(mitral_responses_mean[0,central_glom] - mitral_responses_mean[6,central_glom])
    odorB_Dfrate_mit1 = mean(mitral_responses_mean[0,central_glom+1] - mitral_responses_mean[6,central_glom+1])
    #odorA_Dfrate_mit2 = mean(mitral_responses_mean[5,central_glom+2] - mitral_responses_mean[6,central_glom+2])
    #odorA_Dfrate_mit3 = mean(mitral_responses_mean[5,central_glom+3] - mitral_responses_mean[6,central_glom+3])
    #odorB_Dfrate_mit2 = mean(mitral_responses_mean[0,central_glom+2] - mitral_responses_mean[6,central_glom+2])
    #odorB_Dfrate_mit3 = mean(mitral_responses_mean[0,central_glom+3] - mitral_responses_mean[6,central_glom+3])
    Dfrates = (odorA_Dfrate_mit0,odorA_Dfrate_mit1,odorB_Dfrate_mit0,odorB_Dfrate_mit1)
    #    odorA_Dfrate_mit2,odorA_Dfrate_mit3,odorB_Dfrate_mit2,odorB_Dfrate_mit3)

    ## If one of the arrays is all zeroes, stats.pearsonr gives a one-time warning
    ## BINNED cross correlation
    air_corr = stats.pearsonr(mitral_responses_mean[6,central_glom],\
        mitral_responses_mean[6,central_glom+1])[0]
    odorA_corr = stats.pearsonr(mitral_responses_mean[5,central_glom],\
        mitral_responses_mean[5,central_glom+1])[0]
    odorB_corr = stats.pearsonr(mitral_responses_mean[0,central_glom],\
        mitral_responses_mean[0,central_glom+1])[0]
    if printinfo:
        print "air binned correlation between central sisters = ", air_corr
        print "odor A binned correlation between central sisters = ", odorA_corr
        print "odor B binned correlation between central sisters = ", odorB_corr
        #if isnan(air_corr): print mitral_responses_list, mitral_responses_binned_list

    ## Spike-based cross-correlation
    starttime = REALRUNTIME+SETTLETIME-2*RESPIRATION
    endtime = REALRUNTIME+SETTLETIME
    T = endtime-starttime
    ## Dhawale et al 2010: 5 ms time bin, T=0.5s.
    ## I have T=1.0s as rat resp is 1.0s whereas mouse is 0.5s.
    ## refractory period is typically 1ms, so have that as the bin size:
    ## must ensure that there are never more than one spike per bin per moving window.
    dt = 1e-3
    tcorrlist = arange(-T/4.0,T/4.0+1e-6,dt)
    v1 = get_tvecs(mitral_responses_list,6,central_glom)
    v2 = get_tvecs(mitral_responses_list,6,central_glom+1)
    airxcorrgram = crosscorrgram( v1, v2, dt, T/4.0, starttime, endtime, norm_str )
    v1 = get_tvecs(mitral_responses_list,5,central_glom)
    v2 = get_tvecs(mitral_responses_list,5,central_glom+1)
    odorAxcorrgram = crosscorrgram( v1, v2, dt, T/4.0, starttime, endtime, norm_str )
    v1 = get_tvecs(mitral_responses_list,0,central_glom)
    v2 = get_tvecs(mitral_responses_list,0,central_glom+1)
    odorBxcorrgram = crosscorrgram( v1, v2, dt, T/4.0, starttime, endtime, norm_str )
    #print "unbinned air correlation between central sisters max =",\
    #    max(airxcorrgram)
    #print "unbinned odor A correlation between central sisters max =",\
    #    max(odorAxcorrgram)
    #print "unbinned odor B correlation between central sisters max =",\
    #    max(odorBxcorrgram)

    return((air_corr,odorA_corr,odorB_corr),\
        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
        Dfrates)

def plot_corrs(tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram):
    fig = figure(facecolor='none')
    ax = fig.add_subplot(111)
    plot(tcorrlist, airxcorrgram, color=(0,0,0), label='air')
    plot(tcorrlist, odorAxcorrgram, color=(1,0,0), label='odor A')
    plot(tcorrlist, odorBxcorrgram, color=(0,1,0), label='odor B')
    biglegend(ax=ax)
    axes_labels(ax,'time (s)','spike probability')
    title('Crosscorrelogram between sisters', fontsize=24)
    show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print "Specify data file containing pickled list."
        sys.exit(1)

    (air_corr,odorA_corr,odorB_corr),\
        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
        Dfrates = \
        calc_corrs(filename, "overall")
    plot_corrs(tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram)
