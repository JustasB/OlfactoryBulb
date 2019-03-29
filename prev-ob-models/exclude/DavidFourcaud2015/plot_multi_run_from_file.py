# -*- coding: utf-8 -*-
"""
This script provides analysis and plot functions for LFP recordings issued from multiple simulations while varying a single parameter.
Note that a large part of the oscillation analysis (time frequency map, oscillation detection) is based on OpenElectrophy 0.2

The "main" section provides the way to plot results from a simulation output file. 

Finally, it also allow to plot few spectograms of LFPs in a simulation output file.
"""

from brian import *
from oscillation_analysis import beta_gamma_detection, spectrum_analysis
import pylab as pl
import pickle
from scipy.stats import sem, nanmean, nanstd



def multi_LFP_analysis(results,osc_th=0.06,freq_cut=40.,burn_time=0.8,distinct_figures=False):
    """
    result : list of (param,LFPsig,LFPsr)
    From a list of LFPs with repetition of a set of param:
    * detect oscillations in each LFP
    * compute some oscillation properties 
    *compute and plot stats of properties across LFPs with the same param
    """

    # Analysis of data in results
    beta_stats=[]
    gamma_stats=[]
    params=[]
    all_betas=[]
    all_gammas=[]
    all_LFPlengths=[]
    for param,LFPsig,LFPsr in results: # Check each element recorded
        params.append(param) # useful to get all simulated parameters
        
        # For that element compute beta and gamma oscillations
        list_beta,list_gamma=beta_gamma_detection(LFPsig,LFPsr/Hz,freq_cut=freq_cut,osc_th=osc_th,burn_time=burn_time/second)
        all_betas.append(list_beta)
        all_gammas.append(list_gamma)
        all_LFPlengths.append(1.*LFPsig.size/LFPsr-burn_time)

    # Compute oscillation properties from all networks with the same parameter of interest, then compute mean and std of properties
    final_list_params=unique(params)
    params=array(params)
    beta_stats=[]
    gamma_stats=[]
    for par in final_list_params:
        
        valid_par=where(params==par)[0]
        
        # Prepare to record stats from each recording
        list_beta_stats=[]
        list_gamma_stats=[]
        
        # Useful function to gather parameters of interest (amp and freq) from all oscillations in a list
        def osc_stats(list_osc):
            amps,freqs,timeinoscs=[],[],[]
            for osc in list_osc:
                amps.append(osc.amplitude_max)
                freqs.append(osc.freq_max)
                dur = (osc.time_stop-osc.time_start)/all_LFPlengths[ind] # proportion of total possible time
                timeinoscs.append(dur)
            return mean(amps),mean(freqs),sum(timeinoscs),len(amps)
        
        for ind in valid_par:
            list_beta_stats.append(osc_stats(all_betas[ind]))
            list_gamma_stats.append(osc_stats(all_gammas[ind]))
            
        list_beta_stats=array(list_beta_stats)
        list_gamma_stats=array(list_gamma_stats)
        
        # Now compute stats on osc properties for all networks using the same param (use std or sem as wanted)
        beta_stats.append([[nanmean(list_beta_stats[:,ind]),nanstd(list_beta_stats[:,ind])] for ind in range(list_beta_stats.shape[1])])
        gamma_stats.append([[nanmean(list_gamma_stats[:,ind]),nanstd(list_gamma_stats[:,ind])] for ind in range(list_gamma_stats.shape[1])])
        # beta_stats.append([[nanmean(list_beta_stats[:,ind]),sem(list_beta_stats[~isnan(list_beta_stats[:,ind]),ind])] for ind in range(list_beta_stats.shape[1])])
        # gamma_stats.append([[nanmean(list_gamma_stats[:,ind]),sem(list_gamma_stats[~isnan(list_gamma_stats[:,ind]),ind])] for ind in range(list_gamma_stats.shape[1])])
        
        
    # Tranform into arrays for simpler manipulation (shape : param, osc prop, mean/std)
    beta_stats=array(beta_stats)
    gamma_stats=array(gamma_stats)

    lw=2 # linewidth for all curves

    if distinct_figures:
        figure()
    else:
        figure(10)
    ax=subplot(411) # axis for amp vs param
    errorbar(final_list_params,gamma_stats[:,0,0],yerr=gamma_stats[:,0,1],color='k',linewidth=lw)
    errorbar(final_list_params,beta_stats[:,0,0],yerr=beta_stats[:,0,1],color='r',linewidth=lw)
    ax.set_ylabel("Amplitude")
    pl.xlim(final_list_params[0],final_list_params[-1])
    pl.ylim(ymin=0.2)

    ax=subplot(412) # axis for freq vs param
    errorbar(final_list_params,gamma_stats[:,1,0],yerr=gamma_stats[:,1,1],color='k',linewidth=lw)
    errorbar(final_list_params,beta_stats[:,1,0],yerr=beta_stats[:,1,1],color='r',linewidth=lw)
    pl.ylim(15,90)
    pl.xlim(final_list_params[0],final_list_params[-1])
    ax.set_ylabel("Frequency")

    ax=subplot(413) # total time spent in gamma and in beta
    errorbar(final_list_params,gamma_stats[:,2,0],yerr=gamma_stats[:,2,1],color='k',linewidth=lw)
    errorbar(final_list_params,beta_stats[:,2,0],yerr=beta_stats[:,2,1],color='r',linewidth=lw)
    ax.set_ylabel("Length")
    pl.xlim(final_list_params[0],final_list_params[-1])
    pl.ylim(ymin=0)

    ax=subplot(414) # number of oscillations
    errorbar(final_list_params,gamma_stats[:,3,0],yerr=gamma_stats[:,3,1],color='k',linewidth=lw)
    errorbar(final_list_params,beta_stats[:,3,0],yerr=beta_stats[:,3,1],color='r',linewidth=lw)
    ax.set_ylabel("Number")
    pl.xlim(final_list_params[0],final_list_params[-1])
    pl.ylim(ymin=-1)
    return
    
def LFP_frequency_content(results,burn_time=0.8*second,freq_min=15.,freq_max=90.,verbose=False):
    """
    Compute wavelet map (or FFT ?) of each LFP, to get an average spectral content for the whole LFP
    Then compute mean and std (sem ?) of spectral profiles across LFP with the same param
    Make a plot with all average profile on the same plot
    """

    params=[]
    all_spectrums=[]
    
    for param,LFPsig,LFPsr in results: # Check each element recorded
        params.append(param) # useful to gather all simulated params
        
        Pxx,freqs=spectrum_analysis(LFPsig,LFPsr/Hz,burn_time=burn_time/second,plot_fig=False,verbose=verbose,return_full=True)
        all_spectrums.append(Pxx)
        
    all_spectrums=array(all_spectrums)
    params=array(params)
    fig=figure()
    ax=fig.add_subplot(111)
    mask_freq=(freqs>=freq_min)&(freqs<=freq_max)
    for ind,param in enumerate(unique(params)):
        mask=params==param
        ax.plot(freqs[mask_freq],all_spectrums[mask,:][:,mask_freq].mean(axis=0),label=str(param-4.),color=pl.cm.autumn(1.*(ind+1)/(len(unique(params)))),linewidth=2)
    ax.legend()
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("LFP average spectrum")

def plot_few_spectrograms(results,simul_to_plot=[0],osc_th=0.06,freq_cut=40.,burn_time=0.8*second):
    """
    Detect oscillations and plot spectrogram for all simul indices in simul_to_plot
    """
            
    for ind in simul_to_plot:
        param,LFPsig,LFPsr=results[ind] # Check each element recorded
        print "param: ",param
    
        # For that element compute beta and gamma oscillations
        list_beta,list_gamma=beta_gamma_detection(LFPsig,LFPsr/Hz,freq_cut=freq_cut,osc_th=osc_th,burn_time=burn_time/second,plot_fig=True)
        
    return

if __name__=="__main__":

    distinct_figures=False

    # Files to plot
    list_filenames=["put_filename_root_here"]
    list_filenames=["test_dict_connectproba"]
    
    # Detection osc params
    osc_th=[0.2, 0.2] # oscillation threshold
    freq_cut=40.
    burn_time=0.8*second
    
    # Spectrum params
    freq_min=15.
    freq_max=90.

    for filename in list_filenames:

        # Load data
        fid=open(filename+"_multi.out","rb")
        param_dict,results=pickle.load(fid)
        fid.close()

        print "Simulation params of plotted file: ",filename 
        for key,val in param_dict.items():
            print key,': ',val
        print
        
        tmp_results=[res[:3] for res in results]  # To keep only LFP and param
        multi_LFP_analysis(tmp_results,osc_th=osc_th,freq_cut=freq_cut,burn_time=burn_time,distinct_figures=distinct_figures)
        LFP_frequency_content(tmp_results,burn_time=burn_time,freq_min=freq_min,freq_max=freq_max)
        #~ plot_few_spectrograms(tmp_results,simul_to_plot=range(6),osc_th=osc_th,freq_cut=freq_cut,burn_time=burn_time)
    
show()
