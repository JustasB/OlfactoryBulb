# -*- coding: utf-8 -*-
"""
This script provides analysis and plot functions for detailed recordings issued from a single simulation.
Note that a large part of the oscillation analysis (time frequency map, oscillation detection) is based on OpenElectrophy 0.2

The "main" section provides the way to plot results from a simulation output file 
"""

from brian import *
from oscillation_analysis import beta_gamma_detection, spectrum_analysis
import pylab as pl
import pickle


def single_LFP_analysis(result_dict,osc_th=0.06,freq_cut=40.,burn_time=0.8*second,xinf=None,xsup=None,clim_max=None):

    # Oscillation analysis
    LFPsig,LFPsr,SpikesM_spikes=result_dict['LFP_values'][0,:],result_dict['LFP_sr'],result_dict['SpikesM_spikes']
    list_beta,list_gamma=beta_gamma_detection(LFPsig,LFPsr/Hz,freq_cut=freq_cut,verbose=False,plot_fig=True,osc_th=osc_th,burn_time=burn_time/second)

    fig=gcf()
    tfr_ax=fig.get_axes()
    
    if len(tfr_ax)>0: # test if time-frequency axis has been created (requires OpenElectrophy)
        tfr_ax[0].set_xlim(xinf,xsup)
        tfr_ax[1].set_ylim(-0.5,0.5)
        img=tfr_ax[0].get_images()[0]
        if clim_max is not None:
            img.set_clim(0,clim_max) # usually 0.5 for gamma only and 0.15 competition

    if 1:

        spectrum_analysis(LFPsig,LFPsr,burn_time=burn_time,plot_fig=True,verbose=True)
        beta_stats=[]
        gamma_stats=[]
        # Useful function to gather parameters of interest (amp and freq) from all oscillations in a list
        def osc_stats(list_osc):
            amps,freqs,timeinoscs=[],[],[]
            for osc in list_osc:
                amps.append(osc.amplitude_max)
                freqs.append(osc.freq_max)
                dur = osc.time_stop-osc.time_start
                timeinoscs.append(dur)
                print 'timeinosc=', sum(timeinoscs)
            return mean(amps),mean(freqs),std(amps),std(freqs),sum(timeinoscs)

        beta_stats.append(osc_stats(list_beta))
        gamma_stats.append(osc_stats(list_gamma))
        
        # Plot data

        SpikesM_spikes=result_dict['SpikesM_spikes']
        N_mitral=result_dict['N_mitral']
        LFP_times=result_dict['LFP_times']
        LFP_values=result_dict['LFP_values']
        if result_dict.has_key('avgIsyn_times'):
            avgIsyn_times=result_dict['avgIsyn_times']
            avgIsyn_values=result_dict['avgIsyn_values']
        M_gEinj=result_dict['M_gEinj']
        Mpot_times=result_dict['Mpot_times']
        Mpot_values=result_dict['Mpot_values']
        Mvars_times=result_dict['Mvars_times']
        Isyn_values=result_dict['Mvars_Isyn_values']
        Gsyn_values=result_dict['Mvars_Gsyn_values']
        use_granule=result_dict.has_key('SpikesG_spikes')
        if use_granule:
            SpikesG_spikes=result_dict['SpikesG_spikes']
            Gpot_times=result_dict['Gpot_times']
            Gpot_values=result_dict['Gpot_values']

        # Big figure with rasters, LFP, potentials and inhibitory input 
        fig=figure(figsize=(8,10))

        ax=fig.add_subplot(611)
        print "number of mitral spikes : ",len(SpikesM_spikes)
        if len(SpikesM_spikes)>0:
            pl.scatter(array(SpikesM_spikes)[:,1],array(SpikesM_spikes)[:,0],s=2)
            sps=array(SpikesM_spikes)
            rates=[]
            for i in range(N_mitral):
                mitral_mask=(sps[:,0]==i)&(sps[:,1]*second>burn_time)
                rates.append([i,mitral_mask.sum()/(LFP_times.max()-burn_time/second)])
            rates=array(rates)
        if use_granule:
            print "number of granule spikes : ",len(SpikesG_spikes)
            if len(SpikesG_spikes)>0:
                pl.scatter(array(SpikesG_spikes)[:,1],array(SpikesG_spikes)[:,0]+N_mitral,color='r',s=2)
        ax.set_ylabel("raster")
        ax.set_xlim((xinf,xsup))

        ax=fig.add_subplot(612,sharex=ax)
        plot(LFP_times,LFP_values[0,:],linewidth=2)
        ax.set_ylabel("LFP")
        ax.set_xlim((xinf,xsup))

        ax=fig.add_subplot(613,sharex=ax)
        for i in range(Mpot_values.shape[0]):
            ax.plot(Mpot_times,Mpot_values[i,:]+0.04*i,linewidth=2)
        ax.set_ylabel("Mitral potential")
        ax.set_xlim((xinf,xsup))

        ax=fig.add_subplot(614,sharex=ax)
        if result_dict.has_key('avgIsyn_times'):        
            ax.plot(avgIsyn_times,avgIsyn_values[0,:])
        else:
            for i in range(Isyn_values.shape[0]):
                ax.plot(Mvars_times,Isyn_values[i,:])
        ax.set_ylabel("Syn Inh I")
        ax.set_xlim((xinf,xsup))

        ax=fig.add_subplot(615,sharex=ax)
        for i in range(Gsyn_values.shape[0]):
            ax.plot(Mvars_times,Gsyn_values[i,:])
        ax.set_ylabel("Syn Inh conduct.")
        ax.set_xlim((xinf,xsup))


        if use_granule:
            ax=fig.add_subplot(616,sharex=ax)
            for i in range(Gpot_values.shape[0]):
                ax.plot(Gpot_times,Gpot_values[i,:])
            ax.set_ylabel("Granule pot")
            ax.set_xlim((xinf,xsup))
            
        checked_cells=[0,1,2,3,4]    
        # Mitral internal variables
        figure()
        title("Mitral internal variables")
        list_colors=['r','b','m']
        for i in checked_cells:
            for name,col in zip(['X','Y','W'],list_colors):
                plot(Mvars_times, result_dict['Mvars_'+name+'_values'][i,:], label=name,color=col)
                if name=='W':
                    nn=result_dict['Mvars_'+name+'_values'][i,:].size
                    print "Mean ",name," : ",result_dict['Mvars_'+name+'_values'][i,nn//2:].mean()
            legend()
            
        # Mitral excitatory inputs
        figure()
        subplot(211)
        bar(arange(M_gEinj.size),M_gEinj)
        ylim(6,8)
        title("Neuron excitatory inputs")
        subplot(212)
        bar(rates[:,0],rates[:,1])
        ylabel("Firing rates (Hz)")
        ylim(0,80)


if __name__=="__main__":

    # File to plot
    list_filenames=["put_filename_root_here"]
    list_filenames=["test_dict_simplegamma"]
    
    # Detection osc params
    osc_th=[0.1, 0.1]  # Threshold for beta and gamma
    freq_cut=40. 
    burn_time=0.8*second

    for filename in list_filenames:

        # Load data
        fid=open(filename+"_single.out","rb")
        param_dict,result_dict=pickle.load(fid)
        fid.close()

        print "Simulation params of plotted file: ",filename
        for key,val in param_dict.items():
            print key,': ',val
        print

        single_LFP_analysis(result_dict,osc_th=osc_th,freq_cut=freq_cut,burn_time=burn_time)

show()