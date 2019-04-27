# -*- coding: utf-8 -*-
"""
This script contains the main fucntion to simulate the network and the associated dictionary of default parameters.

The "main" section provide a way to simulate and record 
either a single run with lot of details 
or multiple runs (with multiprocessing) with only LFPs
"""

#~ import brian_no_units
from brian import *

from numpy.random import seed

from model_mitral_clean import *
from model_granule_clean import *
from populationstatemonitor import *

from oscillation_analysis import beta_gamma_detection, spectrum_analysis
from plot_single_run_from_file import single_LFP_analysis
from plot_multi_run_from_file import multi_LFP_analysis,LFP_frequency_content


# Dictionary of default parameters for simulations, used as a basis to generate dictionary of specific simulations
param_dict_default=dict(
                        # General simulation params
                        dt=0.05*ms,
                        sim_time=4.*second,
                        t_init_M=0.*ms, # time start of mitral sensory excitation (before, no input)
                        t_init_G=0.*ms, # time start of granule centrifugal modulation (before, base input is used)
                        
                        # Mitral intrinsic and input params
                        N_mitral = 100,
                        num_activated_mitral=100, # introduction of this parameter simplify the random selection of activated mitral cells at each run
                        M_gEinj= (linspace,(6.1,7.6,100),siemens*meter**-2), # (linspace(6.1,7.6,N_mitral))*siemens*meter**-2, # MC max sensory input conductances
                        M_gEinj_shift=0., # Used only to test gamma for distinct average M_gEinj
                        M_gEinj_base=None, # use None for no slow respiratory modulation of MC input, otherwise generally 4.0*siemens*meter**-2 
                        M_taumKs = 7.*ms, 
                        M_tauI=7.*ms, # weak inh decay time
                        M_tauI_G=7.*ms, # strong inh decay time
                        M_gI_cst= 20.*siemens*meter**-2, # mitral constant inhibitory input
                        recorded_mitrals=[10,20,30,40,50,60,70,80,90,99], # if return_details=True only !
                        
                        # Granule input params
                        use_granule=False,
                        N_granule=100,
                        G_input='constant', # choose between "constant", "ramp", "sinusoid"
                        G_I_base=-4.*nA, # before t_init or for "constant" input
                        G_I_min=-4.*nA, # min of sinusoid
                        G_I_max=-.1*nA, # max of sinusoid
                        recorded_granules=[10,20,30,40,50,60,70,80,90,99], # if return_details=True only !
                        
                        # Respiratory rhythm parameters
                        freq_modul=2.*Hz, 
                        MC_phase_dispersion=1.5, # SD of Gaussian distribution of phase shifts
                        GC_phase_dispersion=0.2, # SD of Gaussian distribution of phase shifts
                        GC_phase_shift=-pi/2, # Phase shift of granule centrifugal input vs mitral sensory input
                        
                        # Connectivity parameters
                        weakinh_connectproba=1.,
                        weakinh_weight=1.,
                        weakinh_mean_delay=9.*ms, # final delays: uniform in the range "mean +/- dispersion"
                        weakinh_delay_dispersion=4.*ms, # final delays: uniform in the range "mean +/- dispersion"
                        stronginh_connectproba=0.5,
                        stronginh_weight=1.,
                        exc_weight=1., # Automatically multiplied by 2.5 if use_AMPA_STP is set to True
                        use_AMPA_STP=False, # MC to GC short term plasticity (calibrated for a pure "depression")
                        )

def reseau_mitral_granule(param_dict,return_details=True,report='text'):
    """
    Main function to launch a single run of the network with a set of parameters given by param_dict
    Check param_dict_default for the list of parameters and their roles
    
    If return_details is set to True, a result_dict is produced with many details recorded during the simulations,
    otherwise only the LFP is returned
    
    "report" corresponds to the BRIAN "brian.run" function parameter
    """

    # Clear Brian objects in memory
    reinit_default_clock()
    clear(True)
    seed()

    # Parameters
    defaultclock.dt = param_dict['dt']
    sim_time = param_dict['sim_time']
    t_init_M = param_dict['t_init_M']
    t_init_G = param_dict['t_init_G']
    use_granule = param_dict['use_granule']
    
    # Mitral model and initialization
    N_mitral= param_dict['N_mitral']
    M = NeuronGroup(N_mitral, Mitral_eqs, threshold=-30*mV, reset=Mitral_reset,freeze=False,compile=True)
    M.V=(-60.+rand(N_mitral)*15.)*mV
    M.W=0.1*(1.+0.5*rand(N_mitral))
    M.X=0.3*(1.+0.5*rand(N_mitral))
    M.Y=0.
    func_name,func_param,func_unit=param_dict['M_gEinj']
    M.gEinj= func_name(*func_param)*func_unit+param_dict['M_gEinj_shift']*func_unit
    
    if param_dict['num_activated_mitral']<N_mitral:
        # Randomly put to 0 the activation of N_mitral-param_dict['num_activated_mitral']
        M.gEinj[np.random.permutation(N_mitral)[:100-int(param_dict['num_activated_mitral'])]]=0
    
    if param_dict['M_gEinj_base'] is None:
        gEinj_base=M.gEinj # no oscillation of MC sensory inputs
    else:
        gEinj_base=param_dict['M_gEinj_base']
        
    M.taumKS=param_dict['M_taumKs']
    M.tauI=param_dict['M_tauI']
    M.tauI_G=param_dict['M_tauI_G']
    M.gI_cst=param_dict['M_gI_cst']
    freq_modul = param_dict['freq_modul']
    phase_shifts = param_dict['MC_phase_dispersion']*randn(N_mitral) #3.0 # 0.3 // 2 // 5  (fig5C) #(fig3 and fig5 use 1.5) # fig4 uses 2.5 
    
    # M-M Weak inhibition
    dm,dd=param_dict['weakinh_mean_delay'],param_dict['weakinh_delay_dispersion']
    Cmm_weak_inh = Connection(M,M,'rI',
                                                  weight=param_dict['weakinh_weight'],
                                                  delay=(dm-dd,dm+dd), # mean delay +/- delay dispersion
                                                  sparseness=param_dict['weakinh_connectproba'])
    
    # Only for LFP recording
    Cmm_weak_inh2 = Connection(M,M,'rI2',weight=1.0)
    
    gEinj_max=M.gEinj.copy()
    @network_operation()
    def mitral_activation(clock):
        if clock.t>t_init_M:
            M.gEinj=(gEinj_max+gEinj_base)/2.+(gEinj_max-gEinj_base)/2.*cos(2*pi*freq_modul*Hz*clock.t+phase_shifts)
        else :
            M.gEinj=0.*siemens*meter**-2

    #  Granule params
    if use_granule:
        # Granule model and initialization
        N_granule=param_dict['N_granule']
        G_input=param_dict['G_input']
        phase_shifts_g=param_dict['GC_phase_dispersion']*randn(N_granule) # shift of phase to avoid artifact synchrony 0.2 base / 5 for vigile conditions
        G_I_base=param_dict['G_I_base']
        G_I_min=param_dict['G_I_min']
        G_I_max=param_dict['G_I_max']

        # Granule model and initialization
        G=NeuronGroup(N_granule,QIF_eqs,threshold=0.*mV,reset=-70.*mV,freeze=True,compile=True)
        G.V=-70.*mV
        G.Iinj=G_I_base

        # Dendro-dendritic connections
        if not(param_dict['use_AMPA_STP']):
            Cmg_AMPA = Connection(M,G,'sE',
                                                    weight=param_dict['exc_weight'],
                                                    sparseness=param_dict['stronginh_connectproba'],
                                                    delay=1.*ms)
        else: #withSTP
            Cmg_AMPA = Connection(M,G,'sE',
                                                    weight=2.5*param_dict['exc_weight'], # 2.5 factor to compensate for the depressed strength of AMPA in the beta regime
                                                    sparseness=param_dict['stronginh_connectproba'],
                                                    delay=1.*ms) #
            mystp=STP(Cmg_AMPA,taud=150*ms,tauf=1*ms,U=1.)
        
        Cgm_GABA = Connection(G,M,'sI_G',delay=1.*ms) # synaptic weight is given later
        connected=Cmg_AMPA.W.transpose().copy() # copy and transpose M->G connection array
        mask_connect=(connected>0)
        Cgm_GABA.connect(G,M,mask_connect*param_dict['stronginh_weight']) # finally create symmetrical connections
        
        # Time dependent activation of granule cells
        @network_operation()
        def granule_activation(clock):
            if clock.t<t_init_G:
                G.Iinj=G_I_base
            else:               
                if G_input=='sinusoid':
                    G.Iinj=(G_I_max+G_I_min)/2.+(G_I_max-G_I_min)/2.*cos(2*pi*freq_modul*Hz*clock.t+phase_shifts_g+param_dict['GC_phase_shift'])
                elif G_input=='ramp': # standard ramp : from -3.5nA to 0.5nA
                        G.Iinj=G_I_min* (1.-(clock.t-t_init_G)/(sim_time-t_init_G)) + G_I_max*(clock.t-t_init_G)/(sim_time-t_init_G) # current ramp, linear from first to second current value 
                else:
                    G.Iinj=G_I_base
                

    # Monitors
    LFP=PopulationStateMonitor(M,'LFP') # actually it is "fake" Isyn which omits random delays (only spikes convolved with an IPSC waveform, averaged across all cells)
    
    if return_details: # Additional detailed recordings
        recorded_mitrals=param_dict['recorded_mitrals']
        SpikesM = SpikeMonitor(M)
        Mpot=StateMonitor(M,'V',record=recorded_mitrals)
        avgIsyn=PopulationStateMonitor(M,'Isyn')
        Mvars=MultiStateMonitor(M,['Gsyn','Isyn','W','X','Y'],record=recorded_mitrals)
        if use_granule:
            recorded_granules=param_dict['recorded_granules']
            SpikesG = SpikeMonitor(G)
            Gpot=StateMonitor(G,'V',record=recorded_granules)
            Gvars=MultiStateMonitor(G,['sE'],record=recorded_granules)
            
    # Simulation
    run(sim_time,report=report)

    if return_details:
                
        result_dict={}
        result_dict['SpikesM_spikes']=SpikesM.spikes
        if use_granule:
            result_dict['SpikesG_spikes']=SpikesG.spikes
            result_dict['N_granule']=N_granule
            result_dict['Gpot_times']=Gpot.times
            result_dict['Gpot_values']=Gpot.values
        result_dict['N_mitral']=N_mitral
        result_dict['LFP_times']=LFP.times
        result_dict['LFP_values']=LFP.values
        result_dict['Mpot_times']=Mpot.times
        result_dict['Mpot_values']=Mpot.values
        result_dict['avgIsyn_times']=avgIsyn.times
        result_dict['avgIsyn_values']=avgIsyn.values
        result_dict['Mvars_times']=Mvars.times
        for name,mon in Mvars.items():
            result_dict['Mvars_'+name+'_values']=Mvars[name].values        
        result_dict['LFP_sr']=1./LFP.clock.dt
        result_dict['M_gEinj']=M.gEinj
        
        return result_dict
    else:
        print "One simulation done..."
        return LFP[0],1./LFP.clock.dt

if __name__=="__main__":
    
    
    from params_for_article_fig import *
    
    # Parameters to save simulation output
    save_output=False
    filename="test_dict_simplegamma" # _multi.out ou _single.out sont automatiquement ajouté après le nom suivant le type de simul     
    
    # Detection osc params
    osc_th=[0.1,0.1] # 0.2 for constant gamma or beta, 0.1 in presence of slow external modulations
    freq_cut=40.
    burn_time=1.*second
    
    if 1:  # To launch one network (with figures)

        # ## Choose params
        param_dict=param_dict_default.copy()
        # ## Or choose one of the predefind parameter set (defined for article figures)
        # param_dict=gamma_single_step_dict
        # param_dict=beta_single_ramp_dict
        # param_dict=competition_single_base_dict
        # param_dict=competitionSTP_single_lowintensity_dict
        # param_dict=competitionSTP_single_mediumintensity_dict
        # param_dict=competitionSTP_single_highintensity_dict
        # param_dict=competitionSTP_single_lowcentrifugal_dict
        # param_dict=competitionSTP_single_highcentrifugal_dict
        # param_dict=competitionSTP_single_awake
        # param_dict=competitionSTP_single_awake_lowcentrifugal
        
        # Run single network 
        result_dict = reseau_mitral_granule(param_dict)
        if save_output:
            print "Saving"
            import pickle
            fid=open(filename+"_single.out","wb")
            pickle.dump([param_dict,result_dict],fid)
            fid.close()
            print "Saved"
            
        # Plot detailed results
        single_LFP_analysis(result_dict,osc_th=osc_th,freq_cut=freq_cut,burn_time=burn_time)
    
    else: # To launch multiple network on different CPUs (only a synthesis figure is done)
        
        # Replace the parameter to vary with a numpy array of values to simulate
        param_dict=param_dict_default.copy()
        param_dict['weakinh_connectproba']=linspace(0.05,1.0,10)
        # ## Or choose one of the predefind parameter set (defined for article figures)
        # param_dict=gamma_multi_weakinhconnect_dict
        # param_dict=gamma_multi_numact_dict
        # param_dict=gamma_multi_weakinhweight_dict
        # param_dict=gamma_multi_tauGABA_dict
        # param_dict=gamma_multi_taumKs_dict
        # param_dict=gamma_multi_MgEinj_dict
        # param_dict=beta_multi_tauGABA_dict
        # param_dict=beta_multi_stronginhweight_dict
        # param_dict=beta_multi_GIinj_dict
        # param_dict=beta_multi_excweight_dict
        # param_dict=competition_multi_intensity
        # param_dict=competitionSTP_multi_intensity
        # param_dict=competitionSTP_multi_centrifugal_dict

        number_of_runs = 10 # Number of runs for each paramater (with different inputs and random connectivity)

        # Construct the list of param dictionaries for all simulations
        list_dicts=[]
        list_params=[]
        for key,val in param_dict.items():
            if (type(val)==ndarray)or((type(val)==list)and key[:8]!='recorded'):
                for param in val:
                    tmp_dict=param_dict.copy()
                    tmp_dict[key]=param
                    for i in range(number_of_runs):
                        list_dicts.append(tmp_dict)
                        list_params.append(param)
        
        # Third run multiple networks
        import multiprocessing as mp
        from functools import partial
        pool=mp.Pool(10)
        print "Start simulations (no idea of how long it will be)"
        results=pool.map(partial(reseau_mitral_granule,return_details=False,report=None),list_dicts)
        results=[(par,)+rr for par,rr in zip(list_params,results)] # Complete each result with its parameter
        print "finished : ",results

        if save_output:
            print "Saving"
            import pickle
            fid=open(filename+"_multi.out","wb")
            pickle.dump([param_dict,results],fid)
            fid.close()
            print "Saved"

        # Example of analysis from data in results
        distinct_figures=False
        tmp_results=[res[:3] for res in results]  # To keep only LFP and param
        multi_LFP_analysis(tmp_results,osc_th=osc_th,freq_cut=freq_cut,burn_time=burn_time,distinct_figures=distinct_figures)
        # Spectrum plot
        freq_min, freq_max =15.,90.
        LFP_frequency_content(tmp_results,burn_time=burn_time,freq_min=freq_min,freq_max=freq_max)
        
    show()
