# -*- coding: utf-8 -*-
"""
This script provides parameter dictionaries for all article figures.
Note the dictionaries are based on the default dictionary from reseau_mitral_granule_fig_param_dic.py
"""

from brian import *
from reseau_mitral_granule_fig_param_dic import param_dict_default

# Figure 2B-2C
gamma_single_step_dict=param_dict_default.copy()
gamma_single_step_dict['t_init_M']=1.*second

# Figure 2D
gamma_multi_weakinhconnect_dict=param_dict_default.copy()
gamma_multi_weakinhconnect_dict['weakinh_connectproba']=linspace(0.05,1.,20)

# Figure 2E
gamma_multi_numact_dict=param_dict_default.copy()
gamma_multi_numact_dict['num_activated_mitral']=linspace(5,100,20)

# Figure 2F
gamma_multi_weakinhweight_dict=param_dict_default.copy()
gamma_multi_weakinhweight_dict['weakinh_weight']=[x*ms for x in linspace(0.1,2.,20)]

# Figure 2G
gamma_multi_tauGABA_dict=param_dict_default.copy()
gamma_multi_tauGABA_dict['M_tauI']=linspace(0.1,2.,20)

# Figure 2H
gamma_multi_taumKs_dict=param_dict_default.copy()
gamma_multi_taumKs_dict['M_taumKs']=[x*ms for x in linspace(1.,20.,20)]

# Figure 2I (TODO: improve to give final parameters instead of shifts)
gamma_multi_MgEinj_dict=param_dict_default.copy()
gamma_multi_MgEinj_dict['M_gEinj_shift']=linspace(-1.2,2.5,10)

# Figure 3B-C
beta_single_ramp_dict=param_dict_default.copy()
beta_single_ramp_dict['sim_time']=3.5*second
beta_single_ramp_dict['use_granule']=True
beta_single_ramp_dict['t_init_G']=2.*second
beta_single_ramp_dict['G_input']='ramp'
beta_single_ramp_dict['G_I_base']=-3.5*nA
beta_single_ramp_dict['G_I_min']=-3.5*nA
beta_single_ramp_dict['G_I_max']=0.5*nA

# Figure 3D, also used as base for beta_multi
beta_base_dict=param_dict_default.copy()
beta_base_dict['use_granule']=True
beta_base_dict['G_input']='constant'
beta_base_dict['G_I_base']=-1.5*nA

# Figure 3E
beta_multi_tauGABA_dict=beta_base_dict.copy()
beta_multi_tauGABA_dict['M_tauI_G']=[x*ms for x in linspace(3.,16.,20)]

# Figure 3F
beta_multi_stronginhweight_dict=beta_base_dict.copy()
beta_multi_stronginhweight_dict['stronginh_weight']=linspace(.1,3.,20)

# Figure 3G
beta_multi_GIinj_dict=beta_base_dict.copy()
beta_multi_GIinj_dict['G_I_base']=[x*nA for x in linspace(-2.5,.2,20)]

# Figure 3H
beta_multi_excweight_dict=beta_base_dict.copy()
beta_multi_excweight_dict['exc_weight']=linspace(.1,2.,20)

# Figure 4E-4G
competition_single_base_dict=param_dict_default.copy()
competition_single_base_dict['M_gEinj_base']=4.*siemens*meter**-2
competition_single_base_dict['M_gEinj']=(linspace,(6.6,8.1,100),siemens*meter**-2)
competition_single_base_dict['use_granule']=True
competition_single_base_dict['G_input']='sinusoid'

# Figure 5C1-5C2 low intensity
competitionSTP_single_lowintensity_dict=competition_single_base_dict.copy()
competitionSTP_single_lowintensity_dict['use_AMPA_STP']=True
competitionSTP_single_lowintensity_dict['MC_phase_dispersion']=5.

# Figure 5C1-5C2 medium intensity
competitionSTP_single_mediumintensity_dict=competition_single_base_dict.copy()
competitionSTP_single_mediumintensity_dict['use_AMPA_STP']=True
competitionSTP_single_mediumintensity_dict['MC_phase_dispersion']=1.3

# Figure 5C1-5C2 high intensity
competitionSTP_single_highintensity_dict=competition_single_base_dict.copy()
competitionSTP_single_highintensity_dict['use_AMPA_STP']=True
competitionSTP_single_highintensity_dict['MC_phase_dispersion']=0.2

# Figure 5D1-5D2 no STP
competition_multi_intensity=competition_single_base_dict.copy()
competition_multi_intensity['MC_phase_dispersion']=[5., 4., 3., 2.5, 2.0, 1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.]

# Figure 5D1-5D2 with STP
competitionSTP_multi_intensity=competition_single_base_dict.copy()
competitionSTP_multi_intensity['MC_phase_dispersion']=[5., 4., 3., 2.5, 2.0, 1.5, 1.0, 0.5, 0.2, 0.1, 0.05, 0.]
competitionSTP_multi_intensity['use_AMPA_STP']=True

# Figure 6B-6C low centrifugal input
competitionSTP_single_lowcentrifugal_dict=competition_single_base_dict.copy()
competitionSTP_single_lowcentrifugal_dict['use_AMPA_STP']=True
competitionSTP_single_lowcentrifugal_dict['G_I_max']=-0.5*nA

# Figure 6B-6C high centrifugal input
competitionSTP_single_highcentrifugal_dict=competition_single_base_dict.copy()
competitionSTP_single_highcentrifugal_dict['use_AMPA_STP']=True
competitionSTP_single_highcentrifugal_dict['G_I_max']=0.*nA

# Figure 6D-6E high centrifugal input
competitionSTP_multi_centrifugal_dict=competition_single_base_dict.copy()
competitionSTP_multi_centrifugal_dict['use_AMPA_STP']=True
competitionSTP_multi_centrifugal_dict['G_I_max']=[x*nA for x in [-0.6, -0.5, -0.4, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.]]

# Figure 7 high centrifugal, used as a basis dict for the low centrifugal
competitionSTP_single_awake=competition_single_base_dict.copy()
competitionSTP_single_awake['sim_time']=2.*second
competitionSTP_single_awake['freq_modul']=8.*Hz
competitionSTP_single_awake['M_gEinj']=(linspace,(6.1,7.6,100),siemens*meter**-2)
competitionSTP_single_awake['M_gEinj_base']=5.9*siemens*meter**-2
competitionSTP_single_awake['use_AMPA_STP']=True
competitionSTP_single_awake['GC_phase_dispersion']=5.
# Figure 7 low centrifugal
competitionSTP_single_awake_lowcentrifugal=competitionSTP_single_awake.copy()
competitionSTP_single_awake_lowcentrifugal['G_I_max']=-1.5*nA