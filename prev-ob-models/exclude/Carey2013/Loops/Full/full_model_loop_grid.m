function full_model_loop(input)
    
    paroverrides = {...
        'ORNGain', [0.8 1 1.2], ...
        'PGMCS_gSyn', 0.5, ...% FF inhib
        ...'PG1MCS_gSyn', 0.5, ...%(ET-PG-MC)
        ...'PG2MCS_gSyn', 0.5, ...%(ORN-PG-MC)
        'MCGC_g_syn', 0.5, ...
        'ORNPG_gSyn', 0.35*100, ...
        'ORNET_gSyn', 1.5*21, ...%from ORN to ET (original ORN-ET calib)  - divided by 21 in vfield
        'ORNMC_gSyn', 0.30*100/3.9, ...%from ORN to MC (along 50-50 line @ ~140 Hz) (doubled to account for RI)
        'ETMC_gSyn', 0.30, ...%from ET to MC (along 50-50 line @ ~140 Hz) (doubled to account for RI)
        'ETPG_gSyn', 0.35*3.9... %following same "effectiveness" scaling as calculated for pexcite
        };
    
     doloop_grid('orn_inputs_depr_shortepsc_1sniff', @ET_MCRI_pPGslow_pexcite, input, paroverrides{:})
    
end
