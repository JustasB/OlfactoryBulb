function etmcp_loop(input)
    
    paroverrides = {'ORNGain',[1], 'ET_gORN', 1.5*21, 'ES_gSyn', 0:0.025:0.5, 'MC_gORN', (0:0.025:0.5).*100./3.9 };

    paroverrides{6} = paroverrides{6}(input);
    
    doloop_grid('orn_inputs_depr_shortepsc_1sniff', @ET_MCRI_pexcite, input, 'calc_pexcite_charges', true, 'MCGC_g_syn',0, paroverrides{:})
    
end
