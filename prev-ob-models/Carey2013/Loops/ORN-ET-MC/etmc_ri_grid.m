function etmc_ri_grid(input)
    
    paroverrides = {'ORNGain',100,'ES_gSyn',[0.1 0.3 0.5],'MCGC_g_syn',[0], 'save_traces', true};
    
    doloop_grid('orn_glom17_depr', @ET_MCRI, input, paroverrides{:})
    
end
