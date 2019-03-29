function etmc_ri_loop(input)

    
    paroverrides = {'ORNGain',[1:5].*21,'ES_gSyn',[0.1 0.2:0.2:2.0],'MCGC_g_syn',[0]};
    paroverrides{4} = paroverrides{4}(input);
    
    doloop_grid('orn_inputs_depr_shortepsc_1sniff', @ET_MCRI, input, paroverrides{:})
    
end
