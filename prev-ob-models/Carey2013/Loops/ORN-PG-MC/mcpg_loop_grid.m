function mcpg_loop_grid(idx)
    
    paroverrides = {'ORNGain', [10:2:22], 'MCGC_g_syn', 0, 'PGMCS_gSyn', [0.2 0.5 1.0 1.5 2.0 2.5 3.0], 'ORNPG_gain', [20 25 30 35 40], 'PGMCS_tc', {'170'}};
    
    %PG_tc: {'140','140b','170','170b','200','200b'}
    
    [g1, g2] = ind2sub([7, 5], idx);
    paroverrides{6} = paroverrides{6}(g1);
    paroverrides{8} = paroverrides{8}(g2);
    
    doloop_grid('orn_inputs_depr_shortepsc_1sniff', @MCRI_PGslow, idx, paroverrides{:})
    
end
