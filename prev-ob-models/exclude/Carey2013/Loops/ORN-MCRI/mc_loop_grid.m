function mc_loop(input)

% MC_alone:
%     paroverrides = {'ORNGain',[2.5:2.5:30], 'MCGC_g_syn',[0], 'MC_gL', 0.025:0.025:0.2, 'MC_Eleak', -90:5:-55};
%     [g1, g2] = ind2sub([12, 8], input);
%     paroverrides{2} = paroverrides{2}(g1);
%     paroverrides{6} = paroverrides{6}(g2);

% MCRI:
    paroverrides = {'ORNGain',[5:5:30], 'MCGC_g_syn',[0 0.025 0.05 0.075 0.10 0.15 0.20 0.25], 'MCGC_T_decay', [150 500]};
    [g1, g2] = ind2sub([6, 2], input);
    paroverrides{2} = paroverrides{2}(g1);
    paroverrides{6} = paroverrides{6}(g2);
    
    doloop_grid('orn_inputs_depr_shortepsc_1sniff', @MCRI, input, paroverrides{:})
    
end
