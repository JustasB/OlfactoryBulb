function [traces_depr, traces_no_depr, deconv_traces] = preprocess_orn_inputs_with_depression(rawtraces)
% Converts raw ORN calcium imaging data to mean ORN firing rate output, 
% removing bleaching and including ORN depression 
    N = size(rawtraces,1);
    num_orn_cells = 100;
    orn_firing_rate = 50;
    
    ep = epsc(0:500,1000); %samplingrate = 1000 (1-ms bins)
    
    deconv_traces = preprocess_remove_bleaching(rawtraces);
    
    traces_no_depr = zeros(N,26000);
    traces_depr = zeros(N,26000);
    
    for n = 1:N
        sr = simulateORN(deconv_traces(n,:), ones(num_orn_cells,1).*orn_firing_rate);
        traces_no_depr(n,:) = convolve_and_chop(mean(sr));
        traces_depr(n,:) = traces_no_depr(n,:);
    end
    
    function outtr = convolve_and_chop(tr)
        outtr = conv(tr,ep);
        outtr = outtr(1:26000);
    end
    
end


