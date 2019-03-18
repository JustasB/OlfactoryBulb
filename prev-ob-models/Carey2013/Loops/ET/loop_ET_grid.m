function loop_ET_grid(inputgain, inputtrace)
    S_ET = load('ET4_starts_subset');
    
    S_orn = load(inputtrace);
    
    nICS = length(S_ET.gL);
    
    data = struct('spikes', cell(1,nICS), 'which', cell(1,nICS), 'gL', cell(1,nICS), 'vL', cell(1,nICS), 'gain', cell(1,nICS));
    
    total_start = tic;
    
    for i = 1:nICS
        disp(['Integrating ic ' num2str(i) '...'])
        data(i) = ET_with_ics(S_orn.trace, S_orn.samplingrate, inputgain, S_ET.gL(i), S_ET.vL(i), S_ET.all_ics(i,:));
        disp(['... finished ic ' num2str(i) '.'])
    end
    
    disp('Total time ...');
    toc(total_start)
    
    save(['ET4_grid_' num2str(inputgain) '_' inputtrace],'data')
    
end