function loop_ET_blank
    S = load('ET4_starts_subset');
    
    %3 seconds of blankness...
    blank_trace = zeros(1, 30);
    blank_sr = 10;
    
    nICS = size(S.gL,1);
    
    data = struct('spikes', cell(1,nICS), 'which', cell(1,nICS), 'gL', cell(1,nICS), 'vL', cell(1,nICS));
    
    total_start = tic;
    
    for i = 1:nICS
        disp(['Integrating ic ' num2str(i) '...'])
        data(i) = ET_with_ics(blank_trace, blank_sr, S.gL(i), S.vL(i), S.all_ics(i,:));
        disp(['... finished ic ' num2str(i) '.'])
    end
    
    disp('Total time ...');
    toc(total_start)
    
    save('ET4_blank','data')
end