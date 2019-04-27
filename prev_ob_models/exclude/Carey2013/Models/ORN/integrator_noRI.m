function data = integrator_noRI(odefun, odeopts, ics, PARS, ORNtrace, ORNsamplingrate)
    
    ics = ics(:);
    
    PARS.intrace = ORNtrace;
    PARS.ORNsamplingfactor = 1000 / ORNsamplingrate;
    
    %% Integration tolerances
    atol = 1e-6;
    rtol = 1e-6;
    options = odeset('reltol',rtol,'abstol',atol,'InitialStep',0.5);
    options = odeset(options, odeopts); %override odeopts
    
    %% time span
    tstart = 0;
    tend = (length(ORNtrace) - 1) * PARS.ORNsamplingfactor;
    
    %% integrate
    tout = tstart;
    xout = ics';
    teout = [];
    ieout = [];
    
    tic
    while tstart < tend
        [t,x,te,xe,ie] = ode15s(odefun,[tstart tend],ics,options,PARS);
        
        nt = length(t);
        tout = [tout; t(2:nt)];
        xout = [xout; x(2:nt,:)];
        tstart = t(nt);
        
        if ~isempty(te)
            teout = [teout; te];
            ieout = [ieout; ie];
            
            ics = xe(end,:);
            
            options = odeset(options,'InitialStep',t(nt)-t(nt-1)); %use most recent timestep
        end
    end
    if isfield(PARS,'save_traces')
        data = struct('T', tout, 'X', xout, 'spikes', teout, 'which', ieout, 'pars', PARS);
    else
        data = struct('spikes', teout, 'which', ieout);
    end
    toc

end
