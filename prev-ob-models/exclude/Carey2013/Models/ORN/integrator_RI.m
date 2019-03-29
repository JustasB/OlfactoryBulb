function data = integrator_RI(odefun, odeopts, ics, PARS, ORNtrace, ORNsamplingrate)
    
    ics = ics(:);
    
    PARS.intrace = ORNtrace;
    PARS.ORNsamplingfactor = 1000 / ORNsamplingrate;
    
    PARS.last_MC_spike_time = -Inf;
    PARS.current_MC_recurrent_inhibition_decay_amplitude = 0;
    
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
            
            tslp = tstart - PARS.last_MC_spike_time; %previous ISI
            current_MC_recurrent_inhib_value = PARS.current_MC_recurrent_inhibition_decay_amplitude*exp(-tslp / PARS.MCGC_T_decay) - exp(-tslp / PARS.MCGC_T_rise);
            
            PARS.last_MC_spike_time = tstart;
            PARS.current_MC_recurrent_inhibition_decay_amplitude = 1 + current_MC_recurrent_inhib_value;
            
            options = odeset(options,'InitialStep',t(nt)-t(nt-1)); %use most recent timestep
        end
    end
    if isfield(PARS,'save_traces')
        data = struct('T', tout, 'X', xout, 'spikes', teout, 'which', ieout, 'pars', PARS);
    else
        data = struct('spikes', teout, 'which', ieout);
    end
    if isfield(PARS,'calc_pexcite_charges')
        MC_V = xout(:,7);
        MC_syn = xout(:,16);
        I_ETMC = (PARS.ES_gSyn.*MC_syn.*(MC_V-PARS.ES_vRev));
        data.ETMCcharge = sum(I_ETMC(2:end).*diff(tout));
        
        I_ORNMC = -PARS.MC_gORN*ORNtrace;
        data.ORNMCcharge = sum(I_ORNMC(2:end).*(1000/ORNsamplingrate));
    end
    toc
end
