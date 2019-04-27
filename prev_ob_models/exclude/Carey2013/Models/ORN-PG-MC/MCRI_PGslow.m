function [data, PARS] = MCRI_PGslow(ORNtrace, ORNsamplingrate, varargin)
    paroverrides = varargin;
        
    %% set up params
    PARS = MC_pars;
	PARS = PG_pars(PARS);
    
    %inhibitory slow synapse parameters are loaded after PGMCS_tc = {'140', '140b', '170', '170b', '200', '200b'} is set
    [tf, p] = ismember('PGMCS_tc',paroverrides(1:2:end));
    if tf
        PARS.PGMCS_tc = paroverrides{2*p};
    else
        PARS.PGMCS_tc = '140';
    end
    
    SS = load(['SS' PARS.PGMCS_tc '_pars.mat']);
    PARS = InhibSlowSyn_pars(SS.SS,PARS,'PGMCS');
    
    MCVinit = -66.5;

    PARS.PGMCS_gSyn = 10;
    PARS.MCGC_g_syn = 0;
    PARS.ORNGain = 15;
    PARS.ORNPG_gain = 15;
    
    %load the rest of the overrides here:
    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    
    ORNtrace = ORNtrace .* PARS.ORNGain;
    PARS.ORNPG_gain = PARS.ORNPG_gain ./ PARS.ORNGain;
    
    %% set initial conditions
    N = 0;
    ics = zeros(15,1);
	
    ics(N+1) = MCVinit;
    ics(N+2:N+3) = MCNaChanInit(ics(N+1));
    ics(N+4:N+5) = MCKfastChanInit(ics(N+1));
    ics(N+6:N+7) = MCKaChanInit(ics(N+1));
    ics(N+8:N+9) = MCKslowChanInit(ics(N+1));
    ics(N+10:N+11) = SlowSynInit(0,0,0);
    
    N = 11;
    ics(N+1) = -69;
    ics(N+2) = PGKspChanInit(ics(N+1));
    ics(N+3:N+4) = PGKaChanInit(ics(N+1));
    
    odeopts = odeset('Events',@spikedetect_RI);
    
    data = integrator_RI(@vfield_MCRI_PGslow, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect_RI(~,y,~)
    value = [y(12) y(1)]; %detect Vm of 0
    isterminal = [0 1]; %don't stop for PG, just MC
    direction = [1 -1]; %for MC, only downward crossings (~1 ms later for MC->GC recurrent inhib.)
end
