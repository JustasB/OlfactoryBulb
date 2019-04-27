function [data, PARS] = MCRI_PG(ORNtrace, ORNsamplingrate, varargin)
    paroverrides = varargin;
    
    %% set up params
    PARS = MC_pars;
	PARS = PG_pars(PARS);
    PARS = InhibSyn_pars(PARS,'PGMC');
    
    MCVinit = -66.5;

    PARS.PGMC_gSyn = 10;
    PARS.MCGC_g_syn = 0;
    PARS.ORNGain = 15;
    PARS.ORNPG_gain = 15;

    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    
    ORNtrace = ORNtrace .* PARS.ORNGain;
    PARS.ORNPG_gain = PARS.ORNPG_gain ./ PARS.ORNGain;
    
    %% set initial conditions
    N = 0;
    ics = zeros(14,1);
	
    ics(N+1) = MCVinit;
    ics(N+2:N+3) = MCNaChanInit(ics(N+1));
    ics(N+4:N+5) = MCKfastChanInit(ics(N+1));
    ics(N+6:N+7) = MCKaChanInit(ics(N+1));
    ics(N+8:N+9) = MCKslowChanInit(ics(N+1));
    ics(N+10) = ActiveSynInit(ics(N+1), PARS.PGMC_vHalf, PARS.PGMC_kAct, ...
        PARS.PGMC_alpha, PARS.PGMC_beta);
    
    N = 10;
    ics(N+1) = -69;
    ics(N+2) = PGKspChanInit(ics(N+1));
    ics(N+3:N+4) = PGKaChanInit(ics(N+1));
    
    odeopts = odeset('Events',@spikedetect_RI);
    
    data = integrator_RI(@vfield_MCRI_PG, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect_RI(~,y,~)
    value = [y(11) y(1)]; %detect Vm of 0
    isterminal = [0 1]; %don't stop for PG, just MC
    direction = [1 -1]; %for MC, only downward crossings (~1 ms later for MC->GC recurrent inhib.)
end