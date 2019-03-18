function data = ET_MCRI_pexcite(ORNtrace, ORNsamplingrate, varargin)
    
    %% set up params
    PARS = MC_pars;
	PARS = ET_pars_long(PARS);
	PARS = ExciteSyn_pars(PARS);
	
    MCVinit = -66.5;

    PARS.MCGC_g_syn = 0;
    PARS.ET_gORN = 0.5;
    PARS.ES_gSyn = 0;
    PARS.MC_gORN = 0.15;

    paroverrides = varargin;
    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    if isfield(PARS,'ORNGain')
        ORNtrace = ORNtrace .* PARS.ORNGain;
    end
    
    %% set initial conditions
    N = 6;
    ics = zeros(16,1);
    ics(1:6) = [-64.1084    0.0143    0.8435    0.0044    0.5863    0.2299]; %for gL = 2.6, vL = -85
	
    ics(N+1) = MCVinit;
    ics(N+2:N+3) = MCNaChanInit(ics(N+1));
    ics(N+4:N+5) = MCKfastChanInit(ics(N+1));
    ics(N+6:N+7) = MCKaChanInit(ics(N+1));
    ics(N+8:N+9) = MCKslowChanInit(ics(N+1));
    ics(N+10) = ActiveSynInit(ics(N+1), PARS.ES_vHalf, PARS.ES_kAct, ...
        PARS.ES_alpha, PARS.ES_beta);
    
    odeopts = odeset('Events',@spikedetect_RI);
    
    data = integrator_RI(@vfield_ET_MCRI_pexcite, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect_RI(~,y,~)
    value = [y(1) y(7)]; %detect Vm of 0
    isterminal = [0 1]; %don't stop for ET, just MC
    direction = [1 -1]; %for MC, only downward crossings (~1 ms later for MC->GC recurrent inhib.)
end