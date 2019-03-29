function data = MCRI(ORNtrace, ORNsamplingrate, varargin)
    
    %% set up params
    PARS = MC_pars;
    PARS.MCGC_g_syn = 0.05;
    PARS.ORNGain = 12.5;

    paroverrides = varargin;
    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    if isfield(PARS,'ORNGain')
        ORNtrace = ORNtrace .* PARS.ORNGain;
    end
    
    
    %% set initial conditions
    ics = zeros(9,1);
    ics(1) = -66.5;
    ics(2:3) = MCNaChanInit(ics(1));
    ics(4:5) = MCKfastChanInit(ics(1));
    ics(6:7) = MCKaChanInit(ics(1));
    ics(8:9) = MCKslowChanInit(ics(1));
    
    odeopts = odeset('Events',@spikedetect_RI);
    
    data = integrator_RI(@vfield_MCRI, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect_RI(~,y,~)
    value = y(1); %detect Vm of 0
    isterminal = 1; %stop
    direction = -1; %only downward crossings (~1 ms later for MC->GC recurrent inhib.)
end