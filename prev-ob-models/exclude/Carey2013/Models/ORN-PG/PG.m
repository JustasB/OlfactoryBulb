function data = PG(ORNtrace, ORNsamplingrate, varargin)
    
    %% set up params
    PARS = PG_pars;
    PARS.PG.Iext = 17.5;
    
    paroverrides = varargin;
    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    
    %% set initial conditions
    %ics = [-69 0.29 0.18 0.49];
    ics = zeros(1,4);
    ics(1) = -69;
    ics(2) = PGKspChanInit(ics(1));
    ics(3:4) = PGKaChanInit(ics(1));
    
    odeopts = odeset('Events',@spikedetect);
    
    data = integrator_noRI(@vfield_PG, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect(~,y,~)
    value = y(1); %detect Vm of 0
    isterminal = 0; %don't stop for ET
    direction = 1;
end