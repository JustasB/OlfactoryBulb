function data = ET(ORNtrace, ORNsamplingrate, varargin)
    
    %% set up params
	PARS = ET_pars_long;
	
    PARS.ES_gSyn = 1;
    
    paroverrides = varargin;
    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    if isfield(PARS,'ORNGain')
        ORNtrace = ORNtrace .* PARS.ORNGain;
    end
    
    %% set initial conditions
    %{
    ics = zeros(6,1);
    ETVinit = -50.6890;
    ics(1) = ETVinit;
    ics(2) = ETKChanInit(ics(1), PARS);
    ics(3) = ETNaPChanInit(ics(1), PARS);
    ics(4:5) = ETCaTChanInit(ics(1), PARS);
    ics(6) = ETHChanInit(ics(1), PARS);
	%}
    
    ics = [-64.1084    0.0143    0.8435    0.2299    0.0044    0.5863]'; %for gL = 2.6, vL = -85
    
    odeopts = odeset('Events',@spikedetect);
    
    data = integrator_noRI(@vfield_ET, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect(~,y,~)
    value = y(1); %detect Vm of 0
    isterminal = 0; %don't stop for ET
    direction = 1;
end