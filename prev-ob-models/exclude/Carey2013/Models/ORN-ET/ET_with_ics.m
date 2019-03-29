function data = ET_with_ics(ORNtrace, ORNsamplingrate, inputgain, gL, vL, ics)
    
    %% set up params
	PARS = ET_pars_short;
	
    PARS.ET.gL = gL;
    PARS.ET.vL = vL;
    
    %% set initial conditions
    odeopts = odeset('Events',@spikedetect);
    
    data = integrator_noRI(@vfield_ET, odeopts, ics, PARS, ORNtrace .* inputgain, ORNsamplingrate);
    data.gL = gL;
    data.vL = vL;
    data.gain = inputgain;
end

function [value,isterminal,direction] = spikedetect(~,y,~)
    value = y(1); %detect Vm of 0
    isterminal = 0; %don't stop for ET
    direction = 1;
end