function [data, PARS] = ET_MCRI_pPGslow_pexcite(ORNtrace, ORNsamplingrate, varargin)
    paroverrides = varargin;
    
    %% set up params
    PARS = MC_pars;
	PARS = ET_pars_long(PARS);
	PARS = PG_pars(PARS);
	PARS = ExciteSyn_pars(PARS,'ETMC');
    PARS = ExciteSyn_pars(PARS,'ETPG');
    
    %inhibitory slow synapse parameters are loaded after PGMCS_tc = {'140', '140b', '170', '170b', '200', '200b'} is set
    [tf, p] = ismember('PGMCS_tc',paroverrides(1:2:end));
    if tf
        PARS.PGMCS_tc = paroverrides{2*p}{1};
    else
        PARS.PGMCS_tc = '140';
    end
    
    SS = load(['SS' PARS.PGMCS_tc '_pars.mat']);
    PARS = InhibSlowSyn_pars(SS.SS,PARS,'PGMCS');
    
    %PARS = InhibSyn_pars(PARS,'PGMC');
    
    MCVinit = -66.5;

    PARS.ORNGain = 100;

    % Synaptic INPUT strengths -- multiply the input signal
    PARS.ORNET_gSyn = 0.8; %from ORN to ET (original ORN-ET calib)  - divided by 21 in vfield
    PARS.ORNPG_gSyn = 1; %this value makes the PG fire earlier - makes sense in order to suppress weak inputs)
    PARS.ORNMC_gSyn = 0.02; %from ORN to MC (along 50-50 line @ ~140 Hz) (doubled to account for RI)
    
    % FF excitatory:
    PARS.ETMC_gSyn = 2; %from ET to MC (along 50-50 line @ ~140 Hz) (doubled to account for RI)
    PARS.ETPG_gSyn = 2; %from ET to PG (original thesis value)
    
    % FF inhib:
    PARS.PG1MCS_gSyn = 0.5; %(ET-PG-MC)
    PARS.PG2MCS_gSyn = 0.5; %(ORN-PG-MC)
    
    % RI:
    PARS.MCGC_g_syn = 0.8;

    for p = 1:2:length(paroverrides)
        PARS.(paroverrides{p}) = paroverrides{p+1};
    end
    
    [tf, p] = ismember('PGMCS_gSyn',paroverrides(1:2:end));
    if tf
        % set both FF inhib:
        PARS.PG1MCS_gSyn = paroverrides{2*p}; %(ET-PG-MC)
        PARS.PG2MCS_gSyn = paroverrides{2*p}; %(ORN-PG-MC)
    end
    

    if isfield(PARS,'ORNGain')
        ORNtrace = ORNtrace .* PARS.ORNGain;
    end
    
    %% set initial conditions
    ics = zeros(29,1);
    
    % ET Cell
    ics(1:6) = [-64.1084    0.0143    0.8435    0.0044    0.5863    0.2299]; %for gL = 2.6, vL = -85
	
    N = 6;
    % MC cell
    ics(N+1) = MCVinit;
    ics(N+2:N+3) = MCNaChanInit(ics(N+1));
    ics(N+4:N+5) = MCKfastChanInit(ics(N+1));
    ics(N+6:N+7) = MCKaChanInit(ics(N+1));
    ics(N+8:N+9) = MCKslowChanInit(ics(N+1));
    
    % synapses onto MC:
    ics(N+10) = ActiveSynInit(ics(N+1), PARS.ETMC_vHalf, PARS.ETMC_kAct, ...
        PARS.ETMC_alpha, PARS.ETMC_beta); %excitatory from ET
    ics(N+11:N+12) = SlowSynInit(0,0,0); %inhib from PG#1
    ics(N+13:N+14) = SlowSynInit(0,0,0); %inhib from PG#2
    
    N = 20;
    % PG Cell #1, between ET and MC
    ics(N+1) = -69;
    ics(N+2) = PGKspChanInit(ics(N+1));
    ics(N+3:N+4) = PGKaChanInit(ics(N+1));
    
    % excitatory synapse between ET and PG#1:
    ics(N+5) = ActiveSynInit(ics(N+1), PARS.ETPG_vHalf, PARS.ETPG_kAct, ...
        PARS.ETPG_alpha, PARS.ETPG_beta);

    
    N = 25;
    % PG Cell #2, between ORN and MC
    ics(N+1) = -69;
    ics(N+2) = PGKspChanInit(ics(N+1));
    ics(N+3:N+4) = PGKaChanInit(ics(N+1));
    

    odeopts = odeset('Events',@spikedetect_RI);
    
    data = integrator_RI(@vfield_ET_MCRI_pPGslow_pexcite, odeopts, ics, PARS, ORNtrace, ORNsamplingrate);
end

function [value,isterminal,direction] = spikedetect_RI(~,y,~)
    value = [y(1) y(7) y(21) y(26)]; %detect Vm of 0
    isterminal = [0 1 0 0]; %don't stop for ET or either PG, just MC
    direction = [1 -1 1 1]; %for MC, only downward crossings (~1 ms later for MC->GC recurrent inhib.)
end