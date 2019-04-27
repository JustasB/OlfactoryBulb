function pars = MC_pars(pars)
%%% Creates and returns struct with parameters for MC model

if nargin < 1
    pars = struct;
end

%% Conductances/capacitance (S/mF)
Cm = 10;  % membrane capacitance, in mF so time unit is ms
pars.MC_gKa = 100 ./ Cm;
pars.MC_gNa = 500 ./ Cm;
pars.MC_gKfast = 500 ./ Cm;
pars.MC_gKslow = 310 ./ Cm;
pars.MC_gNaP = 1.1 ./ Cm;
pars.MC_gL = 0.1; % originally: 1 ./ (10 * Cm); %1/(Rm*Cm); Rm = 10 ohm*m^2

%% Reversal Potentials (mV)
pars.MC_ENa = 45;
pars.MC_EK = -70;
pars.MC_Eleak = -65; %originally: 66.5;	% "leak" potential (mV)

%% External current and capacitance
pars.MC_Iext = 0;

%% Recurrent inhibition
pars.MCGC_T_decay = 50;
pars.MCGC_T_rise = 1;
%p.MCGC_T_latency = 1; %taken care of in MC_spikedetect
pars.MCGC_g_syn = 32 ./ Cm;
pars.MCGC_V_reverse = -70;

%v0 = -66.5;	% initial membrane potential (mV)

