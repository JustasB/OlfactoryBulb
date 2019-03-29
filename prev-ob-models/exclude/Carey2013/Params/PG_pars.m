function pars = PG_pars(pars)
%%% Creates and returns struct with parameters for PG model

if nargin < 1
    pars = struct;
end

%% Conductances (mS/uF) (mS/cm^2, already divided by Cm = 1 uF/cm^2)
pars.PG_gK = 36;
pars.PG_gNa = 120;
pars.PG_gL = 0.3;

pars.PG_gKa = 10;

%% Reversal Potentials (mV)
pars.PG_ENa = 50; %orig 50
pars.PG_EK = -77; %orig -77
pars.PG_EL = -54.4; %leak, orig -54.4


%%
pars.PG_Iext = 0;

%v0 = -66.5;	% initial membrane potential (mV)
