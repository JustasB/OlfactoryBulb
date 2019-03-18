function pars = ET_pars_long(pars)
%%% Creates and returns struct with parameters for ET model, long period
%%% (1000 ms) bursting

if nargin < 1
    pars = struct;
end

%% Reversal potentials
pars.ET_vL = -85; %originally -62.5;
pars.ET_vNa = 45;
pars.ET_vK = -105;
pars.ET_vH = -35;
pars.ET_vCa = 120;

%% Conductances
pars.ET_gL = 2.5; %%%
pars.ET_gNa = 29.17;
pars.ET_gK = 12.96;
pars.ET_gH = exp(2); %%%
pars.ET_gCaT = 11.59; %%%
pars.ET_gNaP = 6.79; %%%

%% Na
pars.ET_theta_mNa = -20;
pars.ET_sigma_mNa = -5;

%% K
pars.ET_theta_nK = -26;
pars.ET_sigma_nK = -9;
pars.ET_tau_nK = 10;

%% CaT
pars.ET_theta_mCaT = -37.1;
pars.ET_sigma_mCaT = -4.9805; %%%
pars.ET_tau_mCaT = 11;
pars.ET_theta_hCaT = -59.2;
pars.ET_sigma_hCaT = 14.0684; %%%
pars.ET_tau_hCaT = 82;

%% NaP
pars.ET_theta_mNaP = -40; %%%
pars.ET_sigma_mNaP = -6;
pars.ET_theta_hNaP = -54;
pars.ET_sigma_hNaP = 6; %%%
pars.ET_tau_hNaP = 500;

%% H
pars.ET_theta_mH = -89.32;
pars.ET_sigma_mH = 20.855;
pars.ET_tau_mH_T = 773.1;
pars.ET_delta_mH_T = 0.205;
pars.ET_theta_mH_T = -65.95;
pars.ET_sigma_mH_T = 4.44;

%% External current and capacitance
pars.ET_C = 21;
pars.ET_Iext = 0;
