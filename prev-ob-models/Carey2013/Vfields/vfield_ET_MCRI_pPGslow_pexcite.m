function xdot = vfield_ET_MCRI_pPGslow_pexcite(t,x,p)
% Vector field for the ET-MCRI-PG model with two parallel slow PG synapses, parallel excitation.
%
%   INPUTS:
%   t -- current time
%   x -- (29,1) vector of current vector values
%   p -- struct containing parameter values in p.param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time



%% get ORN input for t
Input = p.intrace(floor(t / p.ORNsamplingfactor)+1);

%% ET cell

ET_V = x(1);
ET_nK = x(2);
ET_hNaP = x(3);
ET_mH = x(4);
ET_mCaT = x(5);
ET_hCaT = x(6);

% compute values for the currents
ET_mNa_inf = calc_xinf(ET_V, p.ET_theta_mNa, p.ET_sigma_mNa);
ET_INa = p.ET_gNa.*(1-ET_nK)*ET_mNa_inf.^3*(ET_V-p.ET_vNa);
ET_IK = p.ET_gK.*ET_nK.^4.*(ET_V-p.ET_vK);
ET_IL = p.ET_gL.*(ET_V-p.ET_vL);
ET_IH = p.ET_gH.*ET_mH.*(ET_V-p.ET_vH);
ET_mNaP_inf = calc_xinf(ET_V, p.ET_theta_mNaP, p.ET_sigma_mNaP);
ET_INaP = p.ET_gNaP.*ET_mNaP_inf*ET_hNaP*(ET_V-p.ET_vNa);
ET_ICaT = p.ET_gCaT.*ET_mCaT.^2.*ET_hCaT*(ET_V-p.ET_vCa);


xdot = zeros(size(x));
xdot(1) = -(ET_INa + ET_IK + ET_ICaT + ET_IH + ET_INaP + ET_IL + p.ET_Iext - p.ORNET_gSyn*Input)./p.ET_C;
xdot(2) = calc_xdot(ET_V,ET_nK,p.ET_theta_nK, p.ET_sigma_nK,p.ET_tau_nK);
xdot(3) = calc_xdot(ET_V, ET_hNaP, p.ET_theta_hNaP, p.ET_sigma_hNaP, p.ET_tau_hNaP);

ET_mH_pars = [p.ET_theta_mH, p.ET_sigma_mH, p.ET_theta_mH_T, p.ET_sigma_mH_T, ...
    p.ET_tau_mH_T, p.ET_delta_mH_T];
xdot(4) = calc_H_mdot(ET_V, ET_mH, ET_mH_pars);

xdot(5) = calc_xdot(ET_V,ET_mCaT,p.ET_theta_mCaT, p.ET_sigma_mCaT,p.ET_tau_mCaT);
xdot(6) = calc_xdot(ET_V,ET_hCaT,p.ET_theta_hCaT, p.ET_sigma_hCaT,p.ET_tau_hCaT);

%% PG cell #1 between ET and MC

N = 20;
PG_v = x(N+1);
PG_N = x(N+2);
PG_mKa = x(N+3);
PG_hKa = x(N+4);
PG_syn = x(N+5);

PG_minf = 0.1 .* (PG_v + 40) ./ (0.1 .* (PG_v + 40) + 4 .* exp(-(PG_v + 65) ./ 18) .* (1 - exp(-(PG_v + 40) ./ 10)));

PG_INa = p.PG_gNa .* PG_minf.^3 .* (0.8 - PG_N) .* (PG_v - p.PG_ENa);
PG_IK = p.PG_gK .* PG_N.^4 .* (PG_v - p.PG_EK);

PG_IKa = p.PG_gK .* PG_mKa.^3 .* PG_hKa .* (PG_v - p.PG_EK);

PG_IL = p.PG_gL .* (PG_v - p.PG_EL);

PG_Isyn = p.ETPG_gSyn*PG_syn*(PG_v-p.ETPG_vRev);

%dV/dt:
xdot(N+1) = p.PG_Iext - PG_INa - PG_IK - PG_IKa - PG_IL - PG_Isyn;
xdot(N+2) = PGKspChan(PG_v, PG_N);
xdot(N+3:N+4) = PGKaChan(PG_v, PG_mKa, PG_hKa);

xdot(N+5) = ActiveSyn(PG_v, ET_V, PG_syn, [p.ETPG_vHalf, p.ETPG_kAct, ...
    p.ETPG_alpha, p.ETPG_beta]);


%% PG cell #2 between ORN and MC

N = 25;
PGb_v = x(N+1);
PGb_N = x(N+2);
PGb_mKa = x(N+3);
PGb_hKa = x(N+4);

PGb_minf = 0.1 .* (PGb_v + 40) ./ (0.1 .* (PGb_v + 40) + 4 .* exp(-(PGb_v + 65) ./ 18) .* (1 - exp(-(PGb_v + 40) ./ 10)));

PGb_INa = p.PG_gNa .* PGb_minf.^3 .* (0.8 - PGb_N) .* (PGb_v - p.PG_ENa);
PGb_IK = p.PG_gK .* PGb_N.^4 .* (PGb_v - p.PG_EK);

PGb_IKa = p.PG_gK .* PGb_mKa.^3 .* PGb_hKa .* (PGb_v - p.PG_EK);

PGb_IL = p.PG_gL .* (PGb_v - p.PG_EL);


%dV/dt:
xdot(N+1) = p.PG_Iext - PGb_INa - PGb_IK - PGb_IKa - PGb_IL + p.ORNPG_gSyn*Input; % Note that this is ORN input to PG scaled by a factor!
xdot(N+2) = PGKspChan(PGb_v, PGb_N);
xdot(N+3:N+4) = PGKaChan(PGb_v, PGb_mKa, PGb_hKa);


%% MC

N = 6;
MC_V = x(N+1);
MC_mNa = x(N+2);
MC_hNa = x(N+3);
MC_mKfast = x(N+4);
MC_hKfast = x(N+5);
MC_mKa = x(N+6);
MC_hKa = x(N+7);
MC_mKslow = x(N+8);
MC_hKslow = x(N+9);
%MC_dep = x(N+10);
MC_synET = x(N+10);
MC_synPGr = x(N+11); % slow
MC_synPGg = x(N+12); % slow
MC_synPGrB = x(N+13);% slow  
MC_synPGgB = x(N+14);% slow


MC_INa = p.MC_gNa * MC_mNa^3 * MC_hNa * (MC_V - p.MC_ENa);
MC_INaP = p.MC_gNaP * MCNaPChan(MC_V) * (MC_V - p.MC_ENa);
MC_IKa = p.MC_gKa * MC_mKa * MC_hKa * (MC_V - p.MC_EK);
MC_IKfast = p.MC_gKfast * MC_mKfast^2 * MC_hKfast * (MC_V - p.MC_EK);
MC_IKslow = p.MC_gKslow * MC_mKslow * MC_hKslow * (MC_V - p.MC_EK);
MC_IL = p.MC_gL*(MC_V - p.MC_Eleak);

MC_IsynET = p.ETMC_gSyn*MC_synET*(MC_V-p.ETMC_vRev);
%MC_IsynPG = p.PGMC_gSyn*MC_synPG*(MC_V-p.PGMC_vRev);
MC_IsynPG = p.PGMCS_gSyn*SlowSynConduct(MC_synPGg, p.PGMCS_sp, 4)*(MC_V-p.PGMCS_vRev);
MC_IsynPGB = p.PGMCS_gSynB*SlowSynConduct(MC_synPGgB, p.PGMCS_sp, 4)*(MC_V-p.PGMCS_vRev);

% Recurrent inhibition for MC cell depends on time since last MC spike
tslp = t - p.last_MC_spike_time; %time since last spike
MC_GC_recurrent = (p.current_MC_recurrent_inhibition_decay_amplitude*exp(-tslp / p.MCGC_T_decay) - exp(-tslp / p.MCGC_T_rise)) .* p.MCGC_g_syn .* (MC_V - p.MCGC_V_reverse);

xdot(N+1) = -(MC_INa + MC_INaP + MC_IKa + MC_IKfast + MC_IKslow + MC_IL + p.MC_Iext + MC_IsynET + MC_IsynPG + MC_IsynPGB + MC_GC_recurrent - p.ORNMC_gSyn*Input);
xdot(N+2:N+3) = MCNaChan(MC_V, MC_mNa, MC_hNa);
xdot(N+4:N+5) = MCKfastChan(MC_V, MC_mKfast, MC_hKfast);
xdot(N+6:N+7) = MCKaChan(MC_V, MC_mKa, MC_hKa);
xdot(N+8:N+9) = MCKslowChan(MC_V, MC_mKslow, MC_hKslow);

xdot(N+10) = ActiveSyn(MC_V, ET_V, MC_synET, [p.ETMC_vHalf, p.ETMC_kAct, ...
    p.ETMC_alpha, p.ETMC_beta]);
%xdot(N+11) = ActiveSyn(MC_V, PG_v, MC_synPG, [p.PGMC_vHalf, p.PGMC_kAct, ...
%    p.PGMC_alpha, p.PGMC_beta]);

% Slow synapse variables
% Slow synapse also depends on time since last spike?
xdot(N+11:N+12) = SlowSyn(PG_v, [MC_synPGr MC_synPGg], p.PGMCS_sp);
xdot(N+13:N+14) = SlowSyn(PGb_v, [MC_synPGrB MC_synPGgB], p.PGMCS_sp);
