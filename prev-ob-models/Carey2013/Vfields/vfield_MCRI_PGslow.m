function xdot = vfield_MCRI_PGslow(t,x,p)
% Vector field for the MCRI-PG model with slow synapse.
%
%   INPUTS:
%   t -- current time
%   x -- (14,1) vector of current vector values
%   p -- struct containing parameter values in p.param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time

xdot = zeros(size(x));

%% get ORN input for t
Input = p.intrace(floor(t / p.ORNsamplingfactor)+1);

%% PG cell

N = 11;
PG_v = x(N+1);
PG_N = x(N+2);
PG_mKa = x(N+3);
PG_hKa = x(N+4);

PG_minf = 0.1 .* (PG_v + 40) ./ (0.1 .* (PG_v + 40) + 4 .* exp(-(PG_v + 65) ./ 18) .* (1 - exp(-(PG_v + 40) ./ 10)));

PG_INa = p.PG_gNa .* PG_minf.^3 .* (0.8 - PG_N) .* (PG_v - p.PG_ENa);
PG_IK = p.PG_gK .* PG_N.^4 .* (PG_v - p.PG_EK);

PG_IKa = p.PG_gK .* PG_mKa.^3 .* PG_hKa .* (PG_v - p.PG_EK);

PG_IL = p.PG_gL .* (PG_v - p.PG_EL);

PG_Isyn = -p.ORNPG_gain*Input;

%dV/dt:
xdot(N+1) = p.PG_Iext - PG_INa - PG_IK - PG_IKa - PG_IL - PG_Isyn;
xdot(N+2) = PGKspChan(PG_v, PG_N);
xdot(N+3:N+4) = PGKaChan(PG_v, PG_mKa, PG_hKa);


%% MC




N = 0;
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
MC_synPGr = x(N+10);
MC_synPGg = x(N+11);


MC_INa = p.MC_gNa * MC_mNa^3 * MC_hNa * (MC_V - p.MC_ENa);
MC_INaP = p.MC_gNaP * MCNaPChan(MC_V) * (MC_V - p.MC_ENa);
MC_IKa = p.MC_gKa * MC_mKa * MC_hKa * (MC_V - p.MC_EK);
MC_IKfast = p.MC_gKfast * MC_mKfast^2 * MC_hKfast * (MC_V - p.MC_EK);
MC_IKslow = p.MC_gKslow * MC_mKslow * MC_hKslow * (MC_V - p.MC_EK);
MC_IL = p.MC_gL*(MC_V - p.MC_Eleak);

MC_IsynPG = p.PGMCS_gSyn*SlowSynConduct(MC_synPGg, p.PGMCS_sp, 4)*(MC_V-p.PGMCS_vRev);

% Recurrent inhibition for MC cell depends on time since last MC spike
tslp = t - p.last_MC_spike_time; %time since last spike
MC_GC_recurrent = (p.current_MC_recurrent_inhibition_decay_amplitude*exp(-tslp / p.MCGC_T_decay) - exp(-tslp / p.MCGC_T_rise)) .* p.MCGC_g_syn .* (MC_V - p.MCGC_V_reverse);




xdot(N+1) = -(MC_INa + MC_INaP + MC_IKa + MC_IKfast + MC_IKslow + MC_IL + p.MC_Iext - Input + MC_IsynPG + MC_GC_recurrent);
xdot(N+2:N+3) = MCNaChan(MC_V, MC_mNa, MC_hNa);
xdot(N+4:N+5) = MCKfastChan(MC_V, MC_mKfast, MC_hKfast);
xdot(N+6:N+7) = MCKaChan(MC_V, MC_mKa, MC_hKa);
xdot(N+8:N+9) = MCKslowChan(MC_V, MC_mKslow, MC_hKslow);

% Slow synapse variables
% Slow synapse also depends on time since last spike?
xdot(N+10:N+11) = SlowSyn(PG_v, [MC_synPGr MC_synPGg], p.PGMCS_sp);

