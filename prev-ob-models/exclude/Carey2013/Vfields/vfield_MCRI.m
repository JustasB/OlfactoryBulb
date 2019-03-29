function xdot = vfield_MCRI(t,x,p)
% Vector field for the MC model with recurrent inhibition and ORN input.
%
%   INPUTS:
%   t -- current time
%   x -- (9,1) vector of current vector values
%   p -- struct containing parameter values in p.MC_param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time

%% Phase variables

MC_V = x(1);
MC_mNa = x(2);
MC_hNa = x(3);
MC_mKfast = x(4);
MC_hKfast = x(5);
MC_mKa = x(6);
MC_hKa = x(7);
MC_mKslow = x(8);
MC_hKslow = x(9);

MC_INa = p.MC_gNa * MC_mNa^3 * MC_hNa * (MC_V - p.MC_ENa);
MC_INaP = p.MC_gNaP * MCNaPChan(MC_V) * (MC_V - p.MC_ENa);
MC_IKa = p.MC_gKa * MC_mKa * MC_hKa * (MC_V - p.MC_EK);
MC_IKfast = p.MC_gKfast * MC_mKfast^2 * MC_hKfast * (MC_V - p.MC_EK);
MC_IKslow = p.MC_gKslow * MC_mKslow * MC_hKslow * (MC_V - p.MC_EK);
MC_IL = p.MC_gL*(MC_V - p.MC_Eleak);

%% get ORN input for t
Input = p.intrace(floor(t / p.ORNsamplingfactor)+1);

% Recurrent inhibition for MC cell depends on time since last MC spike
tslp = t - p.last_MC_spike_time; %time since last spike
MC_GC_recurrent = (p.current_MC_recurrent_inhibition_decay_amplitude*exp(-tslp / p.MCGC_T_decay) - exp(-tslp / p.MCGC_T_rise)) .* p.MCGC_g_syn .* (MC_V - p.MCGC_V_reverse);

xdot = zeros(size(x));
xdot(1) = -(MC_INa + MC_INaP + MC_IKa + MC_IKfast + MC_IKslow + MC_IL + p.MC_Iext - Input + MC_GC_recurrent);
xdot(2:3) = MCNaChan(MC_V, MC_mNa, MC_hNa);
xdot(4:5) = MCKfastChan(MC_V, MC_mKfast, MC_hKfast);
xdot(6:7) = MCKaChan(MC_V, MC_mKa, MC_hKa);
xdot(8:9) = MCKslowChan(MC_V, MC_mKslow, MC_hKslow);

