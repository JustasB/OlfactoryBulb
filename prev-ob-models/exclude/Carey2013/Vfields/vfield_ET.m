function xdot = vfield_ET(t,x,p)
% Vector field for the minimal ET model.
%
%   INPUTS:
%   t -- current time
%   x -- (6,1) vector of current vector values
%   p -- struct containing parameter values in p.ET.param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time

%% Phase variables

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

%% get ORN input for t
Input = p.intrace(floor(t / p.ORNsamplingfactor)+1);

xdot = zeros(size(x));
xdot(1) = -(ET_INa + ET_IK + ET_ICaT + ET_IH + ET_INaP + ET_IL + p.ET_Iext - Input)./p.ET_C;
xdot(2) = calc_xdot(ET_V,ET_nK,p.ET_theta_nK, p.ET_sigma_nK,p.ET_tau_nK);
xdot(3) = calc_xdot(ET_V, ET_hNaP, p.ET_theta_hNaP, p.ET_sigma_hNaP, p.ET_tau_hNaP);

ET_mH_pars = [p.ET_theta_mH, p.ET_sigma_mH, p.ET_theta_mH_T, p.ET_sigma_mH_T, ...
    p.ET_tau_mH_T, p.ET_delta_mH_T];
xdot(4) = calc_H_mdot(ET_V, ET_mH, ET_mH_pars);

xdot(5) = calc_xdot(ET_V,ET_mCaT,p.ET_theta_mCaT, p.ET_sigma_mCaT,p.ET_tau_mCaT);
xdot(6) = calc_xdot(ET_V,ET_hCaT,p.ET_theta_hCaT, p.ET_sigma_hCaT,p.ET_tau_hCaT);
