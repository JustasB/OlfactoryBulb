function xdot = vfield_PG(t,x,p)
% Vector field for the isolated PG model.
%
%   USAGE: 
%   xdot = vfield_PG(t,x,p)
%
%   INPUTS:
%   t -- current time
%   x -- (9,1) vector of current vector values
%   p -- struct containing parameter values in p.PG.param_name format
%
%   OUTPUT:
%   xdot -- derivative w.r.t. time

%% Phase variables

PG_v = x(1);
PG_N = x(2);
PG_mKa = x(3);
PG_hKa = x(4);

PG_minf = 0.1 .* (PG_v + 40) ./ (0.1 .* (PG_v + 40) + 4 .* exp(-(PG_v + 65) ./ 18) .* (1 - exp(-(PG_v + 40) ./ 10)));

PG_INa = p.PG.gNa .* PG_minf.^3 .* (0.8 - PG_N) .* (PG_v - p.PG.ENa);
PG_IK = p.PG.gK .* PG_N.^4 .* (PG_v - p.PG.EK);

PG_IKa = p.PG.gK .* PG_mKa.^3 .* PG_hKa .* (PG_v - p.PG.EK);

PG_IL = p.PG.gL .* (PG_v - p.PG.EL);

xdot = x;

%dV/dt:
xdot(1) = p.PG.Iext - PG_INa - PG_IK - PG_IKa - PG_IL;

%dn/dt:

% alphaN = 0.01 .* (PG_v+55) ./ (1 - exp(-(PG_v+55) ./ 10));
% betaN = 0.125 .* exp(-(PG_v + 65) ./ 80);
% 
% xdot(2) = alphaN .* (1 - PG_N) - betaN .* PG_N;

xdot(2) = PGKspChan(PG_v, PG_N);

% % %dPG_mKa:
% PG_mKainf = 1 ./ (1+exp((-47.1-PG_v) ./ 13.1));
% PG_mKatau = 0.04 .* exp(-PG_v/20.3) + 0.2;
% xdot(3) = (PG_mKainf - PG_mKa) ./ PG_mKatau;
% 
% %dPG_hKa:
% PG_hKainf = 1 ./ (1+exp((-(-PG_v-67.3)) ./ 6.3));
% PG_hKatau = 12; %[ish]
% xdot(4) = (PG_hKainf - PG_hKa) ./ PG_hKatau;

xdot(3:4) = PGKaChan(PG_v, PG_mKa, PG_hKa);

