function [P Q r] = FullMixedModel_fixed_point(JJ, II, UUd, DDd, FFd, UUf, DDf, FFf, pd)
%
% One population of excitatory neurons, with mixed recurrent 
% depressing/facilitating synapses
% Threshold-linear response function.
%
% tau dh/dt = - h + J * u * x * E(h)             (eq. 1)
% dxD / dt   = (1-xD)/DD - uD * xD * E(h)        (eq. 2)
% duD / dt   = (UD-uD)/FD + UD * (1-uD) * E(h)   (eq. 3)
% dxF / dt   = (1-xF)/DF - uF * xF * E(h)        (eq. 4)
% duF / dt   = (UF-uF)/FF + UF * (1-uF) * E(h)   (eq. 5)
% 
% E(h) = alpha * (h - theta) .* (h > theta)
%
P = [];
Q = [];

 global theta alpha Ud Dd Fd Uf Df Ff J I ;
 
 % Parameters initialization
 Ud     = UUd;
 Dd     = DDd;
 Fd     = FFd;
 Uf     = UUf;
 Df     = DDf;
 Ff     = FFf;
 J      = JJ;
 I      = II;
 theta  = 3.;
 alpha  = 1.;
 
  % Fixed point search as intersections between two functions:
 h     = -5:0.1:100;
 ud     = Ud * (1 + Fd * (alpha .* (h-theta) .* (h>theta))) ./ (1 + Ud * Fd * (alpha .* (h-theta) .* (h>theta)));
 xd     = 1. ./ (1 + ud .* Dd .* (alpha .* (h-theta) .* (h>theta)));
 uf     = Uf * (1 + Ff * (alpha .* (h-theta) .* (h>theta))) ./ (1 + Uf * Ff * (alpha .* (h-theta) .* (h>theta)));
 xf     = 1. ./ (1 + uf .* Df .* (alpha .* (h-theta) .* (h>theta)));

 g1    = I + J * (pd .* ud .* xd   +  (1. - pd) .* uf .* xf).* (alpha .* (h-theta) .* (h>theta));
 

Ie = I;

C = [];
C(1)=-Dd*Df*Fd*Ff*Ud*Uf*alpha^4;
C(2)=Dd*Fd*Ff*J*Ud*Uf*alpha^4 - Dd*Fd*Ud*alpha^2*(Df*Uf*alpha + Ff*Uf*alpha - 2*Df*Ff*Uf*alpha^2*theta) - Df*Ff*Uf*alpha^2*(Dd*Ud*alpha + Fd*Ud*alpha - 2*Dd*Fd*Ud*alpha^2*theta) + Dd*Df*Fd*Ff*Ie*Ud*Uf*alpha^4 - Dd*Fd*Ff*J*Ud*Uf*alpha^4*pd + Df*Fd*Ff*J*Ud*Uf*alpha^4*pd;
C(3)=Df*Ff*Uf*alpha^2*(Dd*Ud*alpha*theta + Fd*Ud*alpha*theta - Dd*Fd*Ud*alpha^2*theta^2 - 1) - (Dd*Ud*alpha + Fd*Ud*alpha - 2*Dd*Fd*Ud*alpha^2*theta)*(Df*Uf*alpha + Ff*Uf*alpha - 2*Df*Ff*Uf*alpha^2*theta) + Dd*Fd*Ud*alpha^2*(Df*Uf*alpha*theta + Ff*Uf*alpha*theta - Df*Ff*Uf*alpha^2*theta^2 - 1) + Dd*Fd*J*Ud*Uf*alpha^3 + Dd*Ff*J*Ud*Uf*alpha^3 + Fd*Ff*J*Ud*Uf*alpha^3 - Dd*Fd*J*Ud*Uf*alpha^3*pd - Dd*Ff*J*Ud*Uf*alpha^3*pd + Df*Fd*J*Ud*Uf*alpha^3*pd + Df*Ff*J*Ud*Uf*alpha^3*pd + Dd*Df*Fd*Ie*Ud*Uf*alpha^3 + Dd*Df*Ff*Ie*Ud*Uf*alpha^3 + Dd*Fd*Ff*Ie*Ud*Uf*alpha^3 + Df*Fd*Ff*Ie*Ud*Uf*alpha^3 - 4*Dd*Fd*Ff*J*Ud*Uf*alpha^4*theta - 4*Dd*Df*Fd*Ff*Ie*Ud*Uf*alpha^4*theta + 4*Dd*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta - 4*Df*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta;
C(4)=(Dd*Ud*alpha + Fd*Ud*alpha - 2*Dd*Fd*Ud*alpha^2*theta)*(Df*Uf*alpha*theta + Ff*Uf*alpha*theta - Df*Ff*Uf*alpha^2*theta^2 - 1) + (Df*Uf*alpha + Ff*Uf*alpha - 2*Df*Ff*Uf*alpha^2*theta)*(Dd*Ud*alpha*theta + Fd*Ud*alpha*theta - Dd*Fd*Ud*alpha^2*theta^2 - 1) + Ff*J*Uf*alpha^2 + Dd*Fd*Ie*Ud*alpha^2 + Df*Ff*Ie*Uf*alpha^2 + Dd*J*Ud*Uf*alpha^2 + Fd*J*Ud*Uf*alpha^2 + Fd*J*Ud*alpha^2*pd - Ff*J*Uf*alpha^2*pd + Dd*Df*Ie*Ud*Uf*alpha^2 + Dd*Ff*Ie*Ud*Uf*alpha^2 + Df*Fd*Ie*Ud*Uf*alpha^2 + Fd*Ff*Ie*Ud*Uf*alpha^2 - Dd*J*Ud*Uf*alpha^2*pd + Df*J*Ud*Uf*alpha^2*pd - Fd*J*Ud*Uf*alpha^2*pd + Ff*J*Ud*Uf*alpha^2*pd - 3*Dd*Fd*J*Ud*Uf*alpha^3*theta - 3*Dd*Ff*J*Ud*Uf*alpha^3*theta - 3*Fd*Ff*J*Ud*Uf*alpha^3*theta - 3*Dd*Df*Fd*Ie*Ud*Uf*alpha^3*theta - 3*Dd*Df*Ff*Ie*Ud*Uf*alpha^3*theta - 3*Dd*Fd*Ff*Ie*Ud*Uf*alpha^3*theta - 3*Df*Fd*Ff*Ie*Ud*Uf*alpha^3*theta + 3*Dd*Fd*J*Ud*Uf*alpha^3*pd*theta + 3*Dd*Ff*J*Ud*Uf*alpha^3*pd*theta - 3*Df*Fd*J*Ud*Uf*alpha^3*pd*theta - 3*Df*Ff*J*Ud*Uf*alpha^3*pd*theta + 6*Dd*Fd*Ff*J*Ud*Uf*alpha^4*theta^2 - 6*Dd*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^2 + 6*Df*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^2 + 6*Dd*Df*Fd*Ff*Ie*Ud*Uf*alpha^4*theta^2;
C(5)=J*Uf*alpha - (Dd*Ud*alpha*theta + Fd*Ud*alpha*theta - Dd*Fd*Ud*alpha^2*theta^2 - 1)*(Df*Uf*alpha*theta + Ff*Uf*alpha*theta - Df*Ff*Uf*alpha^2*theta^2 - 1) + Dd*Ie*Ud*alpha + Df*Ie*Uf*alpha + Fd*Ie*Ud*alpha + Ff*Ie*Uf*alpha + J*Ud*alpha*pd - J*Uf*alpha*pd - 2*Ff*J*Uf*alpha^2*theta - 2*Dd*Fd*Ie*Ud*alpha^2*theta - 2*Df*Ff*Ie*Uf*alpha^2*theta - 2*Dd*J*Ud*Uf*alpha^2*theta - 2*Fd*J*Ud*Uf*alpha^2*theta - 2*Fd*J*Ud*alpha^2*pd*theta + 2*Ff*J*Uf*alpha^2*pd*theta - 2*Dd*Df*Ie*Ud*Uf*alpha^2*theta - 2*Dd*Ff*Ie*Ud*Uf*alpha^2*theta - 2*Df*Fd*Ie*Ud*Uf*alpha^2*theta - 2*Fd*Ff*Ie*Ud*Uf*alpha^2*theta + 2*Dd*J*Ud*Uf*alpha^2*pd*theta - 2*Df*J*Ud*Uf*alpha^2*pd*theta + 2*Fd*J*Ud*Uf*alpha^2*pd*theta - 2*Ff*J*Ud*Uf*alpha^2*pd*theta + 3*Dd*Fd*J*Ud*Uf*alpha^3*theta^2 + 3*Dd*Ff*J*Ud*Uf*alpha^3*theta^2 + 3*Fd*Ff*J*Ud*Uf*alpha^3*theta^2 + 3*Dd*Df*Fd*Ie*Ud*Uf*alpha^3*theta^2 + 3*Dd*Df*Ff*Ie*Ud*Uf*alpha^3*theta^2 + 3*Dd*Fd*Ff*Ie*Ud*Uf*alpha^3*theta^2 + 3*Df*Fd*Ff*Ie*Ud*Uf*alpha^3*theta^2 - 4*Dd*Fd*Ff*J*Ud*Uf*alpha^4*theta^3 - 3*Dd*Fd*J*Ud*Uf*alpha^3*pd*theta^2 - 3*Dd*Ff*J*Ud*Uf*alpha^3*pd*theta^2 + 3*Df*Fd*J*Ud*Uf*alpha^3*pd*theta^2 + 3*Df*Ff*J*Ud*Uf*alpha^3*pd*theta^2 + 4*Dd*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^3 - 4*Df*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^3 - 4*Dd*Df*Fd*Ff*Ie*Ud*Uf*alpha^4*theta^3;
C(6)=Ie - J*Uf*alpha*theta - Dd*Ie*Ud*alpha*theta - Df*Ie*Uf*alpha*theta - Fd*Ie*Ud*alpha*theta - Ff*Ie*Uf*alpha*theta - J*Ud*alpha*pd*theta + J*Uf*alpha*pd*theta + Ff*J*Uf*alpha^2*theta^2 + Dd*Fd*Ie*Ud*alpha^2*theta^2 + Df*Ff*Ie*Uf*alpha^2*theta^2 + Dd*J*Ud*Uf*alpha^2*theta^2 + Fd*J*Ud*Uf*alpha^2*theta^2 + Fd*J*Ud*alpha^2*pd*theta^2 - Ff*J*Uf*alpha^2*pd*theta^2 + Dd*Df*Ie*Ud*Uf*alpha^2*theta^2 + Dd*Ff*Ie*Ud*Uf*alpha^2*theta^2 + Df*Fd*Ie*Ud*Uf*alpha^2*theta^2 - Dd*Fd*J*Ud*Uf*alpha^3*theta^3 - Dd*Ff*J*Ud*Uf*alpha^3*theta^3 + Fd*Ff*Ie*Ud*Uf*alpha^2*theta^2 - Fd*Ff*J*Ud*Uf*alpha^3*theta^3 - Dd*J*Ud*Uf*alpha^2*pd*theta^2 + Df*J*Ud*Uf*alpha^2*pd*theta^2 - Fd*J*Ud*Uf*alpha^2*pd*theta^2 + Ff*J*Ud*Uf*alpha^2*pd*theta^2 - Dd*Df*Fd*Ie*Ud*Uf*alpha^3*theta^3 - Dd*Df*Ff*Ie*Ud*Uf*alpha^3*theta^3 - Dd*Fd*Ff*Ie*Ud*Uf*alpha^3*theta^3 - Df*Fd*Ff*Ie*Ud*Uf*alpha^3*theta^3 + Dd*Fd*Ff*J*Ud*Uf*alpha^4*theta^4 + Dd*Fd*J*Ud*Uf*alpha^3*pd*theta^3 + Dd*Ff*J*Ud*Uf*alpha^3*pd*theta^3 - Df*Fd*J*Ud*Uf*alpha^3*pd*theta^3 - Df*Ff*J*Ud*Uf*alpha^3*pd*theta^3 - Dd*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^4 + Df*Fd*Ff*J*Ud*Uf*alpha^4*pd*theta^4 + Dd*Df*Fd*Ff*Ie*Ud*Uf*alpha^4*theta^4;

r =[];
r = roots(C);
r = r(find( (abs(imag(r))<1e-8) & (r>theta)) )
if (isempty(r)), r = nan; end; 
end











