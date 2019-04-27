function [P Q ] = FullModel_fixed_point(JJ, II, UU, DD, FF)
%
% One population of excitatory neurons, with recurrent depressing/facilitating
% synapses
% Threshold-linear response function.
%
% tau dh/dt = - h + J * u * x * E(h)        (eq. 1)
% dx / dt   = (1-x)/D - u * x * E(h)        (eq. 2)
% du / dt   = (U-u)/F + U * (1-u) * E(h)    (eq. 3)
% 
% E(h) = alpha * (h - theta) .* (h > theta)
%
 global theta alpha U D F J I h0;
 % Parameters initialization
 U     = UU;
 D     = DD;
 F     = FF;
 J     = JJ;
 I     = II;
 theta = 3.;
 alpha = 1.;
 
  % Fixed point search as intersections between two functions:
 h     = -5:0.1:100;
 u     = U * (1 + F * (alpha .* (h-theta) .* (h>theta))) ./ (1 + U * F * (alpha .* (h-theta) .* (h>theta)));
 x     = 1. ./ (1 + u .* D .* (alpha .* (h-theta) .* (h>theta)));
 g1    = I + J * u .* x .* (alpha .* (h-theta) .* (h>theta));
 
 figure(1);  hold on;
 % Fixed point search as intersections between two functions
 P     = plot(h, g1, h, h, '--');
 set(P, 'Color', [0 0 0], 'LineWidth', 2);
% 
%  x = [];
% for hh=1:max(h)
%     x = [x, fzero(@myfun,hh)];
% end
% x = unique(x)
 
 % Numerical evaluation of the fixed point of the system
 h0 = [];
 a = -alpha^2 * D * F * U;
 b =  alpha * (F * (-1 + alpha * J) + D * (-1 + alpha * F * (I + 2 * theta))) * U;
 c = -1. + alpha^2 * F * theta * (-2. * D * I - 2. * J - 1. * D * theta) * U + alpha * (J + D * (I + theta) + F * (I + theta)) * U;
 d = -alpha * J * theta * (1 - alpha * F * theta) * U + I * (1 - alpha * (D + F) * theta * U + alpha^2 * D * F * theta^2 * U);
 
 sol  = [];
 sol(1) = -(b/(3 * a)) - (2^(1/3) * (-b^2 + 3 * a * c))/(3 * a * (-2 * b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3)) + (-2 * b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3)/(3 * 2^(1/3) * a);
 sol(2) = -(b/(3 * a)) + ((1 + sqrt(-1) * sqrt(3)) * (-b^2 + 3 * a * c))/(3*  2^(2/3)*  a*  (-2*  b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3)) - ((1 - sqrt(-1) * sqrt(3)) * (-2 * b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3))/(6 * 2^(1/3) * a); 
 sol(3) = -(b/(3 * a)) + ((1 - sqrt(-1) * sqrt(3)) * (-b^2 + 3 * a * c))/(3*  2^(2/3)*  a*  (-2*  b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3)) - ((1 + sqrt(-1) * sqrt(3)) * (-2 * b^3 + 9 * a * b * c - 27 * a^2 * d + sqrt(-4 * (b^2 - 3 * a * c)^3 + (2 * b^3 - 9 * a * b * c + 27 * a^2 * d)^2))^(1/3))/(6 * 2^(1/3) * a); 

 sol'
 h0 = [];
 for i=1:3,
    if (abs(imag(sol(i))) < 1e-13)
     h0 = [h0, real(sol(i))];
    end
  end 
  
   h0 = h0(find(h0 > theta)); % Consistency check, since E(h) could not be implemented
   h0
 if (I <= theta)
      h0 = [h0(:)', I];
 end
 if (isempty(h0)), x0 = [];
 else
     u0 = U * (1 + F * (alpha .* (h0-theta) .* (h0>theta))) ./ (1 + U * F * (alpha .* (h0-theta) .* (h0>theta)));
     x0 = 1. ./ (1 + u0 .* D .* (alpha .* (h0-theta) .* (h0>theta)));
 end
 
 disp(sprintf('The system has %d equilibrium point(s)', length((h0))));
 figure(1);  hold on;
 % Fixed point search as intersections between two functions
 P     = plot(h, h, '--',h, g1);
 set(P, 'Color', [0 0 0], 'LineWidth', 2);

 % [ho, xo] contains all the equilibrium points of the system to beplotted.
 Q = [];
 for i=1:length(h0),
      out = stability(h0(i), x0(i), u0(i));    
      if (out==1), 
          Q(i) = plot(h0(i), h0(i), 'o');
          set(Q(i), 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 15);
      elseif (out ==-1)
          Q(i) = plot(h0(i), h0(i), 's');
          set(Q(i), 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 15);
      end
 end
 
 xlabel('h', 'FontName', 'Arial', 'FontSize', 15);
 set(gca, 'Box', 'on', 'FontName', 'Arial', 'FontSize', 15);
 xlim([-1 100]);
 ylim([-1 100]);
 hold off;
  axis square
end


function out = myfun(h) 
 global theta alpha U D F J I;
 
 %out = alpha * J (1+alpha * F * (h-theta)) * (h-theta) * U+(I-h) * (1+alpha * (D+F) * (h-theta) * U+alpha^2 * D * F * (h-theta)^2 * U);
 
 
 u     = U * (1 + F * (alpha .* (h-theta) .* (h>theta))) ./ (1 + U * F * (alpha .* (h-theta) .* (h>theta)));
 x     = 1. ./ (1 + u .* D .* (alpha .* (h-theta) .* (h>theta)));
 g1    = I + J * u .* x .* (alpha .* (h-theta) .* (h>theta));
 out = g1 - h;
end




function out = stability(h0, x0, u0)
%
% LOCAL STABILITY ANALYSIS
%
% Computed the Jacobian matrix of the nonlinear system and then it
% calculates its value for each equilibrium point. The eigenvalues of
% the resulting matrix will tell whether the point is stable or not.
%
 global theta alpha U D F J I;

 Jacob = zeros(3,3);
 if (h0 > theta)
  Jacob(1,1) =-1+J*u0*x0*alpha;
  Jacob(1,2) =J*u0*alpha*(h0-theta);
  Jacob(1,3) =J*x0*alpha*(h0-theta);
  Jacob(2,1) =-u0*x0*alpha;
  Jacob(2,2) =-1/D-u0*alpha*(h0-theta);
  Jacob(2,3) =-x0*alpha*(h0-theta);
  Jacob(3,1) =alpha*U*(1-u0);
  Jacob(3,2) =0;
  Jacob(3,3) =-1/F-U*alpha*(h0-theta);
 else
  Jacob(1,1) = -1;
  Jacob(1,2) = 0;
  Jacob(1,3) = 0;
  Jacob(2,1) = 0;
  Jacob(2,2) = -1/D;
  Jacob(2,3) = 0;
  Jacob(3,1) = 0;
  Jacob(3,2) = 0;
  Jacob(3,3) =-1/F;
 end
 
 EE = eig(Jacob);
 if (real(EE(1)) < 0) & (real(EE(2)) < 0) & (real(EE(3)) < 0)
     disp(sprintf('h0 = %f, x0 = %f, u0 = %f: stable.', h0, x0, u0));
     out = 1;
 elseif ((real(EE(1)) * real(EE(2)) * real(EE(3))) ~= 0) & abs(sign(real(EE(1))) + sign(real(EE(2))) + sign(real(EE(3))))
     disp(sprintf('h0 = %f, x0 = %f: unstable.', h0, x0));
     out = -1;
 elseif ((real(EE(1)) * real(EE(2)) * real(EE(3))) == 0)
     disp(sprintf('h0 = %f, x0 = %f, u0 = %f: ???.', h0, x0, u0));
     out = -2;
 end % if
     
end % stability