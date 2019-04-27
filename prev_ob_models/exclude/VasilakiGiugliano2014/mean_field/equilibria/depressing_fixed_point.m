function [P Q ] = depressing_fixed_point(JJ, II)
%
% One population of excitatory neurons, with recurrent depressing synapses
% Threshold-linear response function.
%
% tau dh/dt = - h + J * U * x * E(h)        (eq. 1)
% dx / dt   = (1-x)/D - U * x * E(h)        (eq. 2)
%
% E(h) = alpha * (h - theta) .* (h > theta)
%
 global theta alpha U D J I;
 % Parameters initialization
 U     = 0.8;
 D     = 0.5;
 J     = JJ;
 I     = II;
 theta = 3.;
 alpha = 1.;
 
 if (I<theta) 
     Jmin  = (1 + sqrt(U * D * (theta - I)))^2 / (alpha * U);
     disp(sprintf('Critical J is %f', Jmin));
 end
 % Fixed point search as intersections between two functions:
 h     = -5:0.1:25;
 g1    = J * U * (alpha * (h - theta) .* (h > theta)) ./ (1 + U * D * (alpha * (h - theta) .* (h > theta))) + I;
 
 % Numerical evaluation of the fixed point of the system
 h0 = [];
 a = U * D * alpha;
 b = 1 - alpha * U * (D * I + J + D * theta);
 c = alpha * U * theta * (I * D + J ) - I;
 if ((b^2-4*a*c)>0)
  sol1 = (-b + sqrt(b^2-4*a*c))/(2*a);
  sol2 = (-b - sqrt(b^2-4*a*c))/(2*a);
  h0 = [sol1 , sol2];
 else
   h0 = [];
 end
 h0 = h0(find(h0 > theta)); % Consistency check, since E(h) could not be implemented
 if (I <= theta)
      h0 = [h0(:)', I];
 end
 if (isempty(h0)), x0 = [];
 else
     x0 = 1. ./ (1. + U * D * (alpha * (h0 - theta) .* (h0 > theta)));
 end
 
 disp(sprintf('The system has %d equilibrium point(s)', length((h0))));
 figure(1);  hold on;
 % Fixed point search as intersections between two functions
 P     = plot(h, g1, h, h, '--');
 set(P, 'Color', [0 0 0], 'LineWidth', 2);

 % [ho, xo] contains all the equilibrium points of the system to beplotted.
 Q = [];
 for i=1:length(h0),
      out = stability(h0(i), x0(i));    
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
 xLim([-1 21]);
 yLim([-1 21]);
 hold off;
end




function out = stability(h0, x0)
%
% LOCAL STABILITY ANALYSIS
%
% Computed the Jacobian matrix of the nonlinear system and then it
% calculates its value for each equilibrium point. The eigenvalues of
% the resulting matrix will tell whether the point is stable or not.
%
 global theta alpha U D J I;

 Jacob = zeros(2,2);
 if (h0 > theta)
  Jacob(1,1) = -1 + alpha * J * U * x0;
  Jacob(1,2) = J * alpha * U * (h0 - theta);
  Jacob(2,1) = -U * alpha * x0;
  Jacob(2,2) = -(1./D + alpha * U * (h0 - theta));
 else
  Jacob(1,1) = -1;
  Jacob(1,2) = 0;
  Jacob(2,1) = 0;
  Jacob(2,2) = - 1./D;      
 end % if
 
 EE = eig(Jacob);
 if (real(EE(1)) < 0) & (real(EE(2)) < 0)
     disp(sprintf('h0 = %f, x0 = %f: stable.', h0, x0));
     out = 1;
 elseif sign(real(EE(1)) * real(EE(2))) < 0
     disp(sprintf('h0 = %f, x0 = %f: unstable.', h0, x0));
     out = -1;
 elseif ((real(EE(1)) * real(EE(2))) == 0)
     disp(sprintf('h0 = %f, x0 = %f: ???', h0, x0));
     out = -2;
 end % if
     
end % stability