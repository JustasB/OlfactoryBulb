function C = generate_connectivity(N, P1, P2)
%
% This function returns a binary connectivity matrix [C], whose generic
% elements C(i,j) are non-zero according to the following scheme
%
% the probability of having both connections C(i,j) and C(j,i) is P2;
% the probability of having only one connections C(i,j) OR C(j,i) is P1;
% the probability of having no connection is 1-(P1+P2);
%
%
% If one wants to generate a usual random connectivity with p = 10% prob.
% and independent probability of (i.e. 10% * 10% = p^2 = 1%) joint 
% connectivity, then please choose 
%
% P1 = 2*p*(1-p)
% P2 = p*p
%


% NOTE:
% 3 events are possible and must be generated...
%
% P_bidirectional  = P2
% P_unidirectional = P1
% P_noconnection   = Po
%
% Let's then take these values, Po, P1, and P2 and place them on the x-axis
% so that they cover entirely the range [0 ; 1] of a uniform R.V.
%
 if ((P1 + P2) > 1) disp('Error! P1 + P2 > 1!'); C = []; return; end; 
 Po = 1 - (P1 + P2); 
 r1 = Po;
 r2 = r1 + P1 * 0.5;
 r3 = r2 + P1 * 0.5;
 r4 = r3 + P2;
 
 for i=1:N,
    for j=(i+1):N,
        r = rand;
        if (r < r1)
         C(i,j) = 0;      C(j,i) = 0;
        elseif (r < r2)
         C(i,j) = 1;      C(j,i) = 0;
        elseif (r < r3)
         C(i,j) = 0;      C(j,i) = 1;
        else % r > r3 and r < r4
         C(i,j) = 1;      C(j,i) = 1;
        end
    end
 end
         
    
end