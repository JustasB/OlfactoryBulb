%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the ligand affinity for a given receptor (Hill
% eq.)
%
% Licurgo de Almeida
% 04/26/2012
%
%
% L = CalculateAffinity(K,A,r,b)
% Inputs:
% K = the odor ligand-receptor affinity;
% A = activity of the afferent network
% r = concentration range.
% b = defines how steep is the slope in L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = CalculateAffinity(K,A,r,b)
if nargin < 4
    b = 1;
end

L = 1 ./ (1 + (10^K ./ 10.^(A * r)).^b);
