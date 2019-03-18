%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates modulation based on Ac cell activation 
%
% Licurgo de Almeida
% 12/10/2013
%
%
% L = ModulationActivation(K,A,b)
% Inputs:
% K = 1/2 modulation;
% A = activity of the afferent network
% b = defines how steep is the slope in L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = ModulationActivation(K,A,b)

L = 1 ./ (1 + (K ./ A).^b);
