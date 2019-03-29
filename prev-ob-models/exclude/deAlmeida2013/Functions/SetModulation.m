function Vlimit = SetModulation(Theta,Mod)
% Vlimit = SetModulation(Theta,Mod)
% This function calculates the limits (Vrest and FThres) for the output
% function
% The input parameters are:
% Theta: ranges for theta min or theta max (it's a 1x2 array)
% Mod: level of modulation
% In output is:
% Vlimit: this variable can be either the current Vrest or FThres for the
% function SetSpikeOutput

% Normalize Mod
I = Mod > 1;
Mod(I) = 1;

% Calculate Vlimit
Vlimit = Theta(1) + (Theta(2) - Theta(1)) * Mod;