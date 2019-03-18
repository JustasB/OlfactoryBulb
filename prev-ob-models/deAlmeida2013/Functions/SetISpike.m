function Iext = SetISpike(V,tau1,tau2,TFirePre,W,E)
% Iext = SetISpike(V,tau1,tau2,TFirePre,W,E)
% This function calculates the external current for a given
% channel (spiking presynaptic neuron)
% V: membrane voltage
% tau1: channel's rising time
% tau2: channel's falling time
% TFirePre: time since last presynaptic spike
% W: synaptic weights
% E: Nernst potential
g = -(exp(-(TFirePre) / tau1) - exp(-(TFirePre) / tau2));
Iext = (W * g) .* (E - V);