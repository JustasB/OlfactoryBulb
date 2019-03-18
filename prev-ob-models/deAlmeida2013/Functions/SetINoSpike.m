function Iext = SetINoSpike(V,PFirePre,W,E)
% Iext = SetINoSpike(V,PFirePre,W,E)
% This function calculates the external current for a given
% channel (non-spiking presynaptic neuron)
% V: membrane voltage
% PFirePre: firing probability presynaptic neuron
% W: synaptic weights
% E: Nernst potential
Iext = (W * PFirePre) .* (E - V);