% Mi = SetBulbSimpleParam(Mi)
% This function set bulb parameters for the simple model
% Change network parameters here

function Mi = SetBulbSimpleParam(Mi)
Mi.f = 30; % frequency of modulation
Mi.r = 60; % sets the average firing rate in response to the odor
% stimulation
Mi.delta = 16; % controls the duration of the positive excursion of
% the oscillation. higher delta = higher synchronization.
Mi.A = 4; % controls the sparseness of the odor representation.
% higher A = low sparseness and vice-versa.
Mi.refrac = 20; % Period o inactivity to make the activation patter of
% Mitral cells from this model similar to the biophysical model. This
% should not be seeing as something like a refractory period of the
% integrate and fire neurons.
end