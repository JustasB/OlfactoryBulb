function [V,TFire] = SpikeOutput(V,TFire,Rest,Thres,Trefrac,Beta,dt,SpikeV,Vhyper)
% [V,TFire] = SpikeOutput(V,TFire,VLimit,Rest,Thres,Trefrac,Beta,dt,SpikeV,Vhyper)
% This function calculates the spiking chance and changes the
% membrane voltage
% TFire: time since last spike
% V: membrane voltage
% Rest: resting potential
% Thres: firing threshold
% Trefrac: inactivity period
% Beta: membrane non-linearity
% dt: timestep
% SpikeV: spike amplitude
% Vhyper: hyperpolarization potential
P_fire = (V - Rest) / (Thres -  Rest);
J = P_fire < 0;
P_fire(J) = 0;
J = P_fire > 1;
P_fire(J) = 1;
J = TFire <= Trefrac;
P_fire(J) = 0;
P_fire = P_fire.^Beta;
Out = rand(length(P_fire),1) <= P_fire * dt;
TFire(Out) = 0; % after firing, TFire returns to 0
TFire(~Out) = TFire(~Out) + dt; %neurons that didn't fire, have their TFire
% incremented by dt
V(Out) = SpikeV;
J = TFire == dt; % neurons that fired in the previous step
V(J) = Vhyper;