function PFire = NoSpikeOutput(V,Rest,Thres,Beta)
% PFire = NoSpikeOutput(V,Rest,Thres,Beta)
% This function calculates the output changes in non-spiking
% neurons
% PFire: output probability
% V: membrane voltage
% Rest: resting potential
% Thres: firing threshold
% Beta: membrane non-linearity
PFire = (V - Rest) / (Thres -  Rest);
J = PFire < 0;
PFire(J) = 0;
J = PFire > 1;
PFire(J) = 1;
PFire = PFire.^Beta;