function [OSN,Pg,M1,M2,Gr] = SetBulbParam(OSN,Pg,M1,M2,Gr,Mod)
% [OSN,Pg,M1,M2,Gr] = SetBulbParam(OSN,Pg,M1,M2,Gr,AChMod)
% This function set bulb parameters for the noradrenaline project, where:
% OSN,Pg,M1,M2,Gr are objects;
% Mod is either between 1 (for Mod fully ON), 0 (for Mod fully OFF) or -1 
% (for Mod dynamic)

% OSN
OSN.AMPAFf.G = 200e-3;
OSN.ModValue = Mod;
OSN.Tmax = [1e-3,1e-3];
OSN.Tmin = [0e-3,0e-3];
OSN.Cmax = 5;

%Pg cells
Pg.AMPAFf.G = 166e-3;
Pg.Tmax = [4e-3,4e-3];
Pg.Tmin = [0e-3,0e-3];
Pg.ModValue = Mod;
Pg.Cmax = 5;

% M1 (Mitral apical)
M1.AMPAFf.G = 160e-3;
M1.GABAFf.G = 380e-3;
M1.tau = 4;
M1.Rsom = 1;
M1.Tmax = [8e-3,8e-3];
M1.Cmax = 5;

% M2 (Mitral soma)
M2.Tmax = [9e-3,1e-3]; % the order here is [Mod OFF, Mod ON]
M2.Tmin = [-1.4e-3,-1.4e-3]; % the order here is [Mod OFF, Mod ON]
M2.ModValue = Mod;
M2.Beta = 2;
M2.K = 2;
M2.Cmax = 5; %Max neuromodulator concentration
M2.GABAFb.G = 180e-3;

% Gr cells
Gr.AMPAFf.G = 20e-3;
Gr.Tmax = [6e-3,6e-3];
Gr.Tmin = [-1e-3,-2.4e-3];
Gr.b = 1.5;
Gr.K = 3;
Gr.Cmax = 5;
Gr.ModValue = Mod;
Gr.Beta = 3;
Gr.Tminne = [0e-3,1e-3];
Gr.bne = 1;
Gr.Kne = 1;
Gr.Cmaxne = 5;
