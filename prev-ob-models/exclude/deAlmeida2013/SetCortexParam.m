function [Ff,Py,Fb] = SetCortexParam(Ff,Py,Fb,Mod)
% [Ff,Py,Fb] = SetBulbParam(Ff,Py,Fb,Mod)
% This function set the cortical parameters, where:
% Ff,Py,Fb are objects;
% Mod is either between 1 (for Mod fully ON), 0 (for Mod fully OFF) or -1 
% (for Mod dynamic)


% Feedforward neurons (Ff)
Ff.AMPAFf.G = 200e-3;
Ff.tau = 5;
Ff.Tmin = [-0.3e-3,-0.3e-3];
Ff.Cmax = 5;
Ff.ModValue = Mod;


%Pyramidal cells
Py.AMPAFf.G = 759e-3;
Py.AMPAFb.G = [510,260];
Py.AMPAFb.tau1 = 1;
Py.AMPAFb.tau2 = 2;
Py.GABAFb.G = [550e-3,550e-3];
Py.GABAFf.G = 55e-3;
Py.AAHP = [40,0];
Py.Cmax = 5;
Py.ModValue = Mod;
Py.tauAHP = 100;
Py.tau = 10;
Py.Beta = 2; 
Py.Tau11 = 800;
Py.Tau01 = 400;
Py.Tau10 = 400;
Py.LearnUnlearn = false;
Py.MAMPAFf = Py.SetConnections(Py.ncells,Py.ConnAMPAFf,'normal');
Py.WAMPAFf = Py.MAMPAFf;


% Feedback cells (Fb)
Fb.AMPAFf.G = [250e-3,60e-3];
Fb.tau = 5;
Fb.Tmin = [0e-3,-0.1e-3];
Fb.Cmax = 5;
Fb.ModValue = Mod;