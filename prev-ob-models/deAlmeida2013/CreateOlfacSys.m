function [OSN,Pg,M1,M2,Gr,Ff,Py,Fb] = CreateOlfacSys(tsim,ncells,Mod,OdorConc,OdorTime)
% [OSN,Pg,M1,M2,Gr,Ff,Py,Fb] = CreateOlfacSysNE(tsim,ncells,Mod,OdorConc,OdorTime)
% This function creates the cell groups and runs the simulation:
% tsim: simulation time
% ncells: number of cells in each cell group
% Mod is a value between 1 (for Mod ON) or 0 (for Mod OFF)
% OdorConc: Odor concentration (between 0 and 1)
% OdorTime: time of odor presentation [tinit,tfinal]

if nargin < 5
    OdorTime = [0, tsim];
end
if nargin < 4
    OdorConc = 0;
end
% Creating cell groups:
% Bulb:

OSN = classOSN(tsim,ncells);
Pg = classPglo(tsim,ncells);
M1 = classMitglo(tsim,ncells);
M2 = classMitsoma(tsim,ncells);
Gr = classGranule(tsim,ncells);

% Cortex:
Ff = classFeedforward(tsim,ncells);
Py = classPyramidal(tsim,ncells);
Fb = classFeedback(tsim,ncells);

% Set network parameters
[OSN,Pg,M1,M2,Gr] = SetBulbParam(OSN,Pg,M1,M2,Gr,Mod);
OSN.InputTimes = OdorTime; % sets the new odor time
OSN.OdorInput = SetOdorInput(OSN);
OSN.OdorInput = OSN.OdorInput .* OdorConc; % changes odor concentration
[Ff,Py,Fb] = SetCortexParam(Ff,Py,Fb,Mod);