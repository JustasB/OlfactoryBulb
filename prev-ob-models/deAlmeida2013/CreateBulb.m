function [OSN,Pg,M1,M2,Gr] = CreateBulb(tsim,ncells,Mod)
% [OSN,Pg,M1,M2,Gr] = CreateBulb(tsim,ncells)
% This function creates the cell groups and runs the simulation:
% tsim: simulation time
% ncells: number of cells in each cell group
% Mod is either true (for Mod ON) or false (for Mod OFF)

% Creating cell groups
OSN = classOSN(tsim,ncells);
Pg = classPglo(tsim,ncells);
M1 = classMitglo(tsim,ncells);
M2 = classMitsoma(tsim,ncells);
Gr = classGranule(tsim,ncells);

[OSN,Pg,M1,M2,Gr] = SetBulbParam(OSN,Pg,M1,M2,Gr,Mod);

[OSN,Pg,M1,M2,Gr] = RunBulb(OSN,Pg,M1,M2,Gr);


end