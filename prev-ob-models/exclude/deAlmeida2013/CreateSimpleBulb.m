function M2 = CreateSimpleBulb(tsim,ncells)
% M2 = CreateSimpleBulb(tsim,ncells)
% This function creates the Mi cells and runs the simulation:
% tsim: simulation time
% ncells: number of cells in each cell group
% Mod is either true (for Mod ON) or false (for Mod OFF)

M2 = classSimpleMitral(tsim,ncells);
M2 = SetBulbSimpleParam(M2);
M2 = CreateOutput(M2);
