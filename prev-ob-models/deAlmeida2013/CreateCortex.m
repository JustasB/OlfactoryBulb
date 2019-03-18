function [Ff,Py,Fb] = CreateCortex(Mi,Mod)
% [Ff,Py,Fb] = CreateCortex(Mi,modulation)
% This function creates the cell groups and runs the simulation:
% Mi: Mitral cell input
% Mod: paramters for a specific modulation configuration.

% Creating cell groups
Ff = classFeedforward(Mi.tsim,Mi.ncells);
Py = classPyramidal(Mi.tsim,Mi.ncells);
Fb = classFeedback(Mi.tsim,Mi.ncells);

[Ff,Py,Fb] = SetCortexParam(Ff,Py,Fb,Mod);

[Ff,Py,Fb] = RunCortex(Mi,Ff,Py,Fb);
end