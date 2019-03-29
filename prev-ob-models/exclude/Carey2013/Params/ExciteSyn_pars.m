function pars = ExciteSyn_pars(pars,prefix)
%%% Creates and returns struct with parameters for excitatory synapse model

if nargin < 1
    pars = struct;
end
if nargin < 2
    prefix = 'ES';
end

%% Kinetics (nS?)
pars.([prefix '_alpha']) = 2.16;
pars.([prefix '_beta']) = 0.216; 
pars.([prefix '_vHalf']) = -20;
pars.([prefix '_kAct']) = 4; 
pars.([prefix '_vRev']) = 0; 
pars.([prefix '_gSyn']) = 1;