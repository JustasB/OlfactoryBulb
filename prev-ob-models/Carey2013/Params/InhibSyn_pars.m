function pars = InhibSyn_pars(pars,prefix)
%%% Creates and returns struct with parameters for synapse model

if nargin < 1
    pars = struct;
end
if nargin < 2
    prefix = 'IS';
end

%% Kinetics (nS?)
pars.([prefix '_alpha']) = 0.1;
pars.([prefix '_beta']) = 0.005; 
pars.([prefix '_vHalf']) = -15;
pars.([prefix '_kAct']) = 3;
pars.([prefix '_vRev']) = -70;
pars.([prefix '_gSyn']) = 1;