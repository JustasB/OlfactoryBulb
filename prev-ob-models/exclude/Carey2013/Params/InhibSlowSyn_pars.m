function pars = InhibSlowSyn_pars(data,pars,prefix)
%%% Creates and returns struct with parameters for 
%%% slow inhibitory synapse model

if nargin < 2
    pars = struct;
end
if nargin < 3
    prefix = 'PGMCS';
end

%% Kinetics (nS?)
pars.([prefix '_svHalf']) = data.svHalf;
pars.([prefix '_sk0']) = data.sk0; 
pars.([prefix '_sTmax']) = data.sTmax;
pars.([prefix '_sk1']) = data.sk1; 
pars.([prefix '_sk2']) = data.sk2; 
pars.([prefix '_sk3']) = data.sk3; 
pars.([prefix '_sk4']) = data.sk4; 
pars.([prefix '_skd']) = data.skd; 
pars.([prefix '_gSyn']) = data.gSyn;
pars.([prefix '_sp']) = data.sp; 
pars.([prefix '_vRev']) = -70; 