
%% Add paths
startup

%% Set up the run parameters
tracefile = 'scaled_orn_input_depr'; % the trace file from the Data/Inputs/ORN_Input_Data directory
modelfnc = @MCRI_PGslow; % which model from Models directory, here MC-PG with slow synapse
ORNscale = 1; % the gain for ORN input


%% Run the model on all of the glomerulus inputs (this will take a very long time)
% In this example, we change the PGMCS.tc parameter from its default value
% to 170 (note the underscore replaces dot), and we turn the special
% parameter 'save_traces' on -- this will save voltage traces, etc. from
% the run. If 'save_traces' is omitted, then only the spike times will be
% saved.
doloop(tracefile, ORNscale, modelfnc, 'PGMCS_tc', '170', 'save_traces', 1)


%% Load the saved output
% The doloop function will save results with a filename following a
% particular pattern. Here we recreate it.
modelname = func2str(modelfnc);
scalingstr = regexprep(num2str(ORNscale),'\.','_');
extrastr = '_PGMCS_tc170_save_traces1'; % concatenate parameter names, values that were changed, '' if none    
loadname = sprintf('%s_%s_gain%s%s',modelname,tracefile,scalingstr,extrastr);
load(loadname);


%% Plot the MC and PG voltage traces together for the first glomerular input
plot(data(1).T, data(1).X(:,1))
