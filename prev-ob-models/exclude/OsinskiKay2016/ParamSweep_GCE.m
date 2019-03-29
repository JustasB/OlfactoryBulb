function [MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file)

% This function runs the OB_network_GCE.m model for every combination of
% parameter values defined in P1(2). Outputs data in cell arrays with
% dimension length(P1.val) x length(P2.val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%
% P1(2)       -  structure containing the following parameter fields:
%                 P1(2).line = row # of parameter in OB_params_GCE.txt file
%                 P1(2).name = parameter name as it appears in OB_params_GCE.txt file
%                 P1(2).val = array containing numerical values
% numtp       -  # time points in simulation
% input_file  -  name of parameter file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%
% MCMAT      -   length(P1.val) x length(P2.val) cell array containing MC data 
% dGCMAT     -   length(P1.val) x length(P2.val) cell array containing GC data 
% ICMAT      -   length(P1.val) x length(P2.val) cell array containing all currents
% MClfpMAT   -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from MC IPSCs 
% dGClfpMAT  -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from GC EPSCs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on running the file:
% - Parameter structures must be created as described in INPUTS
% - InitNetwork_GCE.m must be run before running this file
% - This function runs ILFP_GC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCMAT = cell(length(P1.val),length(P2.val));
dGCMAT = cell(length(P1.val),length(P2.val));
ICMAT = cell(length(P1.val),length(P2.val));

MClfpMAT = cell(length(P1.val),length(P2.val));
dGClfpMAT = cell(length(P1.val),length(P2.val));

numtrials = 1; % only simulate one trial (no averaging)

TextCell = regexp( fileread(input_file), '\n', 'split');

for n1 = 1:length(P1.val)
    
    if isnumeric(P1.val)
        wstring = [P1.name,num2str(P1.val(n1))];
    else
        wstring = [P1.name,P1.val{n1}];
    end
    
    TextCell{P1.line} = sprintf('%s',wstring);
    fid = fopen(input_file, 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);
    for n2 = 1:length(P2.val)
        
        if isnumeric(P2.val)
            wstring = [P2.name,num2str(P2.val(n2))];
        else
            wstring = [P2.name,P2.val{n2}];
        end
        
        TextCell{P2.line} = sprintf('%s',wstring);
        fid = fopen(input_file, 'w');
        fprintf(fid, '%s\n', TextCell{:});
        fclose(fid);
    
        % Run model
        [Mitral GraProximal GraDistal param InputCurrent MitILFPs GraProxILFPs GraDistILFPs] ...
            = ILFP_GCE(numtp, numtrials, input_file);

        MCMAT{n1,n2} = Mitral;
        dGCMAT{n1,n2} = GraDistal;
        ICMAT{n1,n2} = InputCurrent;
        
        MClfpMAT{n1,n2} = MitILFPs;
        dGClfpMAT{n1,n2} = GraDistILFPs;
        
        disp([num2str(length(P1.val)*length(P2.val) - ((n1-1)*length(P2.val) + n2)),' more iterations'])
        
    end
end
