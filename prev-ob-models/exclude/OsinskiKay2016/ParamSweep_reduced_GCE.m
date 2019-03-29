function [MCIMAT, MCVMAT, param, MCilfpMAT, MCvlfpMAT] = ...
    ParamSweep_reduced_GCE(P1,P2,numtp,input_file)

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
% MCIMAT      -   length(P1.val) x length(P2.val) cell array containing only MC currents
% MCVMAT      -   length(P1.val) x length(P2.val) cell array containing only MC V
% MCilfpMAT   -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from MC IPSCs 
% dGCilfpMAT  -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from GC EPSCs
% MCvlfpMAT   -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from MC Vs 
% dGCvlfpMAT  -   length(P1.val) x length(P2.val) cell array containing LFP simulated
%                  from GC Vs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes on running the file:
% - Parameter structures must be created as described in INPUTS
% - InitNetwork_GCE.m must be run before running this file
% - This function runs ILFP_GC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCIMAT = cell(length(P1.val),length(P2.val)); % MC currents
MCVMAT = cell(length(P1.val),length(P2.val)); % MC voltage

MCilfpMAT = cell(length(P1.val),length(P2.val));
% dGCilfpMAT = cell(length(P1.val),length(P2.val));
MCvlfpMAT = cell(length(P1.val),length(P2.val));
% dGCvlfpMAT = cell(length(P1.val),length(P2.val));


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
        [Mitral, ~, ~, param, InputCurrent, MitLFPs, ~] ...
            = IandVLFP_GCE(input_file);
        
        IC.Iext = InputCurrent.Iext;
        IC.Igradistmit = InputCurrent.Igradistmit;
        MCIMAT{n1,n2} = IC;
        
        MCV = zeros(param.nMitral,numtp);
        for n = 1:param.nMitral
            MCV(n,:) = Mitral{n}.V;
        end
        MCVMAT{n1,n2} = MCV;
        
        MCilfpMAT{n1,n2} = MitLFPs.GradistMitGlobal;
        MCvlfpMAT{n1,n2} = MitLFPs.VG;
        
        disp([num2str(length(P1.val)*length(P2.val) - ((n1-1)*length(P2.val) + n2)),' more iterations'])
        
    end
end
