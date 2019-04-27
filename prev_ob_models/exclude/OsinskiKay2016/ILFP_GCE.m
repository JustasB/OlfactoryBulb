function [Mitral GraProximal GraDistal param InputCurrent MitILFPs GraProxILFPs GraDistILFPs] ...
    = ILFP_GCE(numtp, numtrials, input_file)


% Simulate LFP activity using membrane voltage
%
% Boleslaw Osinski (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%
% numtp       -      # time points in simulation
% numtrials   -      # times to run simulation
% input_file  -      name of parameter file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%
% Mtral         -  Structure containing all distal GC dendritic data (except for currents)
% GraProximal   -  Structure containing all proximal GC dendritic data (except for currents, currently unused)
% Gradistal     -  Structure containing all distal GC dendritic data (except for currents)
% param         -  Structure contaioning all network parameters only (not cell parameters)
% InputCurrent  -  Structure containing all currents, including Ca currents
% MitILFPs      -  LFP simulated from MC currents
% GraProxILFPs  -  LFP simulated from proximal GC dendritic currents (currently unused)
% GraDistILFPs  -  LFP simulated from distal GC dendritic currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MitILFPs.GradistMitGlobal = zeros(numtp,numtrials);
MitILFPs.GradistMit1 = zeros(numtp,numtrials);
MitILFPs.GradistMit2 = zeros(numtp,numtrials);
MitILFPs.GradistMit3 = zeros(numtp,numtrials);

MitILFPs.GloMit = zeros(numtp,numtrials);
GraDistILFPs.MitGradist = zeros(numtp,numtrials);
GraProxILFPs.PyrGra = zeros(numtp,numtrials);
GraProxILFPs.GraGra = zeros(numtp,numtrials);

%%% Must sum ALL currents into each neuron type
for ii = 1:numtrials
    [Mitral GraProximal GraDistal param InputCurrent] = OB_network_GCE(input_file);
    % Mitral LFP
    nmit  = length(Mitral);
    % Filter input currrents into MC (ignore respiration)
    Gradistmit_filtered = zeros(nmit,length(Mitral{1}.V));
    for n = 1:nmit
        Gradistmit_filtered(n,:) = InputCurrent.Igradistmit(n,:);
        Gradistmit_filtered(n,:) = smoothts(Gradistmit_filtered(n,:),'b',50);
    end
    MitILFPs.GradistMitGlobal(:,ii) = sum(Gradistmit_filtered,1)/nmit;
    if mod(nmit/3,1) == 0
    MitILFPs.GradistMit1(:,ii) = sum(Gradistmit_filtered(1:(nmit/3),:),1)/(nmit/3);
    MitILFPs.GradistMit2(:,ii) = sum(Gradistmit_filtered((1+nmit/3):2*(nmit/3),:),1)/(nmit/3);
    MitILFPs.GradistMit3(:,ii) = sum(Gradistmit_filtered((1+2*nmit/3):end,:),1)/(nmit/3);
    end
    
    % Distal Granule LFP
    ngra  = length(GraDistal);
    MitGradist_filtered = zeros(ngra,length(GraDistal{1}.V));
    for n = 1:ngra
        MitGradist_filtered(n,:) = GraDistal{n}.V;
        MitGradist_filtered(n,:) = smoothts(MitGradist_filtered(n,:),'b');
    end
    GraDistILFPs.MitGradist(:,ii) = sum(MitGradist_filtered,1)/ngra;
    


end

end
    