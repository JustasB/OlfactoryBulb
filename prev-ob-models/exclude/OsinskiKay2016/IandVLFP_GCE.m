function [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = IandVLFP_GCE(input_file)


% Simulate LFP activity using membrane voltage

% INPUTS
%
% numtp       -      # time points in simulation
% input file  -      name of parameter .txt file 

% OUTPUTS





%%% Must sum ALL currents into each neuron type
    [Mitral GraProximal GraDistal param InputCurrent] = OB_network_GCE(input_file);
    % Mitral LFP
    nmit  = length(Mitral);
    % Filter input currrents into Mit (ignore respiration)
    Gradistmit_filtered = zeros(nmit,length(Mitral{1}.V));
%     Glomit_filtered = zeros(nmit,length(Mitral{1}.V));
    Mitral_filtered = zeros(nmit,length(Mitral{1}.V));
    for n = 1:nmit
%         Gradistmit_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Igradistmit(n,:));
        Gradistmit_filtered(n,:) = InputCurrent.Igradistmit(n,:);
        Gradistmit_filtered(n,:) = smoothts(Gradistmit_filtered(n,:),'b',50);
%         Glomit_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Iglomit(n,:));
%         Glomit_filtered(n,:) = smoothts(Glomit_filtered(n,:),'b');
        Mitral_filtered(n,:) = Mitral{n}.V;
        Mitral_filtered(n,:) = smoothts(Mitral_filtered(n,:),'b',50);
    end
    
    % ILFP
    MitLFPs.GradistMitGlobal = sum(Gradistmit_filtered,1)/nmit;
    if mod(nmit/3,1) == 0
    MitLFPs.GradistMit1 = sum(Gradistmit_filtered(1:(nmit/3),:),1)/(nmit/3);
    MitLFPs.GradistMit2 = sum(Gradistmit_filtered((1+nmit/3):2*(nmit/3),:),1)/(nmit/3);
    MitLFPs.GradistMit3 = sum(Gradistmit_filtered((1+2*nmit/3):end,:),1)/(nmit/3);
    end
    
    % VLFP
    MitLFPs.VG = sum(Mitral_filtered,1)/nmit;
    if mod(nmit/3,1) == 0
    MitLFPs.V1 = sum(Mitral_filtered(1:(nmit/3),:),1)/(nmit/3);
    MitLFPs.V2 = sum(Mitral_filtered((1+nmit/3):2*(nmit/3),:),1)/(nmit/3);
    MitLFPs.V3 = sum(Mitral_filtered((1+2*nmit/3):end,:),1)/(nmit/3);
    end
%     MitILFPs.GloMit = sum(Glomit_filtered,1);
    
    
%     % Proximal Granule LFP
%     ngra  = length(GraProximal);
%     % Filter input currrents into proximal Gra
%     Mitgraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     Pyrgraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     Gragraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     for n = 1:ngra
%         Mitgraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Imitgraprox(n,:));
%         Mitgraprox_filtered(n,:) = smoothts(Mitgraprox_filtered(n,:),'b');
%         Pyrgraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Ipyrgraprox(n,:));
%         Pyrgraprox_filtered(n,:) = smoothts(Pyrgraprox_filtered(n,:),'b');
%         Gragraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Igragra(n,:));
%         Gragraprox_filtered(n,:) = smoothts(Gragraprox_filtered(n,:),'b');
%     end
%     GraILFPs.MitGra(:,ii) = sum(Mitgraprox_filtered,1);
%     GraILFPs.PyrGra(:,ii) = sum(Pyrgraprox_filtered,1);
%     GraILFPs.GraGra(:,ii) = sum(Gragraprox_filtered,1);
    
    % Distal Granule LFP
    ngra  = length(GraDistal);
    Gradist_filtered = zeros(ngra,length(GraDistal{1}.V));
    for n = 1:ngra
%         MitGradist_filtered(n,:) = filtfilt(Hd.Numerator,1,GraDistal{n}.V);
        Gradist_filtered(n,:) = GraDistal{n}.V;
        Gradist_filtered(n,:) = smoothts(Gradist_filtered(n,:),'b');
    end
    GraDistLFPs.VG = sum(Gradist_filtered,1)/ngra;
    


end

    