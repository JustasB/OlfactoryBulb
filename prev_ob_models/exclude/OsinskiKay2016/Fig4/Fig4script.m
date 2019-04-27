% Figure 4 shows the currents, rasters, and LFPs for 3 different values of
% Vrest GC

% params for this fig have been changed after realizing that MCs were not
% being summed in the calculation of IVDCC

% FINAL PARAMS:
% wGABAGR 0.0125
% CCaTh 1.5
% wVDCC 250

clearvars
input_file = 'OB_params_GCE_Fig4.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);


% initialize Fig 4A parameters
% cmat = [0,0,1/ngradist; 0.5/ngradist,0,0.5/ngradist; 1/ngradist,0,0]; % color
cmat = ones(3)*1/ngradist; % bw
fs = 16;% fontsize

P.line = 57;
P.name = 'Vrest ';
VrestGC = -1e-3*[74 68 60];

% initialize Fig 4B params (FFT)
trim = 1000; % trim beginning and end to avoid edge effects
% trim:end-50
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);

faxis = [50 90;20 60;10 40]; % frequency axis for power plots

TextCell = regexp( fileread(input_file), '\n', 'split');

for vv = 1:length(VrestGC)
    
    % write VrestGC to parameter file
    wstring = [P.name,num2str(VrestGC(vv))];
    TextCell{P.line} = sprintf('%s',wstring);
    fid = fopen(input_file, 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);
    
    % simulate model
    [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = ...
        IandVLFP_GCE(input_file);
    
    %%%%%%%%%%%%%%%%%%% Fig 4A %%%%%%%%%%%%%%%%%%%%%
    figure
    set(gcf,'position',[0,400,512,800]);
    xsc = [0 500];% x scale
    
    % plot currents
    ysc = [0 0.013]; % y scale
    % plot ImitgradistAMPA
    subplot(6,1,1)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.ImitgradistAMPA(ii,:)/2,'color',ii*cmat(vv,:)); % color
         plot(timevec,InputCurrent.ImitgradistAMPA(ii,:)/2,'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    subplot(6,1,2)
    hold on
    for ii = 1:ngradist
       %plot(timevec,InputCurrent.ImitgradistNMDA(ii,:),'color',ii*cmat(vv,:)); % color
       plot(timevec,InputCurrent.ImitgradistNMDA(ii,:),'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    % plot ImitgradistVDCC
    subplot(6,1,3)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.ImitgradistVDCC(ii,:),'color',ii*cmat(vv,:)); % color
        plot(timevec,InputCurrent.ImitgradistVDCC(ii,:),'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    % Prelease
    subplot(6,1,4)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.Prelease(ii,:),'color',ii*cmat(vv,:)); % color
        plot(timevec,InputCurrent.Prelease(ii,:),'color',1-ii*cmat(vv,:))
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim([0 1])
    
    % Raster plot
    SpikeV = 65e-3;
        SPIKES = zeros(nmit,ntp);
        for n = 1:nmit
            SPIKES(n,:) = Mitral{n}.S;
        end
    SPIKES = flipud(SPIKES); % order MCs in direction low E (bottom) to high E (top)
    subplot(6,1,5)
%     RasterPlot(SPIKES,dt,ntp*dt,ngradist*cmat(vv,:),fs,0); % color
    RasterPlot(SPIKES,dt,ntp*dt,'k',fs,0)
    xlim([xsc])

    % LFP plots
    subplot(6,1,6)
    hold on
%     plot(timevec,detrend(MitLFPs.GradistMitGlobal),'color',ngradist*cmat(vv,:)); % color
%     plot(timevec,detrend(MitLFPs.VG),'color',[ngradist*cmat(vv,:),0.2],'linewidth',2);
    plot(timevec,detrend(MitLFPs.GradistMitGlobal),'color','k','linewidth',2); % bw
    plot(timevec,detrend(MitLFPs.VG),'color',[0.6,0.6,0.6],'linewidth',2);
    hold off
    set(gca,'fontsize',fs)
    xlim([xsc]);ylim([-.005 0.005])
    xlabel('time (ms)')
%     legend('ILFP','VLFP');legend boxoff
tightfig
    
    
    %%%%%%%%%%%%%%%%%%% Fig 4B %%%%%%%%%%%%%%%%%%%%%
    % Power spectra
    mitFFTGI = fft(detrend(MitLFPs.GradistMitGlobal(trim:end-100),'constant'),NFFT)/L;
    mitFFT1I = fft(detrend(MitLFPs.GradistMit1(trim:end-100),'constant'),NFFT)/L;
    mitFFT2I = fft(detrend(MitLFPs.GradistMit2(trim:end-100),'constant'),NFFT)/L;
    mitFFT3I = fft(detrend(MitLFPs.GradistMit3(trim:end-100),'constant'),NFFT)/L;

    mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;
    mitFFT1V = fft(detrend(MitLFPs.V1(trim:end-100),'constant'),NFFT)/L;
    mitFFT2V = fft(detrend(MitLFPs.V2(trim:end-100),'constant'),NFFT)/L;
    mitFFT3V = fft(detrend(MitLFPs.V3(trim:end-100),'constant'),NFFT)/L;

    scrsz = get(0,'ScreenSize');
    figH=figure;
    set(figH,'position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
    hold on
%     plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color',ngradist*cmat(vv,:)); % color
%     plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[ngradist*cmat(vv,:),0.2],'linewidth',2)
%     plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle','--')
%     plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle',':','linewidth',2);
    % plot(f,2*abs(mitFFT1V(1:NFFT/2+1)),f,2*abs(mitFFT3V(1:NFFT/2+1)));
    plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color','k'); % bw
    plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[0.75,0.75,0.75],'linewidth',2)
    plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color','k','linestyle','-.')
    plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color','k','linestyle',':','linewidth',2);
    
    plot(f,1e-4*ones(1,length(f)),'k--')
    
    hold off
    xlim(faxis(vv,:));ylim([0 2e-3])
    set(gca,'fontsize',fs)
    legend('ILFP - Global','VLFP - Global','ILFP - highest 15 MCs','ILFP - lowest 15 MCs')
    legend boxoff
    xlabel('Frequency (Hz) ')

end


