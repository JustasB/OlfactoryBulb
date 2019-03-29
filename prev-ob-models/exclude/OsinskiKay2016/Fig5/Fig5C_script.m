%% Simulate continuous GC excitability switch of a subpopulation of GCs (Fig 5C)


clear all
input_file = 'OB_params_GCE_Fig5C.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

% set P1
P1.line = 21;
P1.name = 'ExFrac ';
P1.val = 1:-0.1:0.3;


% set P2
P2.line = 57;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-75:1:-55];


% Simulate the model over chosen parameter range
[MCIMAT, MCVMAT, param, MCilfpMAT, MCvlfpMAT] = ...
    ParamSweep_reduced_GCE(P1,P2,numtp,input_file);


% Calculate LFP frequency and power
FmaxIMAT = zeros(length(P1.val),length(P2.val));
maxpwrIMAT = zeros(length(P1.val),length(P2.val));
FmaxVMAT = zeros(length(P1.val),length(P2.val));
maxpwrVMAT = zeros(length(P1.val),length(P2.val));

trim = 1000; % Ignore first 1000 timepoints

L = length(timevec(trim:end-100));  % Length of data that will be analyzed
NFFT = 2^nextpow2(L); % Next power of 2 from L
f = sampf/2*linspace(0,1,NFFT/2+1);
ROI = ceil(8/(f(2)-f(1))):ceil(140/(f(2)-f(1)));
% ROI gives us a window from ~7Hz - 115Hz from which to select peak power.
% This avoids selecting low freuqency components that may dominate the
% spectrum.

% Run FFT on each simulated LFP
for n1 = 1:length(P1.val)
    for n2 = 1:length(P2.val)
        
        MitILFP = MCilfpMAT{n1,n2};
        mitIFFT = fft(detrend(MitILFP(trim:end-100),'constant'),NFFT)/L;

        absmitIFFT = 2*abs(mitIFFT(1:NFFT/2+1));
        maxpwrIMAT(n1,n2) = max(absmitIFFT(ROI));

        maxind = absmitIFFT == maxpwrIMAT(n1,n2);
        FmaxIMAT(n1,n2) = f(maxind);
        
        MitVLFP = MCvlfpMAT{n1,n2};
        mitVFFT = fft(detrend(MitVLFP(trim:end-100),'constant'),NFFT)/L;

        absmitVFFT = 2*abs(mitVFFT(1:NFFT/2+1));
        maxpwrVMAT(n1,n2) = max(absmitVFFT(ROI));

        maxind = absmitVFFT == maxpwrVMAT(n1,n2);
        FmaxVMAT(n1,n2) = f(maxind);
    end
end

FmaxMAT.I = FmaxIMAT;
FmaxMAT.V = FmaxVMAT;
maxpwrMAT.I = maxpwrIMAT;
maxpwrMAT.V = maxpwrVMAT;

save('FmaxMAT_wAMPAp03_wNMDAp04_wN250_wGABAp0135_VrestvsExfrac.mat','FmaxMAT')
save('maxpwrMAT_wAMPAp03_wNMDAp04_wN250_wGABAp0135_VrestvsExfrac.mat','maxpwrMAT')
% save('MCIMAT.mat','MCIMAT')
% save('MCVMAT.mat','MCVMAT')


%% Plot results

clearvars
load FmaxMAT_wAMPAp03_wNMDAp04_wN250_wGABAp0135_VrestvsExfrac
load maxpwrMAT_wAMPAp03_wNMDAp04_wN250_wGABAp0135_VrestvsExfrac

Vrest = 1e-3 * [-75:1:-55];
ExFrac  = 1:-0.1:0.3;
winds = 1:(length(ExFrac));

fs = 14;

figure
set(gcf,'position',[0,400,350,180]);
hold on
for ii = 1:length(winds)
    plot(Vrest*1e3,FmaxMAT.I(winds(ii),:),'.-','color',[0.86*ii/length(winds),0.86*ii/length(winds),0.86*ii/length(winds)])
end
hold off
    set(gca,'fontsize',fs)
     xlabel('V_{rest,GC} (mV)');ylabel('LFP Fq (Hz)')
    xlim([-75 -55])
    ylim([10 80])

figure
set(gcf,'position',[0,400,350,180]);
hold on
for ii = 1:length(winds)
    plot(Vrest*1e3,maxpwrMAT.I(winds(ii),:),'.-','color',[0.86*ii/length(winds),0.86*ii/length(winds),0.86*ii/length(winds)])
end
plot(Vrest*1e3,1e-4*ones(1,length(Vrest)),'k--')
hold off
    set(gca,'fontsize',fs)
     xlabel('V_{rest,GC} (mV)');ylabel('Power')
    xlim([-75 -55])
    ylim([0 2.5e-3])

    
% Fig 5Cii inset
figure
set(gcf,'position',[100,100,130,60])
bar(flipud(maxpwrMAT.I(1:8,19)),'k');set(gca,'XTick',[])

