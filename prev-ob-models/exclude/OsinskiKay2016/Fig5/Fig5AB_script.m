% Figure 5 shows the VLFP and ILFP frequency and power over the full range
% of VrestGC, along with the Spike Synchrony Index



clear variables
input_file = 'OB_params_GCE_Fig5AB.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);


%% Calculate ILFP, VLFP, SFD, & SFC vs VrestGC, 10 trials (Fig. 5A & B)

% set variable parameter
P.line = 57;
P.name = 'Vrest ';
P.val = 1e-3 * (-75:1:-55);
% P.val = 1e-3 * (-75:10:-55);

% initialize params for FFT
trim = 1000; % trim beginning and end to avoid edge effects
% trim:end-50
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);
ROI = ceil(8/(f(2)-f(1))):ceil(140/(f(2)-f(1)));

% Initialize params for SFC
params.tapers = [5 9];
params.pad = 0;
params.Fs = 10000;
params.fpass = [3 100];
params.err = [2 1];
params.trialave = 0;


TextCell = regexp( fileread(input_file), '\n', 'split');

numtrials = 10;


FmaxMAT.I = zeros(length(P.val),numtrials);
maxpwrMAT.I = zeros(length(P.val),numtrials);
SFD.I = zeros(length(P.val),numtrials);
SFC.I = zeros(length(P.val),numtrials);
FmaxMAT.V = zeros(length(P.val),numtrials);
maxpwrMAT.V = zeros(length(P.val),numtrials);
SFD.V = zeros(length(P.val),numtrials);
SFC.V = zeros(length(P.val),numtrials);


for vv = 1:length(P.val)
    for tt = 1:numtrials
        % write VrestGC to parameter file
        wstring = [P.name,num2str(P.val(vv))];
        TextCell{P.line} = sprintf('%s',wstring);
        fid = fopen(input_file, 'w');
        fprintf(fid, '%s\n', TextCell{:});
        fclose(fid);

        % simulate model
        [Mitral, ~, ~, param, ~, MitLFPs, ~] = ...
            IandVLFP_GCE(input_file);

        % Power spectra

        % ILFP
        mitFFTGI = fft(detrend(MitLFPs.GradistMitGlobal(trim:end-100),'constant'),NFFT)/L;
        absmitIFFT = 2*abs(mitFFTGI(1:NFFT/2+1));
        maxpwrMAT.I(vv,tt) = max(absmitIFFT(ROI));
        maxind = absmitIFFT == maxpwrMAT.I(vv,tt);
        FmaxMAT.I(vv,tt) = f(maxind);

        % VLFP
        mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;
        absmitVFFT = 2*abs(mitFFTGV(1:NFFT/2+1));
        maxpwrMAT.V(vv,tt) = max(absmitVFFT(ROI));
        maxind = absmitVFFT == maxpwrMAT.V(vv,tt);
        FmaxMAT.V(vv,tt) = f(maxind);

        % Calculate SFD
        nspikes = 0;
        % NOTE: only sum spikes 100ms after stimulus onset
        for n = 1:nmit
            nspikes = nspikes + sum(Mitral{n}.S(1000:end));
        end
        SpCmaxMATI = ceil(nmit*0.6*FmaxMAT.I(vv,tt));
        SpCmaxMATV = ceil(nmit*0.6*FmaxMAT.V(vv,tt));
        % NOTE: we are only looking at the last 0.6s of the simulation
        SFD.I(vv,tt) = abs(nspikes-SpCmaxMATI);
        SFD.V(vv,tt) = abs(nspikes-SpCmaxMATV);
        
        % Calculate SFC
        spktcell = cell(nmit,1);
        for n = 1:nmit
            spktcell{n} = find(Mitral{n}.S(trim:end) == 1)'/sampf;
        end
        SpkTIMES = struct('f',spktcell);
        
        dataI = repmat(detrend(MitLFPs.GradistMitGlobal(trim:(end-100)))',1,nmit);
        dataV = repmat(detrend(MitLFPs.VG(trim:(end-100)))',1,nmit);
        
        % calculate SFC with finite size corrections
        [CcorrI,~,~,~,~,~,~,~,~,~]=coherencycpt(dataI,SpkTIMES,params,1);
        [CcorrV,~,~,~,~,~,zerosp,~,~,~]=coherencycpt(dataV,SpkTIMES,params,1);
        
        CcorrI(:,logical(zerosp)) = 0;
        CcorrV(:,logical(zerosp)) = 0;

        SFC.I(vv,tt) = mean(max(CcorrI));
        SFC.V(vv,tt) = mean(max(CcorrV));
        
    end
end

save('FmaxMAT.mat','FmaxMAT');
save('maxpwrMAT.mat','maxpwrMAT');
save('SFD.mat','SFD');
save('SFC.mat','SFC');


%% Plot simulated data

clearvars

load FmaxMAT
load maxpwrMAT
load SFD
load SFC

Vrest = 1e-3 * (-75:1:-55);
numtrials = 10;


fs = 12;% fontsize


% Fig 5A.i
figure
set(gcf,'position',[0,400,250,200]);
hold on
shadedErrorBar(Vrest*1e3,mean(FmaxMAT.V,2),std(FmaxMAT.V,0,2),'-k',0.5)
shadedErrorBar(Vrest*1e3,mean(FmaxMAT.I,2),std(FmaxMAT.I,0,2),'-m',0.5)
hold off
set(gca,'fontsize',fs)
xlabel('V_{rest,GC} (mV)');ylabel('LFP Fq (Hz)')
xlim([-75 -55]);ylim([10 90])

% Fig 5A.ii
% plot N-type and NMDA together
E_N = 0;
v = 1e3*((E_N-0.1):0.001:(E_N+0.1)); % in mV
V = -0.1:0.001:0.1; % in V
mbarN = 1./(1 + exp(-(v+45)./7));
Mg_conc = 1;
E_Mg = 0;
gamma = 0.016;
eta = 0.28;
Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma));

figure
set(gcf,'position',[0,400,300,170]);
plot(v,mbarN,'-k',v,Mg_block,':k','LineWidth',2);xlim([-80 0])
set(gca,'fontsize',fs)
legend('N-Type','NMDA','location','southeast')
legend boxoff
xlabel('V_{rest,GC} (mV)');ylabel('Activation')


% Fig 5A.iii
figure
set(gcf,'position',[0,400,250,200]);
hold on
hV = shadedErrorBar(Vrest*1e3,mean(maxpwrMAT.V,2),std(maxpwrMAT.V,0,2),'-k',0.5)
hI = shadedErrorBar(Vrest*1e3,mean(maxpwrMAT.I,2),std(maxpwrMAT.I,0,2),'-m',0.5)
plot(Vrest*1e3,1e-4*ones(1,length(Vrest)),'k--')
hold off
set(gca,'fontsize',fs)
xlabel('V_{rest,GC} (mV)');ylabel('Power')
legend([hI.mainLine,hV.mainLine],{'ILFP','VLFP'});legend boxoff
xlim([-75 -55]);ylim([0 2.3e-3])

% Fig 5B
% Plot SFD and SFC on same plot
figure
set(gcf,'position',[0,400,260,150]);
[hAx, hL1, hL2] = plotyy(Vrest*1e3,mean(SFD.I,2),Vrest*1e3,nanmean(SFC.I,2));
hAx(1).YLim = [0 500];
hAx(1).YTick = [0 100 200 300 400 500];
hAx(2).YLim = [0 1];
hAx(2).YTick = [0 0.2 0.4 0.6 0.8 1];
hL1.Color = 'k'; hL1.LineStyle = '-'; hL1.LineWidth = 2;
hL2.Color = 'k'; hL2.LineStyle = ':'; hL2.LineWidth = 2;
legend('SFD','SFC','location','southeast'); legend boxoff
set(hAx(1),'YColor','k')
set(hAx(2),'YColor','k')
set(hAx(1),'fontsize',fs)
set(hAx(2),'fontsize',fs)
% xlabel('V_{rest,GC} (mV)')
% ylabel('Spike-Freq. Dev.')
xlim([-75 -55])


