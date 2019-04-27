% Figure 7 simulates a fast change in GC excitability

% Note: tsim is set to 1000 for these simulations
% This script simulates the network for 10 different alpha values, but only
% 5 are plotted in Fig 7

% This script uses the tdepVrest family of functions

clearvars
input_file = 'OB_params_GCE_Fig7.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

% define Vrest as sigmoid over simulation time
tt = dt:dt:tsim;
% a = [0.01 0.1 1 10];% sigmoid width parameters
a = logspace(-2,0,10);
VrestGC = zeros(ntp,length(a));
for ii = 1:length(a)
    VrestGC(:,ii) = 1e-3*(-74+14./(1+exp(-a(ii)*(tt-tsim/2))));
end
% plot(tt,VrestGC,'.-','markersize',14)


MCVMAT = cell(length(a),1);
VLFPMAT = zeros(ntp,length(a));
ILFPMAT = zeros(ntp,length(a));

for ii = 1:length(a)
    [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = ...
        IandVLFP_GCE_tdepVrest(input_file, VrestGC(:,ii));
    ILFPMAT(:,ii) = MitLFPs.GradistMitGlobal(:,1);
    VLFPMAT(:,ii) = MitLFPs.VG(:,1);
    MCVMAT{ii} = zeros(ntp,nmit);
    for n = 1:nmit
        MCVMAT{ii}(:,n) = Mitral{n}.V;
    end
end

save('ILFPMAT_tsim1000_Ap01to1_VRGC74to60.mat','ILFPMAT')
save('VLFPMAT_tsim1000_Ap01to1_VRGC74to60.mat','VLFPMAT')
save('MCVMAT_tsim1000_Ap01to1_VRGC74to60.mat','MCVMAT')

%% Plot results

clearvars
load ILFPMAT_tsim1000_Ap01to1_VRGC74to60
load VLFPMAT_tsim1000_Ap01to1_VRGC74to60
load MCVMAT_tsim1000_Ap01to1_VRGC74to60
a = logspace(-2,0,10);
dt = 0.1; %ms
tsim = 1000; %ms
tt = dt:dt:tsim;
ntp = 10000;

% define Vrest
nmit = 45;
VrestGC = zeros(ntp,length(a));
for ii = 1:length(a)
    VrestGC(:,ii) = 1e-3*(-74+14./(1+exp(-a(ii)*(tt-tsim/2))));
end

% Wavelet transform
trim = 100; %trim beginning and end transients from LFP
subsamp = trim:4:(ntp-trim);% don't need full resolution data for wavelet transform
% Wavelet parameters
SF = 1000/dt;%(Hz)
wname='morl'; %select wavelet type; morlet works well
fc=centfrq(wname);
freqrange  = [20 320];
% From matlab page on centfrq: Specifically, scale is inversely proportional
% to frequency with the constant of proportionality being the center frequency
% of the wavelet.
scalerange = fc./(freqrange*(1/SF));
scales = scalerange(end):15:scalerange(1);
psfrq = fliplr(scal2frq(scales,wname,1/SF)); % get fq from scales
% Do continuous wavelet transform (CWT)
coefMAT = cell(length(a),1);
for ii = 1:length(a)
    % from ILFP
    coefs     = cwt(detrend(ILFPMAT(subsamp,ii)),scales,wname);
    coefs     = abs(coefs.*coefs);    
    coefs     = flipud(coefs);    % for some reason, need to flip the coefs mat
    coefMAT{ii} = coefs;
end

% plot LFPs and wavelet transform
fs = 12; % fontsize
YL = [-0.005 0.005]; % ylim for ILFP
CS = [0 3e-4]; % color scale for CWT
XL = cell(1,length(a)/2);XL{end} = '1s';
FH=figure;
set(gcf,'position',[10, 30, 600, 400])
% first plot has labeled axes
subplot(length(a)/2,2,1)
[hAx,hg1,hg2] = plotyy(tt,detrend(ILFPMAT(:,1)),tt,1000*VrestGC(:,1));
set(hg1,'color','k');set(hAx(1),'YColor','k')
ylim(hAx(1),YL);ylim(hAx(2),[-74 -60])
set(hAx(1),'fontsize',fs);set(hAx(2),'fontsize',fs)
% ylabel(hAx(1),'LFP');ylabel(hAx(2),'Vrest_{GC} (mV)')
set(gca,'Xtick',[])
subplot(length(a)/2,2,2)
pcolor(tt(subsamp),psfrq/4, coefMAT{1}); shading interp;
set(gca,'fontsize',fs);caxis(CS);
% colorbar('location','northoutside')
% ylabel('Hz')
set(gca,'Xtick',[]);
% remaining plots unlabeled
for ii = 2:length(a)/2
    subplot(length(a)/2,2,2*ii-1)
    [hAx,hg1,hg2] = plotyy(tt,detrend(ILFPMAT(:,2*ii-1)),tt,1000*VrestGC(:,2*ii-1));
    set(hg1,'color','k');set(hAx(1),'YColor','k');set(hAx(2),'YColor','k')
    ylim(hAx(1),YL);ylim(hAx(2),[-74 -60]);xlabel(XL{ii},'fontsize',14)
    set(hAx(1),'Xtick',[],'Ytick',[]);set(hAx(2),'Xtick',[],'Ytick',[])
    
    subplot(length(a)/2,2,2*ii)
    pcolor(tt(subsamp),psfrq, coefMAT{2*ii-1}); shading interp;
    set(gca,'Xtick',[],'Ytick',[]);caxis(CS);xlabel(XL{ii},'fontsize',16)
end
tightfig(FH)


