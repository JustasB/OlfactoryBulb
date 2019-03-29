% Code for running ParamSeep_GCE on the cronusx server

% Boleslaw Osinski (2015)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
input_file = 'OB_params_GCE.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

%% Run ParamSweep_GCE

% set P1
% P1.line = 45;
% P1.name = 'wGABAGR ';
% P1.val = 0.013:0.001:0.017;

P1.line = 7;
P1.name = 'Wmin ';
P1.val = 0.0104:0.0003:0.0164;




% set P2
% P2.line = 57;
% P2.name = 'Vrest ';
% P2.val = 1e-3 * [-75:1:-55];

P2.line = 45;
P2.name = 'wGABAGR ';
P2.val = 0.003:0.0017:0.026;

% Simulate the model over chosen parameter range
[MCIMAT, MCVMAT, param, MCilfpMAT, MCvlfpMAT] = ...
    ParamSweep_reduced_GCE(P1,P2,numtp,input_file);

%% Calculate LFP frequency and power

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

save('FmaxMAT.mat','FmaxMAT')
save('maxpwrMAT.mat','maxpwrMAT')
% save('MCIMAT.mat','MCIMAT')
% save('MCVMAT.mat','MCVMAT')
