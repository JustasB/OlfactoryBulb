% Code for running ParamSeep_GCE on the cronusx server

% Boleslaw Osinski (2015)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
input_file = 'OB_params_GCE.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

%% Run ParamSweep_GCE

% set P1
% P1.line = 44;
% P1.name = 'wGABAGR ';
% P1.val = 0.013:0.001:0.017;

P1.line = 7;
P1.name = 'Wmin ';
P1.val = 0.0128;


% set P2
P2.line = 57;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-75];

% Simulate the model over chosen parameter range
[MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file);

%% Calculate LFP frequency and power

FmaxMAT = zeros(length(P1.val),length(P2.val));
maxpwrMAT = zeros(length(P1.val),length(P2.val));

trim = 500; % Ignore first 500 timepoints

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
        
        MitLFP = MClfpMAT{n1,n2}.GradistMitGlobal;
        mitFFT = fft(detrend(MitLFP(trim:end-100,1),'constant'),NFFT)/L;

        absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
        maxpwrMAT(n1,n2) = max(absmitFFT(ROI));

        maxind = find(absmitFFT == maxpwrMAT(n1,n2));
        FmaxMAT(n1,n2) = f(maxind);
    end
end


save('FmaxMAT.mat','FmaxMAT')
save('maxpwrMAT.mat','maxpwrMAT')
save('ICMAT.mat','ICMAT')
