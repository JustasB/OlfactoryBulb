
% Perform frequency analysis of the sLFP in the OB network
% Written by Guoshi Li, Cornell University, 2013

clc;
clear all;
close all;

load tt;         
load Vm;     % sLFP, mean membrane somatic voltage of all MCs   
% load Vg;
FILORDER = 1000;

DT = 0.2;          % sampling time: ms
Fs = 1/DT*1000;    % sampling frequency: Hz

T1 = 2000;
T2 = 3000;
n1 = T1/DT+1;
n2 = T2/DT;
maxlags = 2000;    % For auto-correlation!  

Fmax = 100;        % maximal frequency to plot
Fc   = [20 100];   % Cut-off frequency
Wc   = Fc/(Fs/2);  % 

t = tt;
t = t(n1:n2);
y = Vm(n1:n2);
y = y-mean(y);


L = length(y);
NFFT = 2^nextpow2(L);     % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2));
% YY = abs(Y(1:NFFT/2)).^2;

% Y  = fft(y,NFFT);
% YY = 2*abs(Y)/NFFT;
% YY = YY(1:end/2);

f = Fs/2*linspace(0,1,NFFT/2);

m = Fmax/(0.5*Fs)*(0.5*NFFT);
m = floor(m);

%=================================================
xmin = 1000;
xmax = 2000;

% figure;
% plot(t,y,'b');
% title('Original Signal');
% % axis([xmin, xmax, -80, -20]);

% % Plot single-sided amplitude spectrum.
figure;
plot(f(1:m),YY(1:m));
title('Single-Sided Amplitude Spectrum of sLFP')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')


%=================================================

h = fir1(FILORDER, Wc);
x = filtfilt(h,1, y);

% figure;
% freqz(h, 1, 512);


%=========================
% [b, a]=butter(3, Wc);
% x = filtfilt(b, a, y);
% x = filter(b, a, y);
% figure;
% freqz(b, a, 512, Fs);

% hd = dfilt.dffir(h);
% x = filter(hd, y);

%=========================

X = fft(x,NFFT)/L;
XX = 2*abs(X(1:NFFT/2));
[Peak, I] = max(XX);
disp('The oscillation frequency is:');
f(I)
disp('The spectrum peak is:');
Peak

%=======================================
%         Auto-correlation 
%=======================================
u  = mean(x);
yn = x-u;

% [cy, lags] = xcorr(y,'unbiased');
[cy, lags] = xcorr(yn, maxlags,'coeff');

% tau = -(L-1):(L-1);
% cc  = xcov(x, 'coef');
% figure;
% plot(tau, cc);
% axis([-2000,2000,-0.4,1]);

for k=(maxlags+2):length(cy)
%     if cy(k)>cy(k-1) && cy(k)>cy(k+1)&& cy(k)>0
      if cy(k)>cy(k-1) && cy(k)>cy(k+1)  
      break;
    end
end

disp('The oscillation index is:');
Power = cy(k)

% disp('The index is:');
% PI=k-maxlags-1

% disp('The oscillation frequency is:');
% fo = 1/(PI*DT)*1000


xmin1 = 2000;
xmax1 = 3000;


figure;
plot(t, x, 'LineWidth',0.5);
xlabel('ms', 'FontSize',14);
ylabel('mV', 'FontSize',14);
title('Filtered LFP', 'FontSize',14);
set(gca, 'FontSize',12);
axis([xmin1, xmax1, -10, 10]);
box('off');

% Plot single-sided amplitude spectrum.
figure;
% plot(f, 2*abs(Y(1:NFFT/2))) 
plot(f(1:m),XX(1:m));
title('FFT Spectrum', 'FontSize',14)
xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
set(gca, 'FontSize',12);
% axis([0, 150, 0, 2]);
box('off');

% Plot auto-correlation of LFP
figure;
plot(lags, cy);
title('Auto-correlation fo sLFP','FontSize',14);


%=======================================
%     Plot LFP and Auto-Correlation
%=======================================

figure;
subplot(3,1,1);
plot(t-2000, x, 'LineWidth',1);
set(gca, 'FontSize',12);
xlabel('ms', 'FontSize',12,'FontWeight','bold');
ylabel('LFP (mV)', 'FontSize',12,'FontWeight','bold');
axis([0, 1000, -6, 6]);
box('off');

subplot(3,1,2);
plot(lags*DT, cy, 'LineWidth',1);
xlabel('Lags (ms)', 'FontSize',12, 'FontWeight','bold')
set(gca, 'FontSize',12);
axis([-400, 400, -1.0, 1]);
box('off');

subplot(3,1,3);
plot(f(1:m),XX(1:m), 'LineWidth',1);
xlabel('Frequency (Hz)', 'FontSize',12,'FontWeight','bold')
ylabel('Power', 'FontSize',12,'FontWeight','bold')
set(gca, 'FontSize',12);
axis([0, 100, 0, 2.0]);
box('off');




