
% Phase distribution plot and Rasterplot of spike phases
% Written by Guoshi Li, Cornell University, 2013

clc;
clear all;
close all;

load tt;
load Vm;
load OdorA1.dat;
Odor = OdorA1;

nMit  = 25;
nPG   = 25;
nGran = 100;

Nm = 25;
Ng = 100;   

ratio = 0.3;       % factor of the average min-max distance

T1 = 2000;
T2 = 3000;

DT = 0.2;  
n1 = T1/DT+1;
n2 = T2/DT;
 
Fs = 1/DT*1000;    % sampling frequency: Hz
FILORDER = 1000;

Fmax = 100;        % maximal frequency to plot
Fc   = [20 100];   % Cut-off frequency
Wc   = Fc/(Fs/2);  % 

t = tt(n1:n2);
y = Vm(n1:n2);
y = y-mean(y);

L = length(y);
NFFT = 2^nextpow2(L);     % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2));

f = Fs/2*linspace(0,1,NFFT/2);


h = fir1(FILORDER, Wc);
x = filtfilt(h,1, y);

figure;
plot(t,x);
axis([2000, 3000, -8, 8]);
set(gca, 'FontSize',12);
xlabel('ms','FontSize',12);
title('Filtered sLFP', 'FontSize',12);
box('off');

figure;
[p,n]=peakdetect(x);

Np=length(p);
Nn=length(n);

PP=x(p);
PN=x(n);

if (p(1)<n(1))
   flag=1;   % Start with positive peak
else
   flag=0;   % Start with negative peak
end

if (flag==1)
  if(Np>Nn)
     D1=PP(1:end-1)-PN; 
     D2=PN-PP(2:end);
  else
     D1=PP-PN; 
     D2=PN(1:end-1)-PP(2:end);      
  end

else
   if(Np<Nn)
     D1=PP-PN(2:end); 
     D2=PN(1:end-1)-PP;
  else
     D1=PP(1:end-1)-PN(2:end); 
     D2=PN-PP;      
  end
end

DD=[D1; abs(D2)];
disp('The mean is:');
Dm=mean(DD)
threshold = ratio*Dm

% figure;
% hist(DD);

[maxtab, mintab] = peakdet(x, threshold, t);

%===============================================
%       Calcualte the phase of spikes
%===============================================

MT = maxtab(:,1);

TO1 = MT(1);        
TO2 = MT(end);   

Nc = length(MT);
PHASE = [];

for i = 0:1:(nMit-1)
   m=1;
   P=[];
   s = ['load Ms' int2str(i)  ';'];    
   eval(s);   
   ss = ['SpkT = Ms' int2str(i) ';'];    
   eval(ss);  
    
   A = find(SpkT>=TO1 & SpkT<=TO2); 
   L = length(A);  
    
   for j = 1:L
     ST = SpkT(A(j));
     
     for k = 2:Nc
       if (ST < MT(k))
         P(m)=(ST-MT(k-1))/(MT(k)-MT(k-1))*360;  
         if(P(m)>180)
           P(m)=P(m)-360;
         end
         m=m+1;
         break;
       end      
     end       
   end
   
   PHASE = [PHASE; P'];
   ss=['P' int2str(i+1) '=P;'];
   eval(ss);
end    


PHASE=PHASE';

N = length(PHASE);

% Calcualte the phase-locking index

SIN = 0;
COS = 0;

for i = 1:N
  psi = PHASE(i)/360*2*pi;  
  SIN = SIN + sin(psi);
  COS = COS + cos(psi);
end

Ksyn = 1/N*sqrt(SIN^2+COS^2);
disp('The spike phase syn. index is:');
Ksyn



%===============================================

figure;
subplot(2,1,1);
plot(t,x);
set(gca,'FontSize',12);
hold on;
plot(t(p), x(p),'r*');
hold on;
plot(t(n), x(n),'g*');
axis([2500, 3000, -8, 8]);
set(gca, 'XTickLabel',[ ]);
% xlabel('ms');
ylabel('mV', 'FontSize',14);
box('off');


subplot(2,1,2);
plot(t,x);
set(gca,'FontSize',12);
hold on; plot(mintab(:,1), mintab(:,2), 'g*');
plot(maxtab(:,1), maxtab(:,2), 'r*');
axis([2500, 3000, -8, 8]);
xlabel('ms', 'FontSize',14);
ylabel('mV', 'FontSize',14);
box('off');


%========================================

BIN  = -180:30:180;
Pbin = -165:30:165;
Ndist = histc(PHASE, BIN);
Ndist = Ndist(1:end-1);
P_dist= Ndist/sum(Ndist);

figure;
bar(Pbin, P_dist);
set(gca,'FontSize',12);
title('Control', 'FontSize',14);
set(gca, 'XTick',[-180:60:180]);
xlabel('Degree','FontSize',14);
ylabel('Probability','FontSize',14);
axis([-180, 180, 0, 0.5]);
box('off');

hold on;
x = -180:0.5:180;
xx = (2*pi*x)/(360);
y = 0.2*(sin(xx+0.5*pi)+1);
plot(x,y,'r-','LineWidth',2);


%=========================================
dx = 0.2;

figure;
for i = 1:1:nMit
   ss = ['PH = P' int2str(i) ';'];    
   eval(ss);  
   L = length(PH); 
   
   if(L==0)
      Pm(i) = 0;
   else
      Pm(i) = mean(PH);
   end
  
   if(L~=0)
      for (k=1:L)
        x = [i-dx   i+dx ];
        y = [PH(k)  PH(k)];
        plot(x,y,'r','LineWidth',1);
        hold on;
      end 
    end
      
end
plot([0 25.5],[0 0],'k:');

axis([0,25.5, -180,180]);
box('off');
xlabel('MC#', 'FontSize',14);
ylabel('Degree', 'FontSize',14);
title('Control', 'FontSize',14);
set(gca, 'FontSize',12);


Glom = 1:25;
width = 0.6;

figure;
subplot(2,1,1);
bar(Glom, Odor, width);
set(gca, 'FontSize',12);
set(gca,'XTickLabel',[ ]);
ylabel('nA', 'FontSize',14);
title('Control', 'FontSize',14);
axis([0 26 0 1.0]);
box('off');

subplot(2,1,2);
bar(Glom,Pm);
axis([0,25.5, -130,100]);
box('off');
xlabel('MC#', 'FontSize',14);
ylabel('Degree', 'FontSize',14);
set(gca, 'FontSize',12);


