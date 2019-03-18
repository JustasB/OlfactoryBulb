
% Generate raster plots in the OB network
% Written by Guoshi Li, Cornell University, 2013
% Odor is presented at T = 2000 ms

clc;
clear all;
close all;

NTCE = 0;     % NTCE=0: Plot GCs; otherwise, plot MCs and PGs only

load OdorA1.dat;   % load steady-state odor values
Odor = OdorA1;

nMit  = 25;
nPG   = 25;
nGran = 100;

Nm = 25;
Ng = 100;

dm = 0.4;
dp = 0.4;
dg = 0.4;

T_Start = 1000;          % start time of calculation     
T_End   = 3000;          % end time of calculation  
T  = T_End - T_Start;    % Total duration in ms

min_T = T_Start;
max_T = T_End;

TP1 = 1000;     % Start of spontaneous activity
TP2 = 2000;     % End of spontaneous activity
TO1 = 2000;     % Start of odor stimulus
TO2 = 3000;     % End of odor stimulus
TP  = (TP2-TP1)/1000;
TO  = (TO2-TO1)/1000;

%============================================
%        Generate raster plot 
%============================================
% for mitral cells
figure;

for i = 0:1:(nMit-1)
       
    n = i+1;
    s = ['load Ms' int2str(i)  ';'];    
    eval(s);   
   
    ss = ['SpkT = Ms' int2str(i) ';'];    
    eval(ss);  
   
   % Spontanous rate 
    A = find (SpkT>=TP1 & SpkT<TP2); 
    FM_SP(n,1) = length(A)/TP;
    
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    FM(n,1) = length(A)/TO;
   
    L = length(SpkT);
    if (L~=0)  
     for k = 1:L
      if (SpkT(k) > T_Start)
      
       x = [SpkT(k) SpkT(k)];
       y = [n-dm    n+dm ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end
    end
   end
   
end


xlabel('ms', 'FontSize',14);
ylabel('MC #', 'FontSize',14);
set(gca, 'FontSize',12);
axis([1000,max_T,0,26]);
box('off');



% For PG cells
figure;
for i = 0:1:(nPG-1)
    
    n=i+1;  
    
    s = ['load Pd' int2str(i) ';'];    
    eval(s);
    
    ss = ['SpkT = Pd' int2str(i) ';'];    
    eval(ss);  
    
   % Spontanous rate 
    A = find (SpkT>=TP1 & SpkT<TP2); 
    FP_SP(n,1) = length(A)/TP;
    
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    FP(n,1) = length(A)/TO;
    
   L = length(SpkT);
   if (L~=0)
    for k = 1:L
     if (SpkT(k) > T_Start)
      
      x = [SpkT(k) SpkT(k)];
      y = [n-dp    n+dp ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end  
    end
   end 
    
end

  
box('off');
xlabel('ms', 'FontSize',14);
ylabel('PG #', 'FontSize',14);
set(gca, 'FontSize',12);
axis([1000,max_T,0,26]);


% % For granule cells

if(NTCE==0)
    
figure;
for i = 0:1:(nGran-1)
    
    n=i+1;  
    
    s = ['load Gd' int2str(i) ';'];    
    eval(s);
    
    ss = ['SpkT = Gd' int2str(i) ';'];    
    eval(ss);  
   
   % Spontanous rate 
    A = find (SpkT>=TP1 & SpkT<TP2); 
    FG_SP(n,1) = length(A)/TP;
    
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    FG(n,1) = length(A)/TO;
    
   L = length(SpkT);
   if (L~=0)
    for k = 1:L
     if (SpkT(k) > T_Start)
      
      x = [SpkT(k) SpkT(k)];
      y = [n-dg    n+dg ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end  
    end
   end 
    
  end


axis([1000,max_T, 0,101]);
box('off');
xlabel('ms', 'FontSize',14);
ylabel('GRAN #', 'FontSize',14);
set(gca, 'FontSize',12);

end

%=======================================

fM_SP = mean(FM_SP);
fP_SP = mean(FP_SP);

if (NTCE==0)
fG_SP = mean(FG_SP);
fG = mean(FG);
end

fM = mean(FM);
fP = mean(FP);

disp('The spontaneous firing rates are:');
fM_SP
fP_SP

if (NTCE==0)
  fG_SP
end

disp('The firing rates during odor presentation is:');
fM
fP

if (NTCE==0)
  fG
end


%===================================================
%               Plot Input-Output
%===================================================
width = 0.6;
glom = 1:25;
     
Fm = FM-FM_SP;   % Odor coding rate

figure;
subplot(3,1,1);
bar(glom, Odor, width);
box('off');
axis([0,26,0,1.05]);
ylabel('Input (nA)', 'FontSize',12);

subplot(3,1,2);
bar(glom, FM, width);
box('off');
axis([0,26,0,30]);
ylabel('Odor-evoked Rate', 'FontSize',10);

subplot(3,1,3);
bar(glom, Fm, width);
axis([0,26,-10,30]);
xlabel('Glom #');
ylabel('Odor-coding Rate', 'FontSize',10);
box('off');




