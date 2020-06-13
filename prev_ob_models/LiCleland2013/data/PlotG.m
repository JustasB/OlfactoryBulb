
% Plot PG-->MC and GC-->MC GABAA conductance in the OB network
% Written by Guoshi Li, Cornell University, 2013

clc;
clear all;
close all;

NTCE = 0;    % 0: With GCs, for the full model; generate figures similar to Fig. 7E-F
             % 1: no GCs, for the glomerular model, generate figures
             % similar to Fig. 5E-F
             

DT = 0.2;    % sampling time: ms
T1 = 1000;
T2 = 3000;
n1 = T1/DT+2;
n2 = T2/DT;

nMit = 25;
nPG  = 25;

load tt;
t  = tt(n1:n2);
t = (t-2000)/1000;

% Load PG-->MC conductances
for i = 0:1:(nMit-1)
     s = ['load GABApm' int2str(i) ';'];    
     eval(s);
     s = ['U=GABApm' int2str(i) ';'];
     eval(s);
     U = U(n1:n2);
     s = ['GABApm' int2str(i) '=U*1e3'  ';'];
     eval(s);

end


% Load GC-->MC conductances
if (NTCE==0)
  for i = 0:1:(nMit-1) 
     s = ['load Ggm' int2str(i) ';'];    
     eval(s);
     s = ['U=Ggm' int2str(i) ';'];
     eval(s);
     U = U(n1:n2);
     s = ['Ggm' int2str(i) '=U*1e3'  ';'];
     eval(s);
  end
end

%================================================
%               For PG-->MC 
%================================================
if (NTCE==1)
    figure;
    subplot(7,1,1);
    plot(t,GABApm0,'b-');
    title('G_P_G_-_>_M_C', 'FontSize',14);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    % axis([-1,1,0,15]);
    box('off');

    subplot(7,1,2);
    plot(t,GABApm3,'b');
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,3);
    plot(t,GABApm8,'b');
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,4);
    plot(t,GABApm12,'b');
    set(gca, 'XTickLabel',[ ]);
    ylabel('nS', 'FontSize',14);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,5);
    plot(t,GABApm16,'b');
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,6);
    plot(t,GABApm20,'b');
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,7);
    plot(t,GABApm24,'b');
    set(gca, 'FontSize',12);
    xlabel('Sec', 'FontSize',14);
    set(gca, 'YTickLabel',[ ]);
    box('off');

end

%=====================================
%           For GC-->MC 
%=====================================
if (NTCE==0)
    xmax = 1.001;
    ymax = 20;
    
    figure;
    subplot(3,1,1);
    plot(t,Ggm0,'b','LineWidth',2);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'FontSize',12);
    axis([-0.2,xmax,0,ymax]);
    % set(gca, 'YTick',[0:15:30]);
    set(gca, 'YTickLabel',[ ]);
    title('G_G_C_-_>_M_C', 'FontSize',14);
    legend('MC1');
    box('off');

    subplot(3,1,2);
    plot(t,Ggm12,'b','LineWidth',2);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'FontSize',12);
    axis([-0.2,xmax,0,ymax]);
    % set(gca, 'YTick',[0:15:30]);
    ylabel('nS', 'FontSize',14);
    legend('MC13');
    box('off');

    subplot(3,1,3);
    plot(t,Ggm23,'b','LineWidth',2);
    set(gca, 'FontSize',12);
    xlabel('Sec', 'FontSize',14);
    axis([-0.2, xmax, 0, ymax]);
    % set(gca, 'YTick',[0:15:30]);
    set(gca, 'YTickLabel',[ ]);
    legend('MC24');
    box('off');

end





