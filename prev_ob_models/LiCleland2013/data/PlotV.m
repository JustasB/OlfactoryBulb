
% Plot cell voltages in the OB network
% Written by Guoshi Li, Cornell University, 2013

clc;
clear all;
close all;

NTCE = 0;    % 1: For the glomerular model, produce plots similar to Fig. 5A-D in the paper
             % 0: For the full model, produce plots similar to Fig. 7A-D in the paper

Nd = 10;
dt = 0.02;   % simulation step: ms
DT = 0.2;    % sampling time: ms

T1 = 1000;
T2 = 3000;
n1 = T1/DT+2;
n2 = T2/DT;

nMit  = 25;
nPG   = 25;
nGran = 100;

load tt;

t  = tt(n1:n2);

for i = 0:1:(nMit-1)
     s = ['load Vms' int2str(i) ';'];  % Soma voltage   
     eval(s);
     s = ['U=Vms' int2str(i) ';'];
     eval(s);
     U = U(n1:n2);
     s = ['ms' int2str(i) '=U'  ';'];
     eval(s);    
end

for i = 0:1:(nPG-1)
    s = ['load Vpb' int2str(i) ';'];    
    eval(s);    
    s = ['U=Vpb' int2str(i) ';'];
    eval(s);
    U = U(n1:n2);
    s = ['pb' int2str(i)  '=U'  ';'];
    eval(s);     
end

for i = 0:1:(nGran-1)
    s = ['load Vgb' int2str(i) ';'];    
    eval(s);    
    s = ['U=Vgb' int2str(i) ';'];
    eval(s);
    U = U(n1:n2);
    s = ['gb' int2str(i)  '=U'  ';'];
    eval(s);     
end


xmax = 1.001;

if (NTCE==1)    
% Plot MC somatic voltages
    figure;
   subplot(7,1,1);
   plot(t,ms0,'b', 'LineWidth',1);
   axis([T1,T2,-80,40]);
   title('MC', 'FontSize',14);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   set(gca, 'FontSize',12);
   box('off');

   subplot(7,1,2);
   plot(t,ms3,'b', 'LineWidth',1);
   axis([T1,T2,-80,40]);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   set(gca, 'FontSize',12);
   box('off');

   subplot(7,1,3);
   plot(t,ms8,'b', 'LineWidth',1);
   axis([T1,T2,-80,40]);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   set(gca, 'FontSize',12);
   box('off');

   subplot(7,1,4);
   plot(t,ms12,'b');
   axis([T1,T2,-80,40]);
   ylabel('mV', 'FontSize',14);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'FontSize',12);
   box('off');

   subplot(7,1,5);
   plot(t,ms16,'b');
   axis([T1,T2,-80,40]);
   set(gca, 'FontSize',12);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   box('off');

   subplot(7,1,6);
   plot(t,ms20,'b');
   axis([T1,T2,-80,40]);
   set(gca, 'FontSize',12);
   set(gca, 'XTickLabel',[ ]);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   box('off');

   subplot(7,1,7);
   plot(t,ms24,'b');
   axis([T1,T2,-80,40]);
   set(gca, 'FontSize',12);
   set(gca, 'YTick',[-80:40:40]);
   set(gca, 'YTickLabel',[ ]);
   xlabel('Sec', 'FontSize',14);
   box('off');


% Plot PG spine voltage
    figure;
    subplot(7,1,1);
    plot(t,pb0,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    title('PGC', 'FontSize',14);
    box('off');

    subplot(7,1,2);
    plot(t,pb3,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,3);
    plot(t,pb8,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'YTickLabel',[ ]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,4);
    plot(t,pb12,'b');
    axis([T1,T2,-80,40]);
    ylabel('mV', 'FontSize',14);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'FontSize',12);
    box('off');

    subplot(7,1,5);
    plot(t,pb16,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'FontSize',12);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'YTickLabel',[ ]);
    box('off');

    subplot(7,1,6);
    plot(t,pb20,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'FontSize',12);
    set(gca, 'XTickLabel',[ ]);
    set(gca, 'YTick',[-80:40:40]);
    set(gca, 'YTickLabel',[ ]);
    box('off');

    subplot(7,1,7);
    plot(t,pb24,'b');
    axis([T1,T2,-80,40]);
    set(gca, 'FontSize',12);
    set(gca, 'YTick',[-80:40:40]);
    xlabel('Sec', 'FontSize',14);
    set(gca, 'YTickLabel',[ ]);
    box('off');

end


if (NTCE==0)
  tv = (t-2000)/1000;
% Plot MC cells
  figure;
  subplot(2,1,1);
  plot(tv,ms0,'b', tv,ms23,'r','LineWidth',2);   % ms0 & ms23 || ms1 & ms23
  set(gca, 'FontSize',12);
  set(gca, 'YTick',[-80:40:40]);
  ylabel('mV', 'FontSize',14);
  title('MC', 'FontSize',14);
  legend('MC1','MC24');
  axis([-0.2,xmax,-80,40]);
  box('off');  
  
  subplot(2,1,2);
  plot(tv,ms10,'b', tv,ms12,'r','LineWidth',2);  % ms10 & 12||ms10 & ms14
  axis([-0.2,xmax,-80,40]);
  set(gca, 'FontSize',12);
  set(gca, 'YTick',[-80:40:40]);
  xlabel('Sec', 'FontSize',14);
  ylabel('mV', 'FontSize',14);
  legend('MC11','MC13');  
  box('off');


% Plot GC spine voltage
  figure;
  subplot(2,1,1);
  plot(tv,gb13,'b', tv,gb43,'r','LineWidth',2);  % gb0|gb13 & gb43
  axis([-0.2,xmax,-80,40]);
  set(gca, 'FontSize',12);
  title('GC', 'FontSize',14);
  set(gca, 'YTick',[-80:40:40]);
  ylabel('mV', 'FontSize',14);
  legend('GC13','GC43');
  box('off');  
  
  subplot(2,1,2);
  plot(tv,gb66,'b', tv,gb92,'r','LineWidth',2);  % gb66 & gb92
  axis([-0.2,xmax,-80,40]);
  set(gca, 'FontSize',12);
  set(gca, 'YTick',[-80:40:40]);
  xlabel('Sec', 'FontSize',14);
  ylabel('mV', 'FontSize',14);
  legend('GC66','GC92');  
  box('off');

end





