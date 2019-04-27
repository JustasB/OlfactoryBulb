% Code for generating Figure 3

% This script simulates the currents of a single dGC connected to 14 MCs
% Parameters are in OB_params_GCE_Fig3.txt

% Boleslaw Osinski (2015)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
input_file = 'OB_params_GCE_Fig3.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

% P1 turns Calcium dependent incativation (CDI) on and off
P1.line = 19;
P1.name = 'hCaflag ';
P1.val = {'true';'false'};

% P2 controls Vrest, the GC excitability
P2.line = 57;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-70,-60];

% Simulate the model over chosen parameter range
[MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file);


%% Part A

ysc = [0 0.005]; % y scale
xsc = [0 200];% x scale
% plot ImitgradistAMPA
scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.4*scrsz(4)]);
subplot(3,1,1)
plot(timevec,ICMAT{1,1}.ImitgradistAMPA/2,'b',timevec,ICMAT{1,2}.ImitgradistAMPA/2,'r')
set(gca,'fontsize',22)
legend('Low Excitability  (-70 mV)','High Excitability (-60 mV)','location','EastOutside')
legend boxoff
xlim([40 65])
ylim([0 0.01])
title(['AMPA current'])

% plot ImitgradistNMDA;
subplot(3,1,2)
plot(timevec,ICMAT{1,1}.ImitgradistNMDA*1.5,'b',timevec,ICMAT{1,2}.ImitgradistNMDA*1.5,'r')
set(gca,'fontsize',22)
xlim(xsc)
ylim(ysc)
title(['NMDA current'])

% plot ImitgradistVDCC
subplot(3,1,3)
hold on
plot(timevec,ICMAT{1,1}.ImitgradistVDCC,'b',timevec,ICMAT{1,2}.ImitgradistVDCC,'r')
plot(timevec,ICMAT{2,1}.ImitgradistVDCC,'b--',timevec,ICMAT{2,2}.ImitgradistVDCC,'r--')
set(gca,'fontsize',22)
xlim(xsc)
ylim(ysc)
title(['N-type current '])
xlabel('time (ms) ');


%% Parts C and D

% plot [Ca] and P_release together
scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
hold on
h1 = plot(timevec,ICMAT{1,1}.CCa,'b')
h2 = plot(timevec,ICMAT{2,1}.CCa,'b--')
plot(timevec,ICMAT{1,2}.CCa,'r')
plot(timevec,ICMAT{2,2}.CCa,'r--')
hold off
set(gca,'fontsize',22)
xlim(xsc)
ylim([0 1.4])
legend([h1,h2],{'With CDI ','Without CDI '})
legend boxoff
xlabel('time (ms) ');ylabel('[Ca] (\muM) ')

scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
hold on
plot(timevec,ICMAT{1,1}.Prelease,'b')
plot(timevec,ICMAT{2,1}.Prelease,'b--')
plot(timevec,ICMAT{1,2}.Prelease,'r')
plot(timevec,ICMAT{2,2}.Prelease,'r--')
hold off
set(gca,'fontsize',22)
xlim(xsc)
ylim([0 1])
xlabel('time (ms) ');ylabel('P_{release} ')


%% Part B

% plot m
scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.4*scrsz(4)]);
subplot(3,1,1)
hold on
plot(timevec,ICMAT{1,1}.mgradist,'b');
plot(timevec,ICMAT{2,1}.mgradist,'b--')
plot(timevec,ICMAT{1,2}.mgradist,'r');
plot(timevec,ICMAT{2,2}.mgradist,'r--')
hold off
set(gca,'fontsize',22)
xlim(xsc)
title('N-type activation (m, V-dependent) ')

% plot h
subplot(3,1,2)
hold on
h1 = plot(timevec,ICMAT{1,1}.hgradist,'b')
plot(timevec,ICMAT{2,1}.hgradist,'b--')
h2 = plot(timevec,ICMAT{1,2}.hgradist,'r')
plot(timevec,ICMAT{2,2}.hgradist,'r--')
hold off
legend([h1,h2],{'Low Excitability  (-70 mV)','High Excitability (-60 mV)'})
legend boxoff
set(gca,'fontsize',22)
set(gca,'YTickLabel',[0 0.0005 0.001])
xlim(xsc)
ylim([0 0.001])
title('N-type inactivation (h, [Ca]-dependent) ')

% plot m*h
subplot(3,1,3)
hold on
plot(timevec,ICMAT{1,1}.hgradist.*ICMAT{1,1}.mgradist,'b')
plot(timevec,ICMAT{2,1}.hgradist.*ICMAT{2,1}.mgradist,'b--')
plot(timevec,ICMAT{1,2}.hgradist.*ICMAT{1,2}.mgradist,'r')
plot(timevec,ICMAT{2,2}.hgradist.*ICMAT{2,2}.mgradist,'r--')
hold off
set(gca,'fontsize',22)
xlim(xsc)
title('m * h ')
xlabel('time (ms) ');
