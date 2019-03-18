function figure1cde

tipsp = 150;
tsim  = 200;
n_fig = 1;
dt    = 0.01;
vms   = load('vmitfile.dat');

figure(n_fig);hold on
for i=1:20;
    plot([0:(tsim+10)/dt],vms([(tipsp-10)/dt:(tipsp+tsim)/dt],i+1),'k-');
end
xlabel('T (ms)')
ylabel('V (mV)')
%axis([0 tsim/dt -71 -65.5]);
set(gca,'XTick',0:50/dt:tsim/dt)
set(gca,'XTickLabel',{'0','50','100','150','200'})

figure(n_fig+1);hold on
for i=1:19
    Rsdeq(i)=100*(1/((0.2+i*0.1)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
    peakamp(i) = abs(min(vms(:,i+1))-vms(tipsp/dt-1,i+1));
end
Rsdeq(20)=100*(1/((0.5+100*0.3)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
peakamp(20) = abs(min(vms(:,20+1))-vms(tipsp/dt-1,20+1));

plot((Rsdeq),peakamp,'ks')
xlabel('Rsd')
ylabel('IPSP amplitude (mv')

figure(n_fig+2);hold on
for i=1:19
    Rsdeq(i)=100*(1/((0.2+i*0.1)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
    peakdelay(i) = vms(min(find(vms(:,i+1)==min(vms(:,i+1)))),1)-tipsp;
end
Rsdeq(20)=100*(1/((0.5+100*0.3)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
peakdelay(20) = vms(min(find(vms(:,20+1)==min(vms(:,20+1)))),1)-tipsp;

plot(Rsdeq,peakdelay,'ks')
xlabel('Rsd')
ylabel('Peak delay (ms)')

