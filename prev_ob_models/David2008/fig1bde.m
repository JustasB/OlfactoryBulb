function fig1bde

tipsp = 150;
tsim = 200;
n_fig = 1;
dt=0.01;
vms=load('vmitfile.dat');

figure(n_fig);hold on
for i=[2:21];
    plot([0:(tsim+10)/dt],vms([(tipsp-10)/dt:(tipsp+tsim)/dt],i),'k-')
end
xlabel('t (ms)')
ylabel('V (mV)')

set(gca,'XTick',0:50/dt:tsim/dt)
set(gca,'XTickLabel',{'0','50','100','150','200'})

figure(n_fig+1);hold on
for i=1:20
    peakamp(i) = abs(min(vms(:,i+1))-vms(tipsp/dt-1,i+1));
end
plot([1:5:100],peakamp,'k^')
xlabel('IPSP location')
ylabel('IPSP amplitude')
set(gca,'XTick',0:20:100)
set(gca,'XTickLabel',{'0','20','40','60','80','100'})

figure(n_fig+2);hold on
for i=1:20
    peakdelay(i) = vms(min(find(vms(:,i+1)==min(vms(:,i+1)))),1)-tipsp;
end
plot([1:5:100],peakdelay,'k^')
xlabel('IPSP location')
ylabel('Peak delay')
set(gca,'XTick',0:20:100)
set(gca,'XTickLabel',{'0','20','40','60','80','100'})
