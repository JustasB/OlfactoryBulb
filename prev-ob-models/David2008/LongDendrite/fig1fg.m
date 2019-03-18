function fig1fg
%Longdendrite

n_fig = 1;
dt=0.01;
Iinj = 0.1;

vms=load('fig1fg.dat');

figure(n_fig);hold on
for i=1:20
   vref      = vms(490/dt,i+1);
   v2tiers   = 0.66*(vms(end,i+1)-vref)+vref;
   timeok    = vms(max(find(floor(vms(:,i+1)/0.01)==floor(v2tiers/0.01))),1);
   deltaV(i) = vms(end,i+1)-vms(490,i+1);
   Rin(i)    = abs(deltaV(i)/Iinj);
   taum(i)   = vms(floor(timeok/dt),1)-500;
end


figure(n_fig);
hold on
plot([1:5:100],Rin,'k^')
xlabel('Shunt Location')
ylabel('R_{inp}','Fontsize',14)


figure(n_fig+1); hold on
plot([1:5:100],taum,'k^')
xlabel('Shunt Location')
ylabel('\tau_m','Fontsize',14)

%axis([0 100 0 12])