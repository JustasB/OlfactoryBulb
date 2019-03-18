function figure1fg

% Davison cell
n_fig = 1;
dt    = 0.01;
Iinj  = 0.1;

vms=load('fig1fg.dat');

figure(n_fig); hold on
for i=1:20
    vref      = vms(490/dt,i+1);
    v2tiers   = 0.66*(vms(end,i+1)-vref)+vref;
    timeok    = vms(max(find(floor(vms(:,i+1)/0.01)==floor(v2tiers/0.01))),1);
    
    Rsdeq(i) = 100*(1/((0.2+i*0.1)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
    if i==20
       Rsdeq(i) = 100*(1/((0.5+100*0.3)*(1.94e-10)))/(1/((0.2+1*0.1)*(1.94e-10)));
    end
    deltaV(i) = vms(end,i+1) -vms(490,i+1);
    Rin(i)    = abs(deltaV(i)/Iinj);
    taum(i)   = vms(floor(timeok/dt),1)-500;
end

figure(n_fig);
hold on
plot(Rsdeq,Rin,'ks')
xlabel('Intermediate resistance (\Omega.cm)')
ylabel('Input Resistance (\Omega/cm^2)','Fontsize',14)

figure(n_fig+1); hold on
plot(Rsdeq,taum,'ks')
xlabel('Intermediate resistance (\Omega.cm)')
ylabel('\tau_m (ms)','Fontsize',14)
axis([0 100 0 12])