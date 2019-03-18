function figure3

p_ipsc=12;
n_del=60;
n_fig=1;
pa1 = dlmread('fig3_0.500000.dat');
pa2 = dlmread('fig3_0.800000.dat');
pa3 = dlmread('fig3_1.800000.dat');

% pa1
    for p=1:p_ipsc
        for k=1:n_del
             isi(p,k) = pa1(k+(p-1)*n_del,11)-pa1(k+(p-1)*n_del,10);
             retard(p,k) = isi(p,k)-isi(1,1);
             del(p,k)= k-1;
        end
    end
    figure(n_fig);hold on
    for p=1:p_ipsc
        plot(del(p,2:end-5)*2*pi/(isi(p,n_del)),retard(p,2:end-5),'k','Linewidth',1)
    end
    xlabel('IPSC phase','Fontsize',12)
    set(gca,'Fontsize',12)
    set(gca,'XTick',[0,pi,2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    ylabel('\Deltat (ms)','Fontsize',12)
    axis([0 6.29 0 60])

% pa2   
    for p=1:p_ipsc
        for k=1:n_del
             isi(p,k) = pa2(k+(p-1)*n_del,11)-pa2(k+(p-1)*n_del,10);
             retard(p,k) = isi(p,k)-isi(1,1);
             del(p,k)= k-1;
        end
    end
    figure(n_fig+1);hold on
    for p=1:p_ipsc
        plot(del(p,2:end-27)*2*pi/(isi(p,n_del)),retard(p,2:end-27),'k','Linewidth',1)
    end
    xlabel('IPSC phase','Fontsize',12)
    set(gca,'Fontsize',12)
    set(gca,'XTick',[0,pi,2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    ylabel('\Deltat (ms)','Fontsize',12)
    axis([0 6.29 0 60])

% pa3
    for p=1:p_ipsc
        for k=1:n_del
             isi(p,k) = pa3(k+(p-1)*n_del,11)-pa3(k+(p-1)*n_del,10);
             retard(p,k) = isi(p,k)-isi(1,1);
             del(p,k)= k-1;
        end
    end
    figure(n_fig+2);hold on
    for p=1:p_ipsc
        plot(del(p,2:end-34)*2*pi/(isi(p,n_del)),retard(p,2:end-34),'k','Linewidth',1)
    end
    xlabel('IPSC phase','Fontsize',12)
    set(gca,'Fontsize',12)
    set(gca,'XTick',[0,pi,2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    ylabel('\Deltat (ms)','Fontsize',12)
    axis([0 6.29 0 60])