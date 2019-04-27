function fig2ace

n_del = 60 ;    % number of delay step
nmit  = 10;     % number of neurons
n_fig = 10;     % figure index

for ind =1:3
    
    if ind == 1
       pa = dlmread('fig2a.dat');
    elseif ind== 2
       clear('pa')
       pa = dlmread('fig2b.dat');
    elseif ind == 3
       clear('pa')
       pa = dlmread('fig2c.dat');
    end
       
    for k=1:n_del
        for i=1:nmit
            isi(k,i) = pa((k-1)*nmit+i,11)-pa((k-1)*nmit+i,10);
            retard(k,i) = isi(k,i)-isi(1,1);
            del(k,i)= k-1;
        end
    end

    figure(n_fig+ind-1); hold on
    for i=1:nmit
        plot(del(2:end,i)*2*pi/(isi(n_del,i)),retard(2:end,i),'k','Linewidth',1)
    end

    xlabel('IPSC phase','Fontsize',12)
    set(gca,'Fontsize',12)
    set(gca,'XTick',[0,pi,2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    ylabel('\Deltat (ms)','Fontsize',12)
    %title('Influence of IPSC timing & Location','Fontsize',20)
    axis([0 2*pi 0 8])
    
end