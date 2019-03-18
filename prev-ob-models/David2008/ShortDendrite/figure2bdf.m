function figure2bdf

% n_del = 41 pour 40Hz
% n_del = 41 pour 30Hz
% n_del = 61 pour 20 Hz
% inhib_ipsc_del_dx(fich_pa,n_del,n_fig)
% traite les donnees du programme gGABAainfluence
% fich_pa est le fichier des temps de spike des differents essais mis bout
% a bout en colonne.

nmit = 10;
n_del = 60;
n_fig = 10;
for st=1:3
    if st==1
        clear('pa')
        pa = dlmread('fig2bdf_0.500000.dat');
    elseif st==2
        clear('pa')
        pa = dlmread('fig2bdf_0.800000.dat');
    elseif st==3
        clear('pa')
        pa = dlmread('fig2bdf_1.800000.dat');
    end
    % creation de la matrice de dimension (n_dx, n_gGABA, n_del)
    for k=1:n_del
        for i=1:nmit
            isi(k,i) = pa((k-1)*nmit+i,11)-pa((k-1)*nmit+i,10);
            isi(1,i) = pa(i,11)-pa(i,10);
            if isi(k,i)<15
               isi(k,i) = pa((k-1)*nmit+i,12)-pa((k-1)*nmit+i,11);
               isi(1,i) = pa(i,12)-pa(i,11);
               if isi(k,i)<15
                   isi(k,i) = pa((k-1)*nmit+i,13)-pa((k-1)*nmit+i,12);
                   isi(1,i) = pa(i,13)-pa(i,12);
               end
            end
            retard(k,i) = isi(k,i)-isi(1,i);
            del(k,i)= k-1;
        end
    end
    retard;
    del;
    figure(n_fig+st-1); hold on
    for i=1:nmit
        plot(del(3:end,i)*2*pi/(isi(1,i)),retard(3:end,i),'k','Linewidth',1)
    end

    xlabel('IPSC phase','Fontsize',12)
    set(gca,'Fontsize',12)
    set(gca,'XTick',[0,pi,2*pi])
    set(gca,'XTickLabel',{'0','pi','2pi'})
    ylabel('\Deltat (ms)','Fontsize',12)
    %title('Influence of IPSC timing','Fontsize',20)
    axis([0 2*pi 0 8])

end