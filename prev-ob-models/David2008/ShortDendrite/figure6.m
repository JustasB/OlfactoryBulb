function figure6(n_trial)

doublim = 2;
nois_real = 10*n_trial;
pa = dlmread('fig6_200_0.800000.dat');
for i=1:10
    isi{i} = pa(i:10:nois_real,11)-pa(i:10:nois_real,10);
end
isisel{1} = isi{1}(find((isi{1}>doublim) & (isi{i}<100)));
for i=2:10
    isisel{i} = isi{i}(find((isi{i}>doublim+10) & (isi{i}<100)));
end

for i=1:10
    % creation d'une figure automatiquement
    % fit_ML_normal(isisel{i},[0:250]);
    % hold on
    figure;hold on
    a = histc(isisel{i},[0:250],1);
    size(isisel{i});
    a = a/sum(a);
    bar([0:250],a,1,'stack','FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7]);
    axis([0 95 0 0.06])
    set(gca,'YTick',[0,0.03,0.06])
    set(gca,'FontSize',12)
    xlabel('time (ms)','FontSize',12)
    
    [muhat, sigmahat] = normfit(isisel{i});
    x=[0:1:100];
    y=(1/(sigmahat*sqrt(2*pi))) * exp(-((x-muhat).^2)/(2*sigmahat^2));
    plot(x,y,'k','LineWidth',2);
    txt     = sprintf( '\\fontsize{%d}\\bf\\mu = %4.1f \n\\sigma = %4.1f',12,muhat,sigmahat);
    line( muhat + sigmahat + [0 4] , (max(y)/2)*[1 1],'linewidth',3,'color',[0 0 0] );
    text(muhat + sigmahat + 5 , (max(y)/2) , txt);
end

figure; hold on
for i=[3 5 7 9]
    a = histc(isisel{i},[0:250],1);
    a = a/sum(a);
    bar([0:250],a,1,'stack','FaceColor',i*[0.08,0.08,0.08],'EdgeColor','none');
    axis([0 95 0 0.06])
    set(gca,'YTick',[0,0.03,0.06])
    set(gca,'FontSize',12)
    xlabel('time (ms)','FontSize',12)
end

figure
for st = 1:3
    if st == 1
        clear('pa')
        pa  = dlmread('fig6_200_0.500000.dat');
    elseif st == 2
        clear('pa')
        pa  = dlmread('fig6_200_0.800000.dat');
    elseif st == 3
        clear('pa')
        pa  = dlmread('fig6_200_1.800000.dat');
    end
    
    for i=1:10
        isi{i} = pa(i:10:nois_real,11)-pa(i:10:nois_real,10);
    end
    
    isisel{1} = isi{1}(find((isi{1}>doublim) & (isi{i}<100)));
    for i=2:10
        isisel{i} = isi{i}(find((isi{i}>doublim+10) & (isi{i}<100)));
    end

    % calcul de la phase bruitee
    sigi = [0 2.5 5 7.5 10 12.5 15 17.5 20];
    for i=1:10
        delai(i) = mean(isisel{i});
        sigma(i) = std(isisel{i});
    end

    figure(11);hold on
    plot([sigi(1:9)],[delai(2:10)],'ko-')
    set(gca,'XTick',sigi)
    set(gca,'XTickLabel',{'0','2.5','5','7.5','10','12.5','15','17.5','20'})
    xlabel('\sigma_{IPSC} (ms)')
    ylabel('\mu_{ISI} (ms)','Fontsize',14)
    axis([0 20 30 65])

    figure(12);hold on
    plot([sigi(1:9)],[sigma(2:10)],'ko-')
    set(gca,'XTick',sigi)
    set(gca,'XTickLabel',{'0','2.5','5','7.5','10','12.5','15','17.5','20'})
    xlabel('\sigma_{IPSC} (ms)')
    ylabel('\sigma_{ISI} (ms)','Fontsize',14)
    axis([0 20 7 18])
end




