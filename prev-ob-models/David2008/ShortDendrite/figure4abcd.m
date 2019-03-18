function figure4abcd(n_trial,gI)

% n_trial :  number of different noise seed.
% gI      :  gabaergic conductance

doublim = 2; % 5
nois_real = n_trial*7;

s=['fig4_-',num2str(gI,'%2.1f'),'00000.dat'];
pa  = dlmread(s);

for i=1:7
    isi(i,:) = pa(i:7:nois_real,11)-pa(i:7:nois_real,10);
end

for i=1:7
    isisel{i} = isi(i,find(isi(i,:)>doublim+5*(i-1)));
end

% creation d'une figure automatiquement
for i=1:7
    figure
    hold on
    a1 = histc(isisel{i},[0:250],1);
    a1 = a1/sum(a1);
    bar([0:250],a1,1,'stack','FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7]);
    axis([0 95 0 0.06])
    set(gca,'YTick',[0,0.03,0.06])
    xlabel('ISI (ms)')
    set(gca,'FontSize',12)
    
    [muhat, sigmahat] = normfit(isisel{i});
    x=[0:1:100];
    y=(1/(sigmahat*sqrt(2*pi))) * exp(-((x-muhat).^2)/(2*sigmahat^2));
    plot(x,y,'k','LineWidth',2);
    txt     = sprintf( '\\fontsize{%d}\\bf\\mu = %4.1f \n\\sigma = %4.1f',12,muhat,sigmahat);
    line( muhat + sigmahat + [0 4] , (max(y)/2)*[1 1],'linewidth',3,'color',[0 0 0] );
    text(muhat + sigmahat + 5 , (max(y)/2) , txt);
end

