function figure4ef(ntrial)

% n_trial :     number of different noise seed.

doublim = 1;
nois_real = ntrial*7;
ggaba=[-0.1 -0.2 -0.5 -1.0 -2.0];

for g = 1:5
    gI=num2str(ggaba(g),'%2.1f');
    s=['fig4_',gI,'00000.dat'];
    pa  = dlmread(s);

    for i=1:7
        isi(i,:) = pa(i:7:nois_real,11)-pa(i:7:nois_real,10);
    end
    for i=1:7
        isisel{i} = isi(i,find(isi(i,:)>doublim+5*(i-1)));
    end

    % calcul de la phase bruitee
    phase = [0 5 10 15 20 25 30];
    for i=1:7
        delai(g,i) = mean(isisel{i});
        sigma(g,i) = std(isisel{i});
    end
end

figure(1);hold on
for i=1:7
    plot(-ggaba,delai(:,i),'ko-')
end
set(gca,'XScale','log')
%set(gca,'XTick',[0 1 2 3 4])
set(gca,'XTickLabel',[0.1 0.2 0.5 1 2])
xlabel('g_{GABA} (\mu S)','Fontsize',14)
ylabel('\mu_{ISI} (ms)','Fontsize',14)
axis([0 33 30 70])

figure(2);hold on
for i=1:7
    plot(-ggaba,sigma(:,i),'ko-')
end
set(gca,'XScale','log')
set(gca,'XTickLabel',[0.1 0.2 0.5 1 2])
xlabel('g_{GABA} (\mu S)','Fontsize',14)
ylabel('\sigma_{ISI} (ms)','Fontsize',14)
axis([0 33 8.5 12.5])

