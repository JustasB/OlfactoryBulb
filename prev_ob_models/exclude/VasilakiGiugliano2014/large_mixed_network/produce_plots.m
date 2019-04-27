%
% Let's plot the histograms
%
addpath matlab;
load connectivity.mat;

figure(1); clf;

W = load('W_100000.x'); 
[out nil p] = statistics_motifs_pairs_facildepress(C, abs(C), W, 5., 2./3.);
%legend hide;
set(gca, 'FontSize', 12);
print(gcf, '-depsc2', '-loose', 'mixed_population.eps');
print(gcf, '-dpng', 'mixed_population.png');


figure(2); clf
imagesc(W.*abs(C)); colorbar; colormap(gray);
print(gcf, '-depsc2', '-loose', 'mixed_population_W.eps');
print(gcf, '-dpng', 'mixed_population_W.png');


figure(3); clf
%plot_rate_histogram; 
 rates = load('rates.x');  
 delta = 2.;
 edges = 0:delta:120;
 [y1, x] = histc(rates(1:500), edges); 
 [y2, x] = histc(rates(501:1000), edges); 

  y1 = y1 / length(rates);
  y2 = y2 / length(rates);
  hold on;
  B1 = bar(edges, y1, 'histc');
  B2 = bar(edges, y2, 'histc');
  hold off;
  set(B1, 'FaceColor', [0 0 0] + 0.5, 'EdgeColor', [0 0 0] + 0.5);
  set(B2, 'FaceColor', [0 0 0] + 0, 'EdgeColor', [0 0 0] + 0.);
  
  xlabel('Firing Rate [Hz]');
  ylabel('Fraction of neurons');
  
%  xlim([0 max(rates)]);
xlim([0 120]);
%ylim([0 1]);
%  ylim([0 length(rates)/5.]);
set(gca, 'XTick', [0:20:120], 'Box', 'off');
print(gcf, '-depsc2', '-loose', 'mixed_population_rates.eps');
print(gcf, '-dpng', 'mixed_population_rates.png');




maxWeight = 5.;
cropWeight = 2./3.;
N = size(W,1);
NM = floor(N/2.);
CCC = abs(C);
[symmeasure,chancelevel]=symmetry_measure3(W(1:NM,1:NM).*CCC(1:NM,1:NM),maxWeight,cropWeight);
fp = fopen('mixed_population_symmetry.txt', 'w');

fprintf(fp, 'symmetry measure for facilitating subnetwork: %f\n',symmeasure);

[symmeasure,chancelevel]=symmetry_measure3(W(NM+1:N,NM+1:N).*CCC(NM+1:N,NM+1:N),maxWeight,cropWeight);
fprintf(fp, 'symmetry measure for depressing subnetwork: %f',symmeasure);
%display(sprintf('mean connections F->D:%f D->F:%f ',mean(mean(W(NM+1:N,1:NM).*CCC(NM+1:N,1:NM))),mean(mean(W(1:NM,NM+1:N).*CCC(1:NM,NM+1:N)))));
fclose(fp);

