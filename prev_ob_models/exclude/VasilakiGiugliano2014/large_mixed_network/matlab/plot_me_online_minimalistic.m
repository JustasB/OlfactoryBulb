clear all;
clc;

addpath matlab;
%addpath Gen_and_Analize_Connect;


figure(1); 
set(gcf, 'DoubleBuffer', 'on');

load connectivity.mat;

H = mysubplot(3,2,0.1,0.1,0.05);
iter = 1;
mW   = nan*zeros(100,1);
sW   = nan*zeros(100,1);
mmW  = nan*zeros(100,1);
 
sec_old = 0;

while(1)
 ppp = dir('G.x');
% ppp = dir('raster.x');
 sec = str2num(ppp.date(end-1:end));
if (sec~=sec_old)
 sec_old = sec;
 try
 
% % %  axes(H(1));
% % %  W = load('G.x'); weight_distribution(W);
% % %  %plot_rate; draw now;
 
 axes(H(2));
 plot_rate_histogram; %draw now;
 XLIM = get(gca, 'XLim'); YLIM = get(gca, 'YLim');
 hold on; P = plot((1:length(rates))/length(rates)*XLIM(2),rates/max(rates) * YLIM(2), 'r'); hold off;

 
 
 axes(H(3));
 hold on;
 mW(iter) = mean(W(:));
 sW(iter) = std(W(:));
 mmW(iter)= median(W(:));
 plot(1:100, mW, 1:100, sW, 1:100, mmW);
 iter = iter + 1;
 if iter>100, iter = 1; 
     mW   = nan*zeros(100,1);
     sW   = nan*zeros(100,1);
     mmW  = nan*zeros(100,1);
 end;
 hold off
 %plot_raster; draw now;
 
 axes(H(4))
 %imagesc(W>(2./3.*5.)); 
 imagesc(W.*abs(C)); colorbar; colormap(gray);
 

 axes(H(5)); cla
 [out nil p] = statistics_motifs_pairs_normal(abs(C), W, 5., 2./3.);
% %[out nil p] = statistics_motifs_pairs_normal(abs(C), abs(C), 1, 0.001);
 legend hide;
 
 axes(H(6)); cla
 [out nil p] = statistics_motifs_pairs_facildepress(C, abs(C), W, 5., 2./3.);
legend hide;
set(H, 'FontSize', 12);

%[s,p] = symmetry_measure3(W.*abs(C), 5., 2./3.)
maxWeight = 5.;
cropWeight = 2./3.;
N = size(W,1);
NM = floor(N/2.);
CCC = abs(C);
[symmeasure,chancelevel]=symmetry_measure3(W(1:NM,1:NM).*CCC(1:NM,1:NM),maxWeight,cropWeight);
display(sprintf('symmetry measure for facilitating subnetwork: %f two-tailed p-value: %f ',symmeasure,chancelevel));
[symmeasure,chancelevel]=symmetry_measure3(W(NM+1:N,NM+1:N).*CCC(NM+1:N,NM+1:N),maxWeight,cropWeight);
display(sprintf('symmetry measure for depressing subnetwork: %f two-tailed p-value: %f ',symmeasure,chancelevel));
display(sprintf('mean connections F->D:%f D->F:%f ',mean(mean(W(NM+1:N,1:NM).*CCC(NM+1:N,1:NM))),mean(mean(W(1:NM,NM+1:N).*CCC(1:NM,NM+1:N)))));



catch me
;
 end % try-catch
end % if
 pause(1);
end % while