%
%
%

P = dir('W_*.x');
N = length(P);
load connectivity.mat;

time = zeros(length(P),1);
for k=1:length(P),
    fname = P(k).name;
    tmp1  = find(fname=='_') + 1;
    tmp2  = find(fname=='.') - 1;
   
    time(k) = str2num(fname(tmp1:tmp2));
end
[y, I] = sort(time);
P = P(I);

for i=1:N,
 fname = P(i).name
 W     = load(fname);
 
 subplot(N,2,2*i-1);
%subplot(N,1,i);
cla;
 [out nil p] = statistics_motifs_pairs_normal(abs(C), W, 5., 2./3.);
 legend off;
% YLim([0 0.05]);
 %imagesc(W);
 
 subplot(N,2,2*i);
 [out nil p] = statistics_motifs_pairs_facildepress(C, abs(C), W, 5., 2./3.);
  legend off;
% YLim([0 0.05]);

end



% 
% 
% H = [];
% whichones = [1 N];
% for k=1:2
% T1 = title(fname);    
% hh = subplot(2,2,2*k-1); H = [H, hh];
% i =whichones(k);
% fname = P(i).name
% W     = load(fname);
% [out nil p] = statistics_motifs_pairs_normal(abs(C), W, 5., 2./3.);
% legend off;
% %YLim([0 0.05]);
%  
% T2 = title(fname);    
% hh = subplot(2,2,2*k); H = [H, hh];
% [out nil p] = statistics_motifs_pairs_facildepress(C, abs(C), W, 5., 2./3.);
% legend off;
% %YLim([0 0.05]);
% set([T1 T2], 'Interpreter', 'none')
% set(H, 'FontSize', 10);
% end

