function weight_distribution(W)
%
%
%
 N = size(W,1);
 W = W -1000 * eye(N); 
 

maxweight = 5;
delta     = maxweight/40.;
edges     = 0:delta:maxweight;
[y, x] = histc(W(:), edges);

y = y / (N*(N-1));
B = bar(edges, y, 'histc');
set(B, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);

xlabel('W');
ylabel('Fraction');
xlim([0 maxweight])
