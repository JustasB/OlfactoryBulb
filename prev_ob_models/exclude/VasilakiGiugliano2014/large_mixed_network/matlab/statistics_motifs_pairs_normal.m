function [out nil p] = statistics_motifs_pairs_normal(C, A, Amax, h)
%
% out = statistics_motifs_pairs_normal(A, Amax, h)
%
% Note: clipping occur for values below (h * Amax)
%
% A is intended to be the W(:,:) matrix
% C is intended to be the CC(:,:) matrix (i.e. {0,1})
%

%
% In the null-hypothesis that pair-motifs are occurring independently
% (i.e. there is no relationship between A(i,j) and A(j,i)), the following
% is expected:
%
% let's call p = prob{A(i,j) == 1}
%
% motif_1 (i.e. no connection) ---> (1-p) * (1-p)
% motif_2 (i.e. unidirectional) --> 2 * p * (1-p)
% motif_3 (i.e. bidirectional) ---> p * p
%


out = zeros(1,3);

% out(1) counts lack of connections
% out(2) counts unidirectional connections
% out(3) counts bidirectional connections

A = A .* C;
A = (A >= (Amax*h));


n = size(A,1);			   % n is the size of the matrix	
Nmotifs = (n * (n - 1)) / 2.;       % Number of all possible 'pair' motifs. 

for i=1:n,
    for j=i:n,
        if (i~=j)
            index = A(i,j) + A(j,i) + 1;
            out(index) = out(index) + 1; 
        end
    end
end


out = out / Nmotifs;  % Normalization

p   = sum(A(:) == 1) / (n * (n - 1)); % note: A has been already multiplied by C, so diagonal terms are already excluded.
nil = [(1-p)^2, 2 * p * (1-p), p*p];


% n = length(CC(:) > 0);   Which normalization!?

X = [1 2 3];
cla;
hold on;
% B = bar(X, (out - nil) * Nmotifs, 0.5);
% E = errorbar(X, [0 0 0], sqrt(Nmotifs * (1-nil) .* nil / 0.05));
B2 = bar(X+0.2, nil, 0.5);
set(B2, 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
B1 = bar(X, out, 0.5);
set(B1, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
E = errorbar(X+0.2, nil, sqrt((1-nil) .* nil / (0.05 * Nmotifs)));
set(E(1), 'LineStyle', 'none', 'Color', [0.5 0 0]);
hold off;
legend([B1 B2 E(1)], 'Actual', 'Null h.', '95% conf');

%set(B, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'X', '->', '<->'});
%set(gca, 'FontSize', 20, 'Box', 'on');
ylabel('Probability');
YLIM = get(gca, 'YLim');
YLIM(1) = 0;
ylim(YLIM);
end