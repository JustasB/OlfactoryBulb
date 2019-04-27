function [out nil p] = statistics_motifs_pairs_facildepress(Dyn, C, A, Amax, h)
%
% out = statistics_motifs_pairs_facildepress(F, C, A, Amax, h)
%
% Note: clipping occur for values below (h * Amax)
%
% A is intended to be the W(:,:) matrix
% C is intended to be the CC(:,:) matrix (i.e. {0,1})
% Dyn(i,j) is intended to be '1' if the connection is facilitatory, '-1' if
% depressing




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


out = zeros(1,6);
nil = zeros(1,6);
p   = -1;

% out(1) counts lack of connections
% out(2) counts unidirectional connections
% out(3) counts bidirectional connections

A = A .* C;
%A = (A >= (Amax*h)) .* (2 * F - 1);
A = (A >= (Amax*h)) .* Dyn;

% now A(i,j) is '1' for facilitating, '-1' for depressing, '0' for no
% connection.


n = size(A,1);			   % n is the size of the matrix	
Nmotifs = (n * (n - 1)) / 2.;       % Number of all possible 'pair' motifs. 


% classes:
% 1: no connections
% 2: unidirectional facilitating
% 3: unidirectional depressing
% 4: bidirectional facilitating
% 5: bidrectional mixed
% 6: bidirectional depressing
%

for i=1:n,
    for j=i:n,
        if (i~=j)
            %index = A(i,j) + A(j,i) + 1;            
            [x y] = pattern(A(i,j), A(j,i));

            class1 = ((x+y) == 0);
            class2 = ((x+y) == 1);
            class3 = (x == -1);
            class4 = ((x+y) == 3);
            class5 = (y == -1);
            class6 = (x == -2);
            index = 1 * class1 + 2 * class2 + 3 * class3 + 4 * class4 + 5 * class5 + 6 * class6;
            out(index) = out(index) + 1; 
        end
    end
end


out = out / Nmotifs;  % Normalization

 p   = sum(abs(A(:)) == 1) / (n * (n - 1)); % probability of a connection
 % note: A has been already multiplied by C, so diagonal terms are already excluded.

 q   = sum(A(:) == 1) / (n * (n - 1)); % joint probability of a facilitatory connection AND that a connection exist
 q   = q / p;                          % conditional probability of a facilitatory connection (given a connection)[

 
 % nil = [(1-p)^2, 2 * p * (1-p), p*p];
nil(1) = (1-p)^2;
nil(2) = 2* (q * p ) * (1-p);
nil(3) = 2* ((1-q)*p) * (1-p);
nil(4) = (q*p)^2;
nil(5) = 2 *(q*p) * ((1-q)*p);
nil(6) = ((1-q)*p)^2;
% classes:
% 1: no connections
% 2: unidirectional facilitating
% 3: unidirectional depressing
% 4: bidirectional facilitating
% 5: bidrectional mixed
% 6: bidirectional depressing
%

% n = length(CC(:) > 0);   Which normalization!?

X = [1 2 3 4 5 6];
cla;
hold on;
% B = bar(X, (out - nil) * Nmotifs, 0.5);
% E = errorbar(X, [0 0 0], sqrt(Nmotifs * (1-nil) .* nil / 0.05));
B2 = bar(X+0.2, nil, 0.5);
set(B2, 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);
B1 = bar(X, out, 0.5);
set(B1, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
E = errorbar(X+0.3, nil, sqrt((1-nil) .* nil / (0.05 * Nmotifs)));
set(E(1), 'LineStyle', 'none', 'Color', [0.5 0 0]);
hold off;
legend([B1 B2 E(1)], 'Actual', 'Null h.', '95% conf');
%set(B, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
set(gca, 'XTick', [1 2 3 4 5 6], 'XTickLabel', {'X', 'f: ->', 'd: ->', 'f: <->', 'd,f: <->', 'd: <->'});
%set(gca, 'FontSize', 20);
ylabel('Probability');
YLIM = get(gca, 'YLim');
YLIM(1) = 0;
ylim(YLIM);
end







function [x, y] = pattern(a, b)  % Helper function for the motif discrimination
 x = a + b;
 y = a * b;
end