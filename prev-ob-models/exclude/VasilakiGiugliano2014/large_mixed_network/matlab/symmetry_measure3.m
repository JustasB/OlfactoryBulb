function [s,p] = symmetry_measure3(A, Amax, h)
%
% [s,p] = symmetry_measure3(A, Amax, h)
%
% Evaluates the symmetry of an input (square) matrix 
% s = symmetry_measure(A)
%
% Note: clipping occur for values below (h * Amax)
%
% s is a number between 0 and 1
%
% p is the p-value for the Null Hypothesis: s comes from the distribution of 
% symmetry indeces of randomly created connectivity matrices (randomly i.e. 
% connections uniformly distributed between 0 and Amax and clipped to 0 if below h). 
% It is equal to the probablility that a randomly created matrix 
% has a symmetry measure  abs(s) and above and -abs(s) and below (two tailed
% p-value)

if (size(A,1) ~= size(A,2)), % Test whether or not A is square
    disp('Error: the measure is defined only for square matrices!');
    s = -1;	
    p = -1;
    return;
end;

n = size(A,1);			   % n is the size of the matrix	

A = A .* (1 - eye(size(A,1)));   % Elements on the diagonal are set to ?0?
 
tempclip = (A < (h* Amax));	   % Matrix of ?0?-?1? elements..	
A        = A .* (1 - tempclip);  % Hard-clipping for the upper bound

% The ?clipped? matrix [A*] is now composed of two categories of elements:
% those equal to ?0? and those in the range [h* Amax; Amax].

A = A ./ Amax; 			   % Matrix elements are normalized to 1	

% The normalized ?clipped? matrix [A*] is now composed of two categories of 
% elements: those equal to ?0? and those in the range [h ; 1]

uA = triu(A);		% This extracts the upper triangular part of A*
lA = tril(A)';		% This extracts the lower triangular part of A*

% Note: transposing tril(A) makes both uA and lA upper triangular matrices
% [uA] therefore contains A*(i,j), with i<=j, zeros otherwise
% [lA] therefore contains A*(j,i), with i<=j, zeros otherwise

X  = uA(:);			% Converts the matrix uA(i,j) into a vector
Y  = lA(:);			% Converts the matrix ul(i,j) into a vector

% Note: conversion occurs column-wise. Because of the transposition, there is
% the desired correspondence between uAi,j = A*(i,j) and lA(i,j) = A*(j,i), 
% reflected by definition for the corresponding indexes of X and Y. 
 
temp = (X + Y);	     % Because elements of A* were non-negative, temp == 0
			     % if and only if A*(i,j) = A*(j,i) = 0.

% Counting how many times a ?0? occurs in the vector ?temp? identifies (1) 
% all the elements on the diagonal (i.e. that have been zeroed earlier, and that 
% appear once at corresponding positions within X and Y); (2) all those cases in 
% which the elements A*(i,j) = A*(j,i) = 0 because of clipping to zero.
% temp contains many more ?0?, due to (3) the tridiagonal structure of the 
% matrices uA and lA. 
% If the size of uA and lA is n, then the zeros off-diagonal are: (n*n - n)/2 
% (i.e. all the element minus the elements on the diagonal, divided by 2 as we
% count only the lower or upper part of the matrix).
%
% Therefore, since we are only interested in (2) we must account for all the
% uninteresting zeros. In total, the number of trivial ?zeros? in temp are
% (1) n, (2) M, (3) (n*n-n)/2;	note: (1) + (3) = n*(n+1)*0.5.
 
M    = length(find(temp == 0.)) - n*(n+1)*0.5;
% M is the number of (mutually unconnected) pairs that are zero-zero, due to 
% clipping.

% We are now ready to calculate the average absolute difference between the 
% symmetric elements of [A*], |A*(i,j)-A*(j,i)|. If we compute it as 
% ?sum(abs(X-Y))?, we must then normalize it by the number of elements we make 
% our sum on (i.e. n*n, from which we must remove the (mutually unconnected)
% pairs 0-0 that would otherwise bias the measure towards higher symmetry.

K    = n * (n-1) * 0.5 - M;
 
if K>0 % This may not always be the case
    s = 1 - sum(abs(X-Y)) / K;
    
    tmp         =  (1-h) * ((1./3.)* (1-h)^2 + h*(1+h) ) / (1-h^2);
    theory_mean = 1. - tmp;
    theory_var  = ( (1./6.)*(1-h)^4 + (2./3.) * h * (1-h^3) ) / (1-h^2) - tmp^2;
    theory_std   = sqrt((1./(0.5*n*(n-1)*(1-h^2)))*(1 + h^2 / (0.5*n*(n-1)*(1-h^2))) * theory_var);
    
        
    %We calculate the integral between s and +inf and s and -inf
    %erf(x/sigma /sqrt(2)) gives the area of a zero mean gaussian pdf between -x and x (x>0). 
    %The area between -inf and +inf is 1. Subtracting the two gives the area between 
    %x and +inf and -x and -inf (two-tailed p-value)
    
    p=1-erf((abs(s-theory_mean)/theory_std)/sqrt(2));
   
    
else
    s = -1;
    p = -1;
end

end


