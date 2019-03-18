%function StructuredMatrix(N, Pdd, Pff, Pfd, Pdf)
function StructuredMatrix(N, Pdd, Pff, Pfd, Pdf)

% The function constructs a connectivity-like type Matrix F of dimentions NxN, 
% which describes the type of synapses between neurons i and j (facilitating
% or depressing). F[i,j]=1 if connection from j to i is facilitating and 0 
% otherwise. The specific Matrix has spatial structure in the following sense.
% Neurons with index from 1 to floor(N/2) have only facilitating synapses with 
% all other neurons and neurons with index from floor(N/2)+1 to N have
% depressing only synapses. 
%
% In addition, a connectivity-like binary Matrix CC is produced, with CC(i,j)=1 if
% connection from j to i is present and 0 if absent.
% 
% The two matrices are finally combined to Matric FC (functional connectivity), where 
% FC(i,j)=1 if the connection from j to i is present and facilitating, FC(i,j)=-1 if
% the connection is present and depressing and FC(i,j)=0 if no connection is
% present.

F  = zeros(N,N);
CC = zeros(N,N);
M  = floor(N/2); %(~middle)


F(:,1:M) = 1;

CC(1:M,1:M)     = (rand(M,M) < Pff);

CC(M+1:N,M+1:N) = (rand(N-M,N-M) < Pdd);

CC(1:M,M+1:N)   = (rand(M,N-M) < Pfd);

CC(M+1:N,1:M)   = (rand(N-M,M)<Pdf);

FC=CC.*(2*F-1);


MM    = length(find(FC(:)~=0));  % Number of non-zero elements of FC(i,j)
Ctemp = FC';
Ctemp = Ctemp(:);   % A single vector, encoding 'C' row-wise (i.e. lexicographic order)
tmp   = find(Ctemp ~= 0);

fname = 'connectivity.dat';
fp    = fopen(fname, 'wb');
fwrite(fp, N, 'double');  
fwrite(fp, MM, 'double');  
for kk=1:length(tmp)
    fwrite(fp, tmp(kk) * sign(Ctemp(tmp(kk))), 'double');
end
fclose(fp);

C = FC;
save('connectivity.mat', 'C');
