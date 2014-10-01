function [aM, aG, dt] = align(A, B, M)
% align Computes the alignment graph of two input networks based on a node
% to node matching
%
% Input arguments:
% - A, B: Input PPI networks.
% - M: Matching matrix
% Output arguments:
% - aM: Matching matrix with the updated numbering.
% - aG: The alignment graph as a sparse matrix.
% - dt: the time in seconds for the operation.


t0 = clock;

% ensure adjacency, symmetric
A(A > 0) = 1.0;
A = max(A, A');

B(B > 0) = 1.0;
B = max(B, B');

% we assume m <=n 
n = size(A, 1);
m = size(B, 1);
[mb, ma] = find(M);


[s_mb, ind] = sort(mb); % sort matching nodes of smaller graph (B) 
num_matches = length(s_mb);
p_ma = ma(ind); % permute matching nodes of larger graph (A) accordingly
C = A(p_ma, p_ma);
D = B(s_mb, s_mb) + C;
p = spdiags(D, 0);
D = D - spdiags(p, 0, num_matches, num_matches);
[ii, jj] = find(D == 2); % B and the permuted A have elements at these coords

aG = sparse(ii, jj, 1.0, num_matches, num_matches); % alignment graph


aM = sparse(s_mb, p_ma, 1.0, m, n); % matching with the new numbering


dt = etime(clock, t0);

end

