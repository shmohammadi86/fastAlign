function [X, dt] = MAT3_rank(A, B, alpha, iters, H)
% MAT3_rank Computes the similarity of two undirected graphs using 
% IsoRank algorithm with its kernel iteration step implemented 
% as a triple matrix product.
% 
% Input arguments:
% - A, B: the adjacency matrices of the two graphs.
% - alpha: alpha parameter of the algorithm.
% - iters: the number of iterations.
% - H: preferences matrix.
%
% 
% Output arguments:
% - X: the matrix with the similarity scores (similarity matrix).
%     Note that element X(i,j) is the similarity score of node i in B 
%     and node j in A; if B has m nodes and A has n nodes then X is an 
%     m x n matrix.
% - dt: the time in seconds for the operation.
% 
% H is a matrix of m x n nonnegative elements (summing up to 1). It
% contains preknown similarity scores between nodes in the two graphs
% (e.g. sequence similarity scores i.e. sequence data, as opposed to 
% network data already encoded in A, B matrices). Also referred as a 
% preferences matrix. 
%
% If H is not given then it uses a matrix of identical elements 
% ("uniform" initial conditions).
%
% The code is based on full_isorank.m from netalign codes:
% http://www.stanford.edu/~dgleich/publications/2009/netalign/
% (developed by David Gleich).
%
% References:
% Paper where netalign codes were extensively used:
% @inproceedings{bayati_algorithms_2009,
% 	title = {Algorithms for Large, Sparse Network Alignment Problems},
% 	booktitle = {2009 Ninth {IEEE} International Conference on Data Mining},
% 	publisher = {{IEEE}},
% 	author = {M. Bayati and M. Gerritsen and D. F. Gleich and A. Saberi and Y. Wang},
% 	year = {2009},
% 	pages = {705--710}
% }
% 
% Paper about the IsoRank algorithm:
% @article{singh_global_2008,
% 	title = {Global alignment of multiple protein interaction networks with application to functional orthology detection},
% 	volume = {105},
% 	number = {35},
% 	journal = {Proceedings of the National Academy of Sciences},
% 	author = {R. Singh and J. Xu and B. Berger},
% 	year = {2008},
% 	pages = {12763}
% }



B = max(B, B');
P = normout(B);

A = max(A, A');
Q = normout(A);

m = size(P, 1);
n = size(Q, 1);
N = m * n;

if (nargin == 4)
    H = ones(m,n)./N;
end

v = H(:);
v = v ./ csum(v);

t0 = clock;
x = zeros(N,1); 
x = x + v;

iter = 0;
while iter < iters 
    y = x ./ csum(x); 
    X = reshape(x,m,n); 
    Y = reshape(y,m,n);
    for i = 1:n
        X(:,i) = alpha * (P' * (Y * Q(:,i)));
    end
    
    x = reshape(X,N,1);
    y = reshape(Y,N,1);
    
    gamma = csum(y) - csum(x);
    x = x + gamma * v;
    
    iter = iter + 1; 
end

X = reshape(x,m,n);

dt = etime(clock, t0);

end

