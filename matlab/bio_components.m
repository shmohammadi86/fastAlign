function [genes, dt] = bio_components(A, B, M, namesA, namesB, outpath)
% bio_components Computes the connected components of the alignment graph
% of two PPI networks
% 
% Input arguments:
% - A, B: Input PPI networks.
% - M: Matching matrix
% - namesA, namesB: Protein identifiers for the input PPI networks
% - outpath: Path of the file where components are output in a format
%   suitable for further biological processing 
% Output arguments:
% - genes: records with the components in the alignment graph
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

G = sparse(ii, jj, 1.0, num_matches, num_matches); % alignment graph

tG = triu(G);
[ci, sizes] = scomponents(G);
genes = [];
for i=1:length(sizes)
    num = sizes(i);
    if num > 1
        nodesB = find(ci==i);
        mynamesB = {};
        mynamesA = {};
        nodesA = zeros(length(nodesB), 1);
        for j=1:length(nodesB)
            nodesA(j) = p_ma(nodesB(j));
            mynamesB{j} = namesB(nodesB(j));
            mynamesA{j} = namesA(nodesA(j));
        end
        genes(end+1).nodes1 = nodesB;
        genes(end).nodes2 = nodesA;
        genes(end).names1 = mynamesB;
        genes(end).names2 = mynamesA;
        genes(end).size = length(nodesB);
    end
end

f = fopen(outpath, 'w');
fprintf(f, '%d\n', length(genes));
for i=1:length(genes)
    fprintf(f, '%d\n', genes(i).size);
    for j=1:genes(i).size
        fprintf(f, '%s\t', char(genes(i).names1{j}));
    end
    fprintf(f, '\n');
    for j=1:genes(i).size
        fprintf(f, '%s\t', char(genes(i).names2{j}));
    end
    fprintf(f, '\n');
  
end

fclose(f);

dt = etime(clock, t0);

end

