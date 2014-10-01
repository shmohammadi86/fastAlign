This directory contains the source code of the programs used. 

example/
========
A directory with an example of a mat3_greedy run; just run example.m.
Data files and scripts needed are included in data/ and code/ respectively.

Note that the bulk of our experiments additionally parses iso_greedy, iso_hungarian and mat3_auction outputs, saves results in files residing in different folders of the provided datasets archive and iterates over all network pairs; See matlab/bio.m for such an example. 



matlab/
=======

MATLAB files used in the numerical computations in this paper. Additional codes from the netalign project http://www.cs.purdue.edu/homes/dgleich/codes/netalign/ could be included in the search path. For IsoRank computations, native binaries were used. Also look into their documentation part.

MAT3_rank.m
Computes the similarity matrix given the adjacency matrices of two input networks, alpha, number of iterations and the preferences matrix.

greedy_match.m
Computes the matching matrix M (M_ij = 1 iff node i matches with node j) given the similarity matrix of a pair of input networks.

align.m
Also computes the alignment graph of two networks given their adjacency matrices, and the matching.

bio_components.m
Computes and outputs (strongly connected) components in the alignment graph of a pair of input networks.

tab_read.m, dmat_write.m
Respectively for reading the (original) tab-delimited PPI files and outputing dense matrices.

bio.m
The script parsing input and result files in a local file organization mirroring those of actual experiment runs to compute, gather output from different stages (e.g. auction matching) and approaches (e.g. IsoRank binary) and finally produce "uniform" output. It uses all the above routine files.

count_overlap.m, csum.m, normout.m, scomponents.m, sparse_to_csr.m
Files from netalign project (dependencies provided for convenience).
   


auction/
========= 

The source code of the (distributed) auction code including a small example how to use it. The algorithm is detailed in
http://www.sciencedirect.com/science/article/pii/S0167819112000750

