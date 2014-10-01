% folder containing codes this script uses
addpath('code');

% folder with input data (.mat files)
matdir = 'data';

% parameters for similarity matrix computation
alpha = 0.80; 
iters = 20;

% similarity matrix computation method
similarity_m = 'mat3'; 

% matching method
match_m = 'greedy'; 

% species to align
name1 = 'ecoli'; % E. coli
name2 = 'scere'; % S. cerevisiae

% filepaths for input data
mat_path1 = fullfile(matdir, [name1, '.mat']);
mat_path2 = fullfile(matdir, [name2, '.mat']);
h_path = fullfile(matdir, sprintf('%s-%s.mat', name1, name2));

% load files with PPI networks and initial preferences
mdata1 = load(mat_path1);
mdata2 = load(mat_path2);
hdata = load(h_path);

% setup matrices
B = mdata1.A; B = max(B, B');
A = mdata2.A; A = max(A, A');
H = hdata.H;

% similarity matrix computation
[X, dt_similarity] = MAT3_rank(A, B, alpha, iters, H);
         
% matching computation
[M, dt_greedy_match] = greedy_match(X);
[mb, ma] = find(M);

% compute the permutation for sorting mb
[mb_s, perm_s] =  sort(mb); 
ma_perm = ma(perm_s); % ma permuted

% matching names (matching protein identifiers)
names1 = mdata1.names;
names2 = mdata2.names;
matching_names = cell(length(ma_perm), 2);
for q = 1:length(ma_perm)
    ii = mb_s(q);
    jj = ma_perm(q);
    matching_names{q}{1} = names1{ii};
    matching_names{q}{2} = names2{ii};
end

% total computation time (similarity matrix comutation + greedy matching)            
dt_total = dt_similarity + dt_greedy_match;

% compute the alignment graph
[aM, aG, dt] = align(A, B, M);

% number of overlapped edges
edges = nnz(aG)/2;
          
% optionally collect some of the results in a structure
results = [];
results(end+1).alignment = aG;
results(end).match = aM;
results(end).edges = edges;
results(end).dt_similarity = dt_similarity;
results(end).dt_greedy_match = dt_greedy_match;
results(end).dt_total = dt_total;

% Output computation statistics
fprintf(1, '\n\n');
fprintf(1, 'Aligned yeast and bacterium (mat3_greedy)\n');
fprintf(1, '=========================================\n\n');
fprintf(1, 'Construction of the similarity matrix (mat3) : %5.2f secs\n', dt_similarity);  
fprintf(1, 'Matching using greedy mathod (greedy): %5.2f secs\n', dt_greedy_match);
fprintf(1, 'Generated alignment graph consists of %d edges\n', edges);



 







