%% Setup
% Change accordingly
parent = 'experiments';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir, parent)==0
    error('experiment:wrongDir',...
        'experiments must be run from the correct directory, should be %s', ...
        parent);
end

% Change accordingly (or use the provided scripts only)
% our scripts 
addpath('../src/matlab/nsd');

% scripts from netalign project
addpath('../src/matlab/netalign');


%% Initialize
% Change accordingly
workdir = '/home/giorgos/Desktop/bio';

tabdir = fullfile(workdir, 'tab');
matdir = fullfile(workdir, 'mat');
csrdir = fullfile(workdir, 'csr');
evalsdir = fullfile(workdir, 'evals');
matchdir = fullfile(workdir, 'match');
namematchdir = fullfile(workdir, 'namematch');
similaritydir = fullfile(workdir, 'similarity');
resultsdir = fullfile(workdir, 'results');
timedir = fullfile(workdir, 'time');
isorunsdir = fullfile(workdir, 'isoruns');

% matches
m_p = '%s-%s-%s-%03d-%s-match.txt';
mt_p ='%s-%s-%s-%03d-%s-match-time.txt';

% name matches
nm_p = '%s-%s-%s-%03d-%s-namematch.txt';

% similarity
s_p = '%s-%s-%s-%03d-similarity.txt';
st_p ='%s-%s-%s-%03d-similarity-time.txt';

e_p = '%s-%s-%s-%03d-%s-edges.txt'; % number of overlapped edges

r_p = '%s-%s-%s-%03d-%s-results.mat'; % results in matlab

% total time
tt_p = '%s-%s-%s-%03d-%s-total-time.txt';


%% fill up the species structure
species = [];

tab_files = dir([tabdir, '/*.tab']);
for i = 1:length(tab_files)
    fname = tab_files(i).name;
    [pathstr, name, ext, versn] = fileparts(fname);
    tab_path = fullfile(tabdir, [name, ext]);
    [A, names, index, dt] = tab_read(tab_path);
    species(end+1).name = name;
    species(end).nodes = size(A, 1);
end


%% computations
alpha = 0.80; ialpha = int32(alpha * 100);
iters = 20;
similarity_m = 'mat3'; % similarity computation method
match_m = 'greedy'; % matching method


for i = 1:length(species)
    nodes1 = species(i).nodes;
    name1 = species(i).name;
    for k =1:length(species)
        nodes2 = species(k).nodes;
        name2 = species(k).name;
        if nodes1 < nodes2
            mat_path1 = fullfile(matdir, [name1, '.mat']);
            mat_path2 = fullfile(matdir, [name2, '.mat']);
            h_path = fullfile(matdir, sprintf('%s-%s.mat', name1, name2));
            mdata1 = load(mat_path1);
            mdata2 = load(mat_path2);
            hdata = load(h_path);
            B = mdata1.A; B = max(B, B');
            A = mdata2.A; A = max(A, A');
            H = hdata.H;
            
            % similarity matrix computation
            [S, dt_similarity] = MAT3_rank(A, B, alpha, iters, H);
            
            % similarity matrix path
            s_f = sprintf(s_p, name1, name2, similarity_m, ialpha);
            s_path = fullfile(similaritydir, s_f);
            dmat_write(S, s_path);
            
            % similarity matrix time path
            st_f = sprintf(st_p, name1, name2, similarity_m, ialpha);
            st_path = fullfile(timedir, st_f);
            f = fopen(st_path, 'w'); 
            fprintf(f, '%f\n', dt_similarity); 
            fclose(f);
            
            % matching computation
            [M, dt_greedy_match] = greedy_match(S);
            [mb, ma] = find(M);
            % also return the permutation for sorting mb
            % mb_s = mb(p_s)
            [mb_s, perm_s] =  sort(mb); 
            ma_perm = ma(perm_s); % ma permuted, this is what we want (0-based)
            
                        
            % match path
            m_f = sprintf(m_p, name1, name2, similarity_m, ialpha, match_m);
            m_path = fullfile(matchdir, m_f);
            f = fopen(m_path, 'w');
            for q=1:length(ma_perm)
                fprintf(f, '%d\n', ma_perm(q)-1);
            end
            fclose(f);
            
            % name match path
            nm_f = sprintf(nm_p, name1, name2, similarity_m, ialpha, match_m);
            nm_path = fullfile(namematchdir, nm_f);
            f = fopen(nm_path, 'w');
            names1 = mdata1.names;
            names2 = mdata2.names;
            for q = 1:length(ma_perm)
                ii = mb_s(q);
                jj = ma_perm(q);
                fprintf(f, '%s\t%s\n', names1{ii}, names2{jj});
            end
            fclose(f);
                       
            % matching computation time path
            mt_f = sprintf(mt_p, name1, name2, similarity_m, ialpha, match_m);
            mt_path = fullfile(timedir, mt_f);
            f = fopen(mt_path, 'w'); 
            fprintf(f, '%f\n', dt_greedy_match); 
            fclose(f);
            
            
            dt_total = dt_similarity + dt_greedy_match;
            
            % total computation time path
            tt_f = sprintf(tt_p, name1, name2, similarity_m, ialpha, match_m);
            tt_path = fullfile(timedir, tt_f);
            f = fopen(tt_path, 'w'); 
            fprintf(f, '%f\n', dt_total); 
            fclose(f);            
                       
            
            [aM, aG, dt] = align(A, B, M);
            edges = nnz(aG)/2;
            
            % overlapped edges
            e_f = sprintf(e_p, name1, name2, similarity_m, ialpha, match_m);
            e_path = fullfile(resultsdir, e_f);
            f = fopen(e_path, 'w'); 
            fprintf(f, '%d\n', edges); 
            fclose(f);                 
            
            results = [];
            results(end+1).alignment = aG;
            results(end).match = aM;
            results(end).edges = edges;
            results(end).dt_similarity = dt_similarity;
            results(end).dt_greedy_match = dt_greedy_match;
            results(end).dt_total = dt_total;
                       
            r_f = sprintf(r_p, name1, name2, similarity_m, ialpha, match_m);
            r_path = fullfile(resultsdir, r_f);
            save(r_path, 'results');
        end
    end
end


%%
f = fullfile(matdir, 'species-nodes.mat');
nodes = containers.Map;
for i = 1:length(species)
    name = species(i).name;
    nodes(species(i).name) = species(i).nodes;
end

save('species-nodes.mat', 'species', 'nodes');





%% just swap where necessary name columns

alpha = 0.80; ialpha = alpha * 100;
iters = 20;
similarity_m = 'isorank';

ig_m = 'isogreedy';
ih_m = 'isohungarian';

ig_p = 't%03d_%s-%s-final-hsp.txt';
ih_p = 't%03d_%s-%s-final-match.txt';

% isogreedy
ig_files = dir(fullfile(isorunsdir, '*hsp*'));

% isohungarian
ih_files = dir(fullfile(isorunsdir, '*match*'));


% isogreedy case
for i = 1:length(ig_files)
    name = ig_files(i).name;
    ig_path = fullfile(isorunsdir, name);
    parts = regexp(name, '_', 'split');
    name = parts{2};
    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    if nodes(name1) > nodes(name2)
        small_name = name2;
        large_name = name1;
        reversed = 1;
    else
        small_name = name1;
        large_name = name2;
        reversed = 0;
    end
    
    f = fopen(ig_path, 'r');
    triplets = textscan(f, '%s %s %f');
    fclose(f);
    num = length(triplets{1});
    nmig_f = sprintf(nm_p, small_name, large_name, similarity_m, ialpha, ig_m);
    nmig_path = fullfile(namematchdir, nmig_f);
    f = fopen(nmig_path, 'w');
    if reversed == 1
        for k=1:num
            fprintf(f, '%s\t%s\n', triplets{2}{k}, triplets{1}{k});
        end
    else
        for k=1:num
            fprintf(f, '%s\t%s\n', triplets{1}{k}, triplets{2}{k});
        end
    end
    fclose(f);
end
    

% isohungarian case
for i = 1:length(ih_files)
    name = ih_files(i).name;
    ih_path = fullfile(isorunsdir, name);
    parts = regexp(name, '_', 'split');
    name = parts{2};
    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    if nodes(name1) > nodes(name2)
        small_name = name2;
        large_name = name1;
        reversed = 1;
    else
        small_name = name1;
        large_name = name2;
        reversed = 0;
    end
    
    f = fopen(ih_path, 'r');
    triplets = textscan(f, '%s %s %f');
    fclose(f);
    num = length(triplets{1});
    nmih_f = sprintf(nm_p, small_name, large_name, similarity_m, ialpha, ih_m);
    nmih_path = fullfile(namematchdir, nmih_f);
    f = fopen(nmih_path, 'w');
    if reversed == 1
        for k=1:num
            fprintf(f, '%s\t%s\n', triplets{2}{k}, triplets{1}{k});
        end
    else
        for k=1:num
            fprintf(f, '%s\t%s\n', triplets{1}{k}, triplets{2}{k});
        end
    end
    fclose(f);
end
    

%% compute the number of edges from isorank matches
alpha = 0.80; ialpha = alpha * 100;
iters = 20;
similarity_m = 'isorank';


ig_m = 'isogreedy';
ih_m = 'isohungarian';

% isogreedy
ig_files = dir(fullfile(namematchdir, '*isogreedy*'));

% isohungarian
ih_files = dir(fullfile(namematchdir, '*isohungarian*'));

for i = 1:length(ig_files)
    name = ig_files(i).name;
    ig_path = fullfile(namematchdir, name);
    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    mat_path1 = fullfile(matdir, [name1, '.mat']);
    mat_path2 = fullfile(matdir, [name2, '.mat']);
    mdata1 = load(mat_path1);
    mdata2 = load(mat_path2);
    B = mdata1.A; B = max(B, B'); indexB = mdata1.index;
    A = mdata2.A; A = max(A, A'); indexA = mdata2.index;
    
    f = fopen(ig_path, 'r');
    triplets = textscan(f, '%s %s %f');
    fclose(f);
    num = length(triplets{1});
    
    m = size(B, 1);
    n = size(A, 1);
    ii = zeros(num, 1);
    jj = zeros(num, 1);
    
    
    for k=1:num
        bname = triplets{1}{k}; % name of protein from B
        aname = triplets{2}{k}; % name of protein from A
        ii(k) = indexB(bname); 
        jj(k) = indexA(aname);
    end
    
    
    M = sparse(ii, jj, 1, m, n);
    [mb, ma] = find(M);
    [edges, dummyx, dummyy] = count_overlap(B, A, [mb, ma]);
    
    % overlapped edges
    e_f = sprintf(e_p, name1, name2, similarity_m, ialpha, ig_m);
    e_path = fullfile(resultsdir, e_f);
    f = fopen(e_path, 'w'); 
    fprintf(f, '%d\n', edges); 
    fclose(f);                 
    
    results = [];
    % results(end+1).alignment = aG;
    results(end+1).match = M;
    results(end).edges = edges;
    r_f = sprintf(r_p, name1, name2, similarity_m, ialpha, ig_m);
    r_path = fullfile(resultsdir, r_f);
    save(r_path, 'results');
end



for i = 1:length(ih_files)
    name = ih_files(i).name;
    ih_path = fullfile(namematchdir, name);
    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    mat_path1 = fullfile(matdir, [name1, '.mat']);
    mat_path2 = fullfile(matdir, [name2, '.mat']);
    mdata1 = load(mat_path1);
    mdata2 = load(mat_path2);
    B = mdata1.A; B = max(B, B'); indexB = mdata1.index;
    A = mdata2.A; A = max(A, A'); indexA = mdata2.index;
    
    f = fopen(ih_path, 'r');
    triplets = textscan(f, '%s %s %f');
    fclose(f);
    num = length(triplets{1});
    
    m = size(B, 1);
    n = size(A, 1);
    ii = zeros(num, 1);
    jj = zeros(num, 1);
    
    
    for k=1:num
        bname = triplets{1}{k}; % name of protein from B
        aname = triplets{2}{k}; % name of protein from A
        ii(k) = indexB(bname); 
        jj(k) = indexA(aname);
    end
    
    
    M = sparse(ii, jj, 1, m, n);
    [mb, ma] = find(M);
    [edges, dummyx, dummyy] = count_overlap(B, A, [mb, ma]);
    
    % overlapped edges
    e_f = sprintf(e_p, name1, name2, similarity_m, ialpha, ih_m);
    e_path = fullfile(resultsdir, e_f);
    f = fopen(e_path, 'w'); 
    fprintf(f, '%d\n', edges); 
    fclose(f);                 
    
    results = [];
    % results(end+1).alignment = aG;
    results(end+1).match = M;
    results(end).edges = edges;
    r_f = sprintf(r_p, name1, name2, similarity_m, ialpha, ih_m);
    r_path = fullfile(resultsdir, r_f);
    save(r_path, 'results');
end

 
   
%% H:mm.sec to sec    


t_files = dir(fullfile(isorunsdir, '*timing*'));
for i = 1:length(t_files)
    name = t_files(i).name;
    t_path = fullfile(isorunsdir, name);
    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    
    f = fopen(t_path, 'r');
    data = fscanf(f, '%d:%d.%d');
    fclose(f);
    hflag = 0;
    if length(data)==2
        f = fopen(t_path, 'r');
        data = fscanf(f, '%d:%d:%d');
        fclose(f);
        hflag = 1;
    end    
    
    if hflag == 0
        secs = data(1) * 60  + data(2) + data(3)/100;
        fprintf(1, '%s-%s:[%.2f secs] %d min + %d sec + %d/100 sec\n', name1, name2, secs, data(1), data(2), data(3));
    else
       secs = data(1) * 3600 + data(2) * 60 + data(3); 
       fprintf(1, '%s-%s:[%.2f secs] %d hours + %d mins + %d sec\n', name1, name2, secs, data(1), data(2), data(3));
    end
    
    itt_p = '%s-%s-isorank-%03d-total-time.txt';
    itt_f = sprintf(itt_p, name1, name2, ialpha);
    itt_path = fullfile(timedir, itt_f);
    f = fopen(itt_path, 'w');
    fprintf(f, '%f\n', secs);
    fclose(f);
end


%% auction case
% matching (file with the single column) to name pairs file
alpha = 0.80; ialpha = int32(alpha * 100);
iters = 20;

similarity_m = 'mat3';
match_m = 'auction';

m_files = dir(fullfile(matchdir, '*auction*'));

for i = 1:length(m_files)
    name = m_files(i).name;
    m_path = fullfile(matchdir, name);
    ma = load(m_path, 'ascii');
    ma = ma + 1;


    parts = regexp(name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    mat_path1 = fullfile(matdir, [name1, '.mat']);
    mat_path2 = fullfile(matdir, [name2, '.mat']);
    mdata1 = load(mat_path1);
    mdata2 = load(mat_path2);
    
    nm_f = sprintf(nm_p, name1, name2, similarity_m, ialpha, match_m);
    nm_path = fullfile(namematchdir, nm_f);
    f = fopen(nm_path, 'w');
    names1 = mdata1.names;
    names2 = mdata2.names;
    for q = 1:length(ma)
        ii = q;
        jj = ma(q);
        fprintf(f, '%s\t%s\n', names1{ii}, names2{jj});
     end
     fclose(f);
end





alpha = 0.80; ialpha = int32(alpha * 100);

name2nodes = containers.Map;
for i=1:length(species)
    name2nodes(species(i).name) = species(i).nodes;
end

files = dir(fullfile(matchdir, '*auction*'));
for i=1:length(files)
    parts = regexp(files(i).name, '-', 'split');
    name1 = parts{1};
    name2 = parts{2};
    similarity_m = parts{3};
    match_m=parts{5};
    inpath = fullfile(matchdir, files(i).name);
    f = fopen(inpath, 'r');
    ma = fscanf(f, '%d');
    fclose(f);
    ma = ma + 1;
    
    mb = [1:length(ma)];
    M = sparse(mb, ma, 1, name2nodes(name1), name2nodes(name2));
    r_f = sprintf(r_p, name1, name2, similarity_m, ialpha, match_m);
    r_path = fullfile(resultsdir, r_f);
    
    results = [];
    results(end+1).match = M;
    save(r_path, 'results')
end







