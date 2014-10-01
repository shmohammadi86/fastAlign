function [A, names, index, dt] = tab_read(filepath)
% tab_read Parses a .tab file with PPI network data.
%
% A .tab file is a text file with the following format:
% Its fist line is of no interest and it is skipped
% (here: INTERACTOR_A	INTERACTOR_B)
% Each of the following lines contains 2 strings which are the names of 2
% interacting proteins. The .tab files can be downloaded from
% http://groups.csail.mit.edu/cb/mna/
% 
% Input arguments: 
%  - filepath: the .tab file location.
%
% Output arguments:
% - A: the sparse matrix which is the adjacency matrix of PPI data. To 
%     generate the matrix protein names have been mapped to integers and
%     vice versa.
% - names: a cell array with the integer to name mapping, i.e.
%     name{i} is the (protein) name mapped to integer i.
% - index: a Map with the name to integer mapping, i.e.
%     index(name) is the integer mapped to (protein) name. 
% - dt: the time in seconds for the operation.



t0 = clock;

f = fopen(filepath, 'r');
header = textscan(f, '%s', 2);

pairs = textscan(f, '%s %s');
xnames = pairs{1};
ynames = pairs{2};

% build the index
allnames = [xnames; ynames];
num = size(allnames, 1);

names = {};
index = containers.Map;
counter = 1;
for i=1:num
    name = allnames(i);
    name = name{1};
    if ~isKey(index, name)
        names{counter} = name;
        index(name) = counter;
        counter = counter + 1;
    end
end
m = counter - 1;
n = m;
% build the matrix
nnz = size(xnames, 1);
ind_i = zeros(nnz, 1);
ind_j = zeros(nnz, 1);
data = ones(nnz, 1);

for k=1:nnz
    xname = xnames(k);
    yname = ynames(k);
    xname = xname{1};
    yname = yname{1};

    ind_i(k) = index(xname);
    ind_j(k) = index(yname);
end

A = sparse(ind_i, ind_j, data, m, n, nnz);

fclose(f);
dt = etime(clock, t0); 



end

