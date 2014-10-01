function [dt] = dmat_write(A, filepath)
% dmat_write Writes a dense matrix to a file.
% 
% This file will contain the matrix in dmat format:
% First line of this text file (header) contains 2 integers, 
% m (number of rows) and n (number of columns).
% This is followed by m * n lines, each containing a matrix element;
% the matrix is traversed row-wise (row-major order).
% 
% Input arguments:
% - A: the dense matrix.
% - filepath: the .dmat file location.
% Output arguments:
% - dt: the time in seconds for the operation.


t0 = clock;
[m,n] = size(A);
f = fopen(filepath, 'wt');
fprintf(f, '%d %d\n', m, n);
data = A';
fprintf(f, '%.12f\n', data);
fclose(f);
dt = etime(clock, t0); 

end
