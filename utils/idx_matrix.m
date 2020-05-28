function [idx_mat_m, idx_mat_n] = idx_matrix(A)
%IDX_MATRIX calculate the index (subscript) matrix of given matrix (size)
% 
%   This function can generate the row & col index matrixs of given matrix (size).
%   The row/col index matrix's elements can be regard as the row.col subscript 
%   (coordinate) of given matrix's corresponding elements
% 
%   Input:
%   --------
%   - A: 2D matrix or 2 element vecotr(matrix's size)
% 
%   Output:
%   --------
%   - idx_mat_m: row index(subscript) matrix
%   - idx_mat_n: col index(subscript) matrix
% 
% 
%   Note:
%   --------
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-26
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-26   
%   
%   Copyright 2020 Zhihong Zhang

if isvector(A)
    mat_size = A;
elseif ismatrix(A)
    mat_size = size(A);
else
    error("input error");
end

% implement 1
% [idx_mat_m, idx_mat_n] = meshgrid(1:mat_size(1), 1:mat_size(2));
% idx_mat_m = idx_mat_m';
% idx_mat_n = idx_mat_n';

% implement 2
row_idx = 1:mat_size(1);
idx_mat_m = repmat(row_idx',[1, mat_size(2)]);
col_idx = 1:mat_size(2);
idx_mat_n = repmat(col_idx,[mat_size(1), 1]);

end