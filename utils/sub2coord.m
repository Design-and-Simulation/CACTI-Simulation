function coord = sub2coord(element_idx, mat_sz, element_sz, coord_direction)
%SUB2COORD calculate the physical coordinates of given matrix elements
%(assuming the m*n matrix form a m*n blank, each element is at the center
%of corresponding sub-blank, and physical origin point is at the center of the large blank)
% 
%   Input:
%   --------
%   - element_idx: [row, col] of given element, 2x1 int vector
% 
%   - mat_sz: shape of given matrix, 2x1 int vector 
% 
%   - element_sz: physical size of each matrix element, float scalar
% 
%   - coord_direction: str, optional, {'mn', 'xy'}, default = 'xy'; choose to
%       use mn coordinate system direction(down,right) or xy coordinate system
%       direction(right,up) for the output coord
% 
%   Output:
%   --------
%   - coord: physical coordinates of given matrix elements, 2x1 float vector
%   (top(x) left(y) is positive direction)
% 
% 
%   Note:
%   --------
%   This function can aslo deal with row's/col's index matrix's conversion to
%   physical x/y-coordnate matrix, where "element_idx" is a index matrix, mat_sz
%   is the num of rows/columns of this matrix, and element_sz is the 
%   physical height/width of each element of this matrix. Here is an example:
%       [mask_m, mask_n] = meshgrid(1:mat_sz(1), 1:mat_sz(2));
%       mask_m = mask_m';
%       mask_n = mask_n';
%       mask_x = sub2coord(mask_m, mat_sz(1), elem_sz(1), 'mn');
%       mask_y = sub2coord(mask_n, mat_sz(2), elem_sz(2), 'mn');
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang

if nargin < 4
     coord_direction = 'xy';
end

center_idx = (1+mat_sz)/2;
coord = (element_idx - center_idx)*element_sz; 

if strcmp(coord_direction, 'xy')
    % convert from m-n coordinate system direction to x-y coordinate system direction
    coord = [0 1; -1 0]*coord;
end

