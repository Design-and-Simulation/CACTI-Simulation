function element_idx = coord2sub(coord, mat_sz, element_sz, coord_direction, idx_type)
%COORD2SUB alculate the row and column index from given physical coordinates
%(assuming the m*n matrix form a m*n blank, each element is at the center
%of corresponding sub-blank, and physical origin point is at the center of the large blank)
% 
%   Input:
%   --------
%	- coord: physical coordinates of given matrix elements, 2x1 float vector
%           (top(x) left(y) is positive direction)
% 
%   - mat_sz: shape of given matrix, 2x1 int vector
% 
%   - element_sz: physical size of each matrix element, float scalar'
%   - coord_direction: str, optional, {'mn', 'xy'}, default = 'xy'; choose to
%       use mn coordinate system direction(down,right) or xy coordinate system
%       direction(right,up) for the input coord
% 
%   - ind_idx: whether fix output coord to integer, optional, str,
%           {'idx_type', 'float_idx'}, default = 'idx_type')
% 
%   Output:
%   --------
%   - element_idx: [row, col] of given element, 2x1 int vector
% 
%   See also:
%   --------
%   SUB2COORD
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

if nargin < 5
     idx_type = 'int_idx';
end

coord_direction = [lower(char(coord_direction)) ' ']; % Protect against short string
idx_type = [lower(char(idx_type)) ' ']; % Protect against short string

if coord_direction(1)=='x' % xy
    % convert from m-n coordinate system direction to x-y coordinate system direction
    coord = [0 -1; 1 0]*coord;
end

center_idx = (1+mat_sz)/2;
if idx_type(1)=='i' % int_idx
    element_idx = round(coord./element_sz + center_idx);
elseif idx_type(1)=='f' % float_idx
    element_idx = (coord./element_sz + center_idx);
end



end