function distance_map = dist_map(mat_sz, elem_sz, P)
%DIST_MAP	calculating physical distance map of matrix elemetns(physical point set)
% relative to a reference physical point(assuming the physical origin of the coord-sys 
% is at the center of the matrix, and the axis direction is the same with mn-coord-sys's
%  direction, i.e., down and right is the position diretion for first and second dimension)
% 
%	Input:
%   --------
%   - mat_sz: matrix size, int vector
% 
%   - elem_sz: matrix element's physical size, float vector 
% 
%   - P: reference point's physical coordinates, float vector
% 
%   Output:
%   --------
%   - distance_map: distance map between matrix elements and reference point
% col and row idx of mask's all elements (center element (on the optical axis) as origin point)
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang

[mask_m, mask_n] = meshgrid(1:mat_sz(1), 1:mat_sz(2));
mask_m = mask_m';
mask_n = mask_n';

% physical coordinates of mask's elements (center element (on the optical axis) as origin point)
mask_x = sub2coord(mask_m, mat_sz(1), elem_sz(1), 'mn');
mask_y = sub2coord(mask_n, mat_sz(2), elem_sz(2), 'mn');

% distance map between mask element and spot center
distance_map = sqrt((mask_x - P(1)).^2 + (mask_y - P(2)).^2);
end