function  cacti_drawing(sys_params, varargin)
%CACTI_DRAWING illustrating the cacti optical system using given optical system
%parameters
% 
%   Input:
%   --------
%	- sys_params: system parameters
% 
%   - varargin, contain the following optional arguements
%       unit_len[k]     - the unit length of figure coordinate system 
%                         [key, str, '-unit_len_k']
%       unit_len[v]     - the unit length of figure coordinate system 
%                         [value, float, default = smallest element's size]
%       draw_what       - {'draw_sys', 'draw_light', 'all'}
%       dimension       - {'2D', '3D'}
%       fig_handle      - figure handle, optional, default = generate a new axes;
%       
%   Output:
%   --------
%   the optical path figure
% 
%   Note: 
%   --------
%	- physical 3D coordinate system is:x-up, y-out, z-right
% 
%	- the optical element(obj, dmd, mask,..) is placed with their matrix mn 
%   coord-sys oppsite to physical 2D xy coord-sys. (i.e., observing from 
%   anti optical axis direction, the [1 1] element of matrix is at the left
%   up corner of the physical object(obj, dmd or mask)
% 
%   See also:
%   --------
%	LIGHT_DRAWING
% 
%   Info:
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-06-24
%   Last Modified:  Zhihong Zhang, 2020-06-24
%               
%   Copyright 2020 Zhihong Zhang


% set sys params
% obj
drawing_params(1).name = 'obj';
drawing_params(1).pos = sys_params.obj_pos;
drawing_params(1).size = sys_params.obj_size;
drawing_params(1).pattern = sys_params.obj;
drawing_params(1).element_sz = sys_params.obj_pix_size;
drawing_params(1).spot = [[sys_params.obj_pix_size.*sys_params.obj_size/2, sys_params.obj_pos]'; 0];

% lens
drawing_params(2).name = 'lens';
drawing_params(2).pos = sys_params.lens_pos;
drawing_params(2).size = 2*sys_params.lens_radius;
drawing_params(2).spot = [[0,0,sys_params.lens_pos]'; sys_params.lens_radius];


% dmd
drawing_params(3).name = 'dmd';
drawing_params(3).pos = sys_params.dmd_pos;
drawing_params(3).size = sys_params.dmd_size;    
drawing_params(3).pattern = sys_params.dmd;
% drawing_params(3).spot = [squeeze(sys_params.dmd_spot_p(1,1,:)); dmd_spot_r];
drawing_params(3).element_sz = sys_params.dmd_pix_size;

% mask
drawing_params(4).name = 'mask';
drawing_params(4).pos = sys_params.mask_pos;
drawing_params(4).size = sys_params.mask_size;    
drawing_params(4).pattern = sys_params.mask;
% drawing_params(4).spot = [squeeze(sys_params.mask_spot_p(1,1,:)); mask_spot_r];
drawing_params(4).element_sz = sys_params.mask_pix_size;

 % sensor
sensor_point = point_img(drawing_params(1).spot(1:3), sys_params.lens_f);
drawing_params(5).name = 'sensor';
drawing_params(5).pos = sys_params.sensor_pos;
drawing_params(5).size = sys_params.sensor_size;    
% drawing_params(5).pattern = sys_params.sensor_img;
drawing_params(5).spot = [sensor_point; 0]; %in-focus
drawing_params(5).element_sz = sys_params.sensor_pix_size;
fig_drawing_ = figure;
light_drawing(drawing_params,fig_drawing_) %,'draw_sys'
drawnow

end