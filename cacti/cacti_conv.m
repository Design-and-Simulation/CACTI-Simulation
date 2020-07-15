function  image = cacti_conv(obj, dmd, mask, sys_params, ctrl_params)
%CACTI_CONV    Coded Aperture Compressive Temporal Imaging System
%   Simulation of the imaging process of a CACTI system, which contains a  
%   object, a dmd, a lens, a mask and a sensor.
% 
%   Input:
%   --------
%   - obj: object pattern, 2D numeric matrix
% 
%   - dmd: dmd pattern, 2D logical(double) matrix, 0|1
% 
%   - mask: mask pattern, 2D logical matrix (binary mask,0|1) or 2D float
%     matrix(gray mask, [0,1])
% 
%   - sys_params: system parameters, 1*1 struct, with the following fields:
%       obj_pos, obj_size, obj_pix_size;
%       lens_pos, lens_f, lens_radius;
%       dmd_pos, dmd_size, dmd_pix_size;
%       mask_pos, mask_size, mask_pix_size;
%       sensor_pos, sensor_size, sensor_pix_size
% 
%   - ctrl_params: controlling parameters, 1*1 struct, optional, with the 
%     following fields:
%       show_compare_flag   % show obj, dmd, mask and image, default = 1;
%       drawing_flag		% drawing the optical system
%       ideal_sensor_flag   % sensor conducts ideal sampling     
% 
%   Output:
%   --------
%   - image: the image on the sensor, 2D numeric matrix
% 
%   Note: 
%   --------
%   * system description:
%       - physical 3D coordinate system is:x-up, y-out, z-right
% 
%       - the optical element(obj, dmd, mask,..) is placed with their matrix mn 
%       coord-sys oppsite to physical 2D xy coord-sys. (i.e., observing from 
%       anti optical axis direction, the [1 1] element of matrix is at the left
%       up corner of the physical object(obj, dmd or mask)
% 
%   * function note:
% 
%     - ideal_sensor_flag: if true, we assume that the sensor conducts ideal 
%     sampling, i.e.,the sensor image is the same as the image plane, the
%     physical size or AD conversion is not considerated.
% 
% 
%   See also:
%   --------
%   CACTI_IMAGING_CONV, MASK_CONV
% 
%   Log:
%   --------
% 
%   Info:
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-07-15
%   Last Modified:  Zhihong Zhang, 2020-07-15
%               
%   Refs:
%       - [1] T. Bishop, S. Zanetti, and P. Favaro, ¡°Light Field Superresolution,¡± 
%           in Proc. ICCP ¡¯09. IEEE International Conference on Computational Photography, 2009.
% 
%   Copyright (c) 2020 Zhihong Zhang.
% 


%% input setting
% input check
narginchk(4,5);

% input assignment
if nargin < 5
    % parameter
    show_compare_flag = 1;              % show obj, dmd, mask and image
    ctrl_params.ideal_sensor_flag = 0;  % whether to assume sensor is ideal(sensor is the same as image plane)
else
    % convert params to variables for fast calculation
    show_compare_flag = ctrl_params.show_compare_flag;
    ideal_sensor_flag = ctrl_params.ideal_sensor_flag;
end


%% env init
% addpath(genpath('..'));
addpath('./utils');
addpath('./cacti');


%% calculation
% init
% convert params to variables for fast calculation
% obj
obj_pos = sys_params.obj_pos;
obj_size = sys_params.obj_size;
obj_pix_size = sys_params.obj_pix_size;
% lens
lens_pos = sys_params.lens_pos;
lens_f = sys_params.lens_f;
lens_radius = sys_params.lens_radius;
% dmd
dmd_pos = sys_params.dmd_pos;
dmd_size = sys_params.dmd_size;
dmd_pix_size = sys_params.dmd_pix_size;
% mask
mask_pos = sys_params.mask_pos;
mask_size = sys_params.mask_size;
mask_pix_size = sys_params.mask_pix_size;
% sensor
sensor_pos = sys_params.sensor_pos;
sensor_size = sys_params.sensor_size;
sensor_pix_size = sys_params.sensor_pix_size;

obj_row = obj_size(1);
obj_col = obj_size(2);

% judge whether dmd is at the aperture plane
dmd_pos_zero = 0;
if dmd_pos == 0
    dmd_pos_zero = 1;
end

obj = single(obj); % convert from int to single for calculation accuracy
combined_mask_t = zeros(obj_size);    % combined mask's transimttance

% propagating
fprintf('\n* start calculation...\n')

% geometrical measurement
% dmd
if ~dmd_pos_zero
	error(['this code can only deal with the situation where DMD is located at'...
	'the aperture plane, please use "cacti_imaging.m" instead.']); 
end

% spot on the mask
obj_point = sub2coord(obj_size', obj_size', obj_pix_size, 'mn'); 
obj_point = [obj_point; obj_pos];
[mask_spot_p, mask_spot_r] = blur_spot(obj_point, lens_f,...
    2*lens_radius, mask_pos);

% effective mask area
mask_spot_p_mn = ceil(abs(mask_spot_p(1:2))./mask_pix_size);
mask_spot_r_mn = ceil(mask_spot_r./mask_pix_size);
effect_area_r_mn = mask_spot_p_mn + mask_spot_r_mn;

up_left = round(mask_size./2) - effect_area_r_mn;
up_left = max(up_left, [1 1]); % in case of out of boundary
right_down = round(mask_size./2) + effect_area_r_mn - 1;
right_down = min(right_down, mask_size);

croped_mask = mask(up_left(1):right_down(1),up_left(2):right_down(2));

effect_mask = zeros(2*effect_area_r_mn');
croped_mask_size = size(croped_mask);

if all(croped_mask_size == 2*effect_area_r_mn)
	effect_mask = croped_mask;
else
	up_left = effect_area_r_mn - round(croped_mask_size(1)./2);
	up_left = max(up_left, [1 1]); % in case of out of boundary
	right_down = up_left + croped_mask_size;
	right_down = min(right_down, 2*effect_area_r_mn);

	effect_mask(up_left(1):right_down(1),up_left(2):right_down(2)) = croped_mask;
end

% calculate the combined masks' transmittance matrix
% start timer
tic
dmd_shrink_ratio = mask_spot_r/lens_radius;
pix_ratio =(dmd_pix_size.*dmd_shrink_ratio) ./ mask_pix_size;
combined_mask_t = mask_conv(dmd, effect_mask,  pix_ratio);

% crop center area as the sensor's collection
up_left = round((mask_size - obj_size)./2);
combined_mask_t = imcrop(combined_mask_t, [up_left obj_size-1]);
toc

% tmp image
combined_mask_t = imresize(combined_mask_t, obj_size);
tmp_image = combined_mask_t.*obj;

% sensor sampling(assume the sensor is palced at the image plane)
tmp_img_pix_size = -sensor_pos/obj_pos*obj_pix_size; % amplification ratio
if ideal_sensor_flag
    sensor_img = sensor_imaging(tmp_image, 'single');
% 	sensor_img = sensor_imaging(tmp_image);
else
    sensor_img = sensor_imaging(tmp_image, tmp_img_pix_size, sensor_size, sensor_pix_size);
end

% output the sensor image
image = sensor_img;


%% illustration
% show all images
if show_compare_flag
    % show result image
    figure
    imshow(rot90(image,2)), title(['result(convert 180)-' class(image)])
    % show other images
    figure
    subplot(221)
    imshow(obj)
    title('obj')
    subplot(222)
    imshow(dmd)
    title('dmd')
    subplot(223)
    imshow(mask)
    title('mask')
    subplot(224)
    imshow(image)
    title(['sensor image-' class(image)])
end

end

