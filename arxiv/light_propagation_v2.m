%LIGHT_PROPAGATION_V2	Simulation of the imaging process of CACTI system containing 
%a dmd and a  static mask. A plane wave incidentes into the system and propagate forward 
%to form a image on the sensor.
% 
%   Note: 
%   --------
%	- physical 3D coordinate system is:x-up, y-out, z-right
% 
%	- obj, dmd and mask is placed with their matrix mn coord-sys oppsite to
%     physical 2D xy coord-sys. (i.e., observing from anti optical axis direction,
%     the [1 1] element of matrix is at the left up corner of the physical
%     object(obj, dmd or mask)
% 
%     - resampling_factor: For dmd's MASK_PROPAGATE, resampling_factor =
%     50 is recommended, it means the spot on dmd will contain 50*50
%     virtual elements after resampling, which is proper from perspective of
%     memory-consumption and time-consumption. For mask's MASK_PROPAGATE,
%     resampling_factor = default (or = dmd resampling_factor) is recommended, 
%     it means the resampling_factor will be determined by the former mask
%     (i.e. the dmd), so that when combining two mask, they can directly 
%     match well without resize.
% 
%   Version:
%   --------
%   v2.0-beta
% 
%   See also:
%   --------
%   light_propagation_v1
% 
%	Todo:
%   --------
%   * use parallel computation to speed calculating
%   * implement camera sampling model
%   - implement true system data
%   * use matrix representation to substitute iterative for speeding
%     calculation (cellfun() or something)
%   * optimize the params's structure to make them more clear
%  
%	Log:
%   --------
%   - use matrix representation to substitute iterative for speeding
% 
%	Info:
%   --------
%   Created: Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-15
%   Last Modified:  Zhihong Zhang, 2020-04-01
%               
%   Refs:
%       - [1] T. Bishop, S. Zanetti, and P. Favaro, ¡°Light Field Superresolution,¡± 
%           in Proc. ICCP ¡¯09. IEEE International Conference on Computational Photography, 2009.
% 
%   Copyright (c) 2020 Zhihong Zhang.


%% init
clc, clear
close all
addpath('./utils');
addpath('./cacti');


%% env parameters
save_data_flag = 1;     % whether to save workspace data
load_params_flag = 0;   % whether to load sys_params
show_compare_flag = 1;  % whether to show obj, dmd, mask and image
drawing_sys_flag = 1;   % whether to draw the optical system
drawing_light_flag = 1; % whether to draw the light path
ideal_sensor_flag = 1;  % whether to assume sensor is ideal (sensor is the same as image plane)
resampling_factor = 50; % spot contains at least scale_factor.^2 elements

% paths and names
test_name = 'demo1_dmd_up_half';
config_dir = './config/';
result_dir = './result/subpix/';
data_name = [ test_name '_data.mat'];

if ~isfolder(result_dir)
    mkdir(result_dir)
end


%% system parameters (mm)
if load_params_flag
    % laod params
%     load([config_dir 'paparms/params.mat'])
    load([result_dir test_name '_params.mat']);
else
    % set params
    % obj
    obj_pos = -5;               % objective distance -50 mm
    obj_size = [200,200];       % pixel number
    obj_pix_size = 5e-3;        % 0.5um sampling(representing) size
    obj = 255*ones(obj_size);   % obj pattern, white board 
%     load([config_dir 'obj/apple_200.mat'], 'obj')

    % lens
    lens_f = 2.5;               % focal length
    lens_radius = 0.5;          % lens radius

    % DMD
    dmd_pos = 0;                % dmd position
    dmd_size = [100,100];
    dmd_pix_size = 10e-3;       % 10um pixel size
    dmd_t = 0.5;                % dmd's transmission rate [standard]
    dmd = binary_mask(dmd_size, dmd_t);                   % dmd pattern
%     dmd = binary_mask(dmd_size, 'fixed', 'up_half');        % dmd pattern

    % mask
    mask_pos =  4.8;                % dmd position
    mask_size = [200,200];          % pixel number
    mask_pix_size = 5e-3;           % 10um pixel size
    mask_t = 0.5;                   % mask's transmission rate [standard] 
    mask = rand(mask_size) + mask_t-0.5;  % mask pattern, gray
%     mask = binary_mask(mask_size,mask_t);  % mask pattern, binary
%    load([config_dir 'mask/mask_50.mat'], 'mask')

    % sensor
    sensor_pos = [];                % dmd position
    sensor_size = [100, 100];       % pixel number
    sensor_pix_size = 10e-3;        % 5um pixel size
    sensor = zeros(sensor_size);    % image pattern
    % calculate image position(sensor_pos)
    sensor_point = point_img([0,0,obj_pos]', lens_f);
    sensor_pos = sensor_point(3);

    % save params to result file
	params_regexp = join(["obj\w*" "lens\w*" "dmd\w*" "mask\w*" "sensor\w*" "ideal_sensor_flag"], '|');
    save([result_dir test_name '_params.mat'], '-regexp',params_regexp)
    print_params([result_dir test_name '_params_setting.txt'], params_regexp, '-regexp')
end


%% propagating

% tmp image
tmp_image = zeros(obj_size);

obj_row = obj_size(1);
obj_clo = obj_size(2);

fprintf('\n* start calculation...\n')

% get the index matrix of obj and convert to physical coordinates matrix
[obj_M, obj_N] = meshgrid(1:obj_size(1), 1:obj_size(2));
obj_M = obj_M';
obj_N = obj_N';
obj_x = -sub2coord(obj_M, obj_size(1), obj_pix_size, 'mn'); % physical x-coordinate matrix
obj_y = -sub2coord(obj_N, obj_size(2), obj_pix_size, 'mn'); % physical y-coordinate matrix
obj_z = obj_pos*ones(obj_size);                             % physical z-coordinate matrix

% physical coordinates of obj points
obj_points = cat(3, obj_x, obj_y, obj_z);

% geometrical measurement
[dmd_spot_p, dmd_spot_r] = blur_spot(obj_points, lens_f, 2*lens_radius, dmd_pos);
[mask_spot_p, mask_spot_r] = blur_spot(obj_points, lens_f, 2*lens_radius, mask_pos);

% intensity measurement
% start timer
tic
for obj_m = 1:obj_row 
    for obj_n = 1:obj_clo
        
        % intensity measurement, be careful about the coord-sys's direction, mn-coord-sys and current physical xy-coord-sys is opposite
        % dmd  
        [~, dmd_masked_spot] = mask_propagate(dmd, dmd_pix_size,...
            -squeeze(dmd_spot_p(obj_m,obj_n, 1:2)), dmd_spot_r, resampling_factor, 'mn');
        
        % mask
        mask_coverd_t = mask_propagate(mask, mask_pix_size, -squeeze(mask_spot_p(obj_m,obj_n, 1:2)),...
            mask_spot_r, dmd_masked_spot,'mn', 'non-test');

        % image
        tmp_image(obj_m,obj_n) = obj(obj_m,obj_n)*mask_coverd_t; % need to rot 180 degree to form the sensor image
        
    end
     
    % processing bar and light drawing
    if mod(obj_m,10)==0
        % processing bar 
        disp([num2str(obj_m/obj_row*100) '% done£¡']);
        toc
    end
    
end
toc

% sensor sampling
tmp_img_pix_size = -sensor_pos/obj_pos*obj_pix_size; % amplification ratio
if ideal_sensor_flag
    sensor = sensor_imaging(tmp_image);
else
    sensor = sensor_imaging(tmp_image, tmp_img_pix_size, sensor_size, sensor_pix_size);
end


%% save data
if save_data_flag
 save([result_dir data_name], '-regexp', '\w*[^_]\>');
end

%% illustration
% show obj and image
figure
imshow(rot90(sensor,2),[ ])
title('rescaled sensor(rot180)');
figure
imshow(obj,[ ])
title('rescaled obj ')

% show all origin image
if show_compare_flag
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
    imshow(sensor)
    title('sensor')
end




%% show light propagating
if drawing_flag
    % set sys params
    % obj
    sys_params(1).name = 'obj';
    sys_params(1).pos = obj_pos;
    sys_params(1).size = obj_size;
    sys_params(1).pattern = obj;
    sys_params(1).spot = [squeeze(obj_points(1,1,:)); 0];
    sys_params(1).element_sz = obj_pix_size;

    % lens
    sys_params(2).name = 'lens';
    sys_params(2).pos = 0;
    sys_params(2).size = 2*lens_radius;
    sys_params(2).spot = [[0,0,0]'; lens_radius];


    % dmd
    sys_params(3).name = 'dmd';
    sys_params(3).pos = dmd_pos;
    sys_params(3).size = dmd_size;    
    sys_params(3).pattern = dmd;
    sys_params(3).spot = [squeeze(dmd_spot_p(1,1,:)); dmd_spot_r];
    sys_params(3).element_sz = dmd_pix_size;

    % mask
    sys_params(4).name = 'mask';
    sys_params(4).pos = mask_pos;
    sys_params(4).size = mask_size;    
    sys_params(4).pattern = mask;
    sys_params(4).spot = [squeeze(mask_spot_p(1,1,:)); mask_spot_r];
    sys_params(4).element_sz = mask_pix_size;
    
     % sensor
    sensor_point = point_img(squeeze(obj_points(1,1,:)), lens_f);
    sys_params(5).name = 'sensor';
    sys_params(5).pos = sensor_pos;
    sys_params(5).size = sensor_size;    
    sys_params(5).pattern = sensor;
    sys_params(5).spot = [sensor_point; 0]; %in-focus
    sys_params(5).element_sz = sensor_pix_size;
    
    light_drawing(sys_params)
end


%% backup
%{
% real time drawing (when use "parfor", comment the following code for drawing)
if drawing_flag
    sys_params(1).spot = [point; 0];
    sys_params(3).spot = [dmd_spot_p; dmd_spot_r];
    sys_params(4).spot = [mask_spot_p; mask_spot_r];
    sys_params(5).spot = [sensor_point; 0];
    light_drawing(sys_params, 'draw_light', fig_drawing_)
    drawnow
end
%}


