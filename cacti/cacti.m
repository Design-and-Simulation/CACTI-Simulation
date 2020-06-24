function  image = cacti(obj, dmd, mask, sys_params, ctrl_params)
%CACTI    Coded Aperture Compressive Temporal Imaging System
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
%       resampling_factor   % resampling factor for MASK_PROPAGATE
%       ideal_sensor_flag   % sensor conducts ideal sampling
%       TEST_MODE_FLAG_     % test mode for MASK_PROPAGATE, default = 'non-test'
%       parallel_flag       % whether to use parallel calculation in CACTI      
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
%     - resampling_factor: For dmd's MASK_PROPAGATE, resampling_factor =
%     50 is recommended, it means the spot on dmd will contain 50*50
%     virtual elements after resampling, which is proper from perspective of
%     memory-consumption and time-consumption. For mask's MASK_PROPAGATE,
%     resampling_factor = default (or = dmd resampling_factor) is recommended, 
%     it means the resampling_factor will be determined by the former mask
%     (i.e. the dmd), so that when combining two mask, they can directly 
%     match well without resize.
% 
%     - ideal_sensor_flag: if true, we assume that the sensor conducts ideal 
%     sampling, i.e.,the sensor image is the same as the image plane, the
%     physical size or AD conversion is not considerated.
% 
%     - TEST_MODE_FLAG_: the test flag for function MASK_PROPAGATE, it only
%     take effect when parallel_flag=false (in paraller mode, the function 
%     test is not permitted. {'test_dmd', 'test_mask', 'test_all','non_test'}
% 
% 
%   See also:
%   --------
%   POINT_IMAGE, MASK_PROPAGATE£¬ 
% 
%   Log:
%   --------
% 
%   Info:
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-30
%   Last Modified:  Zhihong Zhang, 2020-03-30
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
    dmd_resampling_factor = 50;         % dmd's MASK_PROPAGATE, resampling_factor
    mask_resampling_factor = dmd_resampling_factor;  % mask's MASK_PROPAGATE, resampling_factor
	TEST_MODE_FLAG_ = 'non-test';       % test mode, 'test_dmd', 'test_mask', 'test_all', 'non-test'
    parallel_flag = 1;                  % whether to use parallel calculation in CACTI
    ctrl_params.ideal_sensor_flag = 0;  % whether to assume sensor is ideal(sensor is the same as image plane)
else
    % convert params to variables for fast calculation
    show_compare_flag = ctrl_params.show_compare_flag;
    dmd_resampling_factor = ctrl_params.dmd_resampling_factor;
    mask_resampling_factor = ctrl_params.mask_resampling_factor;
	TEST_MODE_FLAG_ = ctrl_params.TEST_MODE_FLAG_;
    parallel_flag = ctrl_params.parallel_flag;
    ideal_sensor_flag = ctrl_params.ideal_sensor_flag;
end

switch TEST_MODE_FLAG_
    case 'test_dmd'
        DMD_TEST_MODE_FLAG_ = 'test';
        mask_TEST_MODE_FLAG_ = 'non_test';
    case 'test_mask'
        DMD_TEST_MODE_FLAG_ = 'non_test';
        mask_TEST_MODE_FLAG_ = 'test';
    case 'test_all'
        DMD_TEST_MODE_FLAG_ = 'test';
        mask_TEST_MODE_FLAG_ = 'test';
    case 'non_test'
        DMD_TEST_MODE_FLAG_ = 'non_test';
        mask_TEST_MODE_FLAG_ = 'non_test';  
    otherwise
        error('error input - TEST_MODE_FLAG_')
end

if parallel_flag && ~strcmp(TEST_MODE_FLAG_, 'non_test')
    warning('function test is not allowed in parallel mode, omitted - TEST_MODE_FLAG_')
    DMD_TEST_MODE_FLAG_ = 'non_test';
    mask_TEST_MODE_FLAG_ = 'non_test'; 
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

% get the index matrix of obj and convert to physical coordinates matrix
[obj_M, obj_N] = idx_matrix(obj_size);
% physical x-coordinate matrix
obj_x = -sub2coord(obj_M, obj_row, obj_pix_size, 'mn'); 
% physical y-coordinate matrix
obj_y = -sub2coord(obj_N, obj_col, obj_pix_size, 'mn');
% physical z-coordinate matrix
obj_z = obj_pos*ones(obj_size);                                     

% physical coordinates of obj points
obj_points = cat(3, obj_x, obj_y, obj_z);

% geometrical measurement
% dmd
if dmd_pos_zero
    % if dmd_pos = 0, the dmd_masked_spot is the same for all object points, i.e.
    % the spot center is at [0 0 0]', and the spot radius is the lens radius,
    % so it's unnecessary to calculate this info.    
    dmd_spot_p = zeros([obj_size, 3]); dmd_spot_r = lens_radius;
    [~, dmd_masked_spot] = mask_propagate(dmd, dmd_pix_size, -squeeze(dmd_spot_p(1,1, 1:2)),...
        dmd_spot_r, dmd_resampling_factor, 'mn', DMD_TEST_MODE_FLAG_);  
else
    [dmd_spot_p, dmd_spot_r] = blur_spot(obj_points, lens_f,...
        2*lens_radius, dmd_pos);
    % 2d coordinates
    dmd_spot_p_2d = -dmd_spot_p(:,:, 1:2);    
end
% mask
[mask_spot_p, mask_spot_r] = blur_spot(obj_points, lens_f,...
    2*lens_radius, mask_pos);
% 2d coordinates
mask_spot_p_2d = -mask_spot_p(:,:, 1:2);

% calculate the combined masks' transmittance matrix
% start timer
tic
if ~parallel_flag && ~dmd_pos_zero
%   doesn't use parallel calculation && dmd is not located at aperture plane
    for obj_m = 1:obj_row 
        dmd_spot_p_2d_row = squeeze(dmd_spot_p_2d(obj_m,:, :))';
        mask_spot_p_2d_row = squeeze(mask_spot_p_2d(obj_m,:, :))';
        for obj_n = 1:obj_col    
            % dmd_masked_spot
            [~, dmd_masked_spot] = mask_propagate(dmd, dmd_pix_size, dmd_spot_p_2d_row(:,obj_n),...
                dmd_spot_r, dmd_resampling_factor, 'mn', DMD_TEST_MODE_FLAG_);
            % combined mask's transmittance
            mask_coverd_t = mask_propagate(mask, mask_pix_size, mask_spot_p_2d_row(:,obj_n),...
                mask_spot_r, dmd_masked_spot, mask_resampling_factor, 'mn', mask_TEST_MODE_FLAG_);
            % image point gray value: need to rot 180 degree to form the sensor image 
            combined_mask_t(obj_m,obj_n) = mask_coverd_t;        
        end

        % processing bar
        if mod(obj_m,20)==0
            % processing bar 
            disp([num2str(obj_m/obj_row*100) '% done£¡']);
            toc
        end   
    end
    
elseif ~parallel_flag && dmd_pos_zero
%   doesn't use parallel calculation && dmd is located at aperture plane
    for obj_m = 1:obj_row 
        mask_spot_p_2d_row = squeeze(mask_spot_p_2d(obj_m,:, :))';
        for obj_n = 1:obj_col    
            % combined mask's transmittance
            mask_coverd_t = mask_propagate(mask, mask_pix_size, mask_spot_p_2d_row(:,obj_n),...
                mask_spot_r, dmd_masked_spot, mask_resampling_factor, 'mn', mask_TEST_MODE_FLAG_);
            % image point gray value: need to rot 180 degree to form the sensor image 
            combined_mask_t(obj_m,obj_n) = mask_coverd_t;  
%             mask_t_store(obj_m,obj_n) = mask_coverd_t;  % for test
        end

        % processing bar
        if mod(obj_m,20)==0
            % processing bar 
            disp([num2str(obj_m/obj_row*100) '% done£¡']);
            toc
        end   
    end

elseif parallel_flag && ~dmd_pos_zero   
%	use parallel calculation && dmd is not located at aperture plane
    parfor obj_m = 1:obj_row 
        dmd_spot_p_2d_row = squeeze(dmd_spot_p_2d(obj_m,:, :))';
        mask_spot_p_2d_row = squeeze(mask_spot_p_2d(obj_m,:, :))';
        for obj_n = 1:obj_col    
            % dmd_masked_spot
            [~, dmd_masked_spot] = mask_propagate(dmd, dmd_pix_size, dmd_spot_p_2d_row(:,obj_n),...
                dmd_spot_r, dmd_resampling_factor, 'mn');
            % combined mask's transmittance
            mask_coverd_t = mask_propagate(mask, mask_pix_size, mask_spot_p_2d_row(:,obj_n),...
                mask_spot_r, dmd_masked_spot, mask_resampling_factor, 'mn');
            % image point gray value: need to rot 180 degree to form the sensor image 
            combined_mask_t(obj_m,obj_n) = mask_coverd_t;        
        end

        % processing bar
        if mod(obj_m,20)==0
            % processing bar 
            disp([num2str(obj_m/obj_row*100) '% done£¡']);
        end   
    end
    
elseif parallel_flag && dmd_pos_zero   
%	use parallel calculation && dmd is located at aperture plane
    parfor obj_m = 1:obj_row 
        mask_spot_p_2d_row = squeeze(mask_spot_p_2d(obj_m,:, :))';
        for obj_n = 1:obj_col   
            % combined mask's transmittance
            mask_coverd_t = mask_propagate(mask, mask_pix_size, mask_spot_p_2d_row(:,obj_n),...
                mask_spot_r, dmd_masked_spot, mask_resampling_factor, 'mn');
            % image point gray value: need to rot 180 degree to form the sensor image 
            combined_mask_t(obj_m,obj_n) = mask_coverd_t;        
        end

        % processing bar
        if mod(obj_m,20)==0
            % processing bar 
            disp([num2str(obj_m/obj_row*100) '% done£¡']);
        end   
    end   
end
toc

% tmp image
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

