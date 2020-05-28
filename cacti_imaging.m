%CACTI_IMAGING: Imaging with Coded Aperture Compressive Temporal Imaging System
%   Simulation of the imaging process of a CACTI system, which contains an  
%   object, a dmd, a lens, a mask and a sensor.
% 
%   CACTI Setting:
%   --------
%   - obj: object pattern stack(gray image), 3D numeric matrix, [h, w, batch_num]
% 
%   - dmd: dmd pattern stack, 3D logical(double) matrix, 0|1, [h, w, batch_num]
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
%   - image: the result image(gray image) stack, 3D numeric matrix, [h, w, batch_num]
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
%   CACTI, MASK_PROPAGATE£¬ 
% 
%   Log:
%   --------
% 
%   Info:
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-04-01
%   Last Modified:  Zhihong Zhang, 2020-04-01
%               
%   Refs:
%       - [1] T. Bishop, S. Zanetti, and P. Favaro, ¡°Light Field Superresolution,¡± 
%             in Proc. ICCP ¡¯09. IEEE International Conference on Computational Photography, 2009.
%       - [2] dvp, Yang Liu
% 
%   Copyright (c) 2020 Zhihong Zhang.
% 


%% init
clc, clear
% close all
addpath('./utils');
addpath('./cacti');


%% env setting
% flags
save_data_flag = 0;         % save data to the local result file
save_params_flag = 0;       % save data to the local result file
load_ctrl_params_flag = 0;  % 0-setting; 1-load from file; 2-use default
load_sys_params_flag = 0;
show_image_flag = 1;
load_patterns_flag = 1;           % load patterns for obj, dmd and mask
save_dvp_type = 1;                % save simulation dataset for dvp recon algorithm
save_dl_cacti_type = 1;           % save simulation dataset for dl_cacti recon algorithm
whiteboard_flag = 0;			  % whether the obj is a whiteboard

% paths and names
test_name = 'traffic_512_8frame';
config_dir = './config/';
result_dir = './result/dvp_traffic/';
data_name = [ test_name '_data.mat'];
image_name = [ test_name '_image.mat'];


%% control parameter
if load_ctrl_params_flag
    ctrl_params = load([result_dir test_name '_ctrl_params.mat']);
else
    ctrl_params.show_compare_flag = 0;            % whether to show obj, dmd, mask and image
    ctrl_params.drawing_sys_flag = 0;             % whether to draw optical system
    ctrl_params.dmd_resampling_factor = 50;       % dmd's MASK_PROPAGATE, resampling_factor
    ctrl_params.mask_resampling_factor = 50;      % mask's MASK_PROPAGATE, resampling_factor
%     ctrl_params.TEST_MODE_FLAG_ = 'test_mask';     % test mode, 'test_dmd', 'test_mask', 'test_all', 'non_test'
	ctrl_params.TEST_MODE_FLAG_ = 'non_test';     % test mode, 'test_dmd', 'test_mask', 'test_all', 'non_test'
%     ctrl_params.parallel_flag = 0;                % whether to use parallel calculation in CACTI
	ctrl_params.parallel_flag = 1;                % whether to use parallel calculation in CACTI
    ctrl_params.ideal_sensor_flag = 1;            % whether to assume sensor is ideal (sensor is the same as image plane)
end

% params save
if save_params_flag
    save([result_dir test_name '_ctrl_params.mat'], '-struct','ctrl_params')
%     save([config_dir 'ctrl_params/ctrl_params.mat'], '-struct','ctrl_params')
end
% print params settings to txt
print_params([result_dir 'ctrl_params_setting.txt'], ctrl_params)


%% system parameter (mm)
if load_sys_params_flag
    sys_params = load([result_dir test_name '_sys_params.mat']);
else
    % obj
    sys_params.obj_pos = -5;
    sys_params.obj_size = [512,512];
    sys_params.obj_pix_size = 1.15e-3;
    
    % lens
    sys_params.lens_pos = 0;
    sys_params.lens_f = 2.5;
    sys_params.lens_radius = 0.5;
    
    % dmd
    sys_params.dmd_pos = 0;
    sys_params.dmd_size = [100,100];
    sys_params.dmd_pix_size = 10.5e-3;
    
    % mask
    sys_params.mask_pos = 4.98;
%     sys_params.mask_pos = 5;
    sys_params.mask_size = [600,600];
    sys_params.mask_pix_size = 1e-3;
    
    % sensor
    sys_params.sensor_pos = [];
    sys_params.sensor_size = [512, 512];
    sys_params.sensor_pix_size = 1.15e-3;
    
    % calculate image position(sensor_pos)
    sensor_point = point_img([0,0,sys_params.obj_pos]', sys_params.lens_f);
    sys_params.sensor_pos = sensor_point(3);    
end

% params save
if save_params_flag
    save([result_dir test_name '_sys_params.mat'], '-struct','sys_params')
%     save([config_dir 'sys_params/sys_params.mat'], '-struct','sys_params')
end
% print params settings to txt
print_params([result_dir 'sys_params_setting.txt'], sys_params)


%% patterns setting
% load patterns from file
if load_patterns_flag
    config_dataset_dir = '.\config\patterns_generator\';
    
    if whiteboard_flag
        config_dataset_name = '512scale_whiteboard_config.mat';
    else
        config_dataset_name = '512scale_traffic_config.mat';
    end

    config_dataset = load([config_dataset_dir config_dataset_name]);

    obj = config_dataset.obj;
    dmd = config_dataset.dmd;
    mask = config_dataset.mask; 
    
    
else
    
    % obj
    obj = 255*ones(sys_params.obj_size);

    % dmd
    % dmd_t = 0.5;                            % dmd's transmission rate [standard]
    % dmd_t = 1;                            % dmd's transmission rate 
    % dmd = binary_mask(sys_params.dmd_size, dmd_t);     % dmd pattern
    % dmd = binary_mask(sys_params.dmd_size, 'fixed', 'up_half');  % dmd pattern
    dmd = binary_mask(sys_params.dmd_size, 'fixed', 'down_half');  % dmd pattern
    % dmd = load([config_dir 'dmd/dmd.mat']);


    %mask
    % mask_t = 0.5;                           % mask's transmission rate 
    % mask = binary_mask(dmd_size, dmd_t);    % binary mask pattern
    % mu=0.5; sigma=0.5; kernel_size=1;         % gaussian distribution params
    % mask = gray_mask(sys_params.mask_size, [mu,sigma,kernel_size]);  % gray mask pattern
    mask_t = 0.5;
    mask = gray_mask(sys_params.mask_size);      % mask pattern, uniform distribution
    % mask = load([config_dir 'mask/mask.mat']);
end



%% simulation dataset saving setting 
% dataset name
dataset_name = '512scale_traffic_cacti_simu';
cr = size(obj,3); % compressive ratio

if save_dvp_type
	% simu data for dvp reconstruction algorithm
    dvp_simudata_dir = '.\CI algorithm\dvp\dataset\';
	dvp_simudata_name = [dataset_name '.mat'];
    dvp_mask_name = [dvp_simudata_dir dataset_name '_syn_mask.mat'];
% 	syn_mask_name = '512scale_whiteboard_cacti.mat';
end

if save_dl_cacti_type
	% simu data for dl_cacti reconstruction algorithm
    dl_cacti_simudata_dir = '.\CI algorithm\DL-CACTI\dataset\';
	dl_cacti_simudata_name = [dataset_name '.mat'];
	dl_cacti_mask_name = ['mask_' dataset_name '_cr_' num2str(cr) '.mat'];
	dl_cacti_obj_name = ['obj_' dataset_name '_cr_' num2str(cr) '.mat'];
	dl_cacti_meas_name = ['meas_' dataset_name '_cr_' num2str(cr) '.mat'];
end


%% cacti imaging
% number of samples
obj_num = size(obj,3);
dmd_num = size(dmd,3); 
% mask_num = size(mask,3); % mask is fixed

% images preset
if ctrl_params.ideal_sensor_flag == 1
    image = zeros([sys_params.obj_size obj_num]);
else
    image = zeros([sys_params.sensor_size obj_num]);
end

% timer preset
timer = zeros(1,obj_num); % tic/toc is triggered in the cacti function

fprintf('\n* CACTI Imaging...\n')

for k = 1:obj_num
    fprintf('\n---- Image No. %d ----\n', k)
    
    % input
    obj_input = obj(:,:,k);
    dmd_input = dmd(:,:,k);
    
    % cacti imaging
%     image_output = cacti(obj_input, dmd_input, mask, sys_params);
    image_output = cacti(obj_input, dmd_input, mask, sys_params, ctrl_params);

    % output
    image(:,:,k) = image_output;
    
    % time
    timer(k) = toc;
end

% class convert: convert class of "image" to match the class of
% "image_output"
image = feval(class(image_output), image);

% time consumption
total_time = sum(timer, 'all');
fprintf('\n* Total time: %f\n',total_time)


%% save cacti result data
if save_data_flag==1
	% mkdir
	if ~isfolder(result_dir)
		mkdir(result_dir)
    end
    
    % save all data
	save([result_dir data_name]);
    % save images
    save([result_dir data_name], 'image');
end


%% save cacti simulation dataset
% save simu data for dvp reconstruction algorithm
if save_dvp_type==1
    % save cacti simulation dataset for dvp SCI algorithm
    % save and print system parameters to dvp dataset file
    save([dvp_simudata_dir dvp_simudata_name(1:end-4) '_sys_params.mat'], '-struct','sys_params')
    print_params([dvp_simudata_dir dvp_simudata_name(1:end-4) '_sys_params_setting.txt'], sys_params)
    
    % copy image
    dvp_image = image;
    
    % convert image to "single" format and rotate to mask the scene upright
    dvp_image = single(dvp_image);
    dvp_image = rot90(dvp_image,2);
    
    if whiteboard_flag
        % the obj is a whiteboard, used for synthetic mask calibration      
        mask_bayer = dvp_image./single(max(obj, [], 'all'));  % normalized to 0-1
        % save the synthetic mask
        save(dvp_mask_name, 'mask_bayer'); 
    else
        % the obj is the natural scene
        orig_bayer = uint8(obj);
        meas_bayer = single(sum(dvp_image, 3));
    end
    
    % save
    if isfile([dvp_simudata_dir dvp_simudata_name]) 
        if whiteboard_flag
            save([dvp_simudata_dir dvp_simudata_name], 'mask_bayer', '-append');
        else
            save([dvp_simudata_dir dvp_simudata_name], 'orig_bayer', 'meas_bayer', '-append');
        end
    else
        if whiteboard_flag
            save([dvp_simudata_dir dvp_simudata_name], 'mask_bayer', '-v7.3');
        else
            save([dvp_simudata_dir dvp_simudata_name], 'orig_bayer', 'meas_bayer', '-v7.3');
        end        
    end
        
end

% save simu data for dl_cacti reconstruction algorithm
if save_dl_cacti_type
    % save cacti simulation dataset for dl_cacti reconstruction algorithm
    % save and print system parameters to dvp dataset file
    save([dl_cacti_simudata_dir 'sys_params_' dataset_name '.mat'], '-struct','sys_params')
    % print_params([dl_cacti_simudata_dir 'sys_params_setting_' dataset_name '.txt'], sys_params)
    
    % copy image
    dl_cacti_image = image;
    
    % convert image to "single" format and rotate to mask the scene upright
    dl_cacti_image = single(dl_cacti_image);
    dl_cacti_image = rot90(dl_cacti_image,2);
    
    if whiteboard_flag
        % the obj is a whiteboard, used for synthetic mask calibration
        mask = dl_cacti_image./single(max(obj, [], 'all'));  % normalized to 0-1
        % save the synthetic mask
        save([dl_cacti_simudata_dir dl_cacti_mask_name], 'mask'); 
    else
        % the obj is the natural scene
		orig = uint8(obj);
        meas = single(sum(dl_cacti_image, 3));
        % save the object
        save([dl_cacti_simudata_dir dl_cacti_obj_name], 'obj');			
        % save the measurenment
        save([dl_cacti_simudata_dir dl_cacti_meas_name], 'meas');
	end
end


%% illustration
if show_image_flag
    for k = 1:obj_num
        figure
        imshow( image(:,:,k),[])
        title('rescaled result')
    end
end


