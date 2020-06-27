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
save_res_flag = 0;          % save result data to the local directory
save_params_flag = 0;       % save data to the local result file
load_ctrl_params_flag = 0;  % 0-setting; 1-load from file; 2-use default
load_sys_params_flag = 0;
show_image_flag = 1;
load_patterns_flag = 0;           % load patterns for obj, dmd and mask
save_compact_type = 0;            % save simulation dataset into one '.mat' file
save_sep_type = 1;				  % save simulation dataset into separate files (mask.mat, orig.mat, meas.mat)
whiteboard_flag = 1;			  % whether the obj is a whiteboard

% paths and names
% test_name = 'cacti_256_10f';
test_name = 'center_circle_dmd'; % dmd design test
config_dir = './config/';
result_dir = './result/tmp/';
data_name = test_name;

% make dir
if save_res_flag || save_params_flag
	if ~isfolder(result_dir)
		mkdir(result_dir)
	end
end



%% control parameter
if load_ctrl_params_flag
    ctrl_params = load([result_dir test_name '_ctrl_params.mat']);
else
    ctrl_params.show_compare_flag = 0;            % whether to show obj, dmd, mask and image
    ctrl_params.drawing_sys_flag = 1;             % whether to draw optical system
    ctrl_params.dmd_resampling_factor = 50;       % dmd's MASK_PROPAGATE, resampling_factor
    ctrl_params.mask_resampling_factor = 50;      % mask's MASK_PROPAGATE, resampling_factor
%     ctrl_params.TEST_MODE_FLAG_ = 'test_mask';     % test mode, 'test_dmd', 'test_mask', 'test_all', 'non_test'
	ctrl_params.TEST_MODE_FLAG_ = 'non_test';     % test mode, 'test_dmd', 'test_mask', 'test_all', 'non_test'
%     ctrl_params.parallel_flag = 0;                % whether to use parallel calculation in CACTI
	ctrl_params.parallel_flag = 1;                % whether to use parallel calculation in CACTI
    ctrl_params.ideal_sensor_flag = 1;            % whether to assume sensor is ideal (sensor is the same as image plane)
end

% ctrl_params save
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
    sys_params.obj_size = [256 256];
    sys_params.obj_pix_size = 1.15e-3;
    
    % lens
    sys_params.lens_pos = 0;
    sys_params.lens_f = 2.5;
    sys_params.lens_radius = 0.25;
    
    % dmd
    sys_params.dmd_pos = 0;
    sys_params.dmd_size = [256,256];
    sys_params.dmd_pix_size = 2.05e-3;
    
    % mask
    sys_params.mask_pos = 4.95;
%     sys_params.mask_pos = 5;
    sys_params.mask_size = [300,300];
    sys_params.mask_pix_size = 1e-3;
    
    % sensor
    sys_params.sensor_pos = [];
    sys_params.sensor_size = [256, 256];
    sys_params.sensor_pix_size = 1.15e-3;
    
    % calculate image position(sensor_pos)
    sensor_point = point_img([0,0,sys_params.obj_pos]', sys_params.lens_f);
    sys_params.sensor_pos = sensor_point(3);    
end

% sys_params save
if save_params_flag
    save([result_dir test_name '_sys_params.mat'], '-struct','sys_params')
%     save([config_dir 'sys_params/sys_params.mat'], '-struct','sys_params')
end
% print params settings to txt
print_params([result_dir 'sys_params_setting.txt'], sys_params)


%% patterns setting
if load_patterns_flag
	% load patterns from file
    config_dataset_dir = '.\config\patterns_generator\';
    
    if whiteboard_flag
        config_dataset_name = 'whiteboard_config_256_10f.mat';
    else
        config_dataset_name = 'obj_config_256_10f.mat';
    end

    config_dataset = load([config_dataset_dir config_dataset_name]);

    obj = config_dataset.obj;
    dmd = config_dataset.dmd;
    mask = config_dataset.mask; 
    
    
else
	% generate patterns
    % obj
	obj_num = 10;
    obj = 255*ones([sys_params.obj_size obj_num]);

	
    % dmd
	% generate dmd
    % dmd_t = 0.5;                            % dmd's transmission rate [standard]
    % dmd_t = 1;                            % dmd's transmission rate 
    % dmd = binary_mask(sys_params.dmd_size, dmd_t);     % dmd pattern
    % dmd = binary_mask(sys_params.dmd_size, 'fixed', 'up_half');  % dmd pattern
    % dmd = binary_mask(sys_params.dmd_size, 'fixed', 'down_half');  % dmd pattern
	
	% load dmd
	dmd_name = 'center_circle_dmd_256_10f';
    load([config_dir 'dmd/designed_dmd/' dmd_name '.mat']);
	

    % mask
	% generate mask
    % mask_t = 0.5;                           % mask's transmission rate 
    % mask = binary_mask(dmd_size, dmd_t);    % binary mask pattern
    % mu=0.5; sigma=0.5; kernel_size=1;         % gaussian distribution params
    % mask = gray_mask(sys_params.mask_size, [mu,sigma,kernel_size]);  % gray mask pattern
    % mask_t = 0.5;
    % mask = gray_mask(sys_params.mask_size);      % mask pattern, uniform distribution
	
	% load mask
	mask_name = 'mask_300';
    load([config_dir 'mask/' mask_name '.mat']);
end

% append sys_params
sys_params.dmd = dmd;
sys_params.mask = mask;
sys_params.obj = obj;


%% cacti simulation dataset saving setting 
% dataset info
cr = size(obj,3); % compressive ratio
size_scale = size(obj,1); % image size

time_now = datestr(now,'yyyy-mm-dd_HH-MM-SS');

% cacti_dataset_dir = ['.\dataset\data\cacti_' num2str(size_scale) '_' num2str(cr) 'f_' time_now '\']; % cacti dataset saving dir
% cacti_name = ['cacti_%s_' num2str(size_scale) '_' num2str(cr) 'f']; % general name for simulated data

cacti_dataset_dir = ['.\dataset\mask\cacti_' test_name '_'  num2str(size_scale) '_' num2str(cr) 'f' '\']; % cacti mask dataset saving dir
cacti_name = ['cacti_%s_' test_name '_' num2str(size_scale) '_' num2str(cr) 'f']; % general name for simulated mask data


% dataset name
if save_compact_type
	% save simulation dataset into one '.mat' file
	cacti_dataset_name = sprintf([cacti_name '.mat'], 'data');
%     cacti_mask_name = [cacti_name '_syn_mask.mat'];
end

if save_sep_type
	% save simulation dataset into separate files (mask.mat, orig.mat, meas.mat)
	cacti_mask_name = sprintf([cacti_name '.mat'], 'mask');
	cacti_orig_name = sprintf([cacti_name '.mat'], 'orig');
	cacti_meas_name = sprintf([cacti_name '.mat'], 'meas');
end


%%  cacti system drawing
if ctrl_params.drawing_sys_flag == 1
	cacti_drawing(sys_params)
end

%% cacti imaging
% number of samples
obj_num = size(obj,3);
% dmd_num = size(dmd,3); 
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
if save_res_flag==1
	% mkdir
	if ~isfolder(result_dir)
		mkdir(result_dir)
	end
    
    % save all data
	save([result_dir data_name '_data.mat']);
    % save images
    save([result_dir data_name '_image.mat'], 'image');
end


%% save cacti simulation dataset
if save_compact_type==1	|| save_sep_type==1	
	% 1. mkdir
	if ~isfolder(cacti_dataset_dir)
		mkdir(cacti_dataset_dir)
	end

	% 2. save and print parameters to dataset file
	save([cacti_dataset_dir sprintf([cacti_name '.mat'], 'sys_params')], '-struct','sys_params')
	print_params([cacti_dataset_dir sprintf([cacti_name '.txt'], 'sys_params')], sys_params)

	% 3. data arrange
	% copy image
	sensor_image = image;

	% convert image to "single" format and rotate to mask the scene upright
	sensor_image = single(sensor_image);
	sensor_image = rot90(sensor_image,2);

	% data assign
	if whiteboard_flag
		% the obj is a whiteboard, used for synthetic mask calibration
		% ! the variable 'mask' is reused for synthtic mask here
	% 	mask = sensor_image;  % non-normalized
		mask = sensor_image./single(max(obj, [], 'all'));  % normalized to 0-1
	else
		% the obj is the natural scene
		orig = uint8(obj);
		meas = single(sum(sensor_image, 3));
	end
end

% 4. save dataset
% save simulation dataset into one '.mat' file
if save_compact_type==1	
	if isfile([cacti_dataset_dir cacti_dataset_name]) 
		if whiteboard_flag
			save([cacti_dataset_dir cacti_dataset_name], 'mask', '-append');
		else
			save([cacti_dataset_dir cacti_dataset_name], 'orig', 'meas', '-append');
		end
	else
		if whiteboard_flag
			save([cacti_dataset_dir cacti_dataset_name], 'mask', '-v7.3');
		else
			save([cacti_dataset_dir cacti_dataset_name], 'orig', 'meas', '-v7.3');
		end        
	end
end

% the obj is the natural scene
if save_sep_type==1	
	if whiteboard_flag		
		% save the mask
		save([cacti_dataset_dir cacti_mask_name], 'mask');
	else
		% save the object
		save([cacti_dataset_dir cacti_orig_name], 'orig');	
		% save the measurenment
		save([cacti_dataset_dir cacti_meas_name], 'meas');	
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


