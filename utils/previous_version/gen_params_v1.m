function params = gen_params_v1(path_)
%GEN_PARAMS_V1 generate parameters and return them together as a struct type variable, or
%directly load them to the workspace [v1.0]
% 
%   Input:
%   --------
%   - non-input: generate params, but won't save it
%   - path: parameter saving path, str
% 
%   Output:
%   --------
%   - non-output: directly load parameters to the workspace
%   - params: return parameters as a struct type variable, struct
% 
%   Note;
%   --------
%   This version is deprecated, and it can be replaced with gen_params
% 
%   See also;
%   --------
%   GEN_PARAMS
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-15
%   Last Modified:  Zhihong Zhang, 2020-03-21
%               
%   Copyright 2020 Zhihong Zhang


%% inputs setting
% if "path" is not assigned, use tmp data path as default
if nargin == 0
    path_ = './tmp_params_data.mat';
end

%% parameter setting
% objective
obj_pos = -5;           % objective distance -50 mm
obj_size = [200,200]; % pixel number
obj_pix_size = 5e-3;     % 0.5um sampling(representing) size
obj = 255*ones(obj_size);  % obj pattern, white board

% lens
lens_f = 2.5;  % focal length
lens_radius = 0.5; % lens radius

% DMD
% dmd_pos = 0;
dmd_pos = 2.5; % [standard]
dmd_size = [100,100]; % pixel number
dmd_pix_size = 10e-3;  % 10um pixel size
dmd_t = 0.5; % dmd's transmission rate [standard]
% dmd_t = 1; % dmd's transmission rate 
dmd = binary_mask(dmd_size, dmd_t);  % dmd pattern
% load('../config/dmd.mat')

% mask
% mask_pos = 5;
mask_pos = 4.5;% [standard]
mask_size = [200,200]; % pixel number
mask_pix_size = 5e-3;   % 10um pixel size
mask_t = 0.5; % mask's transmission rate [standard]
mask = binary_mask(mask_size,mask_t);  % mask pattern
% load('../config/mask.mat')

% sensor
sensor_pos = [];
sensor_size = [100, 100]; % pixel number
sensor_pix_size = [10e-3, 10e-3]'; % 5um pixel size
sensor = zeros(sensor_size); % image pattern
large_sensor_flag = 1; % assuming sensor is large enough

%% parameter passing
% save params
% save params (save all variables except those ended with "_", which are
% function¡¯s inputs)
save(path_, '-regexp', '\w*[^_]\>');

% load params to workspace
if nargout == 1
    % return a struct type params
    params = load(path_);
else
    % directly load to workspace
    assignin('base', 'params_path',path_);
    evalin('base', 'load(params_path)');
    if nargin == 0
        fprintf('params generated and loaded!\n');
    else
        fprintf('params saved and loaded!\n');
    end
end

% delete tmp data
if nargin == 0
    delete(path_);
end
end