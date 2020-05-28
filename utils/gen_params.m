function params = gen_params(path_, mode_)
%GEN_PARAMS	generate parameters and return them together as a struct type..
% variable,directly load them to the workspace or save them to local [v2.0]
% 
%   Input:
%   --------
%   - path: parameter saving path, str,optional, if not assigned,parameters
%       will not be saved to local
% 
%   - mode: whether load parameters to workspace, take effect only when
%       nargout==0(i.e. function doesn't return params), and path is assigned; 
%       optional str = {'not_load', 'load'}, default = 'load'
% 
%   Output:
%   --------
%   - non-output:
%       * if path is assigned and mode=={'not_load'}, only save params but
%         don't load them to workspace
%       * if path is assigned and mode=={'load'} or is not assigned,
%       load params and save them to local
%       * if path is not assigned, only load parameters to the workspace,
%       but don't save them to local     
% 
%   - params: return parameters as a struct
% 
%   Example£º
%   --------
%   - gen_params('../_trash/gen_params.mat')¡ª¡ªload and save
%   - gen_params('../_trash/gen_params.mat','not_load')¡ª¡ªsave only
%   - gen_params()¡ª¡ªload only
%   - params = gen_params('../_trash/gen_params.mat')¡ª¡ªreturn a struct
%     and save
%   - params = gen_params()¡ª¡ªreturn a struct only
% 
%   See also£º
%   --------
%   GEN_PARAMS_V1
% 
% 
%   Log:
%   --------
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-18
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang


%% inputs setting
% if "path" is not assigned, use tmp data path as default
if nargin == 0
    path_ = './tmp_params_data.mat';
end
% if mode is not assigned, use default value 'load'
if nargin <2
    mode_ = 'load';
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
% save params (save all variables except those ended with "_", which are
% function¡¯s inputs)
save(path_, '-regexp', '\w*[^_]\>');

% load or return params
if nargout == 1
    % return params as a struct
    params = load(path_);
elseif strcmp(mode_, 'load')% nargout == 0
    % load params directly to the workspace, don't save params
    assignin('base', 'params_path', path_);
    evalin('base', 'load(params_path)');
    disp('** params generated and loaded! **');
end
    
% delete tmp data
if nargin == 0
    delete(path_);
else
    disp(['params saved to: ' path_])
end
end