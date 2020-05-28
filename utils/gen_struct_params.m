function parameters = gen_struct_params(path)
%GEN_STRUCT_PARAMS	generate struct (array) type parameters and return it.
%Save is optional.
% 
%   Input:
%   --------
%   - path: parameter saving path, str,optional, if not assigned,parameters
%       will not be saved to local
% 
%   Output:
%   - params: return parameters as a struct (array)
%   - Non-output: save the parameters to local, "path" is needed.
% 
%   Note:
%   --------
%   This function is suitable for generating struct array, with structed
%   data, like a table containging a series of items with the same fields. 
%   For struct, GEN_PARAMS is a better choice.
% 
%   Example£º
%   --------
%   - params = gen_struct_params('../_trash/gen_params.mat'), return and save
%   - gen_params('../_trash/gen_params.mat','not_load'), save only
%   - params = gen_struct_params(), return only
% 
%   See also£º
%   --------
%   GEN_PARAMS
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-22
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-22
%   
%   Copyright 2020 Zhihong Zhang


%% parameter setting
% set sys params
% obj
params(1).name = 'obj';
params(1).pos = [];
params(1).size = [];
params(1).pattern = [];
params(1).spot = [];
params(1).element_sz = [];

% lens
params(2).name = 'lens';
params(2).pos = [];
params(2).size = [];
params(2).spot = [];


% dmd
params(3).name = 'dmd';
params(3).pos = [];
params(3).size = [];    
params(3).pattern = [];
params(3).spot = [];
params(3).element_sz = [];

% mask
params(4).name = 'mask';
params(4).pos = [];
params(4).size = [];    
params(4).pattern = [];
params(4).spot = [];
params(4).element_sz = [];

 % sensor
params(5).name = 'sensor';
params(5).pos = [];
params(5).size = [];    
params(5).pattern = [];
params(5).spot = [];
params(5).element_sz = [];

%% parameter passing
if nargin == 1 && nargout==1
    save(path, 'params');
    disp(['params saved to: ' path]);
    parameters = params;
elseif nargin == 0 && nargout==1
    parameters = params;
elseif nargin == 1 && nargout==0
    save(path, 'params');
    disp(['params saved to: ' path]);
end
end