function print_params(path, params, regexp_flag)
%PRINT_PARAMS export the parameter setting to a txt file(only char, 
% string scalar, and numeric vector are supported)
%   Input:
%   --------
%   - path: parameter setting export path, str(.txt)
% 
%   - params: system parameters, struct, or string vector
% 
%   - regexp_flag: only take effect when params is a string, means params
%     are specified with regular expressionsm, optional str = {'-variable', '-regexp'}, default='-variable'
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang


if nargin<3
    regexp_flag = '-variable';
end

disp(['* print parameter settings to ' char(path) '...'])
% get params name
if isa(params, 'struct') % struct type
    % get field names
    params_names = string(fieldnames(params));
    params_names = params_names';
elseif isa(params, 'string') % string type
    
    if  strcmp(regexp_flag, '-regexp')
        vars = evalin('base','whos');
        vars_name = string({vars.name});
        params_names = regexp(vars_name, join(params, '|'), 'match');
        params_names(cellfun(@isempty,params_names))=[];
        params_names = string(params_names);
    elseif strcmp(regexp_flag, '-variable')
        params_names = params;
    else
        disp('third arugment error');
    end  
else
    disp('secord arugment error');
end

% write to txt
[fp, m] = fopen(path,'w');
if fp==-1
    error(m)
end
formatSpec = '%s: %s\n\n';
for name = params_names
    % get param value
    if isa(params, 'struct')
        name_value  = params.(name);
    elseif isa(params, 'string')
        name_value = evalin('base', name);
    end
    % write to file
    % for string or char type value
    if isa(name_value, 'string') || isa(name_value, 'char')
        fprintf(fp, formatSpec, [name, name_value]);
     % for numeric type scalar or vector
    elseif isa(name_value, 'numeric') && isvector(name_value)
        name_value_str = join(string(name_value), 'x');
        fprintf(fp, formatSpec, [name, name_value_str]);
    else
        disp(['exist not supported data type:' char(name) ',skiped']);
    end
end
    
fclose(fp);
end
            