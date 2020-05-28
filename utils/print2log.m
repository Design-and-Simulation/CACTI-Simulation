function print2log(path, params,mode)
%PRINT2LOG execute commands listed in "params" and print the outputs to a txt file
% 
% 
%	Input:
%   --------
%   - path:     local path to save the printed params
%   - params:   params that need to printed
%   - mode:     write mode, choose to append ('a') or overwrite ('w') the 
%               log file if it's not empty
%               file, str,{'a', 'w'}, default='a'
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang


if nargin < 3
    mode = 'w';
end

% if 'write' mode, clear current log file
if strcmp(mode, 'w') && isfile(path)
    delete(path)
end

diary(path); % assign a file to export command window's contents


% start export
diary on

fprintf(['\n\n----------  ', datestr(now), '  ----------\n\n']);

for command = params
    evalin('base', command);
end

%stop export
diary off

end