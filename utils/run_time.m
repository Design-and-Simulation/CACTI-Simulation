function time = run_time(command_strs, run_times)
%TEST_TIME test command's run time
% 
%   Input:
%   --------
%  command_strs: command that need to test the run time, string vector
%  run_times: command's repeate times
% 
%   Output:
%   --------
%  time: command's run time(s), float

if nargin<2
    run_times = 1;
end

% start timer
tic

% run commands
for k=1:run_times
    for command  = command_strs
        evalin('base', command_strs)
    end
end

% stop timer
time = toc;

fprintf('\n time consumed: %f s\n', time);

end

