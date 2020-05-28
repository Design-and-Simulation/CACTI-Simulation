function b = clamp(a, lower_limit, higher_limit)
%CLAMP limits the values of given array to a certain boundary
% 
%   Input
%   --------
%   a           :    input array, numeric
%   lower_limit :    lower boundary, numeric
%   higher_limit:    lower boundary, numeric
% 
%   Output
%   --------
%   b           :    output array with value clamped, numeric
% 

if nargin<3
    higher_limit = inf;
end

b = max(a, lower_limit);
b = min(b, higher_limit);

end

