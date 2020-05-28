% function [center,radius] = blur_spot(obj, f, D, pos)
function varargout = blur_spot(obj, f, D, pos)
%BLUR_SPOT calculate the center coordinates and radius of the defocus
%circle formed from a thin lens
% 
%   Input:
%   --------
%   - obj: coordinates of object point(or point array) in 3D space, float 
%           3x1 vector(or M*N*3 float array, i.e., M*N points' coord array)
%   - f: focal length of a thin lens, float scalar
% 
%   - D: thin lens diameter
% 
%   - pos: spot's position
% 
%   Output:
%   - center: blur spot center coordinates
% 
%   - radius: blur spot center radius
% 
%   Log:
%   --------
%   * take boundary cases into consideration [20200316]
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang


nout = nargout;

% - calculating spot center
if ndims(obj) == 3
    % for points array, convert from "each channel is a x/y/z coord matrix" to
    % "each row is a x/y/z coord matrix"
    obj = permute(obj, [3 2 1]);
    z = obj(3,1,1); % all points have the same z-coord
    center = pos/z*obj;
    center = permute(center, [3 2 1]); % transpose to origin shape
elseif isvector(obj)
    % P = [-1 0 0; 0 -1 0; 0 0 -1];
    z = obj(3);	  % object posance 
    
    % center = pos/(-z)*P*obj;
    center = pos/z*obj;
end

varargout{1} = center;

% - calculating spot radius
if nout==2
    if pos==0
        radius = D/2;
    else
        radius = D *pos/2*abs(1/f - 1/(-z) - 1/pos);
    end
    varargout{2} = radius;
end
