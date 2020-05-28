function point_img = point_img(obj, f)
%POINT_IMG calculate the image coordinates of an object point propagating through
%a thin lens
% 
%   Input:
%   --------
%   - obj: coordinates of object point in 3D space, float 3x1 vector
%   - f: focal length of a thin lens, float scalar
% 
%   Output:
%   --------
%   - img_point: coordinates of the obj
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Copyright 2020 Zhihong Zhang


P = [-1 0 0; 0 -1 0; 0 0 -1];
z = obj(3);	  % object distance 
point_img = f/(-z-f)*P*obj;

end

