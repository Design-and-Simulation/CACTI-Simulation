function  light_drawing(params, varargin)
%LIGHT_DRAWING illustrating the optical path using given optical system
%parameters and light propagating infomation
% 
%   Input:
%   --------
%   - params: system parameters and light propagating information, nx1
%   struct,with the following fields (For 'draw_light' mode, only spot
%   field is needed)
%       name        - the optical componment's  name, str
%       pos         - the optical componment's  position (z-coordinate) , float
%       size        - the optical componment's size, int 2D vector(matrix's
%                     size) or float scalar(physical radius diameter or length)
%       pattern     - the optical element's  patttern, 2D matrix
%       spot        - information about the light spot on the element,
%                     n¡Á1 float vector, [spot_center; spot_r]
%       element_sz  - the size of the pattern's elements
% 
%   - varargin, contain the following optional arguements
%       unit_len[k]     - the unit length of figure coordinate system 
%                         [key, str, '-unit_len_k']
%       unit_len[v]     - the unit length of figure coordinate system 
%                         [value, float, default = smallest element's size]
%       draw_what       - {'draw_sys', 'draw_light', 'all'}
%       dimension       - {'2D', '3D'}
%       fig_handle      - figure handle, optional, default = generate a new axes;
%       
%   Output:
%   --------
%   the optical path figure
% 
%   Note: 
%   --------
%	- physical 3D coordinate system is:x-up, y-out, z-right
% 
%	- the optical element(obj, dmd, mask,..) is placed with their matrix mn 
%   coord-sys oppsite to physical 2D xy coord-sys. (i.e., observing from 
%   anti optical axis direction, the [1 1] element of matrix is at the left
%   up corner of the physical object(obj, dmd or mask)
% 
%   See also:
%   --------
% 
%   Log:
%   --------
%   * implement 2D drawing, 20200321
%   + implement 3D drawing
% 
%   Info:
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-21
%   Last Modified:  Zhihong Zhang, 2020-03-21
%               
%   Copyright 2020 Zhihong Zhang

%% input setting
var_in = varargin;

% default value
% the minimal element_sz of all optical components
% [note: unused parameter, left for possible extension in the future
unit_len = min(cell2mat({params.element_sz})); 
draw_what = 'draw_all';
dimension = '2D';
fig_handle_ = []; 

% input assignment
while numel(var_in)>0
    if ischar(var_in{1}) || isstring(var_in{1})
        var_item = var_in{1}; 
    elseif isa(var_in{1},'matlab.ui.Figure')
%     elseif isa(var_in{1},'matlab.graphics.axis.Axes')
        fig_handle_ = var_in{1};
        var_in = var_in(2:end);
        continue
    else
        error('input error %s', string(var_in{1}))
    end
    switch true
        case strcmp(var_item,'-unit_len_k')
            if isnumeric(var_in{2}) && isscalar(var_in{2})
                unit_len = var_in{2};
                var_in = var_in(2:end);
            else
                error('input error %s', string(var_in{1}))
            end
        case strcmp(var_item, 'draw_sys') ||  strcmp(var_item, 'draw_light')...
				|| strcmp(var_item, 'draw_all')
            draw_what = var_item;
        case strcmp(var_item, '2D') ||  strcmp(var_item, '3D')
            dimension = var_item;
        otherwise
            error('unknown input parameter - %s', var_item);
    end
    var_in = var_in(2:end);
end
   
%% main function
% get fig handle
if isempty(fig_handle_)
    fig_handle_ = figure();
end
figure(fig_handle_);

hold on
if strcmp(draw_what, 'draw_light')
    draw_light(params);
elseif strcmp(draw_what, 'draw_sys')
    draw_sys2D(params);
elseif strcmp(draw_what, 'draw_all')
    draw_sys2D(params);
    draw_light(params);
end
hold off

% beautify the figure
cur_axes = gca;
if ~(strcmp(cur_axes.XAxisLocation, 'origin') && strcmp(cur_axes.YAxisLocation, 'origin'))
    origin2center(cur_axes);
end

%% sub function
% draw light
function draw_light(params)
%DRAW_LIGHT    draw the light path, containing 3 lights, i.e. two boundary
%lines and one center line
% Input: systems parameters, struct

% get spot info
comp_num = numel(params);
spot = cell2mat({params.spot});
if size(spot,2) ~= comp_num
    error('"spot" paramter is lost for some componments')
end
spot_X = spot(1,:);
spot_Y = spot(2,:);
spot_Z = spot(3,:);
spot_R = spot(4,:);

% calculate line points
% center line is (spot_Z, spot_X)
Line_center_X = spot_X;
Line_center_Z = spot_Z;
% above line
Line_up_X = spot_X + spot_R;
Line_up_Z = spot_Z;
% below line
Line_down_X = spot_X - spot_R;
Line_down_Z = spot_Z;

% draw line: be careful about the axes, we draw the z-x plane
color = [rand,rand,rand];
hold on
plot(Line_up_Z, Line_up_X, '-', Line_center_Z, Line_center_X, '-',...
    Line_down_Z,Line_down_X,'-', 'LineWidth', 1.5, 'color', color)
hold off
end


% draw system
function draw_sys2D(params)
%DRAW_SYS2D    draw the optical system, i.e. all the components.
% Input: systems parameters, struct

% get componment info
comp_num = numel(params);
comp_name = string({params.name});
comp_pos = cell2mat({params.pos});
if numel(comp_pos) ~= comp_num
    error('"pos" paramter is lost for some componments')
end

% calculate physical size
comp_phy_sz = ones(comp_num, 1);
for k = 1:comp_num
    if numel(params(k).size) == 1 && ~isempty(params(k).size)
        comp_phy_sz(k) = params(k).size;
    elseif numel(params(k).size) == 2 && ~isempty(params(k).size)
        if ~isempty(params(k).element_sz)
            % z-x coord-sys, anti-light view, so use size(1), i.e. num of rows for component physical height
            comp_phy_sz(k) =params(k).size(1)*params(k).element_sz; 
        else
            error('"element_sz" paramter is lost for - %s', comp_name(k))
        end
    else
        error('"size" paramter is incorrect for -  %s', comp_name(k))
    end
end

% calculate relative size
% comp_rel_sz = comp_phy_sz./unit_len; 

% calcultate component lines and draw
hold on
for k = 1:comp_num
    endpoints_Z = [comp_pos(k), comp_pos(k)];
    endpoints_X = [-comp_phy_sz(k)/2, comp_phy_sz(k)/2];
    
    % setting plot parameters
    switch comp_name(k)
        case "obj"
            line_str = 'c-';
            line_width = 3;
        case "lens"
            line_str = 'k-';
            line_width = 1;
            endpoints_Z = endpoints_Z + endpoints_X./5;
            endpoints_Z = linspace(endpoints_Z(1),endpoints_Z(2),9);
            endpoints_X = endpoints_X(2).*sqrt(1-(endpoints_Z./endpoints_Z(1)).^2);
            endpoints_Z = [endpoints_Z endpoints_Z];
            endpoints_X = [endpoints_X -endpoints_X];
            
        case "dmd"
             line_str = 'k-.';
            line_width = 3;     
        case "mask"
             line_str = 'k-.';
            line_width = 3;     
        case "sensor"
             line_str = 'c--';
            line_width = 3;
        otherwise
             line_str = 'k-';
            line_width = 3;
    end
    
    % draw line
    plot(endpoints_Z, endpoints_X, line_str, 'LineWidth',line_width);
    text(comp_pos(k),-comp_phy_sz(k)/2*1.1,comp_name(k),'color','k','FontSize',13)
   
end
xlabel('z');ylabel('x');
hold off
end


end