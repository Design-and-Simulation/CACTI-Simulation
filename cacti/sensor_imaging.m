function sensor_img = sensor_imaging(varargin)
%SENSOR_IMAGING simulation of the imaging process of the sensor
%   Assume that, the sensor is located at the image plane.
%   Input
%   --------
%   tmp_image           the image formed by preceding optical components,
%                       2D numeric matrix
% 
%   tmp_img_pixel_size  the physical size of tmp image's pixel (element),
%                       float scalar
% 
%   sensor_size         sensor size (pixel numbers), 2*1 int vector
% 
%   sensor_pixel_size   the physical size of sensor's pixel (element),
%                       float scalar
% 
%   AD                  analog-digital-convert params, [quan_level, sat_value,  
%                       save_class], "quan_level" is quantitative level, int scalar,
%                       default=ceil(sat_value), i.e. default quantitative 
%                       interval is "1". "sat_value" is saturation value, i.e.
%                       the maximum limitation of ADC, float scalar
%                       default=max(tmp_img);
%	save_class			output class of sensor image. 'double'|'single'|'uint8'|
%						'uint16'|'default', optional, default = 'default'				
% 
%   Input
%   --------
%   sensor_img          image captured by sensor
%   
%   Input
%   --------
%   If only 
% 
%   Info:
%   --------
%   Created:    Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-31
%   Last Modified:   Zhihong Zhang, 2020-03-31
%
%   Copyright (c) 2020 Zhihong Zhang


% input setting
narginchk(1,6)

tmp_img = varargin{1};
if nargin==1
	save_class = 'default';
elseif nargin==2
	save_class = varargin{2};
elseif nargin==5 || nargin==6
	tmp_img_pixel_size = varargin{2};
	sensor_size = varargin{3};
	sensor_pixel_size = varargin{4};
	AD = varargin{5};
	if nargin==6
		save_class = varargin{6};
	else
		save_class = 'default';
	end
else
	error('error input')
end




% tmp image size
tmp_img_size = size(tmp_img);
tmp_img = double(tmp_img);
% image brightness record
img_energy_scale = max(tmp_img, [],'all');
% geometrical optics, the image is inverted
if img_energy_scale>0
	tmp_img = rot90(tmp_img./img_energy_scale,2);
else
	% black image (all 0)
	tmp_img = rot90(tmp_img,2);
end

    
% ideal sensor sampling
if nargin<3 || (all(tmp_img_size==sensor_size, 'all') && all(tmp_img_pixel_size==sensor_pixel_size, 'all'))
    % ideal sensor sampling: the sensor pixels are exactly matched with the object
    % points' pixels. And ADC is not considerated.
    sensor_img = tmp_img*img_energy_scale;
    sensor_max_value = max(sensor_img, [],'all');
	if strcmp(save_class, 'default')
		% default save class
		if sensor_max_value < 256
			save_class = 'uint8';
		elseif sensor_max_value < 65536
			save_class = 'uint16';
		else
			save_class = 'double';
		end
	end
    sensor_img = feval(save_class, sensor_img);    
    return
end


% sensor size
if ~isrow(sensor_size)
    sensor_size = sensor_size';
end

%% physical size comparsion and crop 
tmp_img_sz = tmp_img_size*tmp_img_pixel_size;
sensor_sz = sensor_size*sensor_pixel_size;

% tmp_img crop 
tmp_delta_sz =max(tmp_img_sz - sensor_sz, 0);
tmp_delta_pix_num =  round(tmp_delta_sz./tmp_img_pixel_size);
image_captured_pix_num = min(round(sensor_sz./tmp_img_pixel_size), tmp_img_size);

tmp_img_captured_idx = [floor(tmp_delta_pix_num./2)+1, floor(tmp_delta_pix_num./2)+image_captured_pix_num];
tmp_image_captured = tmp_img(tmp_img_captured_idx(1):tmp_img_captured_idx(3),tmp_img_captured_idx(2):tmp_img_captured_idx(4));


% sensor area
sensor_delta_sz =max(sensor_sz - tmp_img_sz, 0);
sensor_delta_pix_num =  round(sensor_delta_sz./sensor_pixel_size);
sensor_used_pix_num = min(round(tmp_img_sz./sensor_pixel_size), sensor_size);

sensor_img_used_idx = [floor(sensor_delta_pix_num./2)+1, floor(sensor_delta_pix_num./2)+sensor_used_pix_num];

%% sensor sampling
sensor_img = zeros(sensor_size);

if all(image_captured_pix_num == sensor_used_pix_num,'all')
    sampled_tmp_img = tmp_image_captured;
elseif all(image_captured_pix_num > sensor_used_pix_num,'all')
    % tmp image has more pixels than sensor, combine tmp image pixels to
    % match one sensor pixel
    block_sz = ceil(image_captured_pix_num./sensor_used_pix_num);
    resize_tmp_img_sz = block_sz.*sensor_used_pix_num; 
    resize_tmp_img = imresize(tmp_image_captured,resize_tmp_img_sz);

    sampled_tmp_img = blkcolfun(resize_tmp_img, block_sz, @sum, 'distinct');
else
    % tmp image has less pixels than sensor, divide tmp image pixel into
    % several sub pixels to match one sensor pixel
    sz_ratio = sensor_used_pix_num./image_captured_pix_num;
    resize__tmp_img = kron(tmp_image_captured, ones(floor(sz_ratio)));
    sampled_tmp_img = imresize(resize__tmp_img, sensor_used_pix_num, 'nearest');
    sampled_tmp_img = sampled_tmp_img/prod(sz_ratio); % each subpixel get partial energy(gray value)
end

sensor_img(sensor_img_used_idx(1):sensor_img_used_idx(3),...
        sensor_img_used_idx(2):sensor_img_used_idx(4)) = sampled_tmp_img;

%% quantification
% recover the image's real brightness(energy)
sensor_img = sensor_img*img_energy_scale;
% ADC
% default AD params
if nargin < 5
    default_sat_value = max(sensor_img, [],'all');
    default_quan_level = ceil(default_sat_value);
    AD = [default_quan_level, default_sat_value];
end

quan_level = AD(1);
sat_value = AD(2);

sensor_img = min(sensor_img, sat_value);
sensor_img = round(sensor_img./(sat_value/quan_level));
if strcmp(save_class, 'default')
	if quan_level < 256
		save_class = 'uint8';
	elseif quan_level < 65536
		save_class = 'uint16';
	else
		save_class = 'double';
	end
end
sensor_img = feval(save_class, sensor_img);

end