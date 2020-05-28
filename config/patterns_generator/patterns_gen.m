%% patterns generate
% generate patterns of object, dmd and mask for cacti_imgaing
% 


%% setting
frame_num = 10;
whiteboard_flag = 0;

obj_size = [512 512 frame_num];
dmd_size = [100 100 frame_num];
mask_size = [600 600];

% for dmd_method==4
dmd_subblank_num = 1;
dmd_subblank_size = [10 10];
dmd_margin = 10; % margin area can not be selected as square center

% for dmd_method==5
% initial size of random dmd before enlarge to given size
% init_dmd_size = [floor(sqrt(frame_num)) floor(sqrt(frame_num)) frame_num]; % default
init_dmd_size = [3 3 frame_num];  % granularity


% for mask_method==2
a=0; b=1;         % uniform distribution params
% for mask_method==3
mu=0.5; sigma=0.5;          % gaussian distribution params
% for mask_method = 2 or 3
kernel_size=1;


% saving date name
% data_dir = '.\CI algorithm\dvp\config\';

if whiteboard_flag
    data_name = '512scale_whiteboard_config.mat';
    obj_save_name = '512scale_whiteboard_512.mat';
else
	data_name = '512scale_traffic_config.mat';
    obj_save_name = '512scale_traffic_512.mat';
    obj_load_dir = 'E:\project\CACTI\simulation\CI algorithm\dvp\\dataset\';
    obj_data_name = 'traffic_bayer.mat';
end

dmd_path = './512scale_dmd_100.mat';
mask_path = './512scale_mask_600.mat';

% obj method: 
% 0-scene (load obj)
% 1-whiteboard; 
if whiteboard_flag
    obj_method = 1;
else
    obj_method = 0;
end


% dmd method:
% 0-load dmd
% 1-N*N single subblanks on; 
% 2-handcraft
% 3-up-half on
% 4-random fixed size blanks
% 5-random shape pattern with large granularity
if isfile(dmd_path)
	dmd_method = 0; % load existing dmd
else
	dmd_method = 1;
	disp('generate dmd');
end  


% mask method:
% 0-load mask
% 1-binary mask
% 2-grayscale mask & uniform distribution
% 3-grayscale mask & gaussian distribution
if isfile(mask_path)
	mask_method = 0; % load existing mask
else
	mask_method = 3;
	disp('generate mask');
end  
    
%% obj
if obj_method == 0
    % load, scene
    loaded_data = load([obj_load_dir, obj_data_name],'orig_bayer');
    obj = loaded_data.orig_bayer(:,:,1:frame_num);    
   
elseif obj_method == 1
    % 1-whiteboard
    obj = uint8(255*ones(obj_size));
end

% save obj
save(obj_save_name, 'obj');


%% dmd
if dmd_method == 0
    dmd = load(dmd_path,'dmd');
    dmd  = dmd.dmd;   
	
elseif dmd_method == 1
    % 1-N*N single subblank on
    dmd = zeros(dmd_size);
    dmd_blk_shape = ceil(sqrt(frame_num));
    dmd_blk_shape = [dmd_blk_shape dmd_blk_shape];
    
	for k = 1:dmd_size(3)
		tmp = zeros(dmd_blk_shape);
		tmp(ind2sub(dmd_blk_shape, k)) = 1;
		dmd_kron = kron(tmp, ones(floor(dmd_size(1:2)./dmd_blk_shape)));
		dmd(:,:,k) = imresize(dmd_kron, dmd_size(1:2));
	end
	
	dmd = double(imbinarize(dmd));
	
elseif dmd_method == 2
	% 3*3 hole£¬ 8frame
% 	dmd = zeros(dmd_size);
% 	dmd(20:22,20:22,1)=1;
% 	dmd(57:59,45:47,2)=1;
% 	dmd(60:62,70:72,3)=1;
% 	dmd(78:80,80:82,4)=1;
% 	dmd(20:22,60:62,5)=1;
% 	dmd(40:42,74:76,6)=1;
% 	dmd(53:55,25:27,7)=1;
% 	dmd(80:82,20:22,8)=1;

	% 10*10 hole£¬ 4frame
	dmd = zeros(dmd_size);
	dmd(20:29,20:29,1)=1;
	dmd(80:89,80:89,2)=1;
	dmd(20:29,80:89,3)=1;
	dmd(80:89,20:29,4)=1;

elseif dmd_method == 3
    % up-half on
    dmd = zeros(dmd_size);
    dmd(1:floor(dmd_size(1)/2),:)=1;
	
elseif dmd_method==4
	% random fixed size subblanks
	% boundary of the left-up corner of subblanks
	
	p_low = [dmd_margin dmd_margin];
	p_up = dmd_size(1:2) - dmd_subblank_size - dmd_margin;
	dmd = zeros(dmd_size);
	for k = 1:frame_num
		% random left-up corner
		left_up_p = [randi([p_low(1), p_up(1)], [1 dmd_subblank_num]);...
			randi([p_low(2), p_up(2)], [1 dmd_subblank_num])];
		right_down_p = left_up_p + dmd_subblank_size'-1;
		tmp = zeros(dmd_size(1:2) );
		
		for m=1:dmd_subblank_num
			tmp(left_up_p(1,m):right_down_p(1,m),left_up_p(2,m):right_down_p(2,m)) = 1;
		end
		
		dmd(:,:,k) = tmp;
	end
	
elseif dmd_method==5
	% generate random shape mask with large white or black  area
	init_dmd = binary_mask(init_dmd_size); % initial dmd patterns
	dmd = zeros(dmd_size);
	
	dmd_size_2h = dmd_size(1:2);           % 1 & 2 channel size of dmd
	kron_scale = floor(dmd_size_2h./init_dmd_size(1:2));
	
	for k=1:frame_num
		tmp_dmd = kron(init_dmd(:,:,k), ones(kron_scale));
		dmd(:,:,k) = imresize(tmp_dmd, dmd_size_2h, 'nearest');
	end

	dmd = double(imbinarize(dmd));
end

% save dmd
if dmd_method ~= 0
	save(dmd_path, 'dmd');
end


%% mask
if mask_method == 0
    mask = load(mask_path,'mask');
    mask = mask.mask;    
elseif mask_method == 1
    % 1-binary mask pattern
    mask = binary_mask(mask_size, mu);    
    
elseif mask_method == 2
    % 2-grayscale mask pattern, uniform distribution
    mask = gray_mask(mask_size, 'rand', [a, b, kernel_size]);      

elseif mask_method == 3
    % 3- grayscale mask pattern, guassian distribution
    mask = gray_mask(mask_size, 'randn', [mu,sigma,kernel_size]);  % gray mask pattern
end

% save mask
if mask_method~=0
	save (mask_path, 'mask');
end


%% save all data
save(data_name, 'obj', 'dmd', 'mask')
