function [trans_ratio, masked_spot_acc] = mask_propagate(mask, mask_pix_size, spot_p, spot_r, varargin)
%MASK_PROPAGATE     calculating the transmittance and masked spot when a bunch
%of light(all-pass or modulated with preceding masks) propagates through a mask
%(assuming the center of the mask is on the optical axis, and the light spot is a circle)
% 
%   Input:
%   --------
%   - mask: given mask, float, 2D matrix
% 
%   - mask_pix_element: physical size of mask's single element, float scalar
% 
%   - spot_p; physical coordinates of the spot center, 2x1 float vector
% 
%   - spot_r: physical radius of the spot
% 
%   - varargin, contain the following optional arguements
%       - resampling_factor: mask's resampling scale factor, float scalar,(when <1, 
%       the mask will be average pooled, when>1, the mask will be resampled to refine)
%       optional, default() = 100, when "masked_spot_acc_pre" is empty; otherwise
%       default() = size(masked_spot_acc_pre,1). See Note below for detailed information.
% 
%       - masked_spot_acc_pre: the accurate(no black margin) masked spot modulated 
%       by preceding masks, which will illuminate on the current mask and be 
%       modulated twice, 2D logical matrix, optional, default = [], i.e. current mask 
%       is not modulated by pre_mask, or the incident spot is a all-pass light spot
% 
%       - coord_direction: str, {'mn','xy'},  optional, default() = 'xy'; choose to
%       use mn coordinate system direction(down,right) or xy coordinate system
%       direction(right,up) for the input coord (spot_p)
% 
%       - TEST_MODE_: str, {'non-test','test'},  optional, default() = 'test', 
%       test current function by showing some figures and data
% 
%   Output:
%   --------
%   - trans_ratio: spot's transmission ratio, float, {[0-1]}
% 
%   - masked_spot_acc: the accurate masked spot modulated by the current and
%   preceding masks, which can be used as next mask-type component's
%   "masked_spot_acc_pre" parameter
% 
%   Note:
%   --------
%   - About mask's resampling:
%	In order to get accurate transmittance,  we take incompletely covered 
%	mask elemetns into consideration by virtually  resampling 
%	the mask element with smaller step, i.e., one elements will 
%	be divided into n*n virtual elements with 1/n length and the same 
%	value, here n is resampling_scale. 
% 
%	Resampling_factor ensures that the spot will contain at least
%	"resampling_factor*resampling_factor" small virtual elements after
%	resampling. The relation between resampling_scale "n" and resampling_factor 
%   "f" is: n = mask_pix_size/(2*spot_r/f). For default, it ensures the spot 
%	circle contains at least 100*100 small virtual elements, which is accurate 
%	enough in   genaral cases. So default is recommended. Use larger 
%	resampling_factor will be slow and memory consumed.
% 
%   Be careful, in order to avoid "resize" (it's time-consumed) when combining 
%   mask, the default value of resampling scale is also set to be related with 
%   "masked_spot_acc_pre", i.e., the  default "resampling_factor" will be set to
%   "size(masked_spot_acc_pre,1)", which makes "masked_spot_acc_pre" and 
%   "masked_spot_acc_cur" has likely size
% 
%   Besides, in order to shrink the mask's size, Resampling_factor < 1 is
%   valid. In this case, the mask will be average pooled to get a new mask with
%   smaller size. And the new mask's elements may be a float between [0 1], which
%   represents the equivanlent transimittance of combined elements of original mask.
% 
%   See also:
%   --------
%   MASK_AVER_POOL
% 
%   Log:
%   --------
%   * speed calculation
%   * take boundary cases into consideration
%   * implement scale resampling
%   * crop the mask to get small mask area that roughly contains the spot,
%     then conduct the resampling on the small mask area to lower memory
%     consumption
%   * change boundary processing method to ensure the cropped mask is a square
%     (by filling the surpassing boundary areas with "0")
%   * change the output, output the square cropped mask
%   * enable <1 resampling_factor
%   * take previous masks into consideration, and calculate their combinated mask's transmittance
% 
%   Info:
%   --------
%   Created:    Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-16
%   Last Modified:   Zhihong Zhang, 2020-03-30
%
%   Copyright 2020 Zhihong Zhang


%% input arguements setting
% input getting
% input num check
narginchk(4, 8);

% default valus
resampling_factor = [];
masked_spot_acc_pre = [];
coord_direction = 'x'; % "xy"
TEST_MODE_ = 'n';      % "non-test"

% input assignment
var_num = numel(varargin);
for k = 1:var_num
    var_v = varargin{k};
    if isnumeric(var_v) && isscalar(var_v)
        if (var_v==fix(var_v)) && (var_v>=1)
            resampling_factor = var_v;
        else
            error('input parameter error - "resampling_factor" ?= %f', var_v);
        end    
    elseif isnumeric(var_v) && ~isscalar(var_v) 
        if isfloat(var_v)
            masked_spot_acc_pre = var_v;
        else
            error('input parameter error - "masked_spot_acc_pre" should be a logical matrix')
        end
        
    elseif ischar(var_v) || isstring(var_v)
        if isstring(var_v)
            var_v = char(var_v); % Protect against short string
        end
        switch true
            case var_v(1)=='t' % "test"
                TEST_MODE_ = 't';
            case var_v(1)=='n' % "non-test"
                TEST_MODE_ = 'n';                
            case var_v(1)=='m' % "mn"
                coord_direction = 'm';
            case var_v(1)=='x' % "xy"
                coord_direction = 'x';                
            otherwise
                error('unknown input parameter - %s', var_item);
        end
    else
        error('input error');
    end
end        

% while numel(varargin)>0  
%     if isnumeric(varargin{1}) && isscalar(varargin{1})
%         if (varargin{1}==fix(varargin{1})) && (varargin{1}>=1)
%             resampling_factor = varargin{1};
%             varargin = varargin(2:end);
%         else
%             error('input parameter error - "resampling_factor" ?= %f', varargin{1});
%         end
%         
%     elseif ~(ischar(varargin{1}) || isstring(varargin{1})) && ismatrix(varargin{1}) && ~isscalar(varargin{1}) 
%         if isfloat(varargin{1})
%             masked_spot_acc_pre = (varargin{1});
%             varargin = varargin(2:end);
%         else
%             error('input parameter error - "masked_spot_acc_pre" should be a logical matrix')
%         end
%     end
%         
%     if ischar(varargin{1}) || isstring(varargin{1})
%         var_item = varargin{1}; 
%         switch true
%             case strcmp(var_item,'test')
%                 TEST_MODE_ = var_item;
%             case strcmp(var_item, 'mn')
%                 coord_direction = var_item;
%             otherwise
%                 error('unknown input parameter - %s', var_item);
%         end
%         varargin = varargin(2:end);
%     end
% end

TEST_MODE_FLAG = (TEST_MODE_=='t');
% masked_spot_acc_pre = double(masked_spot_acc_pre);

if isempty(resampling_factor) % unassigned resampling
    if ~isempty(masked_spot_acc_pre)
        % setting resampling_scale according to pre_mask's size to accelerate
        % calcuation
        resampling_factor = size(masked_spot_acc_pre,1); % [zzh 20200326]   
        % resampling scale: n = mask_pix_size/(2*spot_r/f)
        resampling_scale = max(round(mask_pix_size/(2*spot_r/resampling_factor)), 1);
    else
          % resampling scale: n = mask_pix_size/(2*spot_r/f), f=default=100
        %     resampling_scale = max(ceil(mask_pix_size/(2*spot_r/100)), 1);  
            resampling_scale = ceil(mask_pix_size/(2*spot_r/100));   
    end
else % no mattter to refine elements or combine elements, make sure the resampling_factor is close to given value
        resampling_scale = mask_pix_size/(2*spot_r/resampling_factor);  %[zzh 20200327]
end




% if var_num<1 || isempty(var_in{1})
% 	% use default resampling_scale value when it is unassigned
%     resampling_scale = max(ceil(mask_pix_size/(spot_r/50)), 1);
% else
%     resampling_scale = var_in{1};
% end
% 
% if var_num<2 || isempty(var_in{2})
%     coord_direction = 'xy';
% else
%     coord_direction = var_in{2};
% end
% 
% if var_num<3 || isempty(var_in{3})
%     TEST_MODE_ = false;
% else
%     TEST_MODE_ = var_in{3};
% end

%% preprocessing
% mask's full size
mask_sz = size(mask);

% if the spot is just a point(spot_r==0),directly locate the mask element 
% that containing the spot, the element value is just the transmittance
if spot_r < 1e-30 % equal to 0, consider the spot as a point
    spot_p_mn = coord2sub(spot_p, mask_sz', mask_pix_size, coord_direction);
    trans_ratio = mask(spot_p_mn(1),spot_p_mn(2));
    return
end

% avoid OOM
while resampling_scale >= 2e3
    fprintf('The resampling scale(%d) is too large and may cause OOM.\n',resampling_scale)
    fprintf('Trying %d...\n', resampling_scale/2);  
    resampling_scale = ceil(resampling_scale/2);
    if floor(2*spot_r/(mask_pix_size/resampling_scale)) < 2 
        % spot can cover no more than 4 virtual elements
        warnning_mess = ['The spot radius(%d) is too small\n',...
            'consider to stop the program(Ctrl+C) and change the input\n',...
            'or enlarge the resampling scale boundary in "mask_transmittance" function\n',...
            'or consider the spot as a point and continue the program\n'];
        fprintf(warnning_mess, spot_r);
        fprintf('Press any key to continue, Press Ctrl+C to stop...\n')
        pause
        % consider the spot as a point
        spot_p_mn = coord2sub(spot_p, mask_sz', mask_pix_size, coord_direction, 'ind_idx');
        trans_ratio = mask(spot_p_mn(1),spot_p_mn(2));
        return
    end
end


%% crop the mask and remain the area containing the spot
% suffix '_mn' means matrix coord system

% determine crop area
% idx of spot's center(float£¬ relative to physical origin point)
spot_p_mn = coord2sub(spot_p, mask_sz', mask_pix_size, coord_direction, 'float_idx');
% length of spot radius(float, in matrix coord system, the uint is "1")
spot_r_mn = spot_r/mask_pix_size;
% center elemetn of cropped rectangle
spot_p_mn_int = round(spot_p_mn);
% length of cropped radius
cropped_r_mn = ceil(spot_r_mn) + 1; % margin = 1 element size 
 %[top-left vertex, bottom-right vertex]
crop_rect_mn_tmp = [spot_p_mn_int'-cropped_r_mn, spot_p_mn_int'+cropped_r_mn];
% avoid surpassing the boundary
crop_rect_mn = max(crop_rect_mn_tmp, 1);
crop_rect_mn(3) = min(crop_rect_mn_tmp(3), mask_sz(1));
crop_rect_mn(4) = min(crop_rect_mn_tmp(4), mask_sz(2));

% crop mask (may not be a square at the boundary)
cropped_mask_tmp = mask(crop_rect_mn(1):crop_rect_mn(3), crop_rect_mn(2):crop_rect_mn(4));

% spot is surpassing the mask area, all spot will directly pass through
if isempty(cropped_mask_tmp)
    masked_spot_acc = ones(1);
    trans_ratio = 1;
    return
end    

% a full square shape mask (areas surpassing boundary will be assigned 1)
% use "ones" because we assume that, the spot surpassing the mask will directly pass through
% cropped_mask = ones(2*cropped_r_mn+1); 
cropped_mask = ones(2*cropped_r_mn);
filled_idx  = crop_rect_mn - crop_rect_mn_tmp + [1 1 2*cropped_r_mn+1 2*cropped_r_mn+1];
cropped_mask(filled_idx(1):filled_idx(3),filled_idx(2):filled_idx(4)) = cropped_mask_tmp;

% for test
if TEST_MODE_FLAG
%     crop_idx = mask;
%     crop_idx(crop_rect_mn(1):crop_rect_mn(3), crop_rect_mn(2):crop_rect_mn(4)) = max(max(mask))/2;
    mask_crop_idx = zeros(mask_sz);
    mask_crop_idx(crop_rect_mn(1):crop_rect_mn(3), crop_rect_mn(2):crop_rect_mn(4)) = 0.5;

end

%% change the relative origin point from mask's center to 
%cropped_mask center [matrix coord system]

% idx of cropped_mask's center (float£¬ relative to physical origin point)
% (square shape)
cropped_mask_p_mn = [crop_rect_mn_tmp(1)+crop_rect_mn_tmp(3); crop_rect_mn_tmp(2)+crop_rect_mn_tmp(4)]/2;

% idx of spot's center (float£¬ relative to cropped_mask's center--reference origin point)
spot_p_ref_mn = spot_p_mn - cropped_mask_p_mn; 

%% resampling cropped mask and calculating tranmittance [matrix coord system]
if resampling_scale > 1 
    % need refine to mask element smaller
    resampling_scale = ceil(resampling_scale);
    cropped_mask = kron(cropped_mask, ones(resampling_scale));
elseif resampling_scale < 1 
    % resampling_scale<1, need de-refine to combine small element together and make element larger
    cropped_mask = mask_aver_pool(cropped_mask, resampling_factor);
%     % resize cropped_mask to enable it to be exact divided 
%     crop_mask_sz = 2*cropped_r_mn+1 * [1 1];
%     combine_block_sz = ceil(crop_mask_sz/resampling_factor);
%     crop_mask_reisze_sz = combine_block_sz*resampling_factor;
%     cropped_mask = imresize(cropped_mask,crop_mask_reisze_sz,'nearest');
%     
%     %combine block elements together to get one element, and assign it with
%     %average value of these combined elements. (average transmittance)
%     combine_f =  @(A)(mean(A.data,'all'));
%     cropped_mask = blockproc(cropped_mask, combine_block_sz, combine_f);
end

cropped_mask_sz = size(cropped_mask);
mask_pix_size_mn = 1; %mask_pix_size in matrix coord system is "1"
mask_pix_size_mn = mask_pix_size_mn./resampling_scale;
 

% distance map between cropped_mask elements and spot's center
% distance_map = dist_map(cropped_mask_sz, [mask_pix_size_mn,mask_pix_size_mn], spot_p_ref_mn'); %<distance_map>
% index map
[mask_m, mask_n] = idx_matrix(cropped_mask_sz);
% physical coordinates of mask's elements (center element (on the optical axis) as origin point)
mask_x = sub2coord(mask_m, cropped_mask_sz(1), mask_pix_size_mn, 'mn');
mask_y = sub2coord(mask_n, cropped_mask_sz(2), mask_pix_size_mn, 'mn');
% distance map's square between mask element and spot center
distance_map2 = (mask_x - spot_p_ref_mn(1)).^2 + (mask_y - spot_p_ref_mn(2)).^2;

% index of the spot area on the cropped mask
spot_area_idx = zeros(cropped_mask_sz); 
% spot_area_idx(distance_map<=spot_r_mn) = 1;  % <for distance_map>
spot_area_idx(distance_map2<=spot_r_mn^2) = 1;

% covered mask area (might have black margin, but won't affect calculation of trans_ratio)
masked_spot = cropped_mask.*spot_area_idx; % i.e. masked spot

% get accurate masked spot area (throw the black margin)
row_start = find(sum(spot_area_idx,2),1);
row_end = find(sum(spot_area_idx,2), 1, 'last');
col_start = find(sum(spot_area_idx,1),1);
% col_end = find(sum(spot_area_idx,1), 1, 'last');
col_end = col_start + row_end - row_start; % mask sure the mask is a square

masked_spot_acc_cur = masked_spot(row_start:row_end, col_start:col_end); % accurate masked spot area for current mask
spot_area_idx_acc_cur = spot_area_idx(row_start:row_end, col_start:col_end);
% spot_area_idx_acc_cur = spot_area_idx;
% delete_row_idx = all(spot_area_idx_acc_cur==0,2);
% delete_col_idx = all(spot_area_idx_acc_cur==0,1);
% spot_area_idx_acc_cur(delete_row_idx,:)=[];
% spot_area_idx_acc_cur(:,delete_col_idx)=[];
% masked_spot_acc_cur = masked_spot;
% masked_spot_acc_cur(delete_row_idx,:)=[];
% masked_spot_acc_cur(:,delete_col_idx)=[];

% by accurately calculating, spot is surpassing the mask area, all spot will directly pass through
if isempty(spot_area_idx_acc_cur)
    masked_spot_acc = ones(1);
    trans_ratio = 1;
    return
end   

% calculating pre_mask-cur_masked-spot (i.e. combination of preceding mask & current mask)  area's transmission ratio
if isempty(masked_spot_acc_pre)
    % no pre_mask (i.e. current mask is not modulated by pre_mask,
    % or the incident spot is a all "1" light spot
    masked_spot_acc = masked_spot_acc_cur;
    spot_area_idx_acc = spot_area_idx_acc_cur;
    trans_ratio = sum(sum(masked_spot_acc))/sum(sum(spot_area_idx_acc));
else
    pre_spot_sz = size(masked_spot_acc_pre);
    cur_spot_sz = size(masked_spot_acc_cur);
    sz_ratio = cur_spot_sz./pre_spot_sz;
    
    % pre mask and cur mask has likely size (<=1.1 times), imcrop to match
    if (1/1.1<=sz_ratio(1)) && (sz_ratio(1)<1.0)
        crop_idx = [round((pre_spot_sz-cur_spot_sz)/2)+1,...
            round((pre_spot_sz-cur_spot_sz)/2)+cur_spot_sz];
        masked_spot_acc_pre = masked_spot_acc_pre(crop_idx(1):crop_idx(3),...
        crop_idx(2):crop_idx(4)); 
        
    elseif (1.0<sz_ratio(1)) && (sz_ratio(1)<=1.1)
        crop_idx = [round((cur_spot_sz-pre_spot_sz)/2)+1,...
            round((cur_spot_sz-pre_spot_sz)/2)+pre_spot_sz];
        masked_spot_acc_cur = masked_spot_acc_cur(crop_idx(1):crop_idx(3),...
        crop_idx(2):crop_idx(4));
        spot_area_idx_acc_cur = spot_area_idx_acc_cur(crop_idx(1):crop_idx(3),...
        crop_idx(2):crop_idx(4));
       
    % size ratio of pre mask and cur mask is larger than 2, kron and imresize
    % to match (kron before imresize can improve match performance)
% 	elseif sz_ratio(1) <= 1/2 && sz_ratio(2) <= 1/2
    elseif sz_ratio(1) <= 1/2
        masked_spot_acc_cur = kron(masked_spot_acc_cur, ones(floor(1./sz_ratio)));
        spot_area_idx_acc_cur = kron(spot_area_idx_acc_cur, ones(floor(1./sz_ratio)));
        masked_spot_acc_cur = imresize( masked_spot_acc_cur, pre_spot_sz, 'nearest');
        spot_area_idx_acc_cur = imresize( spot_area_idx_acc_cur, pre_spot_sz, 'nearest'); 
%         masked_spot_acc_cur = imbinarize(masked_spot_acc_cur);     

% 	elseif sz_ratio(1) >= 2 && sz_ratio(2) >= 2r
    elseif sz_ratio(1) >= 2
        masked_spot_acc_pre = kron(masked_spot_acc_pre, ones(floor(sz_ratio)));
        masked_spot_acc_pre = imresize( masked_spot_acc_pre, cur_spot_sz);          

            
    % size ratio of pre mask and cur mask is between 1.1 and 2, imresize to match
    elseif (1.1<=sz_ratio(1)) && (sz_ratio(1)<=2)
        masked_spot_acc_pre = imresize(masked_spot_acc_pre, cur_spot_sz, 'nearest'); 
        
    elseif (1/2<=sz_ratio(1)) && (sz_ratio(1)<=1/1.1)
        masked_spot_acc_cur = imresize(masked_spot_acc_cur, pre_spot_sz); 
%         masked_spot_acc_cur = double(imbinarize(masked_spot_acc_cur));
        spot_area_idx_acc_cur = imresize( spot_area_idx_acc_cur, pre_spot_sz, 'nearest'); 
        
    % pre mask and cur mask is equal, no operation
    end
    
    masked_spot_acc = masked_spot_acc_cur.*masked_spot_acc_pre;
    spot_area_idx_acc = spot_area_idx_acc_cur;
    trans_ratio = sum(sum(masked_spot_acc))/sum(sum(spot_area_idx_acc));
end


% show [for test]
sub_fig_sz = [4,2];
if TEST_MODE_FLAG
    figure(10)
%     figure
    subplot(sub_fig_sz(1), sub_fig_sz(2),1)
    imshow(mask)
    title('mask  cur')
    subplot(sub_fig_sz(1), sub_fig_sz(2),2)
    imshow(mask_crop_idx)
    title('cropped area cur')
    subplot(sub_fig_sz(1), sub_fig_sz(2),3)
    imshow(cropped_mask)
    title('cropped mask cur')
    subplot(sub_fig_sz(1), sub_fig_sz(2),4)
    imshow(spot_area_idx)
    title('spot on cropped mask cur')
    subplot(sub_fig_sz(1), sub_fig_sz(2),5)
    imshow(masked_spot_acc_cur)
    title('masked spot acc cur')
    subplot(sub_fig_sz(1), sub_fig_sz(2),6)
    imshow(masked_spot_acc_pre)
    title('masked spot acc pre')
    subplot(sub_fig_sz(1), sub_fig_sz(2),7)
    imshow(masked_spot_acc)
    title('accurate masked spot combined')    
    subplot(sub_fig_sz(1), sub_fig_sz(2),8)
    imshow(spot_area_idx_acc,[])
    title('accurate spot on cropped mask combined')    
    
end
end

