function dst_mask = mask_aver_pool(mask, dst_size, mode)
%MASK_AVER_POOL average pool the given mask to get a new mask with given size
%   This can combine adjacent mask elements to get a new elements which has
%   equivalent transmittance
% 
%   Input:
%   --------
%   - img:      original image, 2D double matrix
% 
%   - dst_size: size of destination, 2-element int row vector
% 
%   - mode:     image resampling method, str, "fast"| "acc", optional,
%               defaut="fast".
% 
%   Note:
%   --------
%     For "fast" mode,  use "imresize" to enlarge the original image
%     to enable it to be devided exactly, but resize process will
%     involve interpolation, which brings about inaccurate factor.
%     For "acc" mode, use "kron", instead of "imresize" to enlarge
%     the original image, it's accurate, but memory-consumed and 
%     time-consumed.
% 
%   Output:
%   --------
%   - dst_img: image after average pooling, 2D matrix
% 
%   See also:
%   --------
%   MASK_ENLARGE
% 
%   Info:
%   --------
%   Created:    Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-29
%   Last Modified:   Zhihong Zhang, 2020-03-29
%
%   Copyright (c) 2020 Zhihong Zhang


%% imput setting
if nargin<3
    mode = "fast";
end
mode = [lower(char(mode)) ' ']; % Protect against short string

if isvector(dst_size) && numel(dst_size)==2
     if ~isrow(dst_size)
        dst_size = dst_size(:)';
     end
elseif isscalar(dst_size)
    dst_size = [dst_size dst_size];
else
     error("error input - 'dst_size'");
end


%% processing
% input mask size
mask_sz = size(mask);

if mode(1) == 'f'

    % convert mask to double
    if islogical(mask)
        mask = double(mask);
    end    
    
    block_sz = round(mask_sz./dst_size);
    enlarge_mask_sz = block_sz.*dst_size; 
    enlarge_mask = imresize(mask,enlarge_mask_sz);

    dst_mask = blkcolfun(enlarge_mask, block_sz, @mean, 'distinct');   
    
elseif  mode(1) == 'a'
    kron_f = lcm(dst_size, mask_sz(1))./mask_sz(1);
    block_sz = lcm(dst_size, mask_sz(1))./dst_size;
    if prod(mask_sz.*kron_f) > 1e8
        warning_mess = ['The "kron" function will create a large matrix(%d elements), ',...
            'which may cause OOM or Not Responding. Press any key to continue, press Ctrl-C to exit'];
        warning(warning_mess, prod(mask_sz.*kron_f) );
        pause
    end
        
    enlarge_mask = kron(mask,ones(kron_f));

    dst_mask = blkcolfun(enlarge_mask, block_sz, @mean, 'distinct');
else
    error("input error - 'mode'");
end

end