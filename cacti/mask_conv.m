function conv_mask = mask_conv(mask1, mask2, pix_ratio, mode)
%MASK_CONV calculate the convolution of two masks that have fixed pixel size
%ratio
% 
%   Input:
%   --------
%   - mask1: 2D matrix, the first mask
% 
%   - mask1: 2D matrix, the second mask
% 
%   - pix_ratio: mask1's pixel size / mask2's pixel size
% 
%	- mode: char, convolution mode, 'same' | 'valid'
% 
%   Output:
%   - mask_conv: 2D matrix, the convolution result of two mask
% 
%   Note:
%   --------
%	This function can be used to calculate the two-mask optical systems
%	transmitance ratio
% 
%   Info£º
%   --------
%   Created:        Zhihong Zhang <z_zhi_hong@163.com>, 2020-07-15
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-07-15 
%   
%   Copyright 2020 Zhihong Zhang

if nargin<4
	mode = 'same';
end

% output size
out_put_size = size(mask2);

% enlarge small mask use kron for accurate match later
if any(size(mask1) < 8 )
	mask1 = kron(mask1, ones(8));
	pix_ratio = pix_ratio./8;
end

if any(size(mask2) < 8 )
	mask2 = kron(mask2, ones(8));
	pix_ratio = pix_ratio.*8;
end

if pix_ratio<1
	pix_ratio = 1./pix_ratio; % mask2's pixel size / mask1's pixel size
	kron_scale = floor(pix_ratio);
	if kron_scale.^2*numel(mask2) > 1e8
		warning(['The kron operation will create a matrix with more than '...
			'1e8 elements, which may cause OOM. Ctrl-C to stop, or press any key to continue']);
		fprintf('Info: \n file: %s \n line: %d\n', 'mask_conv.m', 66);
		pause
	else
		mask2 = kron(mask2, ones(kron_scale));
	end
	if any(kron_scale~=pix_ratio)
		mask2 = imresize(mask2, pix_ratio/kron_scale, 'nearest'); 
	end
else
	kron_scale = floor(pix_ratio);
	if kron_scale.^2*numel(mask1) > 1e8
		warning(['The kron operation will create a matrix with more than 1e8 elements, '...
			'which may cause OOM. Ctrl-C to stop, or press any key to continue']);
		fprintf('Info: \n file: %s \n line: %d\n', 'mask_conv.m', 66);
		pause
	else
		mask1 = kron(mask1, ones(kron_scale));
	end
	if any(kron_scale~=pix_ratio)
		mask1 = imresize(mask1, pix_ratio/kron_scale, 'nearest'); 
	end
end

conv_mask = conv2(mask2, mask1, mode);

% normalize to the transmitance ratio
conv_mask = conv_mask ./numel(mask1);

% resize
conv_mask = imresize(conv_mask, out_put_size);

end
