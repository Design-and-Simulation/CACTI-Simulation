% function [mask, meas] = dvpdata_normalize(init_mask, init_meas)
function varargout = dvpdata_normalize(varargin)
%DVPDATA_NORMALIZE normalize the mask_bayer and meas_bayer with the
%maximal value of mask_bayer for the dvp algorithm reconstruction.
% 
%	Input
%	------
%	"mask_bayer and meas_bayer" or "dataset path and save path"
% 
%	Output
%	------
%	normalized "mask_bayer and meas_bayer" or directly save normalized
%	dataset
% 

% input setting
if isnumeric(varargin{1})
	init_mask = varargin{1};
	init_meas = varargin{2};
	mask_max = max(init_mask,[],'a');
	mask = init_mask./ mask_max;
	meas = init_meas./ mask_max;
	varargout{1} = mask;
	varargout{2} = meas;
elseif ischar(varargin{1}) || isstring(varargin{1})
	dataset_path = varargin{1};
	if nargin<2
		save_path = [dataset_path(1:end-4) '_normalized.mat']; % default_save_path
	else
		save_path = varargin{2};
	end
	
	% load 
	init_mask = load(dataset_path, 'mask_bayer');
	init_mask = init_mask.mask_bayer;
	init_meas = load(dataset_path, 'meas_bayer');
	init_meas = init_meas.meas_bayer;
	orig_bayer = load(dataset_path, 'orig_bayer');
	orig_bayer = orig_bayer.orig_bayer;
	
	% normalize
	mask_max = max(init_mask,[],'a');
	mask_bayer = init_mask./ mask_max;
	meas_bayer = init_meas./ mask_max;	
	
	% save
	save(save_path,'orig_bayer','mask_bayer','meas_bayer','-v7.3')
	disp(['normalized dataset saved to: ' save_path])
else
	error("error input")
end


end

