function mFiles = FileTraversal(strPath, file_type)
%FileTRAVERSAL get all the file paths in the given directory ( breadth - first 
%traversal)
% 
%   Input:
%   --------
%   - strPath: given directory path, str
% 
%   - file_type: assigned file type, str, e.g., '.jpg'
% 
%   Output:
%   --------
%   - mFiles: given directory path, str
% 
%   Info：
%   --------
%   Created:        guoxiaojie_415, 2014
%   Last Modified:  Zhihong Zhang <z_zhi_hong@163.com>, 2020-03-24    
%   
%   Reference:		https://blog.csdn.net/guoxiaojie_415/article/details/21317323
% 
%	Copyright:
%   版权声明：本文为博主原创文章，遵循 CC 4.0 BY-SA 版权协议，转载请附上原文出处链接和本声明。
%   本文链接：https://blog.csdn.net/guoxiaojie_415/article/details/21317323
% 
%   Copyright 2020 Zhihong Zhang


if nargin<2
	% all file type
	file_type = '.*';
else
	% make sure the type is char
	file_type = char(file_type);
end

strPath = char(strPath);
% make sure Strpath end up with '/' or '\'
if strPath(end)~='\' || strPath(end)~='/'
	strPath = [strPath '/'];
end


mFiles = cell(0,0);
mPath  = cell(0,0);

mPath{1}=strPath;
[~,c] = size(mPath);
while c ~= 0
	strPath = mPath{1};
	Files = dir(fullfile( strPath,'*.*'));
	LengthFiles = length(Files);
	if LengthFiles == 0
		break;
	end
	mPath(1)=[];
	iCount = 1;
	while LengthFiles>0
		if Files(iCount).isdir==1
			if Files(iCount).name ~='.'
				filePath = [strPath  Files(iCount).name '/'];
				[~,c] = size(mPath);
				mPath{c+1}= filePath;
			end
		else
			filePath = [strPath  Files(iCount).name];
			[~,col] = size(mFiles);
			mFiles{col+1}=filePath;
		end

		LengthFiles = LengthFiles-1;
		iCount = iCount+1;
	end
	[~,c] = size(mPath);
end

mFiles = mFiles';

if ~strcmp(file_type, '.*')
	file_num = numel(mFiles);
	k = 1;
	while k <= file_num
		file_name = mFiles{k};
		if ~strcmp(file_name(end-numel(file_type)+1:end), file_type)
			mFiles(k) = [];
			file_num = file_num - 1;
		else
			k = k+1;
		end
		
	end
end

	
end