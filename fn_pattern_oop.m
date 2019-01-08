if(exist('file_path', 'dir')) % check if file_path is a folder
    [file_name,file_path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.png';'*.bmp';'*.jpg'},'open an image of the pattern',file_path);
else
    [file_name,file_path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.png';'*.bmp';'*.jpg'},'open an image of the pattern');
end

image = imread([file_path, file_name]);
image = rgb2gray(image);
image = image > 0;
image = bwmorph(image, 'remove', 1);
%image = image & ~bwmorph(image, 'erode', 1);
%image = bwmorph(image, 'skel', Inf);

figure;
imshow(image);

%% save images as png files
% folder_name = uigetdir
% if(folder_name ~= 0)
%     imwrite(image, strcat(folder_name, '/skeleton_image.png'));
% end

%% For a skeleton image, remove all side brances and only keep the longest one
% CC1 = bwconncomp(image);
% temp_image = zeros(size(image));
% for obj = 1 : CC1.NumObjects
%     temp_image(CC1.PixelIdxList{obj}) = 1;
%     [branch_row, branch_col] = find(bwmorph(temp_image, 'branchpoints'));
%     
%     if(~isempty(branch_row))
%         for point_num = 1 : length(branch_row)
%             min_row = max(branch_row(point_num) - 1, 1);
%             max_row = min(branch_row(point_num) + 1, size(temp_image, 1));
%             min_col = max(branch_col(point_num) - 1, 1);
%             max_col = min(branch_col(point_num) + 1, size(temp_image, 2));
%             
%             % save pixels around the branching point
%             old_pixels = temp_image(min_row : max_row, min_col : max_col);
%             %remove pixels around the branching point
%             temp_image(min_row : max_row, min_col : max_col) = 0;
%             
%             % go through all branches for the current branching point, keep
%             % only 2 longest
%             CC = bwconncomp(temp_image);
%             temp_image(min_row : max_row, min_col : max_col) = old_pixels;
%             
%             if(CC.NumObjects > 2)
%                 size_info = zeros(CC.NumObjects, 2);
%                 for n = 1 : CC.NumObjects
%                     size_info(n, 2) = numel(CC.PixelIdxList{n}(:));
%                     size_info(n, 1) = n;
%                 end
%                 
%                 size_info = sortrows(size_info, -2); % sort rows of the matrix based on the 2nd column in the descending order
%                 for obj_n = 3 : CC.NumObjects
%                     image(CC.PixelIdxList{size_info(obj_n,1)}) = 0;
%                     temp_image(CC.PixelIdxList{size_info(obj_n,1)}) = 0;
%                 end
%             end
%         end
%     end
%         
%     temp_image = temp_image*0;
% end

%% save images as png files
% folder_name = uigetdir
% if(folder_name ~= 0)
%     imwrite(image, strcat(folder_name, '/skeleton_image_filtered.png'));
% end

%% Calculate the total length of the edge
% temp_image = image;
% imSizeRow = size(image, 1);
% imSizeCol = size(image, 2);
% total_length = 0;
% for lin_ind = find(temp_image)
%     [row0, col0] = ind2sub(size(temp_image), lin_ind);
%     temp_image(row0, col0) = 0;
%     for row = max(1, row0 - 1) : min(imSizeRow, row0 + 1)
%         for col = max(1, col0 - 1) : min(imSizeCol, col0 + 1)
%             if temp_image(row, col) > 0
%                 total_length = total_length + sqrt((row0-row)^2+(col0-col)^2);
%             end
%         end
%     end
% end
% 
% interface_density = total_length/numel(image)*size(image, 1)/250*1000; % in 1/mm

blksze = 30;
thresh = 0.5;

[normim, mask] = ridgesegment(image, blksze, thresh);
show(normim,1);
show(mask, 2);

[orientim, reliability] = ridgeorient(normim, 1, 3, 3);
plotridgeorient(orientim, 5, image, 3);

angles = orientim .* (image > 0);
angles = angles(:);
angles = angles(angles > 0);

figure;
histogram(angles*180/pi, 180);

OOP_bm = OOP(angles);