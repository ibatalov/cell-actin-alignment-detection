%load the libraries in the matlab directory, this is fucking stupid
addpath('/Users/ivan/Documents/MATLAB/lab folder/Finger Print Detection - MATLAB')
addpath('/Users/ivan/Documents/MATLAB/lab folder/Quentin Actin Alignment')
addpath('/Users/ivan/Documents/MATLAB/lab folder/bfmatlab')

clearvars 'all';

SampleInfo = inputdlg({'Number of samples:','Excel file title:','nuclei layer number'},'Input',1,{'1','nuclei density','1'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

dapiLayerNumber = str2num(SampleInfo{3});

rootfolder = path;

% get files for each sample
for i=1:SampleCount
    message = sprintf('Select files for sample %i',i);
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
    FileList{i,1}=path;
    rootfolder = path;
    FileList{i,2}= file;
end

xlsname = [FileList{1,1} SampleInfo{2} '.csv'];
    
% Create file for results
fileID = fopen(xlsname, 'a');
fprintf(fileID, 'file name;shapes/mm2;nuclei/mm2;nuclei area fraction\n');

% Loop for each sample
for i=1:SampleCount
    
    disp(strcat('Sample_',num2str(i),'_out of_',num2str(SampleCount)));
    
    sample = FileList{i,2};
    if ~isa(sample,'char')
        PictureCount = length(sample);
    else
        PictureCount = 1;
    end
    
    % Setup
    averageNucleiDensity = 0;
    
    % Loop over each picture for sample #i
    for j = 1:PictureCount
        disp(strcat('Image_',num2str(j),'_out of_',num2str(PictureCount)));
        % Open images and load actin channel in a matrix
        if PictureCount == 1
            picture = sample;
            picturename = [FileList{i,1} sample];
        else
            picture = sample{j};
            picturename = [FileList{i,1} sample{j}];
        end
        
        img = bfopen(picturename);
        
        nuclei = img{1,1}{dapiLayerNumber,1};
        
        info = img{1,4};  % Load OME metadata
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        Size = str2double(info.getPixelsSizeX(0));
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        nuclei = double(nuclei);
        nuclei = nuclei/max(nuclei(:));
        
        [imSizeX,imSizeY] = size(nuclei);
        
        back_filter = fspecial('gaussian', round(60/PixelSize), round(20/PixelSize));
        background = filter2(back_filter, nuclei);
        minimum = -max(max(-background));
        background = background - minimum;
        nuclei = nuclei - background;
        
        sampleN = 256;
        max_map = zeros(imSizeX/sampleN, imSizeY/sampleN);
        for row = 1 : imSizeX/sampleN
            for col = 1 : imSizeY/sampleN
                b_r_1 = (row - 1)*sampleN + 1;
                b_r_2 = row*sampleN;
                b_c_1 = (col - 1)*sampleN + 1;
                b_c_2 = col*sampleN;
                max_map(row,col) = max(max(nuclei(b_r_1:b_r_2,b_c_1:b_c_2)));
            end
        end
        max_map = padarray(max_map, [1 1], 'replicate');
        max_map = interp2(max_map, log2(sampleN), 'spline');
        max_map = max_map(sampleN/2 + 1:size(max_map,1) - sampleN/2 - 1, sampleN/2 + 1:size(max_map,2) - sampleN/2 - 1);
        
        nuclei = nuclei./max_map;
        
        blur_filter = fspecial('disk', round(2/PixelSize));
        nuclei = filter2(blur_filter, nuclei);
        
        upper_thresh = 0.8;
        lower_thresh = 0.05;
        
        nuclei = (nuclei - upper_thresh).*(nuclei < upper_thresh) + upper_thresh;
        nuclei = (nuclei - lower_thresh).*(nuclei > lower_thresh);
        nulei = nuclei/max(max(nuclei));
        
%         [edges, threshold] = edge(nuclei, 'canny');
% %         delta = threshold(2) - threshold(1);
% %         threshold(2) = threshold(2);
% %         threshold(1) = threshold(1) + delta*2/3;
% %         edges = edge(nuclei, 'canny', threshold);
% %         
%         edges = bwmorph(edges, 'dilate', 3);
%         edges = bwmorph(edges, 'skel', Inf);
%         edges = bwmorph(edges, 'hbreak');
% %          while nnz(bwmorph(edges, 'endpoints')) > 0
% %              edges = edges - bwmorph(edges, 'endpoints');
% %          end
%         
%         figure
%         imshow(edges);
% %         
%         features = detectMSERFeatures(nuclei, 'RegionAreaRange', [round((5/PixelSize)^2) round((20/PixelSize)^2)]);
%         figure
%         imshow(nuclei);
%         hold on
%         plot(features);
         
        binary_nuclei = nuclei > 0.1;
        
        se = strel('disk',round(1.5/PixelSize));
        
        binary_nuclei = imdilate(binary_nuclei, se);
        binary_nuclei = imerode(binary_nuclei, se);
        
        binary_nuclei = imerode(binary_nuclei, se);
        binary_nuclei = imdilate(binary_nuclei, se);
        
        %% remove small shapes
        CC = bwconncomp(binary_nuclei);
        pixelLists = CC.PixelIdxList;
        
        max_nucleus_area = 300; % in µm2
        min_nucleus_area = (4/PixelSize)^2;
        nuclei_count = 0;
        shape_count = length(pixelLists);
        total_nuc_area = 0;
        
        for shapeNum = 1 : shape_count
            if numel(pixelLists{shapeNum}) < min_nucleus_area
                binary_nuclei(pixelLists{shapeNum}) = 0;
            else if(numel(pixelLists{shapeNum}) < max_nucleus_area)
                    nuc_area = numel(pixelLists{shapeNum})*PixelSize^2;
                    nuclei_count = nuclei_count + 1;
                    total_nuc_area = total_nuc_area + nuc_area;
                else
                    nuc_area = numel(pixelLists{shapeNum})*PixelSize^2;
                    nuclei_count = nuclei_count + floor(nuc_area/max_nucleus_area) + 1;
                    total_nuc_area = total_nuc_area + nuc_area;
                end
            end
        end
        
%         figure
%         imshow(binary_nuclei);
%         title(picture);
%         
        image_area = imSizeX*imSizeY*PixelSize^2/(1000^2); % in mm2
        total_nuc_area_mm = total_nuc_area/(1000^2);
        % save data in excel file
        fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f\n'], shape_count/image_area, nuclei_count/image_area, total_nuc_area_mm/image_area);
        
    end
    
    %Save data
    fprintf(fileID, '\n\n');
    
end

fclose('all');
   