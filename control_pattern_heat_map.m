
if(~exist('rootfolder', 'var'))
    rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
end

heatMapWidth = 40; %in microns
heatMapHeight = 40; %in microns

SampleInfo = inputdlg({'Number of samples:','Excel file title:','actin layer:','alpha-actinin layer:', 'fibronectin layer:'},'Input',1,{'1','heat_map_data','2','3','4'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

actinLayerNumber = str2num(SampleInfo{3});
actininLayerNumber = str2num(SampleInfo{4});
fibronectinLayerNumber = str2num(SampleInfo{5});

% get files for each sample
for i=1:SampleCount
    message = sprintf('Select files for sample %i',i);
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
    if(exist(path, 'dir'))
        rootfolder = path;
    end
    FileList{i,1}=path;
    rootfolder = path;
    FileList{i,2}= file;
end

% Loop for each sample
for i=1:SampleCount
    
    disp('Sample 1')
    sample = FileList{i,2};
    if ~isa(sample,'char')
        PictureCount = length(sample);
    else
        PictureCount = 1;
    end
    
    % Setup
    All_angles = [];
    fbs_All_angles = [];
    average_PercentActinArea = 0;
    average_fbs_PercentActinArea = 0;
    
    cm_orient_block_data = cell(heatMapWidth/1);
    fb_orient_block_data = cell(heatMapWidth/1);
    
    imageData = cell(PictureCount,6);
    
    % Loop over each picture for sample #i
    for j = 1:PictureCount
        disp(['Processing image ', int2str(j), ' out of ', int2str(PictureCount) '...']);
        
        discardImage = 0;
        
        % Open images and load actin channel in a matrix
        if PictureCount == 1
            picture = sample;
            picturename = [FileList{i,1} sample];
        else
            picture = sample{j};
            picturename = [FileList{i,1} sample{j}];
        end
        
        img = bfopen(picturename);
        
        actin = img{1,1}{actinLayerNumber,1};
        sarcomeres = img{1,1}{actininLayerNumber,1};
        fibronectin = img{1,1}{fibronectinLayerNumber,1};
        
        info = img{1,4};  % Load OME metadata
        PixelSize = double(info.getPixelsPhysicalSizeX(0).value());
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        
        % take into account only actin that is near sarcomeres
        
        if ~isa(sarcomeres, 'double'),
            sarcomeres = double(sarcomeres);
            actin = double(actin);
            fibronectin = double(fibronectin);
        end
        
        [imSizeX,imSizeY] = size(sarcomeres);
        
        fibronectin = fibronectin/max(fibronectin(:));
        fibronectin = medfilt2(fibronectin,[3 3]);
        fibronectin2 = fibronectin > 0.22;
        fibronectin2 = bwmorph(fibronectin2, 'close');
        fibronectin2 = bwmorph(fibronectin2, 'open');
        %fibronectin = fibronectin2 .* fibronectin;
        %fibronectin = fibronectin2;
        
        if j == 1
            %% 
            imshow(fibronectin);
            disp('choose 2 pattern corners along the same side of a line' );
            [movingPointsX, movingPointsY] = getpts;
            close;
            % use only last 2 points
            movingPoints = [movingPointsX(1:2), movingPointsY(1:2)];
            patternOrigin = movingPoints(1,:);
            distance = sqrt((movingPointsX(1) - movingPointsX(2))^2 + (movingPointsY(1) - movingPointsY(2))^2);
           
            fixedPoints = [movingPoints(1,1), movingPoints(1,2);
                movingPoints(1,1), movingPoints(1,2)+distance];
            
            transform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
            [outboundsX, outboundsY] = outputLimits(transform,[1 imSizeX],[1 imSizeY]);
            patternOrigin(1) = patternOrigin(1) - (outboundsX(1) - 1);
            patternOrigin(2) = patternOrigin(2) - (outboundsY(1) - 1);
            newPatternOrigin = patternOrigin;
            reference_fn_image = imwarp(fibronectin, transform);
            figure
            imshow(reference_fn_image);
            
            fn_pattern = reference_fn_image(round(patternOrigin(2)) : round(patternOrigin(2)) + round(heatMapHeight/PixelSize) - 1, round(patternOrigin(1)) : round(patternOrigin(1)) + round(heatMapWidth/PixelSize) - 1);
            
        else
            while 1
                
                imshow(fibronectin);
                disp('choose 2 pattern corners along the same side of a line' );
                [movingPointsX, movingPointsY] = getpts;
                close;
                
                pointsEntered = numel(movingPointsX);
                
                if pointsEntered == 2
                    
                    movingPoints = [movingPointsX, movingPointsY];
                    newPatternOrigin = movingPoints(1,:);
                    
                    distance = sqrt((movingPoints(1,1) - movingPoints(2,1))^2 + (movingPoints(1,2) - movingPoints(2,2))^2);
                    fixedPoints = [movingPoints(1,1), movingPoints(1,2);
                    movingPoints(1,1), movingPoints(1,2)+distance];
                    
                    transform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
                else
                    if numel(movingPointsX) == 1
                        movingPoints = [movingPointsX, movingPointsY];
                        newPatternOrigin = movingPoints(1,:);
                        % transform stays the same, so empty line here
                    else
                        discardImage = 1;
                    end
                end
                
                newPatternOrigin = transformPointsForward(transform, newPatternOrigin);
                
                [outboundsX, outboundsY] = outputLimits(transform,[1 imSizeX],[1 imSizeY]);
                newPatternOrigin(1) = newPatternOrigin(1) - (outboundsX(1) - 1);
                newPatternOrigin(2) = newPatternOrigin(2) - (outboundsY(1) - 1);
                
                new_fn_image = imwarp(fibronectin, transform);
                overlay_patterns(fn_pattern > 0.15, new_fn_image > 0.15, newPatternOrigin);
                
                disp('red: first pattern, green: current pattern. If patterns do not match, choose shift vector (first point - new pattern, second - old)');
                [inputX, inputY] = getpts;
                while numel(inputX) == 2
                    close;
                    newPatternOrigin(1) = newPatternOrigin(1) + (inputX(2) - inputX(1));
                    newPatternOrigin(2) = newPatternOrigin(2) + (inputY(2) - inputY(1));
                    overlay_patterns(fn_pattern > 0.15, new_fn_image > 0.15, newPatternOrigin);
                    [inputX, inputY] = getpts;
                end
                if isempty(input('Type any key to fix the pattern manually or Enter to continue'))
                    break;
                end
            end
            
            
            clearvars ptsCurrent featuresCurrent validPtsCurrent indexPairs matchedReference matchedCurrent
        end
        
        if discardImage == 0
            imageData{j,1} = actin;
            imageData{j,2} = sarcomeres;
            imageData{j,3} = fibronectin;
            imageData{j,4} = PixelSize;
            imageData{j,5} = transform;
            imageData{j,6} = newPatternOrigin;
        end
        
    end
    
    %%
    for j = 1: size(imageData,1)
        if ~isempty(imageData{j,1})
            
            disp(['Analyzing alignment of image # ', int2str(j), ' out of ', int2str(size(imageData,1)) '...']);
            
            actin = imageData{j,1};
            sarcomeres = imageData{j,2};
            fibronectin = imageData{j,3};
            PixelSize = imageData{j,4};
            transform = imageData{j,5};
            newPatternOrigin = imageData{j,6};
            
            blksze = floor(3/PixelSize);
            [imSizeX,imSizeY] = size(sarcomeres);
            
            sarcomereMask = zeros(imSizeX, imSizeY);
            actinMask = zeros(imSizeX, imSizeY);
            
            threshold = (0.22*max(sarcomeres(:))+0.78*min(sarcomeres(:)));
            actinThreshold = (0.05*max(actin(:))+0.95*min(actin(:)));
            d = floor(1.5*blksze);
            
            %disp('making sarcomere mask');
            
            circle = zeros(2*d+1);
            for x = -d:d
                for y = -d:d
                    circle(x+d+1,y+d+1) = x*x+y*y < d*d;
                end
            end
            
            for x = d+1:imSizeX-d
                for y = d+1:imSizeY-d
                    if sarcomeres(x,y) > threshold
                        sarcomereMask(x-d:x+d, y-d:y+d) = sarcomereMask(x-d:x+d, y-d:y+d) | circle;
                    end
                    if actin(x,y) > actinThreshold
                        actinMask(x-d:x+d, y-d:y+d) = actinMask(x-d:x+d, y-d:y+d) | circle;
                    end
                end
            end
            
            % erode 60%, because cells' area is higher than sarcomere coverage
            sarcomereMask = bwmorph(sarcomereMask,'erode',floor(d*0.6));
            actinMask = bwmorph(actinMask,'erode',floor(d*0.6));
  
            if ~exist('thresh','var')
                thresh = actinDetectTest(actin,blksze,PixelSize);
            end
            
            %%
            actin = imwarp(actin, transform);
            sarcomeres = imwarp(sarcomeres, transform);
            sarcomereMask = imwarp(sarcomereMask, transform);
            actinMask = imwarp(actinMask, transform);
            
            %fibronectin = imwarp(fibronectin, transform);
            %fn_pattern = fibronectin(round(newPatternOrigin(2)) : round(newPatternOrigin(2)) + round(heatMapHeight/PixelSize) - 1, round(newPatternOrigin(1)) : round(newPatternOrigin(1)) + round(heatMapWidth/PixelSize) - 1);
            %% splitting orientations between blocks
            block_num = 0;
            block_y = 0;
            
            cm_orientim = create_orientation_image(actin,blksze,thresh,sarcomereMask);
            %fb_orientim = create_orientation_image(actin,blksze,thresh,~sarcomereMask);
            
            disp('starting to sort data');
            for col = 1 : size(cm_orientim, 2)
                    block_num = floor((col - newPatternOrigin(1))/(1/PixelSize));
                    
                    block_num = rem(block_num, heatMapWidth/1);
                    
                    if block_num <= 0
                        block_num = block_num + heatMapWidth/1;
                    end
                    
                    if nnz(cm_orientim(:,col)) > 0
                        cm_orient_block_data{block_num} = [cm_orient_block_data{block_num}; cm_orientim(cm_orientim(:,col) ~= 0 ,col)];
                    end
                    
            end
            disp('data sorted');

            % Total number of actin positive pixels in the skeleton image
            % Sarcomere density = total/(image area)
            PercentActinArea = nnz(sarcomereMask)/((imSizeX-20)*(imSizeY-20));
            %PercentFibroblastActinArea = nnz(actinMask.*(~sarcomereMask))/((imSizeX-20)*(imSizeY-20));
            
            average_PercentActinArea = average_PercentActinArea + PercentActinArea;
            %average_fbs_PercentActinArea = average_fbs_PercentActinArea + PercentFibroblastActinArea;
        else
            disp(['Image # ', int2str(j), ' is discarded!']);
        end
        %}
    end
    
    cm_alignment_data = zeros(heatMapWidth/1, 6);
    all_cm_angles = [];
    %fb_alignment_data = zeros(heatMapWidth/1, heatMapHeight/1, 6);
    for block_num = 1 : heatMapWidth/1
            angles = cm_orient_block_data{block_num};
            all_cm_angles = [all_cm_angles; angles];

            if numel(angles) > 5
                
                % Convert radians to degrees
                angles_deg = rad2deg(angles);
                
                Mean = mean(angles_deg);
                Std = std(angles_deg);
                Median = median(angles_deg);
                
                % Create histogram
                [n,xout] = hist(angles_deg,180);
                dx = xout(2)-xout(1);                   % calc a single bin width
                n = n / sum( n*dx );                    % normalize histogram to have area of 1
                
                % Find mode
                [~,I] = max(n);
                Mode = xout(I);
                
                OrientationOrderParameter = OOP(angles);
                cm_alignment_data(block_num, :) = [numel(angles), Mean, Std, Median, Mode, OrientationOrderParameter];
            else
                cm_alignment_data(block_num, :) = [0, 0, 0, 0, 0, 0];
            end
    end
    
    average_values = [0; 0; 0; 0; 0; 0]; % average for the pattern location
    
    for k = 1 : 6
        average_values(k) = sum(cm_alignment_data(:,k))/nnz(cm_alignment_data(:,k));
    end
    
    all_cm_angles_deg = rad2deg(all_cm_angles);
    Mean = mean(all_cm_angles_deg);
    Std = std(all_cm_angles_deg);
    Median = median(all_cm_angles_deg);
    % Create histogram
    [n,xout] = hist(all_cm_angles_deg,180);
    dx = xout(2)-xout(1);                   % calc a single bin width
    n = n / sum( n*dx );                    % normalize histogram to have area of 1
    % Find mode
    [~,I] = max(n);
    Mode = xout(I); 
    OrientationOrderParameter = OOP(all_cm_angles);
    
    overall_values = [average_values(1), Mean, Std, Median, Mode, OrientationOrderParameter];
    %overall_values = average_values;
    
    average_PercentActinArea = average_PercentActinArea/length(imageData);
    %average_fbs_PercentActinArea = average_fbs_PercentActinArea/PictureCount;
    
    data_label = {'Number of orientation vectors in each block';
        'Mean orientation angle';
        'Standatd deviation of the mean angle';
        'Median orientation angle';
        'Most probable orienation angle';
        'OOP'};
    
    %Save data
    %save(strcat(path,'results.mat'), 'fn_pattern', 'cm_alignment_data', 'fb_alignment_data');
    %save(strcat(path,'results.mat'), 'fn_pattern', 'cm_alignment_data');
    
    %orig_pattern = imread('D:\Ivan\2015\2015-05-07 Matlab heat map stuff\original pattern\BM_pattern_for_heat_map_contour_1.png');
    orig_pattern = imread('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-10-26 20x20 heat maps stuff/orig_pattern.bmp');
    orig_pattern = rgb2gray(orig_pattern);
    orig_pattern = double(orig_pattern);
    orig_pattern = orig_pattern/max(orig_pattern(:));
    %orig_pattern = orig_pattern > 0.99;
    %orig_pattern = imresize(orig_pattern, size(fn_pattern), 'Method', 'bilinear');
    imshow(orig_pattern)
    bin_size = size(orig_pattern, 2)/heatMapWidth;
    
    for mapNum = 1 : 6
        valueHeatMap = zeros(size(orig_pattern,1), size(orig_pattern,2));
        
        for block_num = 1 : heatMapWidth/1
                c1 = max(1,round((block_num - 1)*bin_size));
                c2 = round(min(heatMapWidth*bin_size,block_num*bin_size));
                
                valueHeatMap(:, c1: c2) = cm_alignment_data(block_num, mapNum);
                
        end
        
        low_hue = 240/360;
        high_hue = 0;
        palitra_length = round(size(orig_pattern, 1)*2/3);
        font_size = 10;
        
        max_value = round(max(valueHeatMap(:)),4, 'significant');
        min_value = round(min(valueHeatMap(:)),4, 'significant');
        avg_value = round(average_values(mapNum),4, 'significant');
        overall_value = round(overall_values(mapNum),4, 'significant');
        
        hsv = zeros(size(orig_pattern,1), size(orig_pattern,2),3);
        hsv(:,:,1) = low_hue + (high_hue - low_hue)*(valueHeatMap - min(valueHeatMap(:)))/(max(valueHeatMap(:)) - min(valueHeatMap(:)));
        hsv(:,:,2) = 1;
        hsv(:,:,3) = orig_pattern*0.8 + 0.1;
        %     figure;
        %     imshow(OOP_map);
        %     figure
        %     imshow(hsv2rgb(hsv));
        
        scale_hsv = zeros(palitra_length, 120,3);
        scale_hsv(:,:,3) = 1;
        scale_hsv(:,:,2) = 0;
        scale_hsv(:,5:40,2) = 1;
        
        for row = 1 : palitra_length
            scale_hsv(row,5:40,1) = high_hue - (high_hue - low_hue)*(row/palitra_length);
        end
        
        min_y = 1;
        max_y = palitra_length;
        avg_y = round(1 + (palitra_length - 1)*(max_value - avg_value)/(max_value - min_value));
        overall_y = round(1 + (palitra_length - 1)*(max_value - overall_value)/(max_value - min_value));
        
        min_text_y = 1;
        max_text_y = max_y - font_size - 1;
        avg_text_y = avg_y - round(font_size/2);
        overall_text_y = overall_y - round(font_size/2);
        
        delta = abs(avg_y - overall_y);
        if(avg_y ~= overall_y)
            if(delta < font_size + 2)
                shift = font_size + 2 - delta;
                if(avg_y < overall_y)
                    shift = - shift;
                end
                avg_text_y = avg_text_y + round(shift/2);
                overall_text_y = overall_text_y - round(shift/2);
            end
        end
        
        scale_hsv(max_y,5:43,3) = 0;
        scale_hsv(min_y,5:43,3) = 0;
        scale_hsv(avg_y,5:43,1) = 0;
        scale_hsv(avg_y,5:43,2) = 1;
        if(avg_y ~= overall_y)
            scale_hsv(overall_y,5:43,1) = 170;
            scale_hsv(overall_y,5:43,2) = 1;
            
            textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', font_size);
            strings = uint8([num2str(max_value) 0 num2str(min_value) 0 num2str(avg_value) 0 num2str(overall_value)]);
            labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([45, min_text_y; 45, max_text_y; 45, avg_text_y; 45, overall_text_y]));
            %     figure;
            %     imshow(labeled);
        else
            textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', font_size);
            strings = uint8([num2str(max_value) 0 num2str(min_value) 0 num2str(avg_value)]);
            labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([45, min_text_y; 45, max_text_y; 45, avg_text_y]));
        end
        scale_bar_microns = floor(120*PixelSize/10)*10;
        scale_bar_px = scale_bar_microns/PixelSize;
        left_margin = round((120 - scale_bar_px)/2);
        
        scale_bar = ones(120,120,3);
        scale_bar(80:100, left_margin:round(left_margin + scale_bar_px), :) = 0;
        string = uint8([num2str(scale_bar_microns) ' um']);
        textInserter2 = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', 30);
        scale_bar = step(textInserter2, scale_bar, string, int32([16 40]));
        %     figure;
        %     imshow(scale_bar);
        
        final_image = ones(size(valueHeatMap,1), size(valueHeatMap,2) + size(labeled,2),3);
        final_image(:,1:size(valueHeatMap,2),:) = hsv2rgb(hsv);
        final_image(1:size(labeled,1),size(valueHeatMap,2) + 1:size(valueHeatMap,2) + size(labeled,2),:) = labeled;
        final_image(size(final_image,1) - size(scale_bar,1) + 1:size(final_image,1),size(final_image,2) - size(scale_bar,2) + 1:size(final_image,2),:) = scale_bar;
        figure;
        finalFigure = imshow(final_image);
        title(data_label(mapNum));
    end
    clearvars img actin sarcomeres;
end

fclose('all');
clearvars thresh;