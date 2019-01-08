% actinDetectMulti
% To use on multilayer .lsm images of cells sained for actin and alpha-actinin. Saves everything to excel file.
% Code for automatic detection of actin filament alignment in cells
% stained for F-actin (such as phalloidin conjugated to FITC)
%
% Adapted from Function to demonstrate use of fingerprint code
%
% Argument:   Load image of actin stained cells, file should
%             have .lsm format grayscale at least 8-bit depth
%
% Returns:    *.actinDetect.Settings.mat file containing actin alignment
%               data and parameters needed to run actinDetectMulti on a folder
%               full of image files

% Adapted from:
% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
% January 2005
%
% Created by Adam W. Feinberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Havard University, Cambridge, MA 02138

% Updated January 2014 by Quentin Jallerat
% Updated March 2014 by Ivan Batalov

%ask user for number of samples

%load the libraries in the matlab directory, this is fucking stupid
addpath('/Users/ivan/Documents/MATLAB/lab folder/Finger Print Detection - MATLAB')
addpath('/Users/ivan/Documents/MATLAB/lab folder/Quentin Actin Alignment')
addpath('/Users/ivan/Documents/MATLAB/lab folder/bfmatlab')

clearvars 'all';



SampleInfo = inputdlg({'Number of samples:','Excel file title:','actin layer:','alpha-actinin layer:', 'fibronectin layer:'},'Input',1,{'1','heat_map_data','2','3','4'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

actinLayerNumber = str2num(SampleInfo{3});
actininLayerNumber = str2num(SampleInfo{4});
fibronectinLayerNumber = str2num(SampleInfo{5});

rootfolder = path;

% get files for each sample
for i=1:SampleCount
    message = sprintf('Select files for sample %i',i);
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
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
    
    cm_orient_block_data = cell(155/1, 250/1);
    fb_orient_block_data = cell(155/1, 250/1);
    
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
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        
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
            %% Transforming all layers to make the FN pattern rectangular (250x155 µm)
            imshow(fibronectin);
            disp('choose 2 pattern corners along the vertical axis' );
            [movingPointsX, movingPointsY] = getpts;
            close;
            % use only last 2 points
            movingPoints = [movingPointsX(1:2), movingPointsY(1:2)];
            patternOrigin = movingPoints(1,:);
            
            
            
%             fixedPoints = [movingPoints(1,1), movingPoints(1,2);
%                 movingPoints(1,1), movingPoints(1,2)+250/PixelSize;
%                 movingPoints(1,1)+151/PixelSize, movingPoints(1,2)+250/PixelSize;
%                 movingPoints(1,1)+151/PixelSize, movingPoints(1,2)];
            
            distance = sqrt((movingPointsX(1) - movingPointsX(2))^2 + (movingPointsY(1) - movingPointsY(2))^2);

            fixedPoints = [movingPoints(1,1), movingPoints(1,2);
                movingPoints(1,1), movingPoints(1,2)+distance];
            
            transform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
            %transform = estimateGeometricTransform(movingPoints, fixedPoints, 'similarity');
            [outboundsX, outboundsY] = outputLimits(transform,[1 imSizeX],[1 imSizeY]);
            patternOrigin(1) = patternOrigin(1) - (outboundsX(1) - 1);
            patternOrigin(2) = patternOrigin(2) - (outboundsY(1) - 1);
            newPatternOrigin = patternOrigin;
            reference_fn_image = imwarp(fibronectin, transform);
            figure
            imshow(reference_fn_image);
            
            fn_pattern = reference_fn_image(round(patternOrigin(2)) : round(patternOrigin(2)) + round(250/PixelSize) - 1, round(patternOrigin(1)) : round(patternOrigin(1)) + round(155/PixelSize) - 1);
            
%             disp('Detecting reference features...');
%             ptsReference  = detectSURFFeatures(reference_fn_image,'NumOctaves', 10, 'MetricThreshold', 1000, 'NumScaleLevels', 6, 'ROI', [floor(patternOrigin(1)) floor(patternOrigin(2)) floor(151/PixelSize) floor(250/PixelSize)]);
%             ptsReference = detectMSERFeatures(reference_fn_image, 'ThresholdDelta', 2, 'RegionAreaRange', [500 15000],'ROI', [floor(patternOrigin(1)) floor(patternOrigin(2)) floor(151/PixelSize) floor(250/PixelSize)]);
%             [featuresReference,   validPtsReference]  = extractFeatures(reference_fn_image,  ptsReference);
%             disp('detected!');
        else
%                         figure
%                         imshow(fibronectin);
%                         rect = getrect;
%                         close;
            
%             disp('Detecting features...');
%             ptsCurrent = detectSURFFeatures(fibronectin,'NumOctaves', 10, 'MetricThreshold', 1000, 'NumScaleLevels', 6);
%             ptsCurrent = detectMSERFeatures(fibronectin, 'ThresholdDelta', 2, 'RegionAreaRange', [500 15000]);
%             ptsCurrent = detectMSERFeatures(fibronectin, 'ThresholdDelta', 2, 'RegionAreaRange', [500 15000]);
%             
%             [featuresCurrent, validPtsCurrent]  = extractFeatures(fibronectin, ptsCurrent);
%             disp('detected!');
            
%             figure;
%             imshow(fibronectin);
%             hold on
%             plot(ptsCurrent);
            
%             disp('Matching features...');
%             indexPairs = matchFeatures(featuresReference, featuresCurrent);
%             disp('mathed!');
%             disp(strcat('pairs matched: ', num2str(size(indexPairs,1))));
            
%             points1 = validPtsReference(indexPairs(1:20,1));
%             points2 = validPtsCurrent(indexPairs(1:20,2));
%             figure;
%             showMatchedFeatures(reference_fn_image, fibronectin, points1, points2,'montage');
            
%             
%             matchedReference  = validPtsReference(indexPairs(:,1));
%             matchedCurrent = validPtsCurrent(indexPairs(:,2));
%             disp('Finding Transform...');
%             transform = estimateGeometricTransform(matchedCurrent, matchedReference, 'projective');
%             newPatternOrigin = patternOrigin;
%             disp('found!');
%             [outboundsX, outboundsY] = outputLimits(transform,[1 imSizeX],[1 imSizeY]);
%             newPatternOrigin(1) = newPatternOrigin(1) - (outboundsX(1) - 1);
%             newPatternOrigin(2) = newPatternOrigin(2) - (outboundsY(1) - 1);
%             
%             new_fn_image = imwarp(fibronectin, transform);
%             
%             newPatternOrigin = overlay_patterns(fn_pattern, new_fn_image, newPatternOrigin);
            
            while 1
                %disp('bad match');
                
                imshow(fibronectin);
                disp('choose 2 pattern corners' );
                [movingPointsX, movingPointsY] = getpts;
                close;
                
                pointsEntered = numel(movingPointsX);
                
                if pointsEntered >= 2
                    if(pointsEntered == 3)
                    % use only last 4 points
                    movingPoints = [movingPointsX(numel(movingPointsX) - 1 : numel(movingPointsX)), movingPointsY(numel(movingPointsY) - 1 : numel(movingPointsY))];
                    newPatternOrigin = [movingPointsX(1), movingPointsY(1)];
                    else
                        movingPoints = [movingPointsX, movingPointsY];
                        newPatternOrigin = movingPoints(1,:);
                    end
                    
                    distance = sqrt((movingPointsX(1) - movingPointsX(2))^2 + (movingPointsY(1) - movingPointsY(2))^2);
                    fixedPoints = [movingPoints(1,1), movingPoints(1,2);
                    movingPoints(1,1), movingPoints(1,2)+distance];
                
                    transform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
                    %transform = estimateGeometricTransform(movingPoints, fixedPoints, 'similarity');
                    
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
                overlay_patterns(fn_pattern > 0.1, new_fn_image > 0.1, newPatternOrigin);
                
                disp('red: first pattern, green: current pattern. If patterns do not match, choose shift vector (first point - new pattern, second - old)');
                [inputX, inputY] = getpts;
                while numel(inputX) == 2
                    close;
                    newPatternOrigin(1) = newPatternOrigin(1) + (inputX(2) - inputX(1));
                    newPatternOrigin(2) = newPatternOrigin(2) + (inputY(2) - inputY(1));
                    overlay_patterns(fn_pattern, new_fn_image, newPatternOrigin);
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
    for j = 1:PictureCount
        if ~isempty(imageData{j,1})
            
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
            
            %disp('sarcomere mask is made');
            
            %actin = actin.*sarcomereMask;
            
            % show(sarcomereMask);
            % show(actinMask);
            
            if ~exist('thresh','var')
                thresh = actinDetectTest(actin,blksze,PixelSize);
            end
            
            %%
            actin = imwarp(actin, transform);
            sarcomeres = imwarp(sarcomeres, transform);
            sarcomereMask = imwarp(sarcomereMask, transform);
            actinMask = imwarp(actinMask, transform);
            
            %fibronectin = imwarp(fibronectin, transform);
            %fn_pattern = fibronectin(round(newPatternOrigin(2)) : round(newPatternOrigin(2)) + round(250/PixelSize) - 1, round(newPatternOrigin(1)) : round(newPatternOrigin(1)) + round(155/PixelSize) - 1);
            %% splitting orientations between blocks
            block_x = 0;
            block_y = 0;
            
            cm_orientim = create_orientation_image(actin,blksze,thresh,sarcomereMask);
            %fb_orientim = create_orientation_image(actin,blksze,thresh,~sarcomereMask);
            
            for row = 1 : size(cm_orientim, 1)
                for col = 1 : size(cm_orientim, 2)
                    if cm_orientim(row,col) ~= 0
                        % || fb_orientim(row,col) ~= 0
                        block_x = floor((col - newPatternOrigin(1))/(1/PixelSize));
                        block_y = floor((row - newPatternOrigin(2))/(1/PixelSize));
                        
                        block_x = rem(block_x, 155/1);
                        block_y = rem(block_y, 250/1);
                        
                        if block_x <= 0
                            block_x = block_x + 155/1;
                        end
                        if block_y <= 0
                            block_y = block_y + 250/1;
                        end
                        if cm_orientim(row,col) ~= 0
                            cm_orient_block_data{block_x, block_y} = [cm_orient_block_data{block_x, block_y}; cm_orientim(row,col)];
                        end
%                         if fb_orientim(row,col) ~= 0
%                             fb_orient_block_data{block_x, block_y} = [fb_orient_block_data{block_x, block_y}; fb_orientim(row,col)];
%                         end
                    end
                end
            end
            
            % Total number of actin positive pixels in the skeleton image
            % Sarcomere density = total/(image area)
            PercentActinArea = nnz(sarcomereMask)/((imSizeX-20)*(imSizeY-20));
            %PercentFibroblastActinArea = nnz(actinMask.*(~sarcomereMask))/((imSizeX-20)*(imSizeY-20));
            
            average_PercentActinArea = average_PercentActinArea + PercentActinArea;
            %average_fbs_PercentActinArea = average_fbs_PercentActinArea + PercentFibroblastActinArea;            
        end
%}
    end
    
    cm_alignment_data = zeros(155/1, 250/1, 6);
    %fb_alignment_data = zeros(155/1, 250/1, 6);
    
    for block_x = 1 : 155/1
        for block_y = 1 : 250/1
            if numel(cm_orient_block_data{block_x, block_y}) > 1
                angles = [];
                for deltaX = -2 : 2
                    for deltaY = -2 : 2
                        indX = rem(block_x + deltaX + 155/1,155/1);
                        indY = rem(block_y + deltaY + 250/1,250/1);
                        if indX <= 0
                            indX = indX + 155;
                        end
                        if indY <= 0
                            indY = indY + 250;
                        end
                        angles = [angles; cm_orient_block_data{indX, indY}];
                    end
                end
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
                cm_alignment_data(block_x, block_y, :) = [numel(angles), Mean, Std, Median, Mode, OrientationOrderParameter];
            else
                cm_alignment_data(block_x, block_y, :) = [0, 0, 0, 0, 0, 0];
            end
            
%             if numel(fb_orient_block_data{block_x, block_y}) > 1
%                 angles = [];
%                 for deltaX = -2 : 2
%                     for deltaY = -2 : 2
%                         indX = rem(block_x + deltaX + 155/1,155/1);
%                         indY = rem(block_y + deltaY + 250/1,250/1);
%                         if indX <= 0
%                             indX = indX + 155;
%                         end
%                         if indY <= 0
%                             indY = indY + 250;
%                         end
%                         angles = [angles; fb_orient_block_data{indX, indY}];
%                     end
%                 end
%                 % Convert radians to degrees
%                 angles_deg = rad2deg(angles);
%                 
%                 Mean = mean(angles_deg);
%                 Std = std(angles_deg);
%                 Median = median(angles_deg);
%                 
%                 % Create histogram
%                 [n,xout] = hist(angles_deg,180);
%                 dx = xout(2)-xout(1);                   % calc a single bin width
%                 n = n / sum( n*dx );                    % normalize histogram to have area of 1
%                 
%                 % Find mode
%                 [~,I] = max(n);
%                 Mode = xout(I);
%                 
%                 OrientationOrderParameter = OOP(angles);
%                 fb_alignment_data(block_x, block_y, :) = [numel(angles), Mean, Std, Median, Mode, OrientationOrderParameter];
%             else
%                 fb_alignment_data(block_x, block_y, :) = [0, 0, 0, 0, 0, 0];
%             end
        end
    end
    
    average_PercentActinArea = average_PercentActinArea/PictureCount;
    %average_fbs_PercentActinArea = average_fbs_PercentActinArea/PictureCount;
    
    data_label = {'Number of orientation vectors in each block';
        'Mean orientation angle';
        'Standatd deviation of the mean angle';
        'Median orientation angle';
        'Most probable orienation angle';
        'OOP'};
    
    %Save data
    %save(strcat(path,'results.mat'), 'fn_pattern', 'cm_alignment_data', 'fb_alignment_data');
    save(strcat(path,'results.mat'), 'fn_pattern', 'cm_alignment_data');
    
    OOP_map = zeros(round(250/PixelSize), round(155/PixelSize));
    for block_x = 1 : 155/1
        for block_y = 1 : 250/1
            OOP_map(max(1,round((block_y - 1)*1/PixelSize)) : round(min(250/PixelSize,(block_y)*1/PixelSize)), max(1,round((block_x - 1)*1/PixelSize)) : round(min(155/PixelSize,block_x*1/PixelSize))) = cm_alignment_data(block_x, block_y,6);
        end
    end
    hsv = zeros(size(fn_pattern,1), size(fn_pattern,2),3);
    hsv(:,:,1) = (OOP_map - min(OOP_map(:)))/(max(OOP_map(:)) - min(OOP_map(:)))/4;
    hsv(:,:,2) = 1;
    hsv(:,:,3) = fn_pattern*0.4 + 0.3;
    imshow(OOP_map);
    surf(cm_alignment_data(:, :, 1));
    figure
    imshow(hsv2rgb(hsv));
    
    scale_hsv = zeros(800, 120,3);
    scale_hsv(:,:,3) = 1;
    scale_hsv(:,:,2) = 0;
    scale_hsv(:,5:40,2) = 1;
    
    for row = 1 : 800
        scale_hsv(row,5:40,1) = (1 - row/800)/4;
    end
    
    textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', 25);    
    strings = uint8([num2str(round(max(OOP_map(:))*1000)/1000) 0 num2str(round(min(OOP_map(:))*1000)/1000)]);
    labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([42 1; 42 775]));
    figure;
    imshow(labeled);
    
    scale_bar_microns = floor(120*PixelSize/10)*10;
    scale_bar_px = scale_bar_microns/PixelSize;
    left_margin = round((120 - scale_bar_px)/2);
    
    scale_bar = ones(120,120,3);
    scale_bar(80:100, left_margin:round(left_margin + scale_bar_px), :) = 0;
    string = uint8([num2str(scale_bar_microns) ' um']);
    textInserter2 = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', 30);    
    scale_bar = step(textInserter2, scale_bar, string, int32([16 40]));
    figure;
    imshow(scale_bar);
    
    final_image = ones(size(fn_pattern,1), size(fn_pattern,2) + size(labeled,2),3);
    final_image(:,1:size(fn_pattern,2),:) = hsv2rgb(hsv);
    final_image(1:size(labeled,1),size(fn_pattern,2) + 1:size(fn_pattern,2) + size(labeled,2),:) = labeled;
    final_image(size(final_image,1) - size(scale_bar,1) + 1:size(final_image,1),size(final_image,2) - size(scale_bar,2) + 1:size(final_image,2),:) = scale_bar;
    figure;
    imshow(final_image);
         for index = 1 : 6
    
             HM = HeatMap(cm_alignment_data(:,:,index).')
            addTitle(HM, data_label{index})
         end
    
    clearvars img actin sarcomeres;
end

fclose('all');
clearvars thresh;

% 
% fn_temp = fibronectin > 0.1;
% fn_temp = bwmorph(fn_temp,'erode', 1);
% fn_temp = bwmorph(fn_temp,'dilate', 2);
% fn_temp = bwmorph(fn_temp,'erode', 1);
% figure; imshow(fn_temp);
% 
% CC = bwconncomp(fn_temp)
% sizes = [];
% for index = 1 : length(CC.PixelIdxList)
%             currSize = numel(CC.PixelIdxList{index});
%             if currSize > 1
%             sizes = [sizes, numel(CC.PixelIdxList{index})];   
%             end
% end
%         figure
%         histogram(sizes(sizes > 8000), 100);
