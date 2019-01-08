% Alignment analysis of microscopy images (2D or 3D z-stacks) of cells stained for actin
% cytiskeleton. If the images are z-stacks, MIP is analyzed.

% INPUT: lsm image or 'any' microscopy image.
% Must contain (actin) channel, the alignment of which needs to be analyzed
% Can contain alpha-actinin channel to make cardiomyocyte-specific analysis
% If the image is a z-stack, a maximum intensity projections (MIPs) are analyzed
% Image can contain several positions, each position will be analyzed
% separately, then all positions from all files corresponding to the same
% sample are combined for the overall alignmnet analysis

% OUTPUT: a csv file that can be open in MS Excel or a similar software
% File contains the alignment analysis results, including cell area
% coverage, mean/average/mode orientation angles, OOP

% Adapted from: Adam W. Feinberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Havard University, Cambridge, MA 02138

% Created by Ivan Batalov
% Regenerative Biometerials & Therapeutics Group
% Carnegie Mellon University, Pittsburgh, PA, USA

%% Request info and open image files.
inputTitles = {'number of samples:','ouput file title:','actin layer number:','alpha-actinin layer number (-1 if none):'};
defaultInputValues = {'1','actin alignment','2','3'};

SampleInfo = inputdlg(inputTitles,'Input', 1, defaultInputValues);
SampleCount = str2num(SampleInfo{1});
actinLayerNumber = str2num(SampleInfo{3});
actininLayerNumber = str2num(SampleInfo{4});

%rootfolder = path;
if(~exist('rootfolder', 'var') || length(rootfolder) < 2)
    % rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
    rootfolder = 'F:\Images\Ivan'; % put your default folder path/name here. The path format depends on the OS.
end

if(actininLayerNumber > 0)
    analyzeFibroblasts = input('Analyze fibroblasts? (1=yes, 0=no)');
else
    analyzeFibroblasts = 0;
end

% get files for each sample. Multiple files can be open for each sample.
FileList = cell(SampleCount,2);
for i=1:SampleCount
    message = sprintf('Select files for sample %i',i);
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
    FileList{i,1}=path;
    rootfolder = path;
    FileList{i,2}= file;
end

xlsname = [FileList{1,1} SampleInfo{2} '.csv'];

%% Create file for results
fileID = fopen(xlsname, 'a');
if(analyzeFibroblasts)
    fprintf(fileID, 'cardiomyocytes,position,Sarcomere area fraction,Actin area fraction,Mean,Stdev,median,Mode,OOP,,fibroblasts,position,Percent of Actin area,Mean,Stdev,median,Mode,OOP\n');
else
    fprintf(fileID, 'cardiomyocytes,position,Sarcomere area fraction,Actin area fraction,Mean,Stdev,median,Mode,OOP\n');
end

%% Loop for each sample
for i=1:SampleCount
    
    disp(strcat('Sample_', num2str(i)));
    sample = FileList{i,2};
    if ~isa(sample,'char')
        PictureCount = length(sample);
    else
        PictureCount = 1;
    end
    
    % Setup variables for ouput data
    All_angles = [];
    fbs_All_angles = [];
    average_PercentActinArea = 0;
    average_fbs_PercentActinArea = 0;
    average_sarcomere_area = 0;
    
    %% Loop over each file for sample #i
    for j = 1:PictureCount
        disp(strcat('Image_', num2str(j)));
        
        % Open images and load actin channel in a matrix
        if PictureCount == 1
            picture = sample;
            picturename = [FileList{i,1} sample];
        else
            picture = sample{j};
            picturename = [FileList{i,1} sample{j}];
        end
        
        img = bfopen(picturename); % opens microscopy images. Documentation can be found online.
        
        info = img{1,4};  % Load OME metadata
        PixelSize = double(info.getPixelsPhysicalSizeX(0).value());
        Size = info.getPixelsSizeX(0).getValue();
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        d = floor(1.5*blksze);
        disk = strel('disk', d, 0); % used to merge sarcomeres into mask by expanding and shrinking the image
        actin_disk = strel('disk', floor(d/2), 0); % used to merge actin into mask
        
        numberOfPositions = info.getImageCount();
        
        %% Loop over all positions within the current file.
        for position = 1 : numberOfPositions
            disp(strcat('Position_', num2str(position)));
            numberOfChannels = info.getChannelCount(position - 1); % numbering starts from 0, thus the shift
            
            if(actininLayerNumber >= 0)
                mip = MIP(img{position, 1}(:, 1), numberOfChannels, [actinLayerNumber, actininLayerNumber]);
                actin = mip{1};
                sarcomeres = mip{2};
                
                threshold = (0.1*max(sarcomeres(:))+0.9*min(sarcomeres(:)));
                sarcomereMask = bwmorph(sarcomeres > threshold, 'open');
                sarcomeres_border = bwmorph(sarcomereMask, 'remove');
                % dilate and erode to fill the holes between z-lines of
                % sarcomeres
                sarcomereMask = sarcomereMask | imdilate(sarcomeres_border, disk);
                % erode only 60%, because cells' area is higher than sarcomere coverage
                sarcomereMask = bwmorph(sarcomereMask,'erode',floor(d*0.6));
            else
                mip = MIP(img{position, 1}(:, 1), numberOfChannels, actinLayerNumber);
                actin = mip{1};
            end
            
            [imSizeX,imSizeY] = size(actin);
            
            maxActinArea = (imSizeX - 2*0.6*floor(d/2))*(imSizeY - 2*floor(d/2));
            maxSarcomereArea = (imSizeX - 2*0.6*d)*(imSizeY - 2*d);
            
            actinThreshold = (0.05*max(actin(:))+0.95*min(actin(:)));
            actinMask = bwmorph(actin > actinThreshold, 'open');
            actin_border = bwmorph(actinMask, 'remove');
            
            % dilate and erode to fill the holes between actin filaments
            actinMask = actinMask | imdilate(actin_border, actin_disk);
            actinMask = bwmorph(actinMask,'erode',floor(floor(d/2)*0.6));
            if(actininLayerNumber >= 0)
                actinMask = actinMask | sarcomereMask; % expand actin mask to include sarcomere mask
            end
            
            %% Alignment analysis
            if ~exist('thresh','var')
                % determine the optimal threshold for filament detection
                % the same threshold should be used for all images
                thresh = actinDetectTest(actin,blksze,PixelSize);
            end
            
            % get filament orientation image
            if(actininLayerNumber > 0)
                nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,sarcomereMask);
            else
                nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,1);
            end
            
            if(analyzeFibroblasts == 1)
                fibroblast_nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,(~sarcomereMask).*actinMask);
            end
            
            %% Calculate mode of orientation distribution
            % Convert radians to degrees
            nonzero_orientation_angles = rad2deg(nonzero_orientation);
            
            % hist(nonzero_orientation_angles)
            Mean = mean(nonzero_orientation_angles);
            Std = std(nonzero_orientation_angles);
            Median = median(nonzero_orientation_angles);
            
            % Create histogram
            [n,xout] = hist(nonzero_orientation_angles,180);
            dx = xout(2)-xout(1);                   % calc a single bin width
            n = n / sum( n*dx );                    % normalize histogram to have area of 1
            
            % Find mode
            [~,I] = max(n);
            Mode = xout(I);
            
            %% Plot histogram of raw orientation
%             figure, bar(xout,n,'hist')              % plot normalized histogram
%             xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
%             title('Histogram of Actin Orientation Angles')
%             xlabel('Degrees')
%             ylabel('Normalized Occurance')
            
            %% Calculate cell coverage
            PercentActinArea = nnz(actinMask)/maxActinArea;
            average_PercentActinArea = [average_PercentActinArea; PercentActinArea];
            
            if(actininLayerNumber > 0)
                sarcomere_area_fraction = nnz(sarcomereMask)/maxSarcomereArea;
                average_sarcomere_area = [average_sarcomere_area; sarcomere_area_fraction];
            end
            
            %% calculate OOP
            cms_oop = OOP(nonzero_orientation);
            
            %% save data in excel file
            if(actininLayerNumber > 0)
                fprintf(fileID, [picture ',%3i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,,'], position, sarcomere_area_fraction, PercentActinArea, Mean, Std, Median, Mode, cms_oop);
            else
                fprintf(fileID, [picture ',%3i,1,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,,'], position, PercentActinArea, Mean, Std, Median, Mode, cms_oop);
            end
            
            %% Gather data for the entire sample
            All_angles = [All_angles; nonzero_orientation];
            
            %% The same for the fibroblasts
            if(analyzeFibroblasts == 1)
                % Convert radians to degrees
                fbs_nonzero_orientation_angles = rad2deg(fibroblast_nonzero_orientation);
                
                % hist(nonzero_orientation_angles)
                fbs_Mean = mean(fbs_nonzero_orientation_angles);
                fbs_Std = std(fbs_nonzero_orientation_angles);
                fbs_Median = median(fbs_nonzero_orientation_angles);
                
                % Create histogram
                [fbs_n,fbs_xout] = hist(fbs_nonzero_orientation_angles,180);
                fbs_dx = fbs_xout(2)-fbs_xout(1);                   % calc a single bin width
                fbs_n = fbs_n / sum( fbs_n*fbs_dx );                    % normalize histogram to have area of 1
                
                % Find mode
                [~,fbs_I] = max(fbs_n);
                fbs_Mode = xout(fbs_I);
                
                PercentFibroblastActinArea = nnz(actinMask.*(~sarcomereMask))/((imSizeX)*(imSizeY));
                average_fbs_PercentActinArea = [average_fbs_PercentActinArea; PercentFibroblastActinArea];
                
                fbs_oop = OOP(fibroblast_nonzero_orientation);
                
                % save data in excel file
                fprintf(fileID, [picture ',%3i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n'], position, PercentFibroblastActinArea, fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_oop);
                
                fbs_All_angles = [fbs_All_angles; fibroblast_nonzero_orientation];
            else
                % save data in excel file
                fprintf(fileID, '\n');
            end
        end % end of positions loop
    end % end of PictureCount loop
    
    %% Calculate the alignment of the entire sample.
    All_angles_deg = rad2deg(All_angles);
    % Calculate stats for entire sample
    
    Mean = mean(All_angles_deg);
    Std = std(All_angles_deg);
    Median = median(All_angles_deg);
    
    % Create histogram
    [n,xout] = hist(All_angles_deg,180);
    dx = xout(2)-xout(1);                   % calc a single bin width
    n = n / sum( n*dx );                    % normalize histogram to have area of 1
    
    % Find mode
    [~,I] = max(n);
    Mode = xout(I);
        
    %angles here should be in radians, cause cos() and sin() functions are used
    OrientationOrderParameter = OOP(All_angles);
    
    %% The same for fibroblasts
    if(analyzeFibroblasts == 1)
        % Convert radians to degrees
        fbs_All_angles_deg = rad2deg(fbs_All_angles);
        
        % Calculate stats for entire sample
        fbs_Mean = mean(fbs_All_angles_deg);
        fbs_Std = std(fbs_All_angles_deg);
        fbs_Median = median(fbs_All_angles_deg);
        
        % Create histogram
        [fbs_n,fbs_xout] = hist(fbs_All_angles_deg,180);
        fbs_dx = fbs_xout(2)-fbs_xout(1);                   % calc a single bin width
        fbs_n = fbs_n / sum( fbs_n*fbs_dx );                    % normalize histogram to have area of 1
        
        % Find mode
        [~,fbs_I] = max(fbs_n);
        fbs_Mode = xout(fbs_I);
        
        fbs_OrientationOrderParameter = OOP(fbs_All_angles);
        
        %Save data
        fprintf(fileID, 'Average,all,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,,', mean(average_sarcomere_area(:)), mean(average_PercentActinArea(:)), Mean, Std, Median, Mode, OrientationOrderParameter);
        fprintf(fileID, 'Average_fbs,all,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n', mean(average_fbs_PercentActinArea(:)), fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_OrientationOrderParameter);
    else
        if(actininLayerNumber > 0)
            fprintf(fileID, 'Average,all,%3.3f,%3.3f,%3.3f,;%3.3f,%3.3f,%3.3f,%3.3f\n', mean(average_sarcomere_area(:)), mean(average_PercentActinArea(:)), Mean, Std, Median, Mode, OrientationOrderParameter);
        else
            fprintf(fileID, 'Average,all,no,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n', mean(average_PercentActinArea(:)), Mean, Std, Median, Mode, OrientationOrderParameter);
        end
    end
end
fclose('all');
clearvars thresh;
