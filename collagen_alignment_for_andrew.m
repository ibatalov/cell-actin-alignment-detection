clearvars 'all';

inputTitles = {'number of samples:', 'positions', 'channels', 'slices (max.)', 'ouput file title:','data layer number:'};
defaultInputValues = {'1', '1', '1', '1', 'alignment results','1'};

SampleInfo = inputdlg(inputTitles,'Input', 1, defaultInputValues);
SampleCount = str2num(SampleInfo{1});
numberOfPositions = str2num(SampleInfo{2});
numberOfChannels = str2num(SampleInfo{3});
numberOfSlices = str2num(SampleInfo{4});
dataLayerNumber = str2num(SampleInfo{6});

rootfolder = path;
% if(~exist('rootfolder', 'var') || length(rootfolder) < 2)
    % rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
    % rootfolder = 'F:\Images\Ivan';
% end

shift_slises = input('Detect 1st slice based on actin coverage? (1=yes, 0=no)');

% get files for each sample
FileList = cell(SampleCount,2);
for i=1:SampleCount
    message = sprintf('Select files for sample %i',i);
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
    FileList{i,1}=path;
    rootfolder = path;
    FileList{i,2}= file;
end

xlsname = [FileList{1,1} SampleInfo{5} '.csv'];

% Create file for results
fileID = fopen(xlsname, 'a');
fprintf(fileID, 'file name,position,slice, filled area fraction,Mean,Stdev,median,Mode,OOP\n');
% Loop for each sample
for i=1:SampleCount
    
    disp(strcat('Sample_', num2str(i)));
    sample = FileList{i,2};
    if ~isa(sample,'char')
        PictureCount = length(sample);
    else
        PictureCount = 1;
    end
    
    % Setup
    All_angles = cell(numberOfSlices, 1);
    average_area_fraction = cell(numberOfSlices, 1);
    
    % Loop over each picture for sample #i
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
        
        img = bfopen(picturename);
        
        info = img{1,4};  % Load OME metadata
        PixelSize = double(info.getPixelsPhysicalSizeX(0).value());
        Size = info.getPixelsSizeX(0).getValue();
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        % Define blksze (should be 3/pixelsize based on recommendation)
        input_blk = 3;
        input_blk = input(['characteristic width of the fiber in microns (currently ', num2str(input_blk), ' um): ']);
        blksze = floor(input_blk/PixelSize);
        
        numberOfPositions = size(img, 1);
        
        for position = 1 : numberOfPositions
            disp(strcat('Position_', num2str(position)));
            numberOfSlices = size(img{position, 1}, 1)/numberOfChannels;
            
            
            mip = MIP(img{position, 1}(:, 1), numberOfChannels, dataLayerNumber);
            maxActinPx = max(mip{dataLayerNumber}(:));
            slice_shift = 0;
            
            for sliceNum = 1 : numberOfSlices
                disp(strcat('Slice_', num2str(sliceNum)));
                actin = img{position,1}{dataLayerNumber + (sliceNum - 1)*numberOfChannels,1};
                if ~isa(actin, 'double')
                    actin = double(actin);
                end
                
                actin(1,1) = maxActinPx;
                
                d = floor(1.5*blksze);
                disk = strel('disk', d, 0);
                
                [imSizeX,imSizeY] = size(actin);
                
                actinThreshold = (0.05*max(actin(:))+0.95*min(actin(:)));
                
                actinMask = bwmorph(actin > actinThreshold, 'open');
                actin_border = bwmorph(actinMask, 'remove');
                
                actinMask = actinMask | imdilate(actin_border, disk);
                actinMask = bwmorph(actinMask,'erode',floor(d*0.6));
        
                PercentFilledArea = nnz(actinMask)/imSizeX/imSizeY;
                if(shift_slises == 1 && slice_shift + 1 == sliceNum && PercentFilledArea < 0.333)
                    slice_shift = slice_shift + 1;
                    disp(['New slice shift: ', num2str(slice_shift)]);
                else
                    if ~exist('thresh','var')
                        thresh = actinDetectTest(actin,blksze,PixelSize);
                    end
                    
                    nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,1);
                  
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
                    
                    % Plot histogram of raw orientation
                    %                 figure, bar(xout,n,'hist')              % plot normalized histogram
                    %                 xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
                    %                 title('Histogram of Actin Orientation Angles')
                    %                 xlabel('Degrees')
                    %                 ylabel('Normalized Occurance')
                    %
                    % Total number of actin positive pixels in the skeleton image
                    % Sarcomere density = total/(image area)
                    
                    average_area_fraction{sliceNum} = [average_area_fraction{sliceNum}; PercentFilledArea];
                    
                    cms_oop = OOP(nonzero_orientation);
                    
                    fprintf(fileID, [picture ',%3i,%3i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n'], position, sliceNum - slice_shift, PercentFilledArea, Mean, Std, Median, Mode, cms_oop);
                    
                    % Gather data for the entire sample
                    All_angles{sliceNum - slice_shift} = [All_angles{sliceNum - slice_shift} ; nonzero_orientation];
                end
            end
        end
        
    end
    
    angles_all_slices = [];
    
    for sliceNum = 1 : size(All_angles, 1) + 1
        % Convert radians to degrees
        if(sliceNum <= size(All_angles, 1))
            All_angles_deg = rad2deg(All_angles{sliceNum});
            angles_all_slices = [angles_all_slices; All_angles{sliceNum}];
        else
            All_angles_deg = rad2deg(angles_all_slices);
        end
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
        
        %average_PercentActinArea = average_PercentActinArea/PictureCount;
        
        %angles here should be in radians, cause cos() and sin() functions are
        %used
        if(sliceNum <= size(All_angles, 1))
            OrientationOrderParameter = OOP(All_angles{sliceNum});
            
            if(numel(average_area_fraction{sliceNum}) > 0)
                average_actin = sum(average_area_fraction{sliceNum}(:))/numel(average_area_fraction{sliceNum});
            else
                average_actin = 0;
            end
            %fprintf(fileID, 'Average,all,%3i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n', sliceNum, average_actin, Mean, Std, Median, Mode, OrientationOrderParameter);
        else
            OrientationOrderParameter = OOP(angles_all_slices);
            fprintf(fileID, 'Average,all,all,no,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n', Mean, Std, Median, Mode, OrientationOrderParameter);
        end
    end
    clearvars img actin sarcomeres;
end
fclose('all');
clearvars thresh;
