clearvars 'all';

inputTitles = {'number of samples:', 'positions', 'channels', 'slices (approx.)', 'ouput file title:','actin layer number:','alpha-actinin layer number (-1 if none):'};
defaultInputValues = {'1', '1', '4', '1', 'actin alignment','2','3'};

SampleInfo = inputdlg(inputTitles,'Input', 1, defaultInputValues);
SampleCount = str2num(SampleInfo{1});
numberOfPositions = str2num(SampleInfo{2});
numberOfChannels = str2num(SampleInfo{3});
numberOfSlices = str2num(SampleInfo{4});
actinLayerNumber = str2num(SampleInfo{6});
actininLayerNumber = str2num(SampleInfo{7});

%rootfolder = path;
if(~exist('rootfolder', 'var') || length(rootfolder) < 2)
    % rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
    rootfolder = 'F:\Images\Ivan';
end

if(actininLayerNumber < 0)
    analyzeFibroblasts = 0;
else
    analyzeFibroblasts = input('Analyze fibroblasts? (1=yes, 0=no)');
end
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
if(analyzeFibroblasts)
    fprintf(fileID, 'cardiomyocytes;position;slice;Sarcomere area fraction; Actin area fraction;Mean;Stdev;median;Mode;OOP;;fibroblasts;position;slice;Percent of Actin area;Mean;Stdev;median;Mode;OOP\n');
else
    if(actininLayerNumber < 0)
        fprintf(fileID, 'cells;position;slice; Actin area fraction;Mean;Stdev;median;Mode;OOP\n');
    else
        fprintf(fileID, 'cardiomyocytes;position;slice;Sarcomere area fraction; Actin area fraction;Mean;Stdev;median;Mode;OOP\n');
    end
end
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
    average_PercentActinArea = cell(numberOfSlices, 1);
    
    if(analyzeFibroblasts)
        average_fbs_PercentActinArea = cell(numberOfSlices, 1);
        fbs_All_angles = cell(numberOfSlices, 1);
    end
    if(actininLayerNumber >= 0)
        average_sarcomere_area = cell(numberOfSlices, 1);
    end
    
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
        blksze = floor(3/PixelSize);
        
        numberOfPositions = size(img, 1);
        
        for position = 1 : numberOfPositions
            disp(strcat('Position_', num2str(position)));
            numberOfSlices = size(img{position, 1}, 1)/numberOfChannels;
            
            if(actininLayerNumber >= 0)
                mip = MIP(img{position, 1}(:, 1), numberOfChannels, [actinLayerNumber, actininLayerNumber]);
                maxActinPx = max(mip{actinLayerNumber}(:));
                maxSarcomerePx = max(mip{actininLayerNumber}(:));
            else
                mip = MIP(img{position, 1}(:, 1), numberOfChannels, actinLayerNumber);
                maxActinPx = max(mip{actinLayerNumber}(:));
            end
            slice_shift = 0;
            
            for sliceNum = 1 : numberOfSlices
                disp(strcat('Slice_', num2str(sliceNum)));
                actin = img{position,1}{actinLayerNumber + (sliceNum - 1)*numberOfChannels,1};
                if(actininLayerNumber >= 0)
                    sarcomeres = img{position,1}{actininLayerNumber + (sliceNum - 1)*numberOfChannels,1};
                end
                if ~isa(actin, 'double')
                    if(actininLayerNumber >= 0)
                        sarcomeres = double(sarcomeres);
                    end
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
                
                if(actininLayerNumber >= 0)
                    sarcomeres(1,1) = maxSarcomerePx;
                    threshold = (0.1*max(sarcomeres(:))+0.9*min(sarcomeres(:)));
                    
                    sarcomereMask = bwmorph(sarcomeres > threshold, 'open');
                    sarcomeres_border = bwmorph(sarcomereMask, 'remove');
                    
                    sarcomereMask = sarcomereMask | imdilate(sarcomeres_border, disk);
                    % erode 60%, because cells' area is higher than sarcomere coverage
                    sarcomereMask = bwmorph(sarcomereMask,'erode',floor(d*0.6));
                    
                    actinMask = actinMask | sarcomereMask; % expand actin mask to include sarcomere mask
                end
                
                
                PercentActinArea = nnz(actinMask)/((imSizeX-20)*(imSizeY-20));
                if(shift_slises == 1 && slice_shift + 1 == sliceNum && PercentActinArea < 0.333)
                    slice_shift = slice_shift + 1;
                    disp(['New slice shift: ', num2str(slice_shift)]);
                else
                    if ~exist('thresh','var')
                        thresh = actinDetectTest(actin,blksze,PixelSize);
                    end
                    
                    %disp(['Image ' picturename ' of sample ' num2str(i)]);
                    if(actininLayerNumber >= 0)
                        nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,sarcomereMask);
                    else
                        nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,1);
                    end
                    
                    if(analyzeFibroblasts == 1)
                        fibroblast_nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,(~sarcomereMask).*actinMask);
                    end
                    
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
                    
                    average_PercentActinArea{sliceNum} = [average_PercentActinArea{sliceNum}; PercentActinArea];
                    
                    cms_oop = OOP(nonzero_orientation);
                    
                    if(actininLayerNumber >= 0)
                        sarcomere_area_fraction = nnz(sarcomereMask)/((imSizeX-20)*(imSizeY-20));
                        average_sarcomere_area{sliceNum} = [average_sarcomere_area{sliceNum}; sarcomere_area_fraction];
                        % save data in excel file
                        fprintf(fileID, [picture ';%3i;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;'], position, sliceNum - slice_shift, sarcomere_area_fraction, PercentActinArea, Mean, Std, Median, Mode, cms_oop);
                    else
                        fprintf(fileID, [picture ';%3i;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;'], position, sliceNum - slice_shift, PercentActinArea, Mean, Std, Median, Mode, cms_oop);
                    end
                    
                    % Gather data for the entire sample
                    All_angles{sliceNum - slice_shift} = [All_angles{sliceNum - slice_shift} ; nonzero_orientation];
                    
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
                        
                        PercentFibroblastActinArea = nnz(actinMask.*(~sarcomereMask))/((imSizeX-20)*(imSizeY-20));
                        %average_fbs_PercentActinArea = average_fbs_PercentActinArea + PercentFibroblastActinArea;
                        average_fbs_PercentActinArea{sliceNum} = [average_fbs_PercentActinArea{sliceNum}; PercentFibroblastActinArea];
                        
                        fbs_oop = OOP(fibroblast_nonzero_orientation);
                        
                        % save data in excel file
                        fprintf(fileID, [picture ';%3i;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n'], position, sliceNum - slice_shift, PercentFibroblastActinArea, fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_oop);
                        
                        fbs_All_angles{sliceNum - slice_shift} = [fbs_All_angles{sliceNum - slice_shift}; fibroblast_nonzero_orientation];
                    else
                        % save data in excel file
                        fprintf(fileID, '\n');
                    end
                end
            end
        end
        
    end
    
    angles_all_slices = [];
    fbs_angles_all_slices = [];
    
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
            
            if(numel(average_PercentActinArea{sliceNum}) > 0)
                average_actin = sum(average_PercentActinArea{sliceNum}(:))/numel(average_PercentActinArea{sliceNum});
            else
                average_actin = 0;
            end
            if(actininLayerNumber >= 0 && numel(average_sarcomere_area{sliceNum}) > 0)
                average_sarcomeres = sum(average_sarcomere_area{sliceNum}(:))/numel(average_sarcomere_area{sliceNum});
            else
                average_sarcomeres = 0;
            end
        else
            OrientationOrderParameter = OOP(angles_all_slices);
        end
        
        if(analyzeFibroblasts == 1)
            % Convert radians to degrees
            if(sliceNum <= size(All_angles, 1))
                fbs_All_angles_deg = rad2deg(fbs_All_angles{sliceNum});
                fbs_angles_all_slices = [fbs_angles_all_slices; fbs_All_angles{sliceNum}];
            else
                fbs_All_angles_deg = rad2deg(fbs_angles_all_slices);
            end
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
            
            if(sliceNum <= size(All_angles, 1))
                
                %average_fbs_PercentActinArea = average_fbs_PercentActinArea/PictureCount;
                if(numel(average_fbs_PercentActinArea{sliceNum}) > 0)
                    average_fbs_actin = sum(average_fbs_PercentActinArea{sliceNum}(:))/numel(average_fbs_PercentActinArea{sliceNum});
                else
                    average_fbs_actin = 0;
                end
                
                fbs_OrientationOrderParameter = OOP(fbs_All_angles{sliceNum});
            else
                fbs_OrientationOrderParameter = OOP(fbs_angles_all_slices);
            end
            %Save data
            if(sliceNum <= size(All_angles, 1))
                fprintf(fileID, 'Average;all;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;', sliceNum, average_sarcomeres, average_actin, Mean, Std, Median, Mode, OrientationOrderParameter);
                fprintf(fileID, 'Average_fbs;all;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', sliceNum, average_fbs_actin, fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_OrientationOrderParameter);
            else
                fprintf(fileID, 'Average;all;all;no;no;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;', Mean, Std, Median, Mode, OrientationOrderParameter);
                fprintf(fileID, 'Average_fbs;all;all;no;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_OrientationOrderParameter);
            end
            
        else
            if(actininLayerNumber >= 0)
                if(sliceNum <= size(All_angles, 1))
                    fprintf(fileID, 'Average;all;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', sliceNum, average_sarcomeres, average_actin, Mean, Std, Median, Mode, OrientationOrderParameter);
                else
                    fprintf(fileID, 'Average;all;all;no;no;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', Mean, Std, Median, Mode, OrientationOrderParameter);
                end
            else
                if(sliceNum <= size(All_angles, 1))
                    fprintf(fileID, 'Average;all;%3i;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', sliceNum, average_actin, Mean, Std, Median, Mode, OrientationOrderParameter);
                else
                    fprintf(fileID, 'Average;all;all;no;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n', Mean, Std, Median, Mode, OrientationOrderParameter);
                end
            end
        end
        
    end
    clearvars img actin sarcomeres;
end
fclose('all');
clearvars thresh;
