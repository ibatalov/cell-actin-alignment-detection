clearvars 'all';

rootfolder = path;

SampleInfo = inputdlg({'Number of samples:','Excel file title:','layer number:'},'Input',1,{'1','alignment','3'});
SampleCount = str2num(SampleInfo{1});
actinLayerNumber = str2num(SampleInfo{3});

% get files for each sample
for i=1:SampleCount
    message = sprintf('Select files');
    [file,path]=uigetfile({'*.*';'*.lsm';'*.TIF';'*.tif';'*.bmp';'*.jpg'},message,rootfolder,'Multiselect','on');
    FileList{i,1}=path;
    rootfolder = path;
    FileList{i,2}= file;
end

xlsname = [FileList{1,1} SampleInfo{2} '.csv'];

% Create file for results
fileID = fopen(xlsname, 'a');
fprintf(fileID, 'filename;Percent of area;Mean;Stdev;median;Mode;OOP\n');

% Loop for each sample
for i=1:SampleCount
    
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
    
    % Loop over each picture for sample #i
    for j = 1:PictureCount
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
        
        info = img{1,4};  % Load OME metadata
        PixelSize = double(info.getPixelsPhysicalSizeX(0).value());
        Size = info.getPixelsSizeX(0).getValue();
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        
        % take into account only actin that is near sarcomeres
        
        if ~isa(actin, 'double')
            actin = double(actin);
        end
        
        
        [imSizeX,imSizeY] = size(actin);

        if ~exist('thresh','var')
            thresh = actinDetectTest(actin,blksze,PixelSize);
        end
        
        %disp(['Image ' picturename ' of sample ' num2str(i)]);
        [nonzero_orientation, mask] = actinDetectSlice(actin,blksze,thresh,Size,ones(imSizeX, imSizeY));

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
        figure, bar(xout,n,'hist')              % plot normalized histogram
        xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
        title('Histogram of Actin Orientation Angles')
        xlabel('Degrees')
        ylabel('Normalized Occurance')
        
        % Total number of actin positive pixels in the skeleton image
        % Sarcomere density = total/(image area)
        PercentActinArea = nnz(mask)/numel(mask(:));
        average_PercentActinArea = average_PercentActinArea + PercentActinArea;
        cms_oop = OOP(nonzero_orientation);
        
        % save data in excel file
        fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n'], PercentActinArea, Mean, Std, Median, Mode, cms_oop);
        % Gather data for the entire sample
        All_angles = [All_angles;nonzero_orientation];
    end
    
    % Convert radians to degrees
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
    
    average_PercentActinArea = average_PercentActinArea/PictureCount;
    
    %angles here should be in radians, cause cos() and sin() functions are
    %used
    OrientationOrderParameter = OOP(All_angles);
    fprintf(fileID, 'Average;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n\n', average_PercentActinArea, Mean, Std, Median, Mode, OrientationOrderParameter);
    
    clearvars img actin sarcomeres;
end

fclose('all');
clearvars thresh;
