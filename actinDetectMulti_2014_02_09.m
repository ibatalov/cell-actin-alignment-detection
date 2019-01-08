% actinDetectMulti 
% To use on single tifs. Saves everything to excel file. 
% Code for automatic detection of actin filament alignment in cells
% stained for F-actin (such as phalloidin conjugated to FITC)
%
% Adapted from Function to demonstrate use of fingerprint code
%
% Argument:   Load image of actin stained cells, file should
%             TIF format grayscale at least 8-bit depth
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

% Last updated January 2014 by Quentin Jallerat

%ask user for number of samples
SampleInfo = inputdlg({'Number of samples:','Excel file title:'},'Input',1,{'1','actin alignment'});
SampleCount = str2num(SampleInfo{1});
xlsname = ['\DATA\' SampleInfo{2}]; 
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
    sample
    
    % Create Excel sheet
    sheetname =  sample{1}(1:end-5);
    headlines = {'Percent Actin Area','Mean','Std','Median','Mode','OOP'};
    range = 'B1:G1';
    xlswrite(xlsname,headlines,sheetname,range);
    
    % Setup
    All_angles = [];
    clear OOZ;

    % Loop over each picture for sample #i
    for j = 1:PictureCount
        % Open images and load actin channel in a matrix
        picture = sample{j};
        picturename = [FileList{i,1} picture];
        img = bfopen(picturename);
        actin = img{1,1}{1,1};
        info = img{1,4};  % Load OME metadata
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        Size = str2double(info.getPixelsSizeX(0));
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(5/PixelSize);
        
        if ~exist('thresh','var')
            thresh = actinDetectTest(actin,blksze,PixelSize);
        end
        
        disp(['Image ' picturename ' of sample ' num2str(i)]);
        nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size);
        
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
        
%         % Plot histogram of raw orientation
%         figure, bar(xout,n,'hist')              % plot normalized histogram
%         xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
%         title('Histogram of Actin Orientation Angles')
%         xlabel('Degrees')
%         ylabel('Normalized Occurance')

        % Total number of actin positive pixels in the skeleton image
        % Sarcomere density = total/(image area)
        PercentActinArea = length(nonzero_orientation)/((Size-20)*(Size-20));

        % save data in excel file
        picDATA = {picture,PercentActinArea,Mean,Std,Median,Mode};
        range = sprintf('A%i:F%i',j+1,j+1);
        xlswrite(xlsname,picDATA,sheetname,range);
        
        % Gather data for the entire sample
        All_angles = [All_angles;nonzero_orientation];
        
    end
    
    % Calculate stats for entire sample 
    
    Mean = mean(All_angles);
    Std = std(All_angles);
    Median = median(All_angles);
    
    % Create histogram
    [n,xout] = hist(All_angles,180);
    dx = xout(2)-xout(1);                   % calc a single bin width
    n = n / sum( n*dx );                    % normalize histogram to have area of 1
    
    % Find mode
    [~,I] = max(n);
    Mode = xout(I);
    
    PercentActinArea = length(All_angles)/PictureCount/((Size-20)*(Size-20));
    
    OrientationOrderParameter = OOP(All_angles);
    
    %Save data to xls
    sampleDATA = {'Total',PercentActinArea,Mean,Std,Median,Mode,OrientationOrderParameter};
    range = sprintf('A%i:G%i',PictureCount+2,PictureCount+2);
    xlswrite(xlsname,sampleDATA,sheetname,range);
end

   