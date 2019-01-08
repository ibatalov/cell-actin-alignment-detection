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

SampleInfo = inputdlg({'Number of samples:','Excel file title:','actin layer number:','alpha-actinin layer number:'},'Input',1,{'1','actin alignment','2','3'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

actinLayerNumber = str2num(SampleInfo{3});
actininLayerNumber = str2num(SampleInfo{4});

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
fprintf(fileID, 'picture;Percent of Actin area;Mean;Stdev;median;Mode;OOP\n');

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
    
    % Setup
    All_angles = [];

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
        sarcomeres = img{1,1}{actininLayerNumber,1};
        
        info = img{1,4};  % Load OME metadata
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        Size = str2double(info.getPixelsSizeX(0));
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        
        % take into account only actin that is near sarcomeres
       
        if ~isa(sarcomeres, 'double'),
            sarcomeres = double(sarcomeres);
            actin = double(actin);
        end
        
        
        [imSizeX,imSizeY] = size(sarcomeres);
        sarcomereMask = zeros(imSizeX, imSizeY);
        
        threshold = (0.22*max(sarcomeres(:))+0.78*min(sarcomeres(:)));
        d = floor(1.3*blksze)/blksze;
        
        %disp('making sarcomere mask');
        
        circle = zeros(2*d*blksze+1);
        for x = -d*blksze:d*blksze
            for y = -d*blksze:d*blksze
                circle(x+d*blksze+1,y+d*blksze+1) = x*x+y*y < d*blksze*d*blksze;
            end
        end
        
        for x = d*blksze+1:imSizeX-d*blksze
            for y = d*blksze+1:imSizeY-d*blksze   
                if sarcomeres(x,y) > threshold
                    sarcomereMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) = sarcomereMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) | circle;
                end
            end
        end
        
        
        %THIS LINE INVERTS THE MASK FOR CRDIOMYOCYTES, REMOVING THEM FROM
        %ANALYSIS
        sarcomereMask = sarcomereMask == 0;
        
        %disp('sarcomere mask is made');
        
        %actin = actin.*sarcomereMask;
        
        show(sarcomereMask);
        show(sarcomeres);
        
        if ~exist('thresh','var')
            thresh = actinDetectTest(actin,blksze,PixelSize);
        end
        
        %disp(['Image ' picturename ' of sample ' num2str(i)]);
        nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,sarcomereMask);
        
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
        PercentActinArea = length(nonzero_orientation)/((Size-20)*(Size-20));

        % save data in excel file
        fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n'], PercentActinArea, Mean, Std, Median, Mode);
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
    
    PercentActinArea = length(All_angles_deg)/PictureCount/((Size-20)*(Size-20));
    
    %angles here should be in radians, cause cos() and sin() functions are
    %used
    OrientationOrderParameter = OOP(All_angles);
    
    %Save data
    fprintf(fileID, 'Average;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n\n', PercentActinArea, Mean, Std, Median, Mode, OrientationOrderParameter);
    
    clearvars img actin sarcomeres;
end

fclose('all');
clearvars thresh;
   