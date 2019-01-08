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

SampleInfo = inputdlg({'number of samples:','ouput file title:','actin layer number:','alpha-actinin layer number:'},'Input',1,{'1','actin alignment','2','3'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

actinLayerNumber = str2num(SampleInfo{3});
actininLayerNumber = str2num(SampleInfo{4});

rootfolder = path;
if(~exist('rootfolder', 'var') || length(rootfolder) < 2)
    rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
end

analyzeFibroblasts = input('Analyze fibroblasts? (1=yes, 0=no)');

% get files for each sample
FileList = cell(i,2);
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
if(analyzeFibroblasts)
    fprintf(fileID, 'cardiomyocytes;Percent of Actin area;Mean;Stdev;median;Mode;OOP;;fibroblasts;Percent of Actin area;Mean;Stdev;median;Mode;OOP\n');
else
    fprintf(fileID, 'cardiomyocytes;Percent of Actin area;Mean;Stdev;median;Mode;OOP\n');
end
% Loop for each sample
for i=1:SampleCount
    
    disp(strcat('Sample ', num2str(i)));
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
    
    % Loop over each picture for sample #i
    for j = 1:PictureCount
        disp(strcat('Image ', num2str(j)));

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
        PixelSize = double(info.getPixelsPhysicalSizeX(0).value());
        Size = info.getPixelsSizeX(0).getValue();
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        
        % take into account only actin that is near sarcomeres
        
        if ~isa(sarcomeres, 'double')
            sarcomeres = double(sarcomeres);
            actin = double(actin);
        end
        
        
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
        
        %disp(['Image ' picturename ' of sample ' num2str(i)]);
        nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,sarcomereMask);
        
        if(analyzeFibroblasts == 1)
            fibroblast_nonzero_orientation = actinDetectSlice(actin,blksze,thresh,Size,~sarcomereMask);
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
        figure, bar(xout,n,'hist')              % plot normalized histogram
        xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
        title('Histogram of Actin Orientation Angles')
        xlabel('Degrees')
        ylabel('Normalized Occurance')
        
        % Total number of actin positive pixels in the skeleton image
        % Sarcomere density = total/(image area)
        PercentActinArea = nnz(sarcomereMask)/((imSizeX-20)*(imSizeY-20));
        average_PercentActinArea = average_PercentActinArea + PercentActinArea;
        cms_oop = OOP(nonzero_orientation);
        
        % save data in excel file
        fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;'], PercentActinArea, Mean, Std, Median, Mode, cms_oop);
        % Gather data for the entire sample
        All_angles = [All_angles;nonzero_orientation];
        
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
            average_fbs_PercentActinArea = average_fbs_PercentActinArea + PercentFibroblastActinArea;
            fbs_oop = OOP(fibroblast_nonzero_orientation);
            
            % save data in excel file
            fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n'], PercentFibroblastActinArea, fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_oop);
        
            fbs_All_angles = [fbs_All_angles; fibroblast_nonzero_orientation];
        else
            % save data in excel file
            fprintf(fileID, '\n');
        end
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
        
        average_fbs_PercentActinArea = average_fbs_PercentActinArea/PictureCount;
        
        fbs_OrientationOrderParameter = OOP(fbs_All_angles);
        
        %Save data
        fprintf(fileID, 'Average;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;;', average_PercentActinArea, Mean, Std, Median, Mode, OrientationOrderParameter);
        fprintf(fileID, 'Average_fbs;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n\n', average_fbs_PercentActinArea, fbs_Mean, fbs_Std, fbs_Median, fbs_Mode, fbs_OrientationOrderParameter);
    else
        fprintf(fileID, 'Average;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n\n', average_PercentActinArea, Mean, Std, Median, Mode, OrientationOrderParameter);
    end
     
    clearvars img actin sarcomeres;
end

fclose('all');
clearvars thresh;
