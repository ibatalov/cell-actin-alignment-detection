addpath('/Users/ivan/Documents/MATLAB/lab folder/Finger Print Detection - MATLAB')
addpath('/Users/ivan/Documents/MATLAB/lab folder/Quentin Actin Alignment')
addpath('/Users/ivan/Documents/MATLAB/lab folder/bfmatlab')

%clearvars -except path rootfolder;

SampleInfo = inputdlg({'Number of samples:','Excel file title:','alpha-actinin layer number:'},'Input',1,{'1','sarcomere alignment','3'});
SampleCount = str2num(SampleInfo{1});
%xlsname = ['\DATA\' SampleInfo{2}];

actininLayerNumber = str2num(SampleInfo{3});

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

% Create 2 files for results
fileID = fopen(xlsname, 'w');
fprintf(fileID, 'picture;Percent of sarcomere area;Mean;Stdev;median;Mode;OOP\n');

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
    All_rad = [];
    
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
        sarcomeres = double(img{1,1}{actininLayerNumber,1});
        im = sarcomeres;
        
        info = img{1,4};  % Load OME metadata
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        Size = str2double(info.getPixelsSizeX(0));
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        pix2um = 1/PixelSize;
        SarcSpacing = 2*pix2um;          % convert 2 um sarcomere spacing into pixels
        MinSarcSpacing = round(0.75*SarcSpacing); % minimum length of sarcomere spacing
        MaxSarcSpacing = round(1.25*SarcSpacing); % maximum length of sarcomere spacing
        blksze = MaxSarcSpacing;
        
        if ~exist('thresh', 'var')
            % Identify ridge-like regions and normalise image
            index = 0;
            while index < 1;
                %        % Input blocksize over which the standard deviation is determined
                %        blksze = input('Enter blocksize (16): ');
                % Threshold of standard deviation to decide if a block is a ridge region
                thresh = input('Enter Threshold (0.1 - 0.2): ');
                disp('Normalizing Image and Creating Mask' )
                [normim, mask] = ridgesegment(im, round(blksze/2), thresh);
                show(normim,1);
                show(mask, 2);
                % Determine if normalization and mask look good, click on image to
                % accept or press any key to enter new values
                w = input('Accept Threshold (yes = 0, no = 1): ');
                if w == 0
                    disp('Image threshold accepted' )
                    index = 1;
                else
                    disp('Re-analyze imaging...')
                end
            end
        else
            % Identify ridge-like regions and normalise image
            disp('Normalizing Image and Creating Mask' )
            [normim, mask] = ridgesegment(im, round(blksze/2), thresh);
            show(normim,1);
            show(mask, 2);
        end
        
        % Create skeleton image from normalized image
        normim_bin = normim > 0;
        mask_bin = mask > 0;
        normim_mask = normim_bin.*mask_bin;
        norm_bin_skel = bwmorph(normim_mask,'skel',Inf);
        show(norm_bin_skel, 10);
        
        % Determine ridge orientations
        disp('Calculating Ridge Orientations' )
        [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
        % Add pi/2 for direction orthogonal to
        % sarcomere in the direction of force generation
        orientim_perp = orientim + (pi/2);
        plotridgeorient(orientim_perp, 25, im, 4)
        show(reliability,5)
        
        % Determine ridge frequency values across the image
        disp('Calculating Ridge Frequency Values' )
        [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 3, MinSarcSpacing, MaxSarcSpacing);
        
        % Actually I find the median frequency value used across the whole
        % fingerprint gives a more satisfactory result...
        freq = medfreq.*mask;
        
        if ~exist('kx','var')
            index = 0;
            while index < 1;
                % kx controls the sigma in the x direction which is along the
                % filter, and hence controls the bandwidth of the filter.
                kx = input('Enter Sigma along filter (0.4): ');
                % ky controls the sigma across the filter and hence controls the
                % orientational selectivity of the filter. A value of 0.4 for both
                % kx and ky is a good starting point.
                ky = input('Enter Sigma across filter (0.4): ');
                newim = ridgefilter(normim, orientim, freq, kx, ky, 1);
                show(newim,6);
                binim = newim > 0;
                binim_skel = bwmorph(binim,'skel',Inf);
                binim_skel = bwmorph(binim_skel,'clean',Inf);
                show(binim_skel,7);
                b = input('Accept Sarcomere Detection (yes = 0, no = 1): ');
                if b == 0
                    disp('Image accepted' )
                    index = 1;
                else
                    disp('Re-analyze imaging...')
                end
            end
        else
            newim = ridgefilter(normim, orientim, freq, kx, ky, 1);
            show(newim);
            binim = newim > 0;
            binim_skel = bwmorph(binim,'skel',Inf);
            binim_skel = bwmorph(binim_skel,'clean',Inf);
            show(binim_skel,11);
        end
        
        % save filtered image of sarcomere skeleton before manual removal of
        % non-sarcomeres
        type = '.sarcomere.FULL.tif';
        filename3 = [picturename type];
        imwrite(binim_skel,filename3,'Compression','none');
        
        % My awesome idea #2 - cutting too long (and/or too short)
        % uniaxially aligned lines of skeleton
        
        branchingPoints = bwmorph(binim_skel,'branchpoints');
        skeletonMask2 = binim_skel.*(branchingPoints == 0)*0;
        orientTemp = orientim.*skeletonMask2;
        
        
%         for linInd = transpose(find(orientTemp))
%             [row, column] = ind2sub(size(orientTemp), linInd);
%             [arrayRow, arrayColumn] = getUniaxialPiece(row, column, orientTemp);
%             dist = 0;
%             for x1 = 1:size(arrayColumn)
%                 for x2 = x1:size(arrayColumn)
%                     R = sqrt((arrayRow(x1)-arrayRow(x2))*(arrayRow(x1)-arrayRow(x2)) + (arrayColumn(x1)-arrayColumn(x2))*(arrayColumn(x1)-arrayColumn(x2)));
%                     dist = max([dist, R]);
%                 end
%             end
%             % max length of the sarcomere
%             if dist*PixelSize > 7 || dist*PixelSize < 0.5
%                 for z = 1:size(arrayRow)
%                     skeletonMask2(arrayRow(z),arrayColumn(z)) = 0;
%                 end
%             end
%         end
        
        deltaPhi = double(input('Enter the max angle deviation for filetring (in degrees): '));
        % length in µm
        maxLineLength = 8;
        for phi = 0:pi/180:pi
            coorientedPixels = abs(cos(orientTemp-phi)) > cos(deltaPhi*pi/180);
            segments = bwconncomp(coorientedPixels);
            numPixels = cellfun(@numel,segments.PixelIdxList);
            pixelsToDelete = segments.PixelIdxList(numPixels < maxLineLength/PixelSize);
            pixelsToDelete = cat(1,pixelsToDelete{:});
            skeletonMask2(pixelsToDelete) = 1;
        end
        
        show(binim_skel);
        show(skeletonMask2);
        show(binim_skel > skeletonMask2);
        input('Enjoy the view and then press Enter to continue: ');
        
        % My awesome idea #1 - look for the previous/next sarcomere. If
        % there's none - remove from analysis
        orientimTemp = orientim.*binim_skel;
        dim = (MaxSarcSpacing - MinSarcSpacing)/2;
        skeletonMask = zeros(Size);
        
        if ~exist('cosThresh', 'var')
            cosThresh = double(input('Enter the min(<cos(phi)>) between adjacent sarcomeres (0-1): '));
        end
        
        for row = 1 + floor(dim*sqrt(2)+1+SarcSpacing):Size - floor(dim*sqrt(2)+1+SarcSpacing)
            for column = 1+floor(dim*sqrt(2)+1+SarcSpacing):Size-floor(dim*sqrt(2)+1+SarcSpacing)
                
                sum1 = 0;
                sum2 = 0;
                block1size = 0;
                block2size = 0;
                
                if orientimTemp(row,column) ~= 0
                    phi = orientimTemp(row,column);
                    deltaX = round(SarcSpacing*sin(phi));
                    deltaY = round(SarcSpacing*cos(phi));
                    for xx = -floor(dim*sqrt(2)+1):floor(dim*sqrt(2)+1)
                        for yy = -floor(dim*sqrt(2)+1):floor(dim*sqrt(2)+1)
                            if orientimTemp(row+deltaX+xx,column+deltaY+yy) ~= 0
                                if (abs(xx*sin(phi)+yy*cos(phi))+1 <= dim) && abs(xx*cos(phi)-yy*sin(phi))+1 <= dim
                                    sum1 = sum1 + abs(cos(orientimTemp(row+deltaX+xx,column+deltaY+yy)-phi));
                                    block1size = block1size + 1;
                                end
                            end
                            
                            if orientimTemp(row-deltaX+xx,column-deltaY+yy) ~= 0
                                if (abs(xx*sin(phi)+yy*cos(phi))+1 <= dim) && abs(xx*cos(phi)-yy*sin(phi))+1 <= dim
                                    sum2 = sum2 + abs(cos(orientimTemp(row-deltaX+xx,column-deltaY+yy)-phi));
                                    block2size = block2size + 1;
                                end
                            end
                        end
                    end
                end
                if  block1size > 0 && block2size > 0
                    if sum1/block1size > cosThresh && sum2/block2size > cosThresh
                        skeletonMask(row,column) = 1;
                    end
                end
            end
        end
        
                 show(binim_skel);
                 show(skeletonMask);
                 temporaryVariableWillVeverUseIt = input('Enjoy the view and then press Enter to continue: ');
        
        skeletonMask = skeletonMask.*(binim_skel > skeletonMask2);
        binim_skel = binim_skel.*skeletonMask;
        % Remove false sarcomeres at tissue borders selecting regions to be
        % masked. ROI - region of interest.
        index = 1;
        while index < 1;
            disp('Select ROI to exclude from further analysis');
            BW = roipoly(binim_skel);
            BW2 = ~BW;
            binim_skel = binim_skel.*BW2;
            % Use erode mask to remove false sarcomeres from binary skeleton
            show(binim_skel)
            b = input('Select another ROI to exclude? (yes = 0, no = 1): ');
            if b == 0
                disp('Image accepted' )
            else
                disp('Select another ROI...')
                index = 1;
            end
        end
        
        % Multiply orientation angles by the binary skeleton image
        orientim = orientim.*binim_skel;
        orientim_perp = orientim_perp.*binim_skel;
        % Convert 2D-array to 1D vector
        orientation = orientim(:);
        % Keep non-zero values only
        nonzero_orientation = orientation(find(orientation));
        
        % Correct for the fact that the orientation software counts radians
        % clockwise from dead-right as opposed to counter-clockwise
        %         for ii=1:1:length(nonzero_orientation)
        %             if nonzero_orientation(ii) < pi/2
        %                 shift = nonzero_orientation(i);
        %                 nonzero_orientation(ii) = pi - shift;
        %             elseif nonzero_orientation(ii) > pi/2
        %                 shift = pi - nonzero_orientation(i);
        %                 nonzero_orientation(ii) = shift;
        %             end
        %         end
        
        % Correct for the fact that the orientation is measured clockwise
        % instead of counter-clockwise + Rotate angles by pi/2 because
        % z-disks are orthogonal to the sarcomere orientation. Bring all
        % the angles to the [0;pi) range
        nonzero_orientation = 3/2*pi - nonzero_orientation;
        for ii = 1: length(nonzero_orientation)
            if nonzero_orientation(ii) >= pi
                nonzero_orientation(ii) = nonzero_orientation(ii) - pi;
            elseif nonzero_orientation(ii) < 0
                nonzero_orientation(ii) = nonzero_orientation(ii) + pi;
            end
        end
        
        
        % Convert radians to degrees
        nonzero_orientation_angles = rad2deg(nonzero_orientation);
        
        % hist(nonzero_orientation_angles)
        Mean = mean(nonzero_orientation_angles)
        Std = std(nonzero_orientation_angles)
        Median = median(nonzero_orientation_angles)
        
        % Create histogram
        [n,xout] = hist(nonzero_orientation_angles,180);
        title(picturename)
        dx = xout(2)-xout(1);                   % calc a single bin width
        n = n / sum( n*dx );                    % normalize histogram to have area of 1
        
        % Find mode
        [C,I] = max(n);
        Mode = xout(I)
        
        % Total number of sarcomere positive pixels in the skeleton image
        % Sarcomere density = total/(image area)
        Total = length(nonzero_orientation_angles)
        
        % Add data to SUM file
        All_angles = [All_angles;nonzero_orientation_angles];
        All_rad = [All_rad;nonzero_orientation];
        %close all
        
        PercentSarcomereArea = Total/(Size-20)*(Size-20);
        fprintf(fileID, [picture ';%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n'], PercentSarcomereArea, Mean, Std, Median, Mode);
    end
    
    Mean = mean(All_angles);
    Std = std(All_angles);
    Median = median(All_angles);
    OrientationOrderParameter = OOP(All_rad);
    
    % Create histogram
    [n,xout] = hist(All_angles,180);
    title('All_angles')
    dx = xout(2)-xout(1);                   % calc a single bin width
    n = n / sum( n*dx );                    % normalize histogram to have area of 1
    
    [Value,minInd] = min(smooth(n));
    shift = xout(minInd);
    All_angles = All_angles-xout(minInd);
    for ii = 1:length(All_angles)
        if All_angles(ii) < 0
            All_angles(ii) = All_angles(ii) +180;
        elseif All_angles(ii) >= 180
            All_angles(ii) = All_angles(ii) - 180;
        end
    end
    
    [n,xout] = hist(All_angles,180);
    title(picturename)
    dx = xout(2)-xout(1);                   % calc a single bin width
    n = n / sum( n*dx );                    % normalize histogram to have area of 1
    
    xout = xout + shift;
    
    % Plot histogram of raw orientation
    figure, bar(xout,n,'hist')              % plot normalized histogram
    xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
    title('All_angles')
    xlabel('Degrees')
    ylabel('Normalized Occurance');
    
    figure;
    %fitting 2 gaussians in the signal and subtracting linear function
    %for the noise
    [FitResults,LowestError,BestStart,xi,yi,BootResults] = peakfit(transpose([xout; n]),0,0,2,1,0,1,0,1);
    
    figure, plot(xi,[yi;sum(yi,1)]);
    xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
    title('approximation')
    xlabel('Degrees')
    ylabel('Normalized Occurance');
    
    % Find mode
    [C,I] = max(n);
    Mode = xout(I)
    
    Total = length(All_rad);
    PercentSarcomereArea = Total/(Size-20)*(Size-20);
    
    fprintf(fileID, 'Average;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f;%3.3f\n\n', PercentSarcomereArea, Mean, Std, Median, Mode, OrientationOrderParameter);
    fclose(fileID);
    
end

%clear

