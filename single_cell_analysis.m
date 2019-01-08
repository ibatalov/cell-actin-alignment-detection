%% Request data for analysis
SampleInfo = inputdlg({'Number of samples:','Excel file title:', 'nuclei layer number', 'actin layer number:','alpha-actinin layer number:'},'Input',1,{'1','single cell analysis','1','2','3'});
SampleCount = str2num(SampleInfo{1});

dapiLayerNumber = str2double(SampleInfo{3});
actinLayerNumber = str2double(SampleInfo{4});
actininLayerNumber = str2double(SampleInfo{5});

rootfolder = path;

%% get files for each sample
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
fprintf(fileID, 'picture;cell_type;cell_area µm^2;big_semi_axis µm^2;small_semi_axis µm^2;aspect_ratio;ellipticity;orientation deg.;\n');

%% Loop for each sample
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
    
    %% Loop over each picture for sample #i
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
        
        dapi = img{1,1}{dapiLayerNumber,1};
        actin = img{1,1}{actinLayerNumber,1};
        sarcomeres = img{1,1}{actininLayerNumber,1};
        
        info = img{1,4};  % Load OME metadata
        PixelSize = str2double(info.getPixelsPhysicalSizeX(0));
        Size = str2double(info.getPixelsSizeX(0));
        
        fprintf('Image resolution: %5.5f px/µm\n\n', 1/PixelSize);
        
        % Define blksze (should be 3/pixelsize based on recommendation)
        blksze = floor(3/PixelSize);
        
        % take into account only actin that is near sarcomeres
        
        % this part is useless now
        if ~isa(sarcomeres, 'double'),
            sarcomeres = double(sarcomeres);
            actin = double(actin);
            dapi = double(dapi);
        end
        
        %% making actin and sarcomere masks
        thresSarc = (0.07*max(sarcomeres(:))+0.93*min(sarcomeres(:)));
        threshAct = (0.05*max(actin(:))+0.95*min(actin(:)));
        d = floor(1.3*blksze)/blksze;
        
        actin = padarray(actin, [round(2*d) round(2*d)]);
        sarcomeres = padarray(sarcomeres, [round(2*d) round(2*d)]);
        dapi = padarray(dapi, [round(2*d) round(2*d)]);
        
        [imSizeX,imSizeY] = size(sarcomeres);

        disk = strel('disk', d, 0); % used to merge sarcomeres into mask by expanding and shrinking the image
        actin_disk = strel('disk', floor(d/2), 0); % used to merge actin into mask
        
        maxActinArea = (imSizeX - 2*0.6*floor(d/2))*(imSizeY - 2*floor(d/2));
        maxSarcomereArea = (imSizeX - 2*0.6*d)*(imSizeY - 2*d);
        
        sarcThreshold = (0.1*max(sarcomeres(:))+0.9*min(sarcomeres(:)));
        sarcomereMask = bwmorph(sarcomeres > sarcThreshold, 'open');
        sarcomeres_border = bwmorph(sarcomereMask, 'remove');
        % dilate and erode to fill the holes between z-lines of
        % sarcomeres
        sarcomereMask = sarcomereMask | imdilate(sarcomeres_border, disk);
        % erode only 60%, because cells' area is higher than sarcomere coverage
        sarcomereMask = bwmorph(sarcomereMask,'erode',floor(d*0.6));
        
        actinThreshold = (0.05*max(actin(:))+0.95*min(actin(:)));
        actinMask = bwmorph(actin > actinThreshold, 'open');
        actin_border = bwmorph(actinMask, 'remove');
        
        % dilate and erode to fill the holes between actin filaments
        actinMask = actinMask | imdilate(actin_border, actin_disk);
        actinMask = bwmorph(actinMask,'erode',floor(floor(d/2)*0.6));
        actinMask = actinMask | sarcomereMask; % expand actin mask to include sarcomere mask
        
        %%
        % voodoo magic to show the nice-ish image
        processedDapi = dapi./max(dapi(:));
        processedDapi = medfilt2(processedDapi, [3 3]);
        processedDapi = (processedDapi.^0.3);
        processedDapi = processedDapi.*(processedDapi > 0.1);
        
        dapiMask = bwmorph(processedDapi > 0, 'erode', 3);
        dapiMask = bwmorph(dapiMask > 0, 'dilate', 3);
        processedDapi = processedDapi.*dapiMask;
        
        initialMagnification = 60; % Depends on your monitor resolution
        
        %% Use polygonal tool to split cells. To end the loop, select 100% black region (double click on the image also works).
        f = figure;
        isZeroPxSelected = 1;
        while (sum(isZeroPxSelected(:)) ~= 0);
            processedActin = (actin/max(actin(:)));
            processedActin = processedActin.*(processedActin > 0.05);
            processedActin = processedActin.^0.5;
            processedActin(processedDapi > 0) = 0;
            combinedImage = cat(3, actinMask.*(processedDapi == 0)./2, processedActin, processedDapi);
            imshow(combinedImage,'InitialMagnification', initialMagnification);
            actRemove = impoly();
            actRemove = createMask(actRemove);
            isZeroPxSelected = actRemove.*actinMask;
            actinMask = actinMask.*(~actRemove);
        end
        close(f);
        
        
        %% Removes cells that didn't fit the image area
        %  Aslo removes small areas caused by noise
        
        border = d + 3; % border from the edge of the image within wich I don't want any cells
        
        cc = bwconncomp(actinMask, 4);
        imSize = [imSizeX, imSizeY];
        
        for index = 1 : length(cc.PixelIdxList)
            linInd = cell2mat(cc.PixelIdxList(index));
            [rows, columns] = ind2sub(imSize, linInd);
            if (min([rows; columns]) <= border) || (max([rows; columns]) >= imSizeX - border)
                actinMask(linInd) = 0;
            end
            
            if numel(linInd) < (15/PixelSize)^2; % threshold area 20 µm^2
                actinMask(linInd) = 0;
            end
            
        end
        
        %% Final actin mask clean up (if necessary). To end the loop, select 100% black region (double click on the image also works).
        f = figure;
        isZeroPxSelected = 1;
        
        while (sum(isZeroPxSelected(:)) ~= 0);
            imshow(actinMask);
            actRemove = impoly();
            actRemove = createMask(actRemove);
            isZeroPxSelected = actRemove.*actinMask;
            actinMask = actinMask.*(~actRemove);
        end
        
        close(f);
        
        
        %%
        cc = bwconncomp(actinMask, 4);
        cellCount = numel(cc.PixelIdxList);
        cardioCount = 0;
        fibroCount = 0;
        
        clearvars cardio fibro;
        
        cardio(:,:,1) = zeros(imSizeX,imSizeY);
        fibro(:,:,1) = zeros(imSizeX,imSizeY);
        
        % determine the type of each cell and put them in a corresponding
        % array
        for index = 1 : length(cc.PixelIdxList)
            linInd = cell2mat(cc.PixelIdxList(index));
            cellActin = zeros(imSizeX, imSizeY);
            cellActin(linInd) = 1;
            
            actinArea = numel(linInd);
            intersection = cellActin.*sarcomereMask;
            sarcArea = nnz(intersection);
            
            if sarcArea/actinArea > 0.5
                cardioCount = cardioCount + 1;
                cardio(:,:,cardioCount) = cellActin;
            else
                fibroCount = fibroCount + 1;
                fibro(:,:,fibroCount) = cellActin;
            end
        end
        
        allCellList = cat(3, cardio, fibro);
        
        %% Calculate Parameters for each cell
        show(actinMask);
        hold on;
        for number = 1 : size(allCellList, 3)
            area = nnz(allCellList(:,:,number)); % mask area ~ cell area
            [rows, columns] = find(allCellList(:,:,number));
            
            centerRow = sum(rows)/area; % it is not necessarily integer!
            centerColumn = sum(columns)/area; % it is not necessarily integer!
            
            % shifting to the center of inertia
            rows1 = rows - centerRow;
            cols1 = columns - centerColumn;
            
            Irr = sum(cols1.*cols1); % row = y. Irr ~ Iy
            Icc = sum(rows1.*rows1); % column = x. Icc = Ix
            Irc = -sum(rows1.*cols1); % Ixy = Iyx
            
            inertiaTensor = [Irr, Irc; Irc, Icc];
            [eigenVectors, diagMatrix] = eig(inertiaTensor);
            
            % orientation - angle between the vector and the regular
            % x-axis
            
            if diagMatrix(1,1) > diagMatrix(2,2)
                maxMoment = diagMatrix(1,1);
                minMoment = diagMatrix(2,2);
                mainDirection = eigenVectors(:,2);
            else
                maxMoment = diagMatrix(2,2);
                minMoment = diagMatrix(1,1);
                mainDirection = eigenVectors(:,1);
            end
            
            % a - big semi-axis of the fitted ellipse
            % b - small semi-axis of the fitted ellipse
            a = sqrt(area/pi*sqrt(maxMoment/minMoment));
            b = sqrt(area/pi*sqrt(minMoment/maxMoment));
            
            orientation = atan(-mainDirection(1)/mainDirection(2))*180/pi;
            if orientation < 0
                orientation = orientation + 180;
            end
            ellipticity = (a - b)/a;
            aspectRatio = a/b;
            color = 'g';
            % quiver(x,y,u,v);
            % red - cardiomyocytes, green - fibroblasts
            if(number <= cardioCount)
                quiver(centerColumn, centerRow, mainDirection(2)*a, mainDirection(1)*a, '-r', 'linewidth',2);
                color = 'r';
                cellType = 'cardiomyocyte';
            else
                quiver(centerColumn, centerRow, mainDirection(2)*a, mainDirection(1)*a, '-g', 'linewidth',2);
                cellType = 'fibroblast';
            end;
            ellipse(a,b,atan(mainDirection(1)/mainDirection(2)), centerColumn, centerRow, color);
            
            % converts all the relevant stuff to µm (area, inertia moments)
            micronArea = area*(PixelSize^2);
            micronA = a*PixelSize;
            mirconB = b*PixelSize;
            
            fprintf(fileID, [picture ';' cellType ';%3.0f;%3.0f;%3.0f;%3.2f;%3.3f;%3.1f\n'], micronArea, micronA, mirconB, aspectRatio, ellipticity, orientation);
        end
        hold off;   
    end
end

fclose('all');
clearvars thresh;
