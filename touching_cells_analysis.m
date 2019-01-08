%% Request data for analysis
SampleInfo = inputdlg({'Number of samples:','Excel file title:', 'nuclei layer number', 'actin layer number:','alpha-actinin layer number:'},'Input',1,{'1','touching cell analysis','1','2','3'});
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
fprintf(fileID, 'picture;cell_type;cell_area µm^2;big_semi_axis µm^2;small_semi_axis µm^2;aspect_ratio;ellipticity;orientation deg.;connected angle;connected cell\n');

%% Loop for each sample
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
        
        thresSarc = (0.07*max(sarcomeres(:))+0.93*min(sarcomeres(:)));
        threshAct = (0.05*max(actin(:))+0.95*min(actin(:)));
        d = floor(1.3*blksze)/blksze;
        
        actin = padarray(actin, [round(2*d) round(2*d)]);
        sarcomeres = padarray(sarcomeres, [round(2*d) round(2*d)]);
        dapi = padarray(dapi, [round(2*d) round(2*d)]);
        
        [imSizeX,imSizeY] = size(sarcomeres);
        sarcomereMask = zeros(imSizeX, imSizeY);
        actinMask = zeros(imSizeX, imSizeY);
        
        %disp('making sarcomere mask');
        %%
        circle = zeros(2*d*blksze+1);
        for x = -d*blksze:d*blksze
            for y = -d*blksze:d*blksze
                circle(x+d*blksze+1,y+d*blksze+1) = x*x+y*y < d*blksze*d*blksze;
            end
        end
        
        for x = d*blksze+1:imSizeX-d*blksze
            for y = d*blksze+1:imSizeY-d*blksze
                if sarcomeres(x,y) > thresSarc
                    sarcomereMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) = sarcomereMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) | circle;
                end
                
                if actin(x,y) > threshAct
                    actinMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) = actinMask(x-d*blksze:x+d*blksze, y-d*blksze:y+d*blksze) | circle;
                end
            end
        end
        
        sarcomereMask = bwmorph(sarcomereMask,'erode',d);
        actinMask = bwmorph(actinMask,'erode',d);
        
        %% save the objects from the initial mask to use for cell connectivity detection
        cc = bwconncomp(actinMask, 4);
        savedShapes = cc.PixelIdxList;
        
        
        %%
        % voodoo magic to show the nice-ish image
        processedDapi = dapi./max(dapi(:));
        processedDapi = medfilt2(processedDapi, [3 3]);
        processedDapi = (processedDapi.^0.3);
        %processedDapi = processedDapi.*(processedDapi > 0.1);
        
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
        
        %% remove small areas caused by noise
        
        border = d + 3; % border from the edge of the image within wich I don't want any cells
        
        cc = bwconncomp(actinMask, 4);
        imSize = [imSizeX, imSizeY];
        borderCells = zeros(imSize);
        
        for index = 1 : length(cc.PixelIdxList)
            linInd = cell2mat(cc.PixelIdxList(index));
            [rows, columns] = ind2sub(imSize, linInd);
            
            if (min([rows; columns]) <= border) || (max([rows; columns]) >= imSizeX - border)
                borderCells(linInd) = 1;
            end
            if numel(linInd) < (15/PixelSize)^2; % threshold area 20 µm^2
                actinMask(linInd) = 0;
            end
            
        end
        
        %% Final actin mask clean up (if necessary). To end the loop, select 100% black region (double click on the image also works).
        f = figure;
        isZeroPxSelected = 1;
        
        while (sum(isZeroPxSelected(:)) ~= 0);
            imshow(actinMask,'InitialMagnification', initialMagnification);
            actRemove = impoly();
            actRemove = createMask(actRemove);
            isZeroPxSelected = actRemove.*actinMask;
            actinMask = actinMask.*(~actRemove);
        end
        
        close(f);
        
        
        %%
        cc = bwconncomp(actinMask, 4);
        cellCount = numel(cc.PixelIdxList);
        
        % cellInfo(row,column)
        % row = cell number
        % col 1 - cell type: 0 = firboblast, 1 = cardiomyocyte
        % col 2 - originating object (2 cells with
        % the same orig. object are touching)
        cellInfo = zeros(length(cc.PixelIdxList),2);
        cellList = cc.PixelIdxList;
        
        % determine the type of each cell
        for index = 1 : length(cc.PixelIdxList)
            linInd = cell2mat(cc.PixelIdxList(index));
            cellActin = zeros(imSizeX, imSizeY);
            cellActin(linInd) = 1;
            
            actinArea = numel(linInd);
            intersection = cellActin.*sarcomereMask;
            sarcArea = nnz(intersection);
            
            if sarcArea/actinArea > 0.5
                cellInfo(index,1) = 1;
            else
                cellInfo(index,1) = 0;
            end
            
            for obj = 1 : length(savedShapes)
                
                objInd = cell2mat(savedShapes(obj));
                objPixels = zeros(imSizeX, imSizeY);
                objPixels(objInd) = 1;
                
                intersectionObj = cellActin.*objPixels;
                if nnz(intersectionObj) > 0.5*actinArea
                    cellInfo(index,2) = obj;
                    break;
                else if obj == length(savedShapes)
                        show('the originating cluster was not found. SOMETHING IS SUPER-DUPER WRONG HERE!');
                    end
                end
                
            end
        end
        % this matrix shows what cells are connected to other cells
%         connectionMatrix = zeros(length(cellInfo(:,1)));
%         connectionSum = zeros(length(cellInfo(:,1)), 1);
%         
%         for M = 1: length(connectionMatrix)
%             for N = 1 : length(connectionMatrix)
%                 connectionMatrix(M,N) = cellInfo(M,2) == cellInfo(N,2);
%             end
%             connectionSum(M) = sum(connectionMatrix(M,:));
%         end
%         
%         
%         
%         % we want only the groups of at least 2 cells (2 connections - one with
%         % itself, one with the neighbor)
%         cellList = cellList(connectionSum > 1);
%         cellInfo = cellInfo(connectionSum > 1,:);
%         
        clusters = cell(length(savedShapes), 1); % indices of cells belonging to each cluster
        for cNum = 1 : length(savedShapes)
            clusters{cNum,1} = find(cellInfo(:,2) == cNum);
        end
        
        % shows what cells withing the cluster are close to each other
        proximityTable = -ones(length(cellInfo));
        
        %% Calculate Parameters for each cell
        actinMask = actinMask.*(~borderCells);
        show(actinMask);
        hold on;
        for number = 1 : length(cellList)
            proximityTable(number,number) = 0;
            
            linInd = cell2mat(cellList(number));
            cellActin = zeros(imSizeX, imSizeY);
            cellActin(linInd) = 1;
            
            area = numel(linInd); % mask area ~ cell area
            [rows, columns] = ind2sub(imSize, linInd);
            
            centerRow = sum(rows)/area; % it is not necessarily integer!
            centerColumn = sum(columns)/area; % it is not necessarily integer!
            
            if (min([rows; columns]) <= border) || (max([rows; columns]) >= imSizeX - border)
                actinMask(linInd) = 0;
            else
                
                numberOfConnections = 0;
                connectedCell = -1; % 0 - fibroblast, 1 - cardiomyocyte
                connectedAngle = 0;
                neighbors = clusters{cellInfo(number,2),1};
                neighbors = neighbors(neighbors > 0);
                
                for n = 1 : length(neighbors)
                    
                    if(proximityTable(number, neighbors(n)) == -1)
                        neighborInd = cell2mat(cellList(neighbors(n)));
                        neighborActin = zeros(imSizeX, imSizeY);
                        neighborActin(neighborInd) = 1;
                        
                        dilatedCell = bwmorph(cellActin, 'dilate', 3.5/PixelSize); % 3 µm dilation
                        dilatedNeighbor = bwmorph(neighborActin, 'dilate', 3.5/PixelSize); % 3 µm dilation
                        
                        intersection = dilatedCell.*dilatedNeighbor;
                        
                        if nnz(intersection) > 0
                            intArea = nnz(intersection);
                            [r,c] = ind2sub(imSize, find(intersection));
                            
                            ri = sum(r)/intArea;
                            ci = sum(c)/intArea;
                            
                            vectR = ri - centerRow;
                            vectC = ci - centerColumn;
                            vectLength = sqrt(vectR^2 + vectC^2);
                            vectR = vectR/vectLength;
                            vectC = vectC/vectLength;
                            sinA = -vectR;
                            cosA = vectC;
                            angleA = 0;
                            
                            if sinA >= 0
                                angleA = acos(cosA);
                            else
                                if cosA < 0
                                    angleA = pi - asin(sinA);
                                else
                                    angleA = 2*pi + asin(sinA);
                                end
                            end
                            
                            if(angleA == 0)
                                angleA = 2*pi;
                            end
                            
                            if(angleA < 0)
                                show('Ivan, you messed up the angle');
                            end
                            
                            proximityTable(number,neighbors(n)) = angleA;
                           % proximityTable(neighbors(n),number) = angleA;
                        else
                            proximityTable(number,neighbors(n)) = 0;
                            proximityTable(neighbors(n),number) = 0;
                        end
                    end
                    
                    if(proximityTable(number, neighbors(n)) > 0)
                        numberOfConnections = numberOfConnections + 1;
                        connectedCell = cellInfo(neighbors(n),1);
                        connectedAngle = proximityTable(number,neighbors(n));
                    end
                end
                
                if numberOfConnections == 1 || numberOfConnections == 0
                    
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
                    a = sqrt(5*maxMoment/area);
                    b = sqrt(5*minMoment/area);
                    
                    orientation = atan(-mainDirection(1)/mainDirection(2))*180/pi;
                    if orientation < 0
                        orientation = orientation + 180;
                    end
                    ellipticity = (a - b)/a;
                    aspectRatio = a/b;
                    color = 'g';
                    % quiver(x,y,u,v);
                    % red - cardiomyocytes, green - fibroblasts
                    if(cellInfo(number,1) == 1)
                        quiver(centerColumn, centerRow, mainDirection(2)*a, mainDirection(1)*a, '-r', 'linewidth',2);
                        color = 'r';
                        cellType = 'cardiomyocyte';
                    else if cellInfo(number,1) == 0
                            quiver(centerColumn, centerRow, mainDirection(2)*a, mainDirection(1)*a, '-g', 'linewidth',2);
                            cellType = 'fibroblast';
                        end
                    end;
                    ellipse(a,b,atan(mainDirection(1)/mainDirection(2)), centerColumn, centerRow, color);
                    
                    % converts all the relevant stuff to µm (area, inertia moments)
                    micronArea = area*(PixelSize^2);
                    micronA = a*PixelSize;
                    mirconB = b*PixelSize;
                    
                    if(connectedCell == 0)
                        connectedCell = 'fibroblast';
                    else
                        if(connectedCell == 1)
                            connectedCell = 'cardiomyocyte';
                        else
                            connectedCell = 'none';
                        end
                    end
                    
                    connectedAngle = connectedAngle*180/pi;
                    
                    fprintf(fileID, [picture ';' cellType ';%3.0f;%3.0f;%3.0f;%3.2f;%3.3f;%3.1f;%3.0f;' connectedCell '\n'], micronArea, micronA, mirconB, aspectRatio, ellipticity, orientation, connectedAngle);
                end
            end
        end
        hold off;
        
    end
end

fclose('all');
clearvars thresh;
