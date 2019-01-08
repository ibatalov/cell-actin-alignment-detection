% detects non-zero orientations from an actin image
% input:
%       actin - the image to be processed
%       blsze - block size that is used to detect a single orientation
%       thresh - threshold corresponding to the dragient arcoss the fiber.
%       Should between 0 and 1.
%       Size - image size
%       sarcomereMask - a mask to be used to filter out any orientations
%       outside of it. If no mask needed, input '1'
% output: all actin orientations that passed the threshold and masking
function [ nonzero_orientation, mask ] = actinDetectSlice(actin,blksze,thresh,Size,sarcomereMask)
%actinDetectSlice Returns nonzero_orientation_angles for 1 slice
%   Performs actinDetect on a single Slice and returns the list of non zero
%   angles                    
        
        % Identify ridge-like regions and normalise image
        % disp('Normalizing Image and Creating Mask' )
        [normim, mask] = ridgesegment(actin, blksze, thresh);
        %show(mask,1);
        %show(normim,2);
        
        % Create skeleton image from normalized image
        %normim_bin = normim > 0;
        %mask_bin = mask > 0;
        %normim_mask = normim_bin.*mask_bin;
        %norm_bin_skel = bwmorph(normim_mask,'skel',Inf);
        %show(norm_bin_skel, 3);
        
        % Determine ridge orientations
        % disp('Calculating Ridge Orientations' )
        [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
        %plotridgeorient2(orientim, 20, actin, 1000) % made figure number 1000 so it doesn't overwrite any other open figure
        
        % draw masked orientation image for fibroblasts and cardiomyocytes
%         masks = cell(2,1);
%         masks{1} = (~sarcomereMask).*mask;
%         masks{2} = sarcomereMask.*mask;
%         plotridgeorient_masked(orientim, 20, actin, masks, {'-g', '-r'});
        
        %show(reliability,5)
        
        % Only keep orientation values with a reliability greater than 0.5
        reliability_binary = reliability>0.5;
        
        % Remove 10 pixel wide border where orientation values are not accurate
        reliability_binary(:,1:1:10) = 0;
        reliability_binary(1:1:10,:) = 0;
        reliability_binary(:,Size-10:1:Size) = 0;
        reliability_binary(Size-10:1:Size,:) = 0;
        
        %Combine masks for cardiomyocytes and ridge blocks
        mask = mask.*sarcomereMask;
        
        % Multiply orientation angles by the binary mask image to remove
        % data where there are no cardiomyocytes
        newmask = mask.*reliability_binary;
        %Sshow(newmask,3)
        orientim = orientim.*newmask;
        
        % Convert 2D-array to 1D vector
        orientation = orientim(:);
        
        % Keep non-zero values only
        nonzero_orientation = orientation(orientation ~=0);
                
end

