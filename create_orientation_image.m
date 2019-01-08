function [ orientim ] = create_orientation_image(actin,blksze,thresh,sarcomereMask)
%actinDetectSlice Returns nonzero_orientation_angles for 1 slice
%   Performs actinDetect on a single Slice and returns the list of non zero
%   angles                    
        
        % Identify ridge-like regions and normalise image
        % disp('Normalizing Image and Creating Mask' )
        [normim, mask] = ridgesegment(actin, blksze, thresh);
        %show(mask,1);
        %show(normim,2);
        
        % Create skeleton image from normalized image
        normim_bin = normim > 0;
        mask_bin = mask > 0;
        normim_mask = normim_bin.*mask_bin;
        %norm_bin_skel = bwmorph(normim_mask,'skel',Inf);
        %show(norm_bin_skel, 3);
        
        % Determine ridge orientations
        % disp('Calculating Ridge Orientations' )
        [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
        plotridgeorient2(orientim, 10, actin, 4)
        %show(reliability,5)
        
        % Only keep orientation values with a reliability greater than 0.5
        reliability_binary = reliability>0.5;
        
        Size1 = size(reliability_binary,1);
        Size2 = size(reliability_binary,2);
        
        % Remove 10 pixel wide border where orientation values are not accurate
        reliability_binary(:,1:1:10) = 0;
        reliability_binary(1:1:10,:) = 0;
        reliability_binary(:,Size2-10:1:Size2) = 0;
        reliability_binary(Size1-10:1:Size1,:) = 0;
        
        %Combine masks for cardiomyocytes and ridge blocks
        mask = mask.*sarcomereMask;
        
        size(mask);
        size(reliability_binary);
        % Multiply orientation angles by the binary mask image to remove
        % data where there are no cardiomyocytes
        newmask = mask.*reliability_binary;
        %Sshow(newmask,3)
        orientim = orientim.*newmask;        
end

