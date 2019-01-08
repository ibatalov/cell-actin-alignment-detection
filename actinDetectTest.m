%   Determines 'thresh' parameter for actinDetectSlice 
%   Runs actinDetect on 1 slice to determine optimal threshold for actin
%   filament detection
% input:
%       image - the image to be processed
%       blksze - block size that is used to detect a single orientation
%       PixelSize - size of 1 pixel in microns
function [ thresh ] = actinDetectTest( image ,blksze, PixelSize )

thresh = 0; 

% Identify ridge-like regions and normalise image
   index = 0;
   while index < 1;
%      % Threshold of standard deviation to decide if a block is a ridge region
       thresh = input('Enter Threshold (0.1 - 0.2): ');
       disp('Normalizing Image and Creating Mask' )
       [normim, mask] = ridgesegment(image, blksze, thresh);
       show(normim,1);
       show(mask, 2);
     
       % Determine if normalization and mask look good, click on image to
       % accept or press any key to enter new values
       w = input('Accept Threshold (yes = 1, no = 0): ');
       if w == 1
           disp('Image threshold accepted' )
           index = 1;
       else
           disp('Re-analyze imaging...')
       end
       
   end
   
   return 
   
end

