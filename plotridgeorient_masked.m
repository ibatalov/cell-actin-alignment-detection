% PLOTRIDGEORIENT - plot of ridge orientation data
%
% Usage:   plotridgeorient(orient, spacing, im, figno)
%
%        orientim - Ridge orientation image (obtained from RIDGEORIENT)
%        spacing  - Sub-sampling interval to be used in ploting the
%                   orientation data the (Plotting every point is
%                   typically not feasible) 
%        im       - Optional fingerprint image in which to overlay the
%                   orientation plot.
%        figno    - Optional figure number for plot
%
% A spacing of about 20 is recommended for a 500dpi fingerprint image
%
% See also: RIDGEORIENT, RIDGEFREQ, FREQEST, RIDGESEGMENT

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005

function plotridgeorient_masked(orient, spacing, im, mask, line_specs)

    if fix(spacing) ~= spacing
        error('spacing must be an integer');
    end
    
    [rows, cols] = size(orient);
    
    lw = 2;             % linewidth
    len = 0.8*spacing;  % length of orientation lines

    % Subsample the orientation data according to the specified spacing

    s_orient = orient(spacing:spacing:rows-spacing, ...
		      spacing:spacing:cols-spacing);

    xoff = len/2*cos(s_orient);
    yoff = len/2*sin(s_orient);    
    
	show(im); hold on
    
    % Determine placement of orientation vectors
    [x,y] = meshgrid(spacing:spacing:cols-spacing, ...
		     spacing:spacing:rows-spacing);
    
    x = x-xoff;
    y = y-yoff;
    
    % Orientation vectors
    u = xoff*2;
    v = yoff*2;
    
    % remove points based on the mask
    for ind = 1 : length(mask)
        mask_mesh = mask{ind};
        mask_mesh = mask_mesh(spacing:spacing:rows-spacing, ...
            spacing:spacing:cols-spacing) > 0;
        x_pass = x(mask_mesh);
        y_pass = y(mask_mesh);
        u_pass = u(mask_mesh);
        v_pass = v(mask_mesh);
        
        quiver(x_pass,y_pass,u_pass,v_pass,0,line_specs{ind},'linewidth',2);
    end
    
    axis equal, axis ij,  hold off
    