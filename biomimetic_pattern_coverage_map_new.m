%load('D:\Ivan\2015\2015-05-07 Matlab heat map stuff\first results\combined data\contour pattern\data.mat', 'fn_pattern','PixelSize');
load('/Users/ivan/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/first results/combined data/contour pattern/data.mat', 'fn_pattern','PixelSize');

%orig_pattern = imread('D:\Ivan\2015\2015-05-07 Matlab heat map stuff\original pattern\BM_pattern_for_heat_map_contour_1.png');
orig_pattern = imread('/Users/ivan/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/original pattern/BM_pattern_for_heat_map_contour_3.png');
orig_pattern = rgb2gray(orig_pattern);
orig_pattern = double(orig_pattern);
orig_pattern = orig_pattern/max(orig_pattern(:));
orig_pattern = imresize(orig_pattern, size(fn_pattern));

%bw_pattern = imread('D:\Ivan\2015\2015-05-07 Matlab heat map stuff\original pattern\BM_pattern_for_heat_map_big_pattern_solid.png');
bw_pattern = imread('/Users/ivan/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/original pattern/BM_pattern_for_heat_map_big_pattern_solid.png');
bw_pattern = rgb2gray(bw_pattern) == 0;

while 1
    imshow(bw_pattern);
    disp('choose 2 pattern corners' );
    [movingPointsX, movingPointsY] = getpts;
    close;
    pointsEntered = numel(movingPointsX);
    
    if pointsEntered == 2 || pointsEntered == 3
        if(pointsEntered == 3)
            movingPoints = [movingPointsX(numel(movingPointsX) - 1 : numel(movingPointsX)), movingPointsY(numel(movingPointsY) - 1 : numel(movingPointsY))];
            bw_pattern_origin = [movingPointsX(1), movingPointsY(1)];
        else
            movingPoints = [movingPointsX, movingPointsY];
            bw_pattern_origin = movingPoints(1,:);
        end
        
        distance = sqrt((movingPoints(1,1) - movingPoints(2,1))^2 + (movingPoints(1,2) - movingPoints(2,2))^2);
        scale_factor = size(orig_pattern,1)/distance;
        transform = affine2d([scale_factor 0 0; 0 scale_factor 0; 0 0 1]);
    else
        if numel(movingPointsX) == 1
            movingPoints = [movingPointsX, movingPointsY];
            bw_pattern_origin = movingPoints(1,:);
            transform = affine2d([1 0 0; 0 1 0; 0 0 1]);
        else
            discardImage = 1;
        end
    end
    
    bw_pattern_origin = transformPointsForward(transform, bw_pattern_origin);
    
    [outboundsX, outboundsY] = outputLimits(transform,[1 size(bw_pattern, 2)],[1 size(bw_pattern, 1)]);
    bw_pattern_origin(1) = bw_pattern_origin(1) - (outboundsX(1) - 1);
    bw_pattern_origin(2) = bw_pattern_origin(2) - (outboundsY(1) - 1);
    
    new_bw_pattern = imwarp(bw_pattern, transform);
    overlay_patterns(orig_pattern, new_bw_pattern, bw_pattern_origin);
    
    disp('red: first pattern, green: current pattern. If patterns do not match, choose shift vector (first point - new pattern, second - old)');
    [inputX, inputY] = getpts;
    while numel(inputX) == 2
        close;
        bw_pattern_origin(1) = bw_pattern_origin(1) + (inputX(2) - inputX(1));
        bw_pattern_origin(2) = bw_pattern_origin(2) + (inputY(2) - inputY(1));
        overlay_patterns(orig_pattern, new_bw_pattern, bw_pattern_origin);
        [inputX, inputY] = getpts;
    end
    if isempty(input('Type any key to fix the pattern manually or Enter to continue'))
        break;
    end
end

filterSize = round(size(fn_pattern,2)/15);

fn_coverage_map = filter2(ones(filterSize)/filterSize^2,new_bw_pattern);
fn_coverage_map = fn_coverage_map(round(bw_pattern_origin(2)):round(bw_pattern_origin(2))+size(fn_pattern,1)-1,round(bw_pattern_origin(1)):round(bw_pattern_origin(1))+size(fn_pattern,2)-1);

low_hue = 240/360;
high_hue = 0;
palitra_length = 800;
font_size = 20;

max_value = round(max(fn_coverage_map(:)),4, 'significant');
min_value = round(min(fn_coverage_map(:)),4, 'significant');
avg_value = round(sum(new_bw_pattern(:))/length(new_bw_pattern(:)),4, 'significant');

if(false)
    cut_out_fraction = 0.5; % fraction of data that you want to make outside the color range (over or undersaturated) to increase contrast of what is near the average
    half_range = min(avg_value - min_value, max_value - avg_value)*(1 - cut_out_fraction);
    min_thresh = avg_value - half_range;
    max_thresh = avg_value + half_range;
    fn_coverage_map = (fn_coverage_map - min_thresh).*(fn_coverage_map > min_thresh) + min_thresh; %undersaturate data
    fn_coverage_map = (fn_coverage_map - max_thresh).*(fn_coverage_map < max_thresh) + max_thresh; %oversaturate data
    
    max_value = round(max(fn_coverage_map(:)),4, 'significant');
    min_value = round(min(fn_coverage_map(:)),4, 'significant');
end
hsv = zeros(size(orig_pattern,1), size(orig_pattern,2),3);
hsv(:,:,1) = low_hue + (high_hue - low_hue)*(fn_coverage_map - min(fn_coverage_map(:)))/(max(fn_coverage_map(:)) - min(fn_coverage_map(:)));
hsv(:,:,2) = 1;
hsv(:,:,3) = orig_pattern*0.8 + 0.1;
%     figure;
%     imshow(OOP_map);
%     figure
%     imshow(hsv2rgb(hsv));

scale_hsv = zeros(palitra_length, 120,3);
scale_hsv(:,:,3) = 1;
scale_hsv(:,:,2) = 0;
scale_hsv(:,5:40,2) = 1;

for row = 1 : palitra_length
    scale_hsv(row,5:40,1) = high_hue - (high_hue - low_hue)*(row/palitra_length);
end

min_y = 1;
max_y = palitra_length;
avg_y = round(1 + (palitra_length - 1)*(max_value - avg_value)/(max_value - min_value));

min_text_y = 1;
max_text_y = max_y - font_size - 1;
avg_text_y = avg_y - round(font_size/2);

scale_hsv(max_y,5:43,3) = 0;
scale_hsv(min_y,5:43,3) = 0;
scale_hsv(avg_y,5:43,3) = 0;

textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', font_size);
strings = uint8([num2str(max_value) 0 num2str(min_value) 0 num2str(avg_value)]);
labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([45, min_text_y; 45, max_text_y; 45, avg_text_y]));
%     figure;
%     imshow(labeled);

scale_bar_microns = floor(120*PixelSize/10)*10;
scale_bar_px = scale_bar_microns/PixelSize;
left_margin = round((120 - scale_bar_px)/2);

scale_bar = ones(120,120,3);
scale_bar(80:100, left_margin:round(left_margin + scale_bar_px), :) = 0;
string = uint8([num2str(scale_bar_microns) ' um']);
textInserter2 = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', 30);
scale_bar = step(textInserter2, scale_bar, string, int32([16 40]));
%     figure;
%     imshow(scale_bar);

final_image = ones(size(fn_pattern,1), size(fn_pattern,2) + size(labeled,2),3);
final_image(:,1:size(fn_pattern,2),:) = hsv2rgb(hsv);
final_image(1:size(labeled,1),size(fn_pattern,2) + 1:size(fn_pattern,2) + size(labeled,2),:) = labeled;
final_image(size(final_image,1) - size(scale_bar,1) + 1:size(final_image,1),size(final_image,2) - size(scale_bar,2) + 1:size(final_image,2),:) = scale_bar;
figure;
finalFigure = imshow(final_image);
title('Fibronectin density average');

%%
LASTN = maxNumCompThreads(1)
corrData = zeros(round(size(fn_pattern,2)) - 1,2);
for smoothing = 1 : 10 : round(size(fn_pattern,2));
    filterSize = smoothing;
    fn_coverage_map = filter2(ones(filterSize)/filterSize^2,new_bw_pattern);
    fn_coverage_map = fn_coverage_map(round(bw_pattern_origin(2)):round(bw_pattern_origin(2))+size(fn_pattern,1)-1,round(bw_pattern_origin(1)):round(bw_pattern_origin(1))+size(fn_pattern,2)-1);
    corrData(smoothing, 1) = filterSize*PixelSize;
    corrData(smoothing, 2) = corr2(fn_coverage_map, cell_coverage_map);
end
figure
plot(corrData(sum(corrData,2)>0, 1),corrData(sum(corrData,2)>0,2));
ylabel('Correlation coefficient');
xlabel('averaging diameter, µm');
xlim([0,155]);
ylim([0,1]);
ax = gca;