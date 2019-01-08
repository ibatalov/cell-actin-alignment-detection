function final_image =  draw_heat_map(map_data, orig_pattern, min_value, max_value, avg_value, PixelSize, font_size, picTitle)

low_hue = 240/360;
high_hue = 0;
palitra_length = 800;
%font_size = 10;

dim1 = size(map_data, 1);
dim2 = size(map_data, 2);

hsv = zeros(dim1, dim2, 3);
hsv(:,:,1) = low_hue + (high_hue - low_hue)*(map_data - min_value)/(max_value - min_value);
hsv(:,:,2) = 1;
hsv(:,:,3) = orig_pattern*0.8 + 0.1;
%     figure;
%     imshow(OOP_map);
%     figure
%     imshow(hsv2rgb(hsv));

scale_hsv = zeros(palitra_length, 200,3);
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

scale_hsv(max_y-3:max_y,5:43,3) = 0;
scale_hsv(min_y:min_y+3,5:43,3) = 0;
scale_hsv(avg_y-1:avg_y+1,5:43,1) = 0;
scale_hsv(avg_y-1:avg_y+1,5:43,2) = 1;

textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', font_size);
strings = uint8([num2str(max_value) 0 num2str(min_value) 0 num2str(avg_value)]);
labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([45, min_text_y; 45, max_text_y; 45, avg_text_y]));

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

final_image = ones(dim1, dim2 + size(labeled,2),3);
final_image(:,1:dim2,:) = hsv2rgb(hsv);
final_image(1:size(labeled,1),dim2 + 1:dim2 + size(labeled,2),:) = labeled;
final_image(size(final_image,1) - size(scale_bar,1) + 1:size(final_image,1),size(final_image,2) - size(scale_bar,2) + 1:size(final_image,2),:) = scale_bar;
figure;
imshow(final_image);
title(picTitle, 'FontSize', 15);
