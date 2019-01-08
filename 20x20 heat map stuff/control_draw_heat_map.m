function final_image =  control_draw_heat_map(map_data, orig_pattern, min_value, max_value, avg_value, PixelSize, font_size, picTitle)

low_hue = 240/360;
high_hue = 0;

dim1 = size(map_data, 1);
dim2 = size(map_data, 2);
palitra_length = round(dim1*5/7);
palitra_width = round(dim2/5.5);
%font_size = round(min(dim1, dim2)/30);

hsv = zeros(dim1, dim2, 3);
temp_hue = low_hue + (high_hue - low_hue)*(map_data - min_value)/(max_value - min_value);
temp_hue = temp_hue - min(low_hue, high_hue);
temp_hue = temp_hue.*(temp_hue > 0) + min(low_hue, high_hue);

temp_hue = temp_hue - max(low_hue, high_hue);
temp_hue = temp_hue.*(temp_hue < 0) + max(low_hue, high_hue);

fn_pattern_fill = repmat(orig_pattern, ceil(dim1/size(orig_pattern, 1)), ceil(dim2/size(orig_pattern, 2)));
fn_pattern_fill = fn_pattern_fill(1:dim1, 1:dim2);
hsv(:,:,1) = temp_hue;
hsv(:,:,2) = 1;
hsv(:,:,3) = fn_pattern_fill*0.8 + 0.1;
%     figure;
%     imshow(OOP_map);
%     figure
%     imshow(hsv2rgb(hsv));

scale_hsv = zeros(palitra_length, palitra_width,3);
scale_hsv(:,:,3) = 1;
scale_hsv(:,:,2) = 0;
scale_hsv(:,round(palitra_width/24):round(palitra_width/3),2) = 1;

for row = 1 : palitra_length
    scale_hsv(row,round(palitra_width/24):round(palitra_width/3),1) = high_hue - (high_hue - low_hue)*(row/palitra_length);
end

min_y = 1;
max_y = palitra_length;
avg_y = round(1 + (palitra_length - 1)*(max_value - avg_value)/(max_value - min_value));

min_text_y = 1;
max_text_y = max_y - font_size;
avg_text_y = avg_y - round(font_size/2);

shift = 0;
min_dist = round(font_size*1.5 + 2);
if(abs(avg_y - min_y) < min_dist)
    shift = min_dist - (avg_y - min_y); % this shift is more than 0
else if(abs(avg_y - max_y) < min_dist)
    shift = (max_y - avg_y) - min_dist; % this shift is less than 0
    end
end

avg_text_y = avg_text_y + shift;

scale_hsv(max_y-5:max_y,round(palitra_width/24):round(palitra_width/2.79),3) = 0;
scale_hsv(min_y:min_y+5,round(palitra_width/24):round(palitra_width/2.79),3) = 0;
scale_hsv(avg_y-2:avg_y+2,round(palitra_width/24):round(palitra_width/2.79),1) = 0;
scale_hsv(avg_y-2:avg_y+2,round(palitra_width/24):round(palitra_width/2.79),2) = 1;

textInserter = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', font_size);
strings = uint8([num2str(max_value) 0 num2str(min_value) 0 num2str(avg_value)]);
labeled = step(textInserter, hsv2rgb(scale_hsv), strings, int32([round(palitra_width/2.67), min_text_y; round(palitra_width/2.67), max_text_y; round(palitra_width/2.67), avg_text_y]));

scale_bar_microns = floor((palitra_width - 4)*PixelSize/10)*10;
scale_bar_px = scale_bar_microns/PixelSize;
left_margin = round((palitra_width - scale_bar_px)/2);

scale_bar = ones(min(palitra_width, round(dim1*2/7)),palitra_width,3);
scale_bar(round(size(scale_bar, 1)*0.8):size(scale_bar, 1), left_margin:round(left_margin + scale_bar_px), :) = 0;
string = uint8([num2str(scale_bar_microns) ' um']);
textInserter2 = vision.TextInserter('%s', 'LocationSource', 'Input port', 'Color',  [0, 0, 0], 'FontSize', round(font_size));
scale_bar = step(textInserter2, scale_bar, string, int32([round(left_margin*1.1) round(size(scale_bar, 1)*0.8 - font_size*1.15)]));
%     figure;
%     imshow(scale_bar);

final_image = ones(dim1, dim2 + size(labeled,2),3);
final_image(:,1:dim2,:) = hsv2rgb(hsv);
final_image(1:size(labeled,1),dim2 + 1:dim2 + size(labeled,2),:) = labeled;
final_image(size(final_image,1) - size(scale_bar,1) + 1:size(final_image,1),size(final_image,2) - size(scale_bar,2) + 1:size(final_image,2),:) = scale_bar;
figure;
imshow(final_image);
title(picTitle, 'FontSize', 15);
