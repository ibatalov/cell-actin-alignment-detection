% to use this script, you need to load alignment maps for both low and high
% density cases named cm_alignment_data_low and cm_alignment_data_high. The
% same for the fn_pattern_low, fn_pattern_high

%load('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/first results (low density)/combined data/contour pattern/data.mat', 'fn_pattern','cm_alignment_data');
load('E:\Dropbox (RBG)\Dropbox (RBG)\Lab stuff\2015\2015-05-07 Matlab heat map stuff\first results (low density)\combined data\contour pattern\data.mat', 'fn_pattern','cm_alignment_data');
fn_pattern_low = fn_pattern;
cm_alignment_data_low = cm_alignment_data;

PixelSize = 155/size(fn_pattern_low, 2);
bin_size = 1/PixelSize;

%load('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/high density/results_high_density.mat', 'fn_pattern','cm_alignment_data');
load('E:\Dropbox (RBG)\Dropbox (RBG)\Lab stuff\2015\2015-05-07 Matlab heat map stuff\high density\results_high_density.mat', 'fn_pattern','cm_alignment_data');
fn_pattern_high = fn_pattern;
cm_alignment_data_high = cm_alignment_data;
fn_pattern_high = imresize(fn_pattern_high, [round(250/PixelSize), round(155/PixelSize)], 'Method', 'bilinear');

cm_alignment_data_low(:,:,1) = cm_alignment_data_low(:,:,1)/max(max(cm_alignment_data_low(:,:,1)))*100;
cm_alignment_data_high(:,:,1) = cm_alignment_data_high(:,:,1)/max(max(cm_alignment_data_high(:,:,1)))*100;

average_values_low = [0; 0; 0; 0; 0; 0]; % average for the pattern location
average_values_high = [0; 0; 0; 0; 0; 0]; % average for the pattern location

for k = 1 : 6
    average_values_low(k) = sum(sum(cm_alignment_data_low(:,:,k)))/nnz(cm_alignment_data_low(:,:,k));
    average_values_high(k) = sum(sum(cm_alignment_data_high(:,:,k)))/nnz(cm_alignment_data_high(:,:,k));
end

data_label = {'Normalized Cell Occurrence';
    'Mean Orientation Angle';
    'Standard Deviation of the Mean Angle';
    'Median Orientation Angle';
    'Most Probable Orienation Angle';
    'OOP'};
x_axis_label = { 'Normalized Cell Occurrence';
    'Mean Orientation Angle, Degrees'
    'Std. Dev. of Mean Angle, Degrees'
    'Median Orientation Angle, Degrees'
    'Most Probable Orienation Angle, Degrees';
    'OOP'};

%orig_pattern = imread('D:\Ivan\2015\2015-05-07 Matlab heat map stuff\original pattern\BM_pattern_for_heat_map_contour_1.png');
%orig_pattern = imread('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/original pattern/BM_pattern_for_heat_map_contour_4.png');
orig_pattern = imread('E:\Dropbox (RBG)\Dropbox (RBG)\Lab stuff\2015\2015-05-07 Matlab heat map stuff\original pattern\BM_pattern_for_heat_map_contour_4.png');
orig_pattern = rgb2gray(orig_pattern);
orig_pattern = double(orig_pattern);
orig_pattern = orig_pattern/max(orig_pattern(:));
%orig_pattern = orig_pattern > 0.99;
orig_pattern = imresize(orig_pattern, [round(250/PixelSize), round(155/PixelSize)], 'Method', 'bilinear');
%imshow(orig_pattern)

newParrenOrigin_low = matchPatterns(fn_pattern_low, orig_pattern);
newParrenOrigin_high = matchPatterns(fn_pattern_high,orig_pattern);

shift_vector_low = round(newParrenOrigin_low/bin_size);
shift_vector_high = round(newParrenOrigin_high/bin_size);

correlation = [0 0 0 0 0 0];

final_images = cell(2,6);

for mapNum = 1 : 6
    
    map_low = createHeatMap(cm_alignment_data_low(:,:,mapNum), shift_vector_low, PixelSize);
    map_high = createHeatMap(cm_alignment_data_high(:,:,mapNum), shift_vector_high, PixelSize);
    
    correlation(mapNum) = corr2(map_high, map_low);
    
    max_value = round(max([map_low(:); map_high(:)]),3, 'significant');
    min_value = round(min([map_low(:); map_high(:)]),3, 'significant');
    
    max_value = round(max_value,3, 'decimals');
    min_value = round(min_value,3, 'decimals');
    
    avg_value_low = round(average_values_low(mapNum),3, 'significant');
    avg_value_high = round(average_values_high(mapNum),3, 'significant');
    
    avg_value_low = round(avg_value_low,3, 'decimals');
    avg_value_high = round(avg_value_high,3, 'decimals');
    
    %% under- and oversaturate maps for higher contrast (optional)
    if(mapNum == 2 || mapNum == 3 || mapNum == 4 || mapNum == 5)
        cut_out_fraction = 0.1; % fraction of data that you want to make outside the color range (over or undersaturated) to increase contrast of what is near the average
        avg_value = (avg_value_low + avg_value_high)/2;
        
        sorted_data_low = sort(map_low(:));
        sorted_data_high = sort(map_high(:));
        data_size = length(sorted_data_low);
        
        min_value = min(sorted_data_low(round(data_size*cut_out_fraction/2)), sorted_data_high(round(data_size*cut_out_fraction/2)));
        max_value = max(sorted_data_low(round(data_size*(1-cut_out_fraction/2))), sorted_data_high(round(data_size*(1-cut_out_fraction/2))));
        
        max_value = round(max_value,3, 'significant');
        min_value = round(min_value,3, 'significant');
        
        max_value = round(max_value,3, 'decimals');
        min_value = round(min_value,3, 'decimals');
        
        %half_range = min(avg_value - min_value, max_value - avg_value)*(1 - cut_out_fraction);
        %min_thresh = avg_value - half_range;
        %max_thresh = avg_value + half_range;
        map_low = (map_low - min_value).*(map_low > min_value) + min_value; %undersaturate data
        map_low = (map_low - max_value).*(map_low < max_value) + max_value; %oversaturate data
        
        map_high = (map_high - min_value).*(map_high > min_value) + min_value; %undersaturate data
        map_high = (map_high - max_value).*(map_high < max_value) + max_value; %oversaturate data
    end
    %% show final images
   final_images{1,mapNum} = draw_heat_map(map_low, orig_pattern, min_value, max_value, avg_value_low, PixelSize, 94, strcat(data_label(mapNum), ', Low Cell Density'));
   final_images{2,mapNum} = draw_heat_map(map_high, orig_pattern, min_value, max_value, avg_value_high, PixelSize, 94, strcat(data_label(mapNum), ', High Cell Density'));
    
    %% Show correlations:
    %disp(strcat(data_label(mapNum),' correlation:_', num2str(correlation(mapNum))));
 
    %% Draw histograms
    bin_width = (max_value - min_value)/250;
    
    current_plot = figure('units','normalized','position',[0.1 0.1 0.43 0.6])
    hist_low = histogram(map_low(:), 100, 'BinWidth', bin_width);
%    title(strcat(data_label(mapNum), ' Distribution'), 'FontSize', 30);
    xlim([min_value max_value]);
    xlabel(x_axis_label(mapNum), 'FontSize', 30, 'FontName', 'Arial', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontSize', 30, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca,'FontSize',25, 'FontName', 'Arial');
    
    hold on
    hist_high = histogram(map_high(:), 100, 'BinWidth', bin_width);

    hist_low.EdgeColor = 'none';
    hist_high.EdgeColor = 'none';
    legend('\fontsize{30}\fontname{Arial} low cell density', '\fontsize{30}\fontname{Arial} high cell density', 'Location', 'southoutside');

    %saveas(current_plot, strcat('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/vector graphics histograms/', data_label{mapNum}, '.pdf'));
    saveas(current_plot, strcat('E:\Dropbox (RBG)\Dropbox (RBG)\Lab stuff\2015\2015-05-07 Matlab heat map stuff\vector graphics histograms\new\', data_label{mapNum}, '.pdf'));
    
end
%% save images as png files
folder_name = uigetdir
if(folder_name ~= 0)
    for i = 1 : 6
        imwrite(final_images{1,i}, strcat(folder_name, '/', data_label{i}, '_low_density.png'));
        imwrite(final_images{2,i}, strcat(folder_name, '/', data_label{i}, '_high_density.png'));
    end
end