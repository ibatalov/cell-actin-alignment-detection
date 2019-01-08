% to use this script, you need to load alignment maps for both low and high
% density cases named cm_alignment_data_low and cm_alignment_data_high. The
% same for the fn_pattern_low, fn_pattern_high
if(~exist('rootfolder', 'var'))
    rootfolder = '/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff';
end

[file1, path1] = uigetfile({'*.mat', '*.mat select low density data'}, 'select low density data', rootfolder);
load([path1 file1], 'fn_pattern','cm_alignment_data');
fn_pattern_low = fn_pattern;
cm_alignment_data_low = cm_alignment_data;

PixelSize = 40/size(fn_pattern_low, 2);
bin_size = 1/PixelSize;

[file2, path2] = uigetfile({'*.mat', '*.mat select high density data'}, 'select high density data', path1);
load([path2 file2], 'fn_pattern','cm_alignment_data');
fn_pattern_high = fn_pattern;
cm_alignment_data_high = cm_alignment_data;

cm_alignment_data_low(:,1) = cm_alignment_data_low(:,1)/max(cm_alignment_data_low(:,1))*100;
cm_alignment_data_high(:,1) = cm_alignment_data_high(:,1)/max(cm_alignment_data_high(:,1))*100;

average_values_low = [0; 0; 0; 0; 0; 0]; % average per location
average_values_high = [0; 0; 0; 0; 0; 0]; % average per location

for k = 1 : 6
    nnz_values_low = cm_alignment_data_low(:,k);
    nnz_values_low = nnz_values_low(nnz_values_low ~= 0);
    nnz_values_high = cm_alignment_data_high(:,k);
    nnz_values_high = nnz_values_high(nnz_values_high ~= 0);
    average_values_low(k) = sum(nnz_values_low(:))/numel(nnz_values_low);
    average_values_high(k) = sum(nnz_values_high(:))/numel(nnz_values_high);
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

[file, path] = uigetfile({'*.*', '*.* select original pattern image'}, 'select low density data', rootfolder);
orig_pattern = imread([path file]);
orig_pattern = rgb2gray(orig_pattern);
orig_pattern = double(orig_pattern);
orig_pattern = orig_pattern/max(orig_pattern(:));

correlation = [0 0 0 0 0 0];

final_images = cell(2,6);

for mapNum = 1 : 6
    
    map_low = controlCreateHeatMap2(cm_alignment_data_low(:,mapNum), 11, size(orig_pattern, 1), [80, 155]);
    map_high = controlCreateHeatMap2(cm_alignment_data_high(:,mapNum), 10, size(orig_pattern, 1), [80, 155]);
    
    correlation(mapNum) = corr2(cm_alignment_data_low(:,mapNum), cm_alignment_data_high(:,mapNum));
    
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
        cut_out_fraction = 0.2; % fraction of data that you want to make outside the color range (over or undersaturated) to increase contrast of what is near the average
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
    final_images{1,mapNum} = control_draw_heat_map(map_low, orig_pattern, min_value, max_value, avg_value_low, 40/size(orig_pattern, 1), 60, strcat(data_label(mapNum), ', Low Cell Density'));
    final_images{2,mapNum} = control_draw_heat_map(map_high, orig_pattern, min_value, max_value, avg_value_high, 40/size(orig_pattern, 1), 60, strcat(data_label(mapNum), ', High Cell Density'));
    
    %% Show correlations:
%    disp(strcat(data_label(mapNum),' correlation:_', num2str(correlation(mapNum))));
 
    %% Draw histograms
    bin_width = (max_value - min_value)/50;
    
    current_plot = figure('units','normalized','position',[0.1 0.1 0.43 0.6])
    hist_low = histogram(cm_alignment_data_low(:,mapNum), 10, 'BinWidth', bin_width);
%    title(strcat(data_label(mapNum), ' Distribution'), 'FontSize', 30);
    xlim([min_value max_value]);
    xlabel(x_axis_label(mapNum), 'FontSize', 30, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontSize', 30, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    set(gca,'FontSize',25, 'FontName', 'Times New Roman');
    
    hold on
    hist_high = histogram(cm_alignment_data_high(:,mapNum), 10, 'BinWidth', bin_width);

    hist_low.EdgeColor = 'none';
    hist_high.EdgeColor = 'none';
    legend('\fontsize{30}\fontname{Times New Roman} low cell density', '\fontsize{30}\fontname{Times New Roman} high cell density', 'Location', 'southoutside');
    
    if(~exist('save_dir', 'var'))
        save_dir = uigetdir(rootfolder, 'select low density data');
    end
    saveas(current_plot, strcat('/Volumes/Macintosh HD 2/Drop-boxes/Dropbox (RBG)/Lab stuff/2015/2015-05-07 Matlab heat map stuff/vector graphics histograms/', data_label{mapNum}, '.pdf'));
    
end
%% save images as png files
folder_name = uigetdir
if(folder_name ~= 0)
    for i = 1 : 6
        imwrite(final_images{1,i}, strcat(folder_name, '/', data_label{i}, '_low_density.png'));
        imwrite(final_images{2,i}, strcat(folder_name, '/', data_label{i}, '_high_density.png'));
    end
end