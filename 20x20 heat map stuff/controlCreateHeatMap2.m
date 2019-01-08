function heat_map = controlCreateHeatMap2(cm_alignment_data, shift, mapSize, mapSizeMicrons)

bin_size = mapSize/40;

micronSizePx = bin_size;
widthMicron = mapSizeMicrons(2);
heightMicron = mapSizeMicrons(1);
widthPx = widthMicron*micronSizePx;
heightPx = heightMicron*micronSizePx;

heat_map = zeros(heightPx, widthPx);

for block = 1 : widthMicron
    
    transformed_block = block - shift;
    transformed_block = mod(transformed_block, 40);
    
    if(transformed_block == 0)
        transformed_block = 40;
    end
    
    col1 = max(1,round((block - 1)*bin_size));
    col2 = min(widthPx,round(block*bin_size));
    
    heat_map(:, col1 : col2) = cm_alignment_data(transformed_block);
    
end
end