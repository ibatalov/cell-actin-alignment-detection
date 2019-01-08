function heat_map = controlCreateHeatMap(cm_alignment_data, shift, mapSize)

heat_map = zeros(mapSize, mapSize);
bin_size = mapSize/40;

for block = 1 : size(cm_alignment_data, 1)
    
    transformed_block = block - shift;
    transformed_block = mod(transformed_block, 40);
    
    if(transformed_block == 0)
        transformed_block = 40;
    end
    
    col1 = max(1,round((block - 1)*bin_size));
    col2 = min(mapSize,round(block*bin_size));
    
    heat_map(:, col1 : col2) = cm_alignment_data(transformed_block);
    
end
end