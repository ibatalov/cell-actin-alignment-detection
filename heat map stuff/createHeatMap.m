function heat_map = createHeatMap(cm_alignment_data, shift_vector, PixelSize)

bin_size = 1/PixelSize;
heat_map = zeros(round(250/PixelSize), round(155/PixelSize));
    
    for block_x = 1 : 155/1
        for block_y = 1 : 250/1
            transformed_block_x = 155-block_x - shift_vector(1);
            transformed_block_y = 250-block_y - shift_vector(2);
            
            if transformed_block_x < 1
                transformed_block_x = transformed_block_x + (floor(-transformed_block_x/155) + 1)*155;
            else if transformed_block_x > 155
                    transformed_block_x = transformed_block_x - (floor((transformed_block_x - 155)/155) + 1)*155;
                end
            end
            
            if transformed_block_y <= 0
                transformed_block_y = transformed_block_y + (floor(-transformed_block_y/250) + 1)*250;
            else if transformed_block_y > 250
                    transformed_block_y = transformed_block_y - (floor((transformed_block_y - 250)/250) + 1)*250;
                end
            end
            
            heat_map(max(1,round((block_y - 1)*bin_size)) : round(min(250*bin_size,(block_y)*bin_size)), max(1,round((block_x - 1)*bin_size)) : round(min(155*bin_size,block_x*bin_size))) = cm_alignment_data(transformed_block_x, transformed_block_y);
            %OOP_map(max(1,round((block_y - 1)*1/PixelSize)) : round(min(250/PixelSize,(block_y)*1/PixelSize)), max(1,round((block_x - 1)*1/PixelSize)) : round(min(155/PixelSize,block_x*1/PixelSize))) = cm_alignment_data(block_x, block_y, mapNum);
        end
    end