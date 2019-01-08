

function origin = overlay_patterns(original_pattern, new_pattern, origin)

dimY = size(original_pattern,1); % 250 micron
dimX = size(original_pattern,2); % 155 micron

if origin(1) < 0
    origin(1) = origin(1) + (floor(-origin(1)/(dimX)) + 1)*dimX;
else if origin(1) > size(new_pattern, 2)
        origin(1) = origin(1) - (floor((origin(1) - size(new_pattern, 2))/dimX) + 1)*dimX;
    end
end

if origin(2) <= 0
    origin(2) = origin(2) + (floor(-origin(2)/dimY) + 1)*dimY;
else if origin(2) > size(new_pattern, 1)
        origin(2) = origin(2) - (floor((origin(2) - size(new_pattern, 1))/dimY) + 1)*dimY;
    end
end


if (size(new_pattern, 1) < round(origin(2)) + round(dimY) - 1) || (size(new_pattern, 2) < round(origin(1)) + round(dimX) - 1)
    new_expanded_pattern = zeros(round(origin(2)) + round(dimY) - 1, round(origin(1)) + round(dimX) - 1);
    new_expanded_pattern(1 : size(new_pattern,1), 1 : size(new_pattern,2)) = new_pattern;
    new_cropped_pattern = new_expanded_pattern(round(origin(2)) : round(origin(2)) + round(dimY) - 1, round(origin(1)) : round(origin(1)) + round(dimX) - 1);
else
    new_cropped_pattern = new_pattern(round(origin(2)) : round(origin(2)) + round(dimY) - 1, round(origin(1)) : round(origin(1)) + round(dimX) - 1);
end

new_cropped_pattern = new_cropped_pattern > 0.1;
new_cropped_pattern = new_cropped_pattern.*0.3 + 0.3;

figure
title('pattern overlay')
fn_overlay = zeros(size(original_pattern,1), size(original_pattern,2), 3);
fn_overlay(:,:,1) = original_pattern;
fn_overlay(:,:,2) = new_cropped_pattern;
imshow(fn_overlay);
end