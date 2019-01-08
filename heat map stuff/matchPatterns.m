function newPatternOrigin = matchPatterns(fn_pattern, orig_pattern)
tform = affine2d([-1 0 0; 0 -1 0; 0 0 1]);
fn_pattern1 = imwarp(fn_pattern,tform);
newPatternOrigin = [1 1];

overlay_patterns(orig_pattern, fn_pattern1 > 0.15, newPatternOrigin);

disp('red: first pattern, green: current pattern. If patterns do not match, choose shift vector (first point - new pattern, second - old)');
[inputX, inputY] = getpts;
while numel(inputX) == 2
    close;
    newPatternOrigin(1) = newPatternOrigin(1) + (inputX(2) - inputX(1));
    newPatternOrigin(2) = newPatternOrigin(2) + (inputY(2) - inputY(1));
    newPatternOrigin = overlay_patterns(orig_pattern, fn_pattern1 > 0.15, newPatternOrigin);
    [inputX, inputY] = getpts;
end
close;