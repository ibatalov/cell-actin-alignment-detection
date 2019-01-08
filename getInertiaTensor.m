%% Calculates an inertia tensor of a white object (pixels with value 1) in a binary image
function [tensor, center] = getInertiaTensor(image)

area = nnz(image);
[rows, cols] = find(image);

centerRow = sum(rows)/area; % it is not necessarily integer!
centerCol = sum(cols)/area; % it is not necessarily integer!
center = [centerRow centerCol];

% shifting to the center of inertia
rows1 = rows - centerRow;
cols1 = cols - centerCol;

Irr = sum(cols1.*cols1); % row = y. Irr ~ Iy
Icc = sum(rows1.*rows1); % column = x. Icc = Ix
Irc = -sum(rows1.*cols1); % Ixy = Iyx

tensor = [Irr, Irc; Irc, Icc];