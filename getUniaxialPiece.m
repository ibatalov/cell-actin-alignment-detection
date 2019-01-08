function [ arrayRow,arrayColumn ] = getUniaxialPiece( row, column, orientImage )

arrayRow = row;
arrayColumn = column;
phiZero = orientImage(row,column);
imageSize = length(orientImage);

count = 1;
while count > 0
    currentCount = count;
    currentArrayRow = arrayRow;
    currentArrayColumn = arrayColumn;
    
    for i = 1:currentCount
        if currentArrayRow(end-i+1) > 3 && currentArrayRow(end-i+1) < imageSize-2 && currentArrayColumn(end-i+1) > 3 && currentArrayColumn(end-i+1) < imageSize-2
            
%             linInd = find(abs(cos(orientImage(currentArrayRow(end-i+1)-2:currentArrayRow(end-i+1)+2,currentArrayColumn(end-i+1)-2:currentArrayColumn(end-i+1)+2)-phiZero)) > 0.8);
%             [rawRow, rawColumn] = ind2sub([5,5],linInd);
%             rawRow = rawRow + currentArrayRow(end-i+1)-3;
%             rawColumn = rawColumn + currentArrayColumn(end-i+1)-3;
%             for j = 1:length(rawRow)
%                 TF = (~any(arrayColumn(arrayRow == rawRow(j)) == rawColumn(j)));
%                 rawRow(j) = rawRow(j)*TF;
%                 rawColumn(j) = rawColumn(j)*TF;
%             end
% 
%             if any(rawRow.*rawColumn)
%                 nonzeroRow = rawRow(logical(rawRow));
%                 nonzeroColumn = rawColumn(logical(rawColumn));
%                 arrayRow = [arrayRow; nonzeroRow];
%                 arrayColumn = [arrayColumn; nonzeroColumn];
%                 count = count + length(nonzeroRow);
%             end
            for x = -2:2
                for y = -2:2
                    phi = orientImage(currentArrayRow(end-i+1)+x,currentArrayColumn(end-i+1)+y);
                    if phi ~= 0 && ~(x==0 && y==0)
                        % select only co-oriented pixels that are not yet
                        % selected
                        if abs(cos(phi-phiZero)) > 0.8 && ~any((arrayRow == currentArrayRow(end-i+1)+x).*(arrayColumn == currentArrayColumn(end-i+1)+y))
                            arrayRow = [arrayRow; currentArrayRow(end-i+1)+x];
                            arrayColumn = [arrayColumn; currentArrayColumn(end-i+1)+y];
                            count = count + 1;
                        end
                    end
                    
                end
            end
        end
        count = count - 1;
    end
end

end

