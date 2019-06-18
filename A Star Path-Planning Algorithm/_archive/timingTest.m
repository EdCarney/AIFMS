tic

trimmed_MM = MM(26:50,26:35,:);

[rowNum, colNum, widthNum] = size(trimmed_MM);
nodeNum = rowNum * colNum;
edgeNum = rowNum * colNum * 8 * 2;

nodes = zeros(10000000, 3);
edges = zeros(10000000, 3);

% populate node array
for i = 1:rowNum
    for j = 1:colNum
        nodeEntry = (i - 1) * colNum + j;
        nodes(nodeEntry, :) = [nodeEntry, trimmed_MM(i, j, 6), trimmed_MM(i, j, 5)];
    end
end

% populate edge array; assume an 8-connected bidirectional graph
for i = 1:rowNum
    for j = 1:colNum
        for k = 1:8
            edgeEntry = ((i - 1) * colNum + (j - 1)) * 8 + k;

            % determine row offset value
            rowOffset = 0;
            if ismember(k, [1, 2, 3, 9, 10, 11])
                rowOffset = 1;
            elseif ismember(k, [5, 6, 7, 13, 14, 15])
                rowOffset = -1;
            end

            % determine col offset value
            colOffset = 0;
            if ismember(k, [3, 4, 5, 11, 12, 13])
                colOffset = 1;
            elseif ismember(k, [7, 8, 1, 15, 16, 9])
                colOffset = -1;
            end

            % if value of offset is outside of matrix, set row to all
            % zeros; else, populate the edge entry based on the
            % iteration
            if (i + rowOffset < 1 || i + rowOffset > rowNum || j + colOffset < 1 || j + colOffset > colNum)
                edges(edgeEntry, :) = [0, 0, 0];
            else
                fromNode = (i - 1) * colNum + j;
                toNode = ((i + rowOffset) - 1) * colNum + (j + colOffset);
                edges(edgeEntry, :) = [fromNode, toNode, 1];
            end
        end
    end
end

toc