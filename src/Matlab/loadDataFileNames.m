
function fileNameToReadCell = loadDataFileNames(fileNameContainingDataFileNames)

% load filenames names for txt file
fileID = fopen(fileNameContainingDataFileNames);
formatSpec = '%s';
C = textscan(fileID,formatSpec,...            
                'Delimiter', '\n', ...
                'CollectOutput', true);
% C = textscan(fileID,formatSpec,...            
%                 'Delimiter', ';', ...
%                 'CollectOutput', true);

fclose(fileID);

% pick only uncommented rows
numRows = size(C{1},1);
%fileNameToReadCell = cell(1, numRows);
count = 1;
for i = 1:numRows
    temp = C{1}{i};
    if(strcmp(temp(1), '%'))
    else
        fileNameToReadCell{count} = temp;
        count = count + 1;
    end
end

