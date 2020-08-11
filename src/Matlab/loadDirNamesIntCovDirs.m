
function [dirNamePolyFiles, dirNameIntCovFiles] = loadDirNamesIntCovDirs(fileNameContainingDirPath)

% chose the right type of slash
if(ispc)
    chosenSlash = '\';
    notChosenSlash = '/';
elseif(isunix)
    chosenSlash = '/';
    notChosenSlash = '\';
else
    display('Error: system si not unix and not PC...')
    pause
end

% load dir names for txt file
fileID = fopen(fileNameContainingDirPath);
formatSpec = '%s';
C = textscan(fileID,formatSpec,...            
                'Delimiter', '\n', ...
                'CollectOutput', true);
fclose(fileID);

% pick only uncommented rows
numRows = size(C{1},1);
count = 1;
for i = 1:numRows
    temp = C{1}{i};
    if(strcmp(temp(1), '%'))
    else
        if(count == 1)
            dirNamePolyFiles = temp(18:end);
            slashInd = strfind(dirNamePolyFiles, notChosenSlash);
            dirNamePolyFiles(slashInd) = chosenSlash;
        elseif(count == 2)
            
            dirNameIntCovFiles = temp(20:end);
            
            slashInd = strfind(dirNameIntCovFiles, notChosenSlash);
            dirNameIntCovFiles(slashInd) = chosenSlash;
        end
        count = count + 1;
    end
end
