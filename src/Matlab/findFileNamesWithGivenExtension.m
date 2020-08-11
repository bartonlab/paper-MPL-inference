% function that returns a list of files tht have extension given in
% extensionCell present within directory dirName 

% input: dirName               -- complete path of the directory
%        extensionCell         -- cell of all file extensions to select in
%                                 this directory

% output: fileNamesListThisDir -- list of all files with extensions in cell
%                                 extensionCell present in directory
%                                 dirName

% Written: 02-Feb, 2020
% Author: M Saqib Sohail
function fileNamesListThisDir = findFileNamesWithGivenExtension(dirName, extensionCell)


numExtensions = length(extensionCell);

fileNamesThisDir_temp = getFolderContent(dirName, 'files');
numFiles_temp = length(fileNamesThisDir_temp);
fileNamesThisDir_temp_selc = false(1, numFiles_temp);

% select only files with specified extensions
for f = 1:numFiles_temp
    thisFileName = fileNamesThisDir_temp{f};
    % check if thisFileName has the specified extension
    for e = 1:numExtensions
        thisExtenstion = extensionCell{e};
        temp = strfind(thisFileName,'.');
        if(contains(thisFileName(temp(end)+1:end), thisExtenstion))
            fileNamesThisDir_temp_selc(f) = true;
        end
    end
end
fileNamesListThisDir = fileNamesThisDir_temp(fileNamesThisDir_temp_selc);

