% function that returns a list of files tht have the given text anywhere in
% the file name (but not extension)  present within directory dirName 

% input: dirName               -- complete path of the directory
%        textCell              -- cell of all texts to scan in filenames
%                                 (excluding file extension)

% output: fileNamesListThisDir -- list of all files containing the texts
%                                 specified in the cell textCell present in
%                                 directory dirName

% Written: 02-Feb, 2020
% Author: M Saqib Sohail
function fileNamesListThisDir = findFileNamesWithGivenText(dirName, textCell)


numTexts = length(textCell);

fileNamesThisDir_temp = getFolderContent(dirName, 'files');
numFiles_temp = length(fileNamesThisDir_temp);
fileNamesThisDir_temp_selc = false(1, numFiles_temp);

% select only files with specified texts
for f = 1:numFiles_temp
    thisFileName = fileNamesThisDir_temp{f};
    % check if thisFileName has the specified text
    for e = 1:numTexts
        thisText = textCell{e};
        temp = strfind(thisFileName,'.');
        if(contains(thisFileName(1:temp(end)-1), thisText))
            fileNamesThisDir_temp_selc(f) = true;
        end
    end
end
fileNamesListThisDir = fileNamesThisDir_temp(fileNamesThisDir_temp_selc);

