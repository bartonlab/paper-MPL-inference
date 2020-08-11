% this code :
%      1. returns the names of all subfolders/files in the folder 'dirName'

% Last updated 28-Jan 2020
% Author: M Saqib Sohail

function selectedContent = getFolderContent(dirName, tag)


% List of all files and folders in dirName
files = dir(dirName);

% Logical vector with 1 indicating it is a directory
if(strcmp(tag, 'dir'))
    flags = [files.isdir];
elseif(strcmp(tag, 'files'))
    flags = [files.isdir] == 0;
else
    disp(['Tag must be either ''dir'' or ''files'' '])
    return
end

% Extract ONLY directories
selectedContent_temp = files(flags);

selectedContent = {};
i = 0;
for k = 1 : length(selectedContent_temp)
    if(strcmp(selectedContent_temp(k).name, '.') || strcmp(selectedContent_temp(k).name, '..'))
    else
        i = i + 1;
        selectedContent{i} = selectedContent_temp(k).name;
    end
end

return