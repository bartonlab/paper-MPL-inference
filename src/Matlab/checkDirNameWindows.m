% window filename can not have a space in directory/file name
% Last updated Oct 2016
% Author: M Saqib Sohail

% ---------------------------------------
function hasSpace = checkDirNameWindows(sourceSAMDirFileName)

correctedSourceSAMDirFileName = sourceSAMDirFileName;
indOfSpace = find(sourceSAMDirFileName == ' ');
if(~isempty(indOfSpace))
    hasSpace = true;
else
    hasSpace = false;
end
