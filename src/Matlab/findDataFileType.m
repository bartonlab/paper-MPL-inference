% this function takes thisFieName as inpput, identifies the file extension,
% and outputs the correct loadDataOption control flag to load the data file
function loadDataOption = findDataFileType(thisFileName)

thisFileName = lower(thisFileName);
temp = strfind(thisFileName,'.');


if(contains(thisFileName(temp(end)+1:end), 'mat'))  %  MAT file format
    loadDataOption = 1;
elseif(contains(thisFileName(temp(end)+1:end), 'dat')) % DAT file format
    loadDataOption = 2;
elseif(contains(thisFileName(temp(end)+1:end), 'fasta') || contains(thisFileName(temp(end)+1:end), 'fa')) % FASTA file format
    loadDataOption = 3;
elseif(contains(thisFileName(temp(end)+1:end), 'frq')) % FRQ file format
    loadDataOption = 4;
end
