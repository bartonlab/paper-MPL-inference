
function [dirNameData, dirNameAnalysis, dirNameFigures, dirNameABC, dirNameApproxWF] = loadDirNames(fileNameContainingDirPath)

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
            dirNameData = temp(13:end);
            slashInd = strfind(dirNameData, notChosenSlash);
            dirNameData(slashInd) = chosenSlash;
        elseif(count == 2)
            if(length(temp) > 17)
                dirNameAnalysis = temp(17:end);
            else
                dirNameAnalysis = [dirNameData(1:length(dirNameData)-5) '/Analysis/'];
            end
            slashInd = strfind(dirNameAnalysis, notChosenSlash);
            dirNameAnalysis(slashInd) = chosenSlash;
        elseif(count == 3)
            dirNameFigures = [dirNameData(1:length(dirNameData)-6) '/Figures/'];
            slashInd = strfind(dirNameFigures, notChosenSlash);
            dirNameFigures(slashInd) = chosenSlash;
        elseif(count == 4)
            if(length(temp) > 12)
                dirNameABC = temp(12:end);
            else
                dirNameABC = '';
            end
            slashInd = strfind(dirNameABC, notChosenSlash);
            dirNameABC(slashInd) = chosenSlash;
        elseif(count == 5)
            if(length(temp) > 17)
                dirNameApproxWF = temp(17:end);
            else
                dirNameApproxWF = '';
            end
            slashInd = strfind(dirNameApproxWF, notChosenSlash);
             dirNameApproxWF(slashInd) = chosenSlash;
        elseif(count == 6)
        end
        count = count + 1;
    end
end


%
% function [dirName, dirShareName] = loadDirNames(fileNameContainingDirPath)
% fileID = fopen(fileNameContainingDirPath);
% formatSpec = '%s';
% C = textscan(fileID,formatSpec,...            
%                 'Delimiter', '\n', ...
%                 'CollectOutput', true);
% fclose(fileID);
% 
% str1000 = C{1}{1};
% ind_str1000 = strfind(str1000, 'dirName=');
% dirName = str1000(9:end);
% 
% str1001 = C{1}{2};
% ind_str1001 = strfind(str1001, 'dirShareName=');
% if(length(str1001) > 13)
%     dirShareName = str1001(14:end);
% else
%     dirShareName = dirName;
% end
