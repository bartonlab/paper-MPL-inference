function mutationProbabilities = loadMutProb(fileNameContainingDirPath)


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

for i = 1:numRows
    temp = C{1}{i};
    if(strcmp(temp(1), '%'))
    else
        indOfComma = strfind(temp, ',');
        thisStr = temp(1:indOfComma(1)-1);
        if(strcmp(thisStr, 'muAC'))
            muAC = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muAG'))
            muAG = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muAT'))
            muAT = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muAgap'))
            muAgap = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muCA'))
            muCA = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muCG'))
            muCG = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muCT'))
            muCT = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muCgap'))
            muCgap = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muGA'))
            muGA = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muGC'))
            muGC = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muGT'))
            muGT = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muGgap'))
            muGgap = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muTA'))
            muTA = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muTC'))
            muTC = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muTG'))
            muTG = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'muTgap'))
            muTgap = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));    
        elseif(strcmp(thisStr, 'mugapA'))
            mugapA = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'mugapC'))
            mugapC = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'mugapG'))
            mugapG = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        elseif(strcmp(thisStr, 'mugapT'))
            mugapT = str2double(temp(indOfComma(1)+1:indOfComma(2)-1));
        end
    end
end


mutationProbabilities =    [muAC;
                            muAG;
                            muAT;
                            muAgap;
                            muCA;
                            muCG;
                            muCT;
                            muCgap;
                            muGA;
                            muGC;
                            muGT;
                            muGgap;
                            muTA;
                            muTC;
                            muTG;
                            muTgap;
                            mugapA;
                            mugapC;
                            mugapG;
                            mugapT];