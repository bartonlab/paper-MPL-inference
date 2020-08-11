
clc
clear all
close all

normalized = false;

CategoriesStr = 'CatD';
distBasedCategories = [   1  100;
                        101  700;
                        701 1500;
                       1501 6000];
                   
mainDir = pwd;
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system is not unix and not PC...')
    pause
end

histBinDataDir = [mainDir chosenSlash CategoriesStr chosenSlash];
if(exist(histBinDataDir, 'dir') == 0)
    mkdir(histBinDataDir)
end


fileNameContainingDirPath = 'DataPathFile.txt';
[dirNamePatInfoFiles, dirNameIntCovFiles] = loadDirNamesIntCovDirs(fileNameContainingDirPath);


textCell{1} = 'poly';
fileNamesListThisDir = findFileNamesWithGivenText(dirNamePatInfoFiles, textCell);
numPat = length(fileNamesListThisDir);


numPat
for pat = 1:numPat
    pat
    thisPatInfoFileName = fileNamesListThisDir{pat};
    
    indOfDot = strfind(thisPatInfoFileName, '.');
    strTemp1 = thisPatInfoFileName(1:(indOfDot(1)-1));
    
    indOfDash = strfind(thisPatInfoFileName, '-');
    
    covFileName = ['covariance-' strTemp1 '-seq2state.dat'];
    intCovMtx = dlmread([dirNameIntCovFiles covFileName]);
        
    patInfoTable = readtable([dirNamePatInfoFiles thisPatInfoFileName]);
    [numMutsThisPat, temp20] = size(patInfoTable);
    
    indMtx = table2array(patInfoTable(:,[1 2]));
    indMtx = indMtx + 1;

    ithSiteRemovedVec = indMtx(:,2);
    polySiteIndVec = indMtx(:,1);
    
    distBasedCategories(end,2) = ithSiteRemovedVec(end);
    protLen = ithSiteRemovedVec(end);
    
    
    numCategories = size(distBasedCategories,1);

    approxCatSize = round((numMutsThisPat*numMutsThisPat*numPat)/2);

    % initialize
    if(pat == 1)
       allSortedIntCovMatxValuesCatCell = repmat({-1*ones(1, approxCatSize)}, 1, numCategories);
       catCounter = zeros(1, numCategories);
    end

    for nmut = 1:numMutsThisPat
        ithSiteRemoved = ithSiteRemovedVec(nmut);
        %------------------------------------------------------
        % find indices of the protien that lie within the
        % distance of each category
        tempInd1 = ithSiteRemoved - flip(distBasedCategories')';
        validCategoriesRanges1 = sum(tempInd1 < 1, 2) <= 1;
        tempInd1(tempInd1 < 1) = 1;

        tempInd2 = ithSiteRemoved + distBasedCategories;
        validCategoriesRanges2 = sum(tempInd2 > protLen, 2) <= 1;
        tempInd2(tempInd2 > protLen) = protLen;
        
        validCategoriesRanges = cell(1, numCategories);
        for cat = 1:numCategories
            if(validCategoriesRanges1(cat))
                currRangeTemp = validCategoriesRanges{cat};
                selAbsSiteInd = polySiteIndVec((ithSiteRemovedVec >= tempInd1(cat,1)) & (ithSiteRemovedVec <= tempInd1(cat,2)));
                validCategoriesRanges{cat} = [currRangeTemp selAbsSiteInd'];
            end
            if(validCategoriesRanges2(cat))
                currRangeTemp = validCategoriesRanges{cat};
                selAbsSiteInd = polySiteIndVec((ithSiteRemovedVec >= tempInd2(cat,1)) & (ithSiteRemovedVec <= tempInd2(cat,2)));
                validCategoriesRanges{cat} = [currRangeTemp selAbsSiteInd'];
            end
        end
        
        for cat = 1:numCategories
            thisCatCounter = catCounter(cat);
            numNewEntries = length(validCategoriesRanges{cat});
            allSortedIntCovMatxValuesCatCell{cat}((thisCatCounter+1):(thisCatCounter+numNewEntries)) = intCovMtx(polySiteIndVec(nmut), validCategoriesRanges{cat});
            catCounter(cat) = thisCatCounter + numNewEntries;
        end     
    end
end

for cat = 1:numCategories
   binDataFileName = ['binData_' num2str(cat) '.txt'];
   dlmwrite([histBinDataDir binDataFileName], allSortedIntCovMatxValuesCatCell{cat}(1:catCounter(cat)));
end

