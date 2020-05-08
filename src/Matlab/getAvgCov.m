function [avgLD, allPairInVec] = getAvgCov(covMtxThisTime)
numSites = size(covMtxThisTime, 1);
allPairInVec = -1*ones(1, (numSites*(numSites-1)/2));
count = 0;
for i = 1:(numSites-1)
    count = count + 1;
    numEntriesThisI = (numSites - i);
    allPairInVec(count:(count+numEntriesThisI-1)) = covMtxThisTime(i,i+1:end);
    count = count + numEntriesThisI - 1;
end
avgLD = mean(allPairInVec);
