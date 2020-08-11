% input: a cell array called header
% output: a double vector containing allTimePoints of all sequences

% where header conatins 'days: XX, freq: YY, reads: ZZ,' entries
% see Preprocessing_HIV_Inhost_DataFileFormat.m 


function [allTimePointsVec, thisSeqFreqVec, thisSeqNumReedsVec] = getInfoFromHeader(header)

numSeq = length(header);
allTimePointsVec = zeros(1,numSeq);
thisSeqFreqVec = zeros(1,numSeq);
thisSeqNumReedsVec = zeros(1,numSeq);
for i = 1:numSeq
    thisHeader = header{i};
    commaStrIndAll = strfind(thisHeader, ',');
    
    daysStrInd = strfind(thisHeader, 'days: ');
    temp1 = commaStrIndAll(commaStrIndAll > daysStrInd); 
    dayStrEndInd = temp1(1); % the first entry of temp1 is the ',' right after the string 'freq: XX'
    allTimePointsVec(i) = str2double(thisHeader(daysStrInd+6:dayStrEndInd-1));    
    
    
    freqStrInd = strfind(thisHeader, 'freq: ');
    temp2 = commaStrIndAll(commaStrIndAll > freqStrInd); 
    freqStrEndInd = temp2(1); % the first entry of temp1 is the ',' right after the string 'freq: XX'
    thisSeqFreqVec(i) = str2double(thisHeader(freqStrInd+6:freqStrEndInd-1));
    
    readsStrInd = strfind(thisHeader, 'reads: ');
    if(~isempty(readsStrInd))
        readsStrInd = readsStrInd(end);
        temp3 = commaStrIndAll(commaStrIndAll > readsStrInd); 
        readsStrEndInd = temp3(1); % the first entry of temp1 is the ',' right after the string 'freq: XX'
        thisSeqNumReedsVec(i) = str2double(thisHeader(readsStrInd+6:readsStrEndInd-1));
    else
        thisSeqNumReedsVec(i) = -1;
    end    
end