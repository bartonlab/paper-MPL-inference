function findIndexingAndMakeCSVFile(dirNameAnalysis, thisFileNameHUFiles, chosenSlash, refSeq, thisGenomicSegStopInd, thisGenomicSegStartInd, freqCountNTAllbSampAllTP, numFileNameHUFilesCell, consensusSeqTP1, numNT, synNonSyn)

numRowsSynNonSynMtx = size(synNonSyn, 1);
if(numRowsSynNonSynMtx == 1)
    binaryCoding = true;
else
    binaryCoding = false;
end
%------------------------------------------------------------------
% 1.2 Get the ref seq numbering of the 'aligned' HXB2 seq (as it may
%     have inserted gaps due to alignmnet with haplotype sequences)


refSeqAlignedInt = nt2int(refSeq);
refSeqAlignedExtInd = -1*ones(1, length(refSeqAlignedInt));


NTCount = thisGenomicSegStartInd - 1;  

for i = 1:length(refSeqAlignedInt)
    thisRefSeqNTInt = refSeqAlignedInt(i);

    if(thisRefSeqNTInt <= 4 && thisRefSeqNTInt >= 1)
        NTCount = floor(NTCount) + 1;
        refSeqAlignedExtInd(i) = NTCount;
%                count2 = count2 + 1;
%                 HXB2Ind(count2) = i;
    elseif(thisRefSeqNTInt == 16)
        NTCount = NTCount + 0.001;
        refSeqAlignedExtInd(i) = NTCount;
    else
        disp(['NT at position ' num2str(i) ' is ' int2nt(thisRefSeqNTInt) ' -----> must be ACGT or -'])
    end
end

if(floor(NTCount) ~= thisGenomicSegStopInd) % genomeLength is the length of the HXB2 sequence
    disp('Error: NTCount ~= genomeLength, check numbering of alignment in UGENE.')
    NTCount
    thisGenomicSegStopInd
    fileNameRefSeqFasta
    pause
end

protStartInd = find(floor(refSeqAlignedExtInd) == thisGenomicSegStartInd);
protStopInd = find(floor(refSeqAlignedExtInd) == thisGenomicSegStopInd);
indicesOfGenomicSegment = protStartInd(1):protStopInd(end);


haploSeqHXB2NTInt = refSeqAlignedInt(indicesOfGenomicSegment);
haploSeqRefIndexing = refSeqAlignedExtInd(indicesOfGenomicSegment);


haploSeqHXB2NT = int2nt(haploSeqHXB2NTInt);

%------------------------------------------------------------------
% write ref seq numbering (indexing) of this haplotype
indOfDash = strfind(thisFileNameHUFiles, '_');
haploSeqRefIndexingFileName = [thisFileNameHUFiles(1:indOfDash(2)-1) '_RefIndexing.txt'];
if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash haploSeqRefIndexingFileName], 'file') == 2)
    delete([dirNameAnalysis 'Analysis_Misc' chosenSlash haploSeqRefIndexingFileName])
end

dlmwrite([dirNameAnalysis 'Analysis_Misc' chosenSlash haploSeqRefIndexingFileName], haploSeqRefIndexing, 'precision','%.3f');

%-------------------------------------------------------------------------
% save fileNameCSVData containing
% Serial #, RefSeq #, RefSeq NT, TP1 cons, mut1, mut2, mut3, mut4, rev, fwd. mut (rev, fwd. mut not in this code)
% here, we just save the 1st four columns
fileNameCSVData = [thisFileNameHUFiles(1:indOfDash(2)-1)  '_info.csv'];
if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameCSVData], 'file') == 2)
    delete([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameCSVData])
end

% save csv file with
% Serial #, RefSeq #, RefSeq NT, TP1 cons, mut1, mut2, mut3, mut4, rev, fwd. mut (rev, fwd. mut not in this code)
%old% Serial #, HXB2 #, HXB2 NT, clade cons. NT, TP1 cons, mut1, mut2, mut3, mut4, syn/nonSyn, rev, fwd. mut (rev, fwd. mut not in this code)
fprintf('Saving CSVInfo file..')
fileID = fopen([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameCSVData],'w');
if(binaryCoding == true)
    fprintf(fileID,'%s\n', 'Serial no., RefSeq no, RefSeq NT, TP1 cons, NT1, NT2, NT3, NT4, NT5, Syn/NSyn, '); 
else
    fprintf(fileID,'%s\n', 'Serial no., RefSeq no, RefSeq NT, TP1 cons, NT1, NT2, NT3, NT4, NT5, S/NS NT1, S/NS NT2, S/NS NT3, S/NS NT4, S/NS NT5, '); 
end
for nt = 1:length(indicesOfGenomicSegment)

    % list consensus, mut1, mut2, mut3, mut4
    [val, ind] = sort(freqCountNTAllbSampAllTP(:,nt)/numFileNameHUFilesCell, 'descend');
    nonZeroInd_Logical = val ~= 0;
    numNonZeroInd = sum(nonZeroInd_Logical);
    indTempSelected = ind(nonZeroInd_Logical);
    strTemp10 = [];
    for jj = 1:numNonZeroInd
        strTemp10 = [strTemp10 int2nt(indTempSelected(jj)) ', '];
    end    
    strTemp20 = repmat(', ', 1, numNT-numNonZeroInd);
    if(binaryCoding == true)
        strTemp10 = [strTemp10 strTemp20 num2str(synNonSyn(nt)) ', '];
    else
        strTemp10 = [strTemp10 strTemp20];
        for jt = 1:numNT
            strTemp10 = [strTemp10 num2str(synNonSyn(nt,jt)) ', '];
        end
        
    end
    
    if(abs(haploSeqRefIndexing(nt) - floor(haploSeqRefIndexing(nt))) > 0)
        % insertion in the current seq, HXB2 will have  gap
        strTemp = [num2str(nt) ', , ' haploSeqHXB2NT(nt) ', ' int2nt(consensusSeqTP1(nt)) ', ' strTemp10];
    else
        strTemp = [num2str(nt) ', ' num2str(haploSeqRefIndexing(nt)) ', ' haploSeqHXB2NT(nt) ', ' int2nt(consensusSeqTP1(nt)) ', ' strTemp10];
    end
    
    fprintf(fileID,'%s\n',strTemp);


end
fclose(fileID);
disp('done.')