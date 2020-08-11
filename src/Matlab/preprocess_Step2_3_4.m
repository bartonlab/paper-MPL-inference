%%
% Written: 16-Feb, 2020
% Author: M Saqib Sohail
%function [] = markSyn_Mut_v2()
% suggestion for speed + reducing memory usage.  Track codonCompListCell,
% freqOfCodonListCell, countOfCodonsPerSiteAA only for polymorphic sites
% (codons). Right now, the code tracks all sites (codons)


%  code to find the consensus, freq of mut NTs, ref sequence numbering, syn
%  and non syn mutations, and mutation flux vectors.
function preprocess_Step2_3_4(fileNameContainingDirPath, fileNameContainingMutProb, numNT, thisGenomicSegStartInd, thisGenomicSegStopInd, FLAG_binaryApprox, FLAG_SaveFile, FLAG_Epi)


% --------------------------  initializations -----------------------------
mainDir = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system is not unix and not PC...')
    pause
end


%----------------------------------------------------------------------
% 1.1 load dir names
[dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);

FLAG_runPreProcessingStep2 = true;
if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash], 'dir') == 0)
    mkdir([dirNameAnalysis 'Analysis_Misc' chosenSlash])
else
    temp1 = getFolderContent([dirNameAnalysis 'Analysis_Misc' chosenSlash], 'files');
    if(~isempty(temp1))
        FLAG_runPreProcessingStep2 = false;
        disp('Preprocessing step 2 already run on this data. Empty the folder ')
        disp([dirNameAnalysis 'Analysis_Misc' chosenSlash])
        disp('to run preprocessing step 2 on this data again.')
        disp(' ')
    end
    
end

if(FLAG_runPreProcessingStep2 == true)
    disp('Running preprocessing Step 2 - ')
    disp('-----------------------------------------------------------')
    fileNamesThisPat_inAnalysisDir = getFolderContent(dirNameAnalysis, 'files');

    %----------------------------------------------------------------------
    % 1.2 load the .txt file name that contains list of HeaderUpdated fasta files
    if(isempty(fileNamesThisPat_inAnalysisDir))
        disp('Error: Analysis directory is empty. Run preprocess Step 1 before running Step 2')
        pause
    else
        extensionCell{1} = 'txt';
        fileNamesListDataThisPat = findFileNamesWithGivenExtension(dirNameAnalysis, extensionCell);
        numFilesWithTXTExt = length(fileNamesListDataThisPat);
        % find the filename that has '...HU.txt' This is the 
        FLAG_fileNameHAsHUTXT = false(1,numFilesWithTXTExt);
        for e = 1:numFilesWithTXTExt
            thisFile = fileNamesListDataThisPat{e};
            indOfHUTXT = strfind(thisFile, 'HU.txt');
            if(~isempty(indOfHUTXT))
                FLAG_fileNameHAsHUTXT(e) = true;
            end
        end
        if(sum(FLAG_fileNameHAsHUTXT) > 1)
            disp('Error: More than 1 files in the Analysis folder has the string HU.txt in it')
            pause 
        elseif(sum(FLAG_fileNameHAsHUTXT) == 0)
            disp('Error: No .txt file has string HU.txt in it. Run preprocessing step 1')
            pause 
        else
            fileNameContainingHUFiles = fileNamesListDataThisPat{find(FLAG_fileNameHAsHUTXT)};
        end
    end

    fileNameHUFilesCell = loadDataFileNames([dirNameAnalysis fileNameContainingHUFiles]);
    numFileNameHUFilesCell = length(fileNameHUFilesCell);

    %----------------------------------------------------------------------
    % 1.3 get timePointVecFromFileName and bsampleVecFromFileName from the
    %     filenames in fileNameHUFilesCell

    timePointVecFromFileName = -1*ones(1, numFileNameHUFilesCell);
    bsampleVecFromFileName = -1*ones(1, numFileNameHUFilesCell);
    for f = 1:numFileNameHUFilesCell
        thisFileNameHUFiles = fileNameHUFilesCell{f};
        indOfDash = strfind(thisFileNameHUFiles, '_');

        % find time point from file name (time has to be the last _t)
        indOfDashT = strfind(thisFileNameHUFiles, '_t');
        indOfDashT = indOfDashT(end);
        temp1 = (indOfDash - indOfDashT > 0); % (logical)
        allIndOfDashGreatherThanDashT = indOfDash(temp1);
        thisTimePoint = str2double(thisFileNameHUFiles(indOfDashT+2:allIndOfDashGreatherThanDashT(1)-1));
        timePointVecFromFileName(f) = thisTimePoint;

        % find bsample from filename
        indOfbSample = strfind(thisFileNameHUFiles, '_bsample');
        if(isempty(indOfbSample)) % filenames do not contain the string bsample---no bootstrapsampling done. 1 data file per time point
            bsampleVecFromFileName(f) = 1;
        else
            temp2 = (indOfDash - indOfbSample) > 0; % (logical)
            allIndOfDashGreatherThanbSample = indOfDash(temp2);
            str1 = thisFileNameHUFiles(indOfbSample+8:(allIndOfDashGreatherThanbSample(1)-1));
            indOfOF = strfind(str1, 'of');
            bsampleVecFromFileName(f) = str2double(str1(1:indOfOF-1));
            %str1(indOfOF+2:end)
        end

    end

    %----------------------------------------------------------------------
    % 2. 

    % all codon matrix (required for syn/Nsyn)
    mtx = [reshape(repmat([1 2 3 4], 16, 1), 64, 1) reshape(repmat([1 2 3 4], 4,4), 64, 1) repmat([1 2 3 4]', 16, 1)];
    mtx = [mtx; 16 16 16];
    codon2AAList = nt2aa(int2nt(mtx));

    uniqueTimePoints = unique(timePointVecFromFileName);
    numUniqueTimePoints = length(uniqueTimePoints);

    for tp = 1:numUniqueTimePoints
        thisTP = uniqueTimePoints(tp);
        numbSampleThisTP = sum(timePointVecFromFileName == thisTP);
        for bsamp = 1:numbSampleThisTP
            indOfFileToLoad = find((timePointVecFromFileName == thisTP) & (bsampleVecFromFileName == bsamp));
            thisFileNameHUFiles = fileNameHUFilesCell{indOfFileToLoad};

            % load the aligned haplotypes 
            [Header, seqNT_All] = fastaread([dirNameAnalysis thisFileNameHUFiles]);
            % remove 1st seq as that is ref seq (used for alignment purpose only)
            HeaderRefSeq = Header{1};
            refSeq = seqNT_All{1};
            Header = Header(2:end);
            seqNT_All = seqNT_All(2:end);

            %----------------------------------------------------------------------
            % 2.1 - get information from header
            %[thisSeqTimeVec, thisSeqFreqVec, thisSeqNumReadsVec] = getDayFreqNumReedsFromHeader(Header);
            [thisSeqTimeVec, thisSeqFreqVecTemp, thisSeqNumReadsVec] = getInfoFromHeader(Header);

            timePointVec = unique(thisSeqTimeVec);
            numTimePoints = length(timePointVec);

            %----------------------------------------------------------------------
            % 2.2 use numReeds to find frequency of each seq in the MSA
            if(sum(thisSeqNumReadsVec == -1) > 0)
                % numReeds entry not present, use freq entry
                thisSeqFreqVec = thisSeqFreqVecTemp;
            else
                thisSeqFreqVec = [];
                for i = 1:numTimePoints
                    thisTimePoint = timePointVec(i);
                    readsVecThisTimePoint = (thisSeqNumReadsVec(thisSeqTimeVec == thisTimePoint));
                    totalReadsThisTimePoint = sum(readsVecThisTimePoint);
                    seqFreqThisTimePointTemp = readsVecThisTimePoint/totalReadsThisTimePoint;
                    temp1 = seqFreqThisTimePointTemp/sum(seqFreqThisTimePointTemp);
                    temp2 = round(temp1*10000)/10000; % four digitl precision
                    diff = sum(temp2) - 1;

                    % remove diff from largest entry of temp2
                    [val, ind] = max(temp2);
                    temp2(ind) = val - diff;
                    if(abs(sum(temp2) - 1) > 1e-4)
                        disp('Error: check, freq does not sum to 1')
                        disp(num2str(sum(temp2)))
                        pause
                    end
                    seqFreqThisTimePoint = seqFreqThisTimePointTemp;
                    thisSeqFreqVec = [thisSeqFreqVec seqFreqThisTimePoint];
                end
            end


            %----------------------------------------------------------------------
            % 2.3 Get freq of each NT at each site(used from calculating
            %     entropy belo in step 3.1

            msaNT = repmat(' ', size(seqNT_All, 2), length(seqNT_All{1}));

            for k = 1:size(seqNT_All, 2)
               msaNT(k,:) = seqNT_All{k};
            end
            msaNTInt = nt2int(msaNT);
            numSitesNT = size(msaNTInt,2);
            numSitesAA = floor(numSitesNT/3);
            freqCountNT = ntFreqCountCompMSA(msaNTInt, thisSeqFreqVec);

            if(bsamp == 1 && tp == 1)
                freqCountNTAllbSampAllTP = freqCountNT;
            else
                freqCountNTAllbSampAllTP = freqCountNTAllbSampAllTP + freqCountNT;
            end

            % ----- 2.4 outside this loop ----
            %----------------------------------------------------------------------
            % 2.5 - make codonCompListCell from all bsamples for each TP
            % fprintf('Making codonNum MSA...')
            if(tp == 1 && bsamp == 1)
                VAR_MemoryReserveCodonPerSite = 10;
                temp1 = -1*ones(VAR_MemoryReserveCodonPerSite, 3);
                temp2 = -1*ones(VAR_MemoryReserveCodonPerSite, 1);
                codonCompListCell = repmat({temp1}, numUniqueTimePoints, numSitesAA); % lists unique codons at each AA site
                freqOfCodonListCell = repmat({temp2}, numUniqueTimePoints, numSitesAA); % freqs of entries in codonCompListCell
                countOfCodonsPerSiteAA = -1*ones(numUniqueTimePoints, numSitesAA); % counts how many entries in codonCompListCell at each AA site
            end

            % generate codonCompListCell
            for i = 1:numSitesAA
                uniqueCodonsThisSiteInt = unique(msaNTInt(:,(i-1)*3 + [1 2 3]),'rows');
                numUniqueCodonsThisSiteInt = size(uniqueCodonsThisSiteInt, 1);
                thisSite_CodonColumn = msaNTInt(:,(i-1)*3 + [1 2 3]);
                freqOfUniqueCodonsThisSiteInt = 10*ones(numUniqueCodonsThisSiteInt, 1);
                for j = 1:numUniqueCodonsThisSiteInt
                    thisCodon = uniqueCodonsThisSiteInt(j,:);
                    thisCodonNum = codon2num(thisCodon);
                    indToReplace = sum(thisSite_CodonColumn == thisCodon,2) == 3;
                    %msaCodonNum(indToReplace, i) = thisCodonNum;
                    freqOfUniqueCodonsThisSiteInt(j) = sum(thisSeqFreqVec(indToReplace));
                end
                if(countOfCodonsPerSiteAA(tp,i) == -1)
                    codonCompListCell{tp,i}(1:numUniqueCodonsThisSiteInt,:) = uniqueCodonsThisSiteInt;
                    freqOfCodonListCell{tp,i}(1:numUniqueCodonsThisSiteInt) = freqOfUniqueCodonsThisSiteInt;
                    countOfCodonsPerSiteAA(tp,i) = numUniqueCodonsThisSiteInt;
                else
                    % conbine the unique codon from this (tp,bsamp) run with current
                    % entries of codonCompListCell
                    temp3 = [codonCompListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i),:);
                             uniqueCodonsThisSiteInt];
                    temp4 = [freqOfCodonListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i)); freqOfUniqueCodonsThisSiteInt];

                    [newUniqueCodonsThisSiteInt,a2,a3] = unique(temp3,'rows');
                    newNumUniqueCodonsThisSiteInt = length(a2);
                    newFreqOfUniqueCodonsThisSiteInt = zeros(newNumUniqueCodonsThisSiteInt, 1);
                    for jj = 1:newNumUniqueCodonsThisSiteInt
                        newFreqOfUniqueCodonsThisSiteInt(jj) = sum(temp4(a3 == jj));
                    end

                    codonCompListCell{tp,i}(1:newNumUniqueCodonsThisSiteInt,:) = newUniqueCodonsThisSiteInt;
                    freqOfCodonListCell{tp,i}(1:newNumUniqueCodonsThisSiteInt) = newFreqOfUniqueCodonsThisSiteInt;
                    countOfCodonsPerSiteAA(tp,i) = newNumUniqueCodonsThisSiteInt;
                end
            end
        end
        
        %----------------------------------------------------------------------
        % 2.6 this loop normalizes the frequencies to 1 to account for multiple
        % bsamples
        if(numbSampleThisTP > 1)
            for i = 1:numSitesAA
                freqOfCodonListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i)) = freqOfCodonListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i))/numbSampleThisTP;
            end
        end

        %----------------------------------------------------------------------
        % 2.7 mark REF codon 
        if(tp == 1)
           temp1 = [-1 -1 -1];
           msaREFCodonCell = repmat({temp1}, 1, numSitesAA); % initialize an empty cell of codons
        end
        for i = 1:numSitesAA

            % find REF reference for this TP
            if(countOfCodonsPerSiteAA(tp,i) == 1)
                % if only 1 codon is present at this TP, that is the REF
                msaREFCodonCell{tp, i} = codonCompListCell{tp,i}(1,:);
            elseif(tp == 1)
                % if more than 1 codon, but this is the first TP, then the REF
                % is the most frequent Codon
                [a1, indOfMaxFreqEntryTemp] = max(freqOfCodonListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i)));
                msaREFCodonCell{tp, i} = codonCompListCell{tp,i}(indOfMaxFreqEntryTemp,:);
            else

                % if more than 1 codon and this is NOT the 1st TP, then the REF
                % is the same as last time point
                msaREFCodonCell{tp, i} = msaREFCodonCell{tp-1, i};
                
%                 % the ref codon must be in population at this TP. If not,
%                 % ref codon should be the most frequent codon at the
%                 % previous TP provided that it is in population at this TP
%                 % too
%                 repREFCodonTemp = repmat(msaREFCodonCell{tp,i}, countOfCodonsPerSiteAA(tp,i), 1);
%                 codonCompListTemp = codonCompListCell{tp,i}(1:countOfCodonsPerSiteAA(tp,i),:);
% 
%                 % true indicates a match: REF codon exists in MSA at this TP
%                 tempVec = sum(repREFCodonTemp == codonCompListTemp, 2) == 3;
% 
%                 % condition = true when REF codon exists in MSA at this TP
%                 if(sum(tempVec) ~= 0)
%                     % can keep current REF sequence as it is present in the population
%                     % --- do nothng
%                 else
%                     % 
%                     [tp i]
%                     [a1, indOfMaxFreqEntryTemp] = max(freqOfCodonListCell{tp-1,i}(1:countOfCodonsPerSiteAA(tp-1,i)))
%                     codonCompListCell{tp-1,i}(indOfMaxFreqEntryTemp,:)
%                     pause
%                 end
                
            end        
        end

        
        %----------------------------------------------------------------------
        % 2.4 - find consensus sequence
        if(tp == 1)
            [a1, consensusSeqTP1] = max(freqCountNTAllbSampAllTP);
    %         if(sum(consensusSeqTP1 < 1 | consensusSeqTP1 > 4) > 0)
    %             disp('Warning: sequence contains at least 1 site where consensus is non-ACGT')
    %             consensusSeqTP1(find(consensusSeqTP1 < 1 | consensusSeqTP1 > 4))
    %         end
            consensusSeqNTTP1 = int2nt(consensusSeqTP1);
            numSitesAA = floor(length(consensusSeqTP1)/3);
        end
    end

    %------------------------------------------------------------------------
    % 3.1 calculate entropy 

    numNTSitesCons = 0;
    numNTSites2NTs = 0;
    numNTSites3NTs = 0;
    numNTSites4NTs = 0;
    numNTSites5orMoreNTs = 0;
    counter = 0;
    for nt = 1:numSitesNT
        temp1 = freqCountNTAllbSampAllTP(:,nt)/numFileNameHUFilesCell;
        temp2 = temp1(temp1 ~= 0);
        if(length(temp2) == 1)
            numNTSitesCons = numNTSitesCons + 1;
        else
            counter = counter + 1;
            temp2Sort = sort(temp2, 'descend');
            entropyAllMuts(counter) = -temp2Sort'*log(temp2Sort);

            temp2SortApprox = [temp2Sort(1); sum(temp2Sort(2:end))];
            entropyBinApprox(counter) = -temp2SortApprox'*log(temp2SortApprox);

            temp2Sort_length = length(temp2Sort);
            if(temp2Sort_length == 2)
                numNTSites2NTs = numNTSites2NTs + 1;
            elseif(temp2Sort_length == 3)
                numNTSites3NTs = numNTSites3NTs + 1;
            elseif(temp2Sort_length == 4)
                numNTSites4NTs = numNTSites4NTs + 1;
            elseif(temp2Sort_length > 4)
                numNTSites5orMoreNTs = numNTSites5orMoreNTs + 1;
            end
        end
    end

    figure
    upLimLine = max([entropyAllMuts entropyBinApprox]) + 0.1;
    plot(0:0.01:upLimLine, 0:0.01:upLimLine, '--', 'color', [0 0 0], 'LineWidth', 0.5)
    hold on
    plot(entropyAllMuts, entropyBinApprox, 'x')
    xlabel('Entropy based on all mutants at NT site')
    ylabel('Entropy based on binary approx. at NT site')
close all

    numNTSitesPoly = numSitesNT - numNTSitesCons;
    disp(['Number of NT sites: ' num2str(numSitesNT)])
    disp(['Conserved: ' num2str(numNTSitesCons) ', percentage: ' num2str(numNTSitesCons/numSitesNT, '%1.3f')])
    disp(['Polymorphic: ' num2str(numNTSitesPoly) ', percentage: ' num2str(numNTSitesPoly/numSitesNT, '%1.3f')])
    disp(' ')
    disp(['Number of sites with 2 NT: ' num2str(numNTSites2NTs) ', percentage: ' num2str(numNTSites2NTs/numSitesNT, '%1.3f')])
    disp(['Number of sites with 3 NT: ' num2str(numNTSites3NTs) ', percentage: ' num2str(numNTSites3NTs/numSitesNT, '%1.3f')])
    disp(['Number of sites with 4 NT: ' num2str(numNTSites4NTs) ', percentage: ' num2str(numNTSites4NTs/numSitesNT, '%1.3f')])
    disp(['Number of sites with 5 or more NT: ' num2str(numNTSites5orMoreNTs) ', percentage: ' num2str(numNTSites5orMoreNTs/numSitesNT, '%1.3f')])

    %----------------------------------------------------------------------
    % 4.1 Mark Syn/NonSyn
    synOrNonSynBeforeFilt = markSynNonSyn(numNT, numSitesNT, numSitesAA, numUniqueTimePoints, codonCompListCell, countOfCodonsPerSiteAA, msaREFCodonCell);
    siteNTSynNonSynReshaped = reshape(synOrNonSynBeforeFilt, 5, numSitesNT);
    numNTsAtSite = sum(siteNTSynNonSynReshaped ~= -1);
    numNTsAtSite(numNTsAtSite == 0) = 1;
    

    % 4.2 For BINARY APPROXIMATION
    % find if SC at this site is syn/nonSyn in binary approx.
    synNonSyn_BinaryApprox = -1*ones(1, numSitesNT);
    sitesWithMoreMuts = find(numNTsAtSite > 1);
    sitesWithMoreMutsLen = length(sitesWithMoreMuts);
    for i = 1:sitesWithMoreMutsLen
        thisSite = sitesWithMoreMuts(i);
        thisSiteAllNTs = siteNTSynNonSynReshaped(:, thisSite);
        nonWTInd = setdiff([1 2 3 4 16], consensusSeqTP1(thisSite));
        nonWTInd(nonWTInd == 16) = 5; % because intNT value of '-' is 16, while it is the 5th entry at each NTsite
        thisSiteNonWTNTs = thisSiteAllNTs(nonWTInd);
%         if(sum(thisSiteNonWTNTs == 1) > 0)
%             synNonSyn_BinaryApprox(thisSite) = 1;
%         else
%             synNonSyn_BinaryApprox(thisSite) = 0;
%         end

        if(sum(thisSiteNonWTNTs == 2) > 0 || consensusSeqTP1(thisSite) == 16)
            synNonSyn_BinaryApprox(thisSite) = 2;
        elseif(sum(thisSiteNonWTNTs == 1) > 0)
            synNonSyn_BinaryApprox(thisSite) = 1;
        else
            synNonSyn_BinaryApprox(thisSite) = 0;
        end
    end

    %----------------------------------------------------------------------
    % 5.1 make files containing reference indexing, and CSVInfo file containing
    %     Serial #, RefSeq #, RefSeq NT, TP1 cons, mut1, mut2, mut3, mut4, mut5
    
    findIndexingAndMakeCSVFile(dirNameAnalysis, thisFileNameHUFiles, chosenSlash, refSeq, thisGenomicSegStopInd, thisGenomicSegStartInd, freqCountNTAllbSampAllTP, numFileNameHUFilesCell, consensusSeqTP1, numNT, synNonSyn_BinaryApprox);

    
    %----------------------------------------------------------------------
    % 1.3 load mutation probabilities
    fileNameAndPathContainingMutProbPath = [mainDir chosenSlash 'Data_Misc' chosenSlash 'MutationProbabilities' chosenSlash fileNameContainingMutProb];
    mutationProbabilities = loadMutProb(fileNameAndPathContainingMutProbPath);
    muAC = mutationProbabilities(1);
    muAG = mutationProbabilities(2);
    muAT = mutationProbabilities(3);
    muAgap = mutationProbabilities(4);
    muCA = mutationProbabilities(5);
    muCG = mutationProbabilities(6);
    muCT = mutationProbabilities(7);
    muCgap = mutationProbabilities(8);
    muGA = mutationProbabilities(9);
    muGC = mutationProbabilities(10);
    muGT = mutationProbabilities(11);
    muGgap = mutationProbabilities(12);
    muTA = mutationProbabilities(13);
    muTC = mutationProbabilities(14);
    muTG = mutationProbabilities(15);
    muTgap = mutationProbabilities(16);
    mugapA = mutationProbabilities(17);
    mugapC = mutationProbabilities(18);
    mugapG = mutationProbabilities(19);
    mugapT = mutationProbabilities(20);

    muAX = muAC + muAG + muAT + muAgap;
    muCX = muCA + muCG + muCT + muCgap;
    muGX = muGA + muGC + muGT + muGgap;
    muTX = muTA + muTC + muTG + muTgap;
    mugapX = mugapA + mugapC + mugapG + mugapT;

    mu_X = [muAX; muCX; muGX; muTX; mugapX];
    mu_Xmean = [sum([muAC muAG muAT muAgap])/sum([muAC muAG muAT muAgap] ~= 0);
                sum([muCA muCG muCT muCgap])/sum([muCA muCG muCT muCgap] ~= 0);
                sum([muGA muGC muGT muGgap])/sum([muGA muGC muGT muGgap] ~= 0);
                sum([muTA muTC muTG muTgap])/sum([muTA muTC muTG muTgap] ~= 0);
                sum([mugapA mugapC mugapG mugapT])/sum([mugapA mugapC mugapG mugapT] ~= 0)];
    
    muX_ = [muCA + muGA + muTA + mugapA;
            muAC + muGC + muTC + mugapC;
            muAG + muCG + muTG + mugapG;
            muAT + muCT + muGT + mugapT;
            muAgap + muCgap + muGgap + muTgap];
    muX_mean = [sum([muCA muGA muTA mugapA])/sum([muCA muGA muTA mugapA] ~= 0);
            sum([muAC muGC muTC mugapC])/sum([muAC muGC muTC mugapC] ~= 0);
            sum([muAG muCG muTG mugapG])/sum([muAG muCG muTG mugapG] ~= 0);
            sum([muAT muCT muGT mugapT])/sum([muAT muCT muGT mugapT] ~= 0);
            sum([muAgap muCgap muGgap muTgap])/sum([muAgap muCgap muGgap muTgap] ~= 0)];

    NTMuMat = [1-muAX muAC muAG muAT muAgap;
               muCA 1-muCX muCG muCT muCgap;
               muGA muGC 1-muGX muGT muGgap;
               muTA muTC muTG 1-muTX muTgap;
               mugapA mugapC mugapG mugapT 1-mugapX];
    

    %----------------------------------------------------------------------
    % 5.1 - calculate mutation vectors For BINARY APPROXIMATION

    consensusSeqTP1_16 = consensusSeqTP1;
    consensusSeqTP1_16(consensusSeqTP1_16 == 5) = 16;
    mut1Vec = zeros(1, numSitesNT); % int value of first mut NT at all sites
    mut2Vec = zeros(1, numSitesNT);
    mut3Vec = zeros(1, numSitesNT);
    mut4Vec = zeros(1, numSitesNT);
    for i = 1:numSitesNT
        temp1 = freqCountNT(:,i) > 0;
        if(sum(temp1) == 1)% conserved site
        else
            nonZeroInd = freqCountNT(:,i) > 0;
            nonZeroNT = find(nonZeroInd);
            mutNTs = setdiff(nonZeroNT, consensusSeqTP1_16(i));
            mutNTAndFreq = [mutNTs freqCountNT(mutNTs,i)];
            [temp20, temp21] = sort(mutNTAndFreq(:,2), 'descend');
            mutNTAndFreq_Sorted = mutNTAndFreq(temp21, :);

            if(sum(temp1) > 1) % polysite
                mut1Vec(i) = mutNTAndFreq_Sorted(1,1);   
            end
            if(sum(temp1) > 2) % polysite
                mut2Vec(i) = mutNTAndFreq_Sorted(2,1);   
            end
            if(sum(temp1) > 3) % polysite
                mut3Vec(i) = mutNTAndFreq_Sorted(3,1);   
            end
            if(sum(temp1) > 4) % polysite
                mut4Vec(i) = mutNTAndFreq_Sorted(4,1);
            end
        end
    end

    % make site specific mu vector for binary approximation
    
    % 1st code had mu_X(1) instead of mu_Xmean(1), which resulted in over
    % estimate of mutational force and underestimation of selection because
    % of high mutation probability. Now. we use mean mut prob at conserved
    % sites and actual NT-2-MUT probabaility at poly sites.
    
    % for the conserved sites, mu = sum of all three mutations away
    % from the founder allele
    muSiteWiseWT2Mut = zeros(1, numSitesNT);
    muSiteWiseWT2Mut(numNTsAtSite == 1 & consensusSeqTP1 == 1) = mu_Xmean(1);
    muSiteWiseWT2Mut(numNTsAtSite == 1 & consensusSeqTP1 == 2) = mu_Xmean(2);
    muSiteWiseWT2Mut(numNTsAtSite == 1 & consensusSeqTP1 == 3) = mu_Xmean(3);
    muSiteWiseWT2Mut(numNTsAtSite == 1 & consensusSeqTP1 == 4) = mu_Xmean(4);
    muSiteWiseWT2Mut(numNTsAtSite == 1 & (consensusSeqTP1 == 5 | consensusSeqTP1 == 16)) = mu_Xmean(5);
    
    muSiteWiseMut2WT = zeros(1, numSitesNT);
    muSiteWiseMut2WT(numNTsAtSite == 1 & consensusSeqTP1 == 1) = muX_mean(1);
    muSiteWiseMut2WT(numNTsAtSite == 1 & consensusSeqTP1 == 2) = muX_mean(2);
    muSiteWiseMut2WT(numNTsAtSite == 1 & consensusSeqTP1 == 3) = muX_mean(3);
    muSiteWiseMut2WT(numNTsAtSite == 1 & consensusSeqTP1 == 4) = muX_mean(4);
    muSiteWiseMut2WT(numNTsAtSite == 1 & (consensusSeqTP1 == 5 | consensusSeqTP1 == 16)) = muX_mean(5);
    
    % for sites with 2,3,4 NT
    sitesWithMuts = find(numNTsAtSite > 1);
    
    sitesWithMutsLen = length(sitesWithMuts);
    for i = 1:sitesWithMutsLen
        thisSite = sitesWithMuts(i);
        %allPresentNTs = find(siteNTSynNonSynReshaped(:,thisSite) ~= -1)
        allPresentNTs = find(freqCountNTAllbSampAllTP(:,thisSite) ~= 0); % also works for non multiples of 3
        mutNTs = setdiff(allPresentNTs , consensusSeqTP1(thisSite));
        WTNT = consensusSeqTP1(thisSite);
        WTNT(WTNT == 16) = 5;
        mutNTs(mutNTs == 16) = 5;
        muSiteWiseWT2Mut(thisSite) = sum(NTMuMat(WTNT, mutNTs));
        muSiteWiseMut2WT(thisSite) = sum(NTMuMat(mutNTs, WTNT));
    end

    firstInd = repmat(1:numSitesNT, 1, numSitesNT);
    secondInd = reshape(repmat(1:numSitesNT, numSitesNT,1),1, numSitesNT*numSitesNT);
    twoPointSiteIDs = [secondInd(firstInd > secondInd);
             firstInd(firstInd > secondInd)];
    numTwoPointSiteIDs = size(twoPointSiteIDs, 2);
    
    if(FLAG_Epi == true)
        muSiteWiseWT2Mut_i = muSiteWiseWT2Mut(twoPointSiteIDs(1,:));
        muSiteWiseWT2Mut_j = muSiteWiseWT2Mut(twoPointSiteIDs(2,:));
        muSiteWiseMut2WT_iPlusj = muSiteWiseMut2WT(twoPointSiteIDs(1,:)) + muSiteWiseMut2WT(twoPointSiteIDs(2,:));
    end

    thisFileNameHUFiles = fileNameHUFilesCell{1};
    indOfDash = strfind(thisFileNameHUFiles, '_');
    fileNameMutVecs = [thisFileNameHUFiles(1:indOfDash(2)-1) '_mutVecs.txt'];
    fileNameMutVecsEpi = [thisFileNameHUFiles(1:indOfDash(2)-1) '_mutVecsEpi.txt'];
    fileNameSynNonSyn = [thisFileNameHUFiles(1:indOfDash(2)-1) '_SynNonSyn.txt'];
    fileNameCSVData = [thisFileNameHUFiles(1:indOfDash(2)-1) '_info.csv'];
    
    if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecs], 'file') == 2)
        delete([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecs])
    end
    
    if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecsEpi], 'file') == 2)
        delete([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecsEpi])
    end
    
    if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameSynNonSyn], 'file') == 2)
        delete([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameSynNonSyn])
    end
    
    if(FLAG_SaveFile == true)
        fprintf('Saving mutation rates...')
        if(FLAG_binaryApprox == true)
            dlmwrite([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecs], [muSiteWiseWT2Mut; muSiteWiseMut2WT])
            if(FLAG_Epi == true)
                dlmwrite([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecsEpi], [muSiteWiseWT2Mut_j; muSiteWiseWT2Mut_i; muSiteWiseMut2WT_iPlusj])
            end
        
            dlmwrite([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameSynNonSyn], synNonSyn_BinaryApprox)
            disp('done.')
        
        else % if FLAG_binaryApprox == false (Full Potts model)
        
        end
    end
end
