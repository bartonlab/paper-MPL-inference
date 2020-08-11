% preprocess UNSW HCV haplotype data for MPL analysis

% input is reconstructed haplotypes (From QuasiRecomb based pipe line of
% Umer), aligned to reference sequence and manually checked for codon
% correct aligment. All time point sequences are already aligned to each
% other

% 1. Get timepoint information from filename, order files w.r.t. time 
% 2. Make header compatible with MPL, rewigth frequencies
% 3. Generate new .txt file conatianing names of header updated
%    fasta files with ref seq to be used by the alignment function

% Written: 14-Feb, 2020
% Author: M Saqib Sohail

% update: now does not combine all tp files into a single file.
%         Also combines step 4 (reweighting of freqs so freq of haplo at
%         each time poit sum to 1)

function [] = preprocess_Step1_v2(fileNameContainingDirPath, fileNameFastaFilesWithHU, FLAG_SaveFile, FLAG_firstSeqIsRef, FLAG_useFreqEntry)


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

[dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
if(exist(dirNameAnalysis, 'dir') == 0)
    mkdir(dirNameAnalysis)        
end


disp('Running preprocessing Step 1 - preparing data files for MPL')
disp('-----------------------------------------------------------')
fileNamesThisPat_inAnalysisDir = getFolderContent(dirNameAnalysis, 'files');


% % if analysis folder is empty, preprocess data .fa/.fasta files
%if(isempty(fileNamesThisPat_inAnalysisDir))
% run this ONLY if analysis folder does not contain fastaFilesHU.txt file
if(~(exist([dirNameAnalysis fileNameFastaFilesWithHU], 'file') == 2))

    
    %extensionCell{1} = 'fasta';
    %extensionCell{2} = 'fa';
    %fileNamesListDataThisPat = findFileNamesWithGivenExtension(dirNameData, extensionCell)
    
    fileNamesListDataThisPat = findFileNamesWithGivenText(dirNameData, {'bsample1of1_t'});
    
    %fileNamesListThisDir = findFileNamesWithGivenText(mainDir, textCell);
    %----------------------------------------------------------------------
    % specific to Umer'sNGS pipeline data

    

    %------------------------------------------------------------------
    % 1.1 Get timepoint information from filename
    %     order files w.r.t. time 
    %length(fileNamesListDataThisPat)
    numFiles = length(fileNamesListDataThisPat);
    timeVecTemp = zeros(1, numFiles);
    for f = 1:numFiles
        thisFileNameThisPat = fileNamesListDataThisPat{f};
        indOfDash = strfind(thisFileNameThisPat, '_');
        indOfDot = strfind(thisFileNameThisPat, '.');
        indOfTime = strfind(thisFileNameThisPat, '_t');
        
        if(isempty(indOfTime))
            disp('Error: Can not extract time information from filename.')
            disp('       Filename does not contain the string DPI...')
            pause
            break
        end
        if(length(indOfTime) == 1)
             % '_t' occurs just once in the file name, it is certainly the timepoint string
        else
            % however, if there are more than 1 '_t' strings, handle it
            disp('Error: Multiple _t strings in file name. CASE NOTHANDLED YET')
            pause
            break
        end

        timeVecTemp(f) = str2double(thisFileNameThisPat(indOfTime(end)+2:indOfDot(end)-1));
    end

    [timeVec, timeVecSortInd] = sort(timeVecTemp);
    fileNamesListDataThisPat = fileNamesListDataThisPat(timeVecSortInd);

    %------------------------------------------------------------------
    % 2. Make header compatible with MPL, rewigth frequencies
    
    
    
    for f = 1:numFiles
        disp(['Processing file ' num2str(f) ' of ' num2str(numFiles)])
        thisFileNameThisPat = fileNamesListDataThisPat{f};
        indOfDash = strfind(thisFileNameThisPat, '_');
        indOfDot = strfind(thisFileNameThisPat, '.');
        
        % load data
        [Header,seqNT_All] = fastaread([dirNameData thisFileNameThisPat]);

        if(iscell(Header))
        else
            HeaderTemp = Header;
            Header = cell(1,1);
            Header{1} = HeaderTemp;
            seqNT_AllTemp = seqNT_All;
            seqNT_All = cell(1,1);
            seqNT_All{1} = seqNT_AllTemp;
        end

        % 2.1 check if Header is compatible with MPL or not 
        indOfDays = strfind(Header{end}, 'days:');
        if(isempty(indOfDays))
            FLAG_updateHeader = true;
        else
            FLAG_updateHeader = false;
            HeaderNew = Header;
            seqNT_AllNew = seqNT_All;
        end

        
        % 2.2 Make header compatible with MPL, i.e., 
        %     Header is in format >........, days: X, freq: Y, reads: Z,
        if(FLAG_updateHeader)
            initialNumEntries = size(seqNT_All,2);
            temp2 = repmat(' ', 1, 150);
            temp3 = repmat(' ', 1, length(seqNT_All{1}));
            HeaderNew = repmat({temp2},1,initialNumEntries);
            seqNT_AllNew = repmat({temp3},1,initialNumEntries); 


            numHaploThisTP = length(Header);
            for h = 1:numHaploThisTP
                thisHeader = Header{h};
                indOfDash = strfind(thisHeader, '_');
                freqThisHaplo = str2double(thisHeader(indOfDash+1:end));


                HeaderNew{h} = [thisHeader ', days: ' num2str(timeVec(f)) ', freq: ' num2str(freqThisHaplo) ', reads: ' num2str(-1) ','];
                if(FLAG_firstSeqIsRef && h == 1)
                    HeaderNew{h} = thisHeader;
                end
                seqNT_AllNew{h} = seqNT_All{h};
            end
        end
        
        % -----------------------------------------------------------------
        % 2.3 Reweight frequncies if needed, 
        
        % load frequency info from Header
        [thisSeqTimeVec, thisSeqFreqVecTemp, thisSeqNumReadsVec] = getInfoFromHeader(HeaderNew(2:end));
        
        timePointVec = unique(thisSeqTimeVec);
        numTimePoints = length(timePointVec);

        % use numReeds to find frequency of each seq in the MSA
        %if(sum(thisSeqNumReadsVec == -1) > 0)
        if(FLAG_useFreqEntry == true)
            % numReeds entry not present, use freq entry
            thisSeqFreqVec = thisSeqFreqVecTemp;
            disp('Using freq: entry of header to find frequency...')
        else
            disp('Using reads entry of header to calculate frequency...')
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

        
        % 2.3.1 Check if reweighting of frequencies is needed
        FLAG_normalize = false;
        thisSeqFreqVec_reweighted = zeros(1, length(thisSeqFreqVec));
        for t = 1:numTimePoints
            freqVecThisTP = thisSeqFreqVec(thisSeqTimeVec == timePointVec(t));
            temp1 = sum(freqVecThisTP);
            disp(['Time point ' num2str(timePointVec(t)) ', sum of frequenices of all haplotypes: ' num2str(temp1)])

            % perform normalization if sum of freq at ths TP is less than 0.99

            if(abs(temp1-1) > 0.01)
                FLAG_normalize = true;
            end

            if(FLAG_normalize)
                % rewight only if normalize flag is set
                freqVecThisTP_reweighted = (freqVecThisTP./sum(freqVecThisTP));
                thisSeqFreqVec_reweighted(thisSeqTimeVec == timePointVec(t)) = freqVecThisTP_reweighted;
            else
                thisSeqFreqVec_reweighted(thisSeqTimeVec == timePointVec(t)) = freqVecThisTP;
            end
        end
 

        if(FLAG_normalize)
            disp(' ')
            disp('Reweighting required. After rewighting:')
            disp(' ')
            entryCounter = 0;
            for t = 1:numTimePoints
                freqVecThisTP = thisSeqFreqVec_reweighted(thisSeqTimeVec == timePointVec(t));
                temp1 = sum(freqVecThisTP);
                disp(['Time point ' num2str(timePointVec(t)) ', sum of frequenices of all haplotypes: ' num2str(temp1)])

                % make new headers

                for h = 1:length(freqVecThisTP)

                    thisHeader = HeaderNew{h+1};
                    indTemp = strfind(thisHeader, ', days:');
                    freqThisHaplo = freqVecThisTP(h);
                    entryCounter = entryCounter + 1;

                    % (thisSeqNumReadsVec(entryCounter-1)) cuz entryCounter starts
                    % from 1 in line #104 instead of 0. Reason is Header{1} is
                    % RefSeq
                    HeaderNew{entryCounter} = [thisHeader(1:indTemp(1)-1) ', days: ' num2str(timePointVec(t)) ', freq: ' num2str(freqThisHaplo) ', reads: ' num2str(thisSeqNumReadsVec(entryCounter)) ','];
                end
            end
        end
        
        % save files with new header
        % 2.4 save all seqs at all time points in 1 single fasta file
        
        fileNameFASTAHU = [thisFileNameThisPat(1:indOfDot(end)-1) '_HU' thisFileNameThisPat(indOfDot(end):end)];
        if(FLAG_normalize)
            indOfDot = strfind(fileNameFASTAHU, '.');
            fileNameFASTAHU = [fileNameFASTAHU(1:indOfDot(end)-1) '_NewFreq' fileNameFASTAHU(indOfDot(end):end)];
        end
        if(exist([dirNameAnalysis fileNameFASTAHU], 'file') == 2)
            delete([dirNameAnalysis fileNameFASTAHU])
        end
        
        fastawrite([dirNameAnalysis fileNameFASTAHU], HeaderNew, seqNT_AllNew);        
        fileNameFASTAHUCell{f} = fileNameFASTAHU;
    end
    disp('Modifying headers of FASTA files to contain Frequency information and copying to ANalysis folder...done')
    
    %------------------------------------------------------------------
    % 3. Generate new .txt file conatianing names of header updated
    %    fasta files with ref seq to be used by the alignment function

    if(FLAG_SaveFile == 1)
        fprintf('Saving .txt file containing names of fasta format data files with modified headers...')

        fileID = fopen([dirNameAnalysis fileNameFastaFilesWithHU],'w');
        for f = 1:numFiles
            if(f == 1)
                fprintf(fileID,'%s\n', '% This file specifies the filenames of the Header Updated FASTA files.');
                fprintf(fileID,'%s\n', '% the syntax is:');
                fprintf(fileID,'%s\n', '% filename1.ext');
                fprintf(fileID,'%s\n', '%');
                fprintf(fileID,'%s\n', '% Note: do not place a comma or semicolon at the end of filename');
                fprintf(fileID,'%s\n', '% -------------------------------------------------------------------------');
            end
            fprintf(fileID,'%s\n',fileNameFASTAHUCell{f}); 
        end
        fclose(fileID);
        disp('done.')
        disp(' ')
    elseif(FLAG_SaveFile == 0)
        disp('Warning: .txt data file not saved as FLAG_SaveFile flag not set.')
    else
        disp('Error: case undefined')
        pause
    end        
else
    disp('aborted.')
    disp('Warning: Exiting without running.')
    disp('         Analysis folder already contains output of Preprocesing step 1.')
    disp(' ')
    disp('Empty the analysis folder and rerun preprocessing steps.')
    disp(dirNameAnalysis)
    disp(' ')
end
