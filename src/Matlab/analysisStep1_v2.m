% Written: 02-March, 2020
% Author: M Saqib Sohail

% update 06-May-2020
%         Linear interpolation option added
% update 3-Jun-2020
%         streamlined inputs
% update 16-Jun-2020
%         option to chose reference sequence by FLAG_UserProvidedRefSeq
%         FLAG_UserProvidedRefSeq = true, user supplies refSeq as an ACGT
%         sequence as the 4rth input parameter to this function 
%         FLAG_UserProvidedRefSeq = false, ref seq is the consensus seq of
%         time point 1 

% function that implements MPL (binary approx), with epistasis

 
% 1. Load data
% 2. Calculate frequencies from MSA
% 3. Calculate terms required for MPL estimate
% 4. Set file names to save analysis results
% 5. Calculating MPL estimates
% 6. Save estimates, other data to file
%function analysisStep1_v2(fileNameContainingDirPath, priorConstSC, FLAG_stratonovich, FLAG_MarkAccessibility, FLAG_SaveFile, FLAG_SaveIntCovMtx, FLAG_useFreqEntry, FLAG_troubleShoot, FLAG_linearInt);
function analysisStep1_v2(fileNameContainingDirPath, priorConst, FLAG_vector, varargin)
    
if(length(priorConst) == 1)
    priorConstSC = priorConst(1);
    priorConstEpi = 0;
elseif(length(priorConst) == 2)
    priorConstSC = priorConst(1);
    priorConstEpi = priorConst(2);
    disp('Warning: Recombination probability not provided. Setting it to 0...press anykey to continue...')
    pause
    recombProb = 0;
elseif(length(priorConst) == 3)
    priorConstSC = priorConst(1);
    priorConstEpi = priorConst(2);
    recombProb = priorConst(3);
else
    disp('Wraning: priorConst should either contain 1 or 2 reg values only. using first two values as \gamma_{sc} and \gamma_{epi}.')
    priorConstSC = priorConst(1);
    priorConstEpi = priorConst(2);
end

if(length(FLAG_vector) == 7)
    FLAG_stratonovich = FLAG_vector(1);
    FLAG_MarkAccessibility = FLAG_vector(2);
    FLAG_UserProvidedRefSeq = FLAG_vector(3);
    FLAG_SaveIntCovMtx = FLAG_vector(4);
    FLAG_useFreqEntry = FLAG_vector(5);
    FLAG_troubleShoot = FLAG_vector(6);
    FLAG_linearInt = FLAG_vector(7);
    FLAG_Epi = false;
elseif(length(FLAG_vector) == 8)
    FLAG_stratonovich = FLAG_vector(1);
    FLAG_MarkAccessibility = FLAG_vector(2);
    FLAG_UserProvidedRefSeq = FLAG_vector(3);
    FLAG_SaveIntCovMtx = FLAG_vector(4);
    FLAG_useFreqEntry = FLAG_vector(5);
    FLAG_troubleShoot = FLAG_vector(6);
    FLAG_linearInt = FLAG_vector(7);
    FLAG_Epi = FLAG_vector(8);
end

if(FLAG_UserProvidedRefSeq == false && nargin == 4)
    disp('Warning: FLAG_UserProvidedRefSeq is UNSET. Ignoring user provided reference sequence. SET this flag to use user provided sequence.')
elseif(FLAG_UserProvidedRefSeq == true && nargin == 3)
    disp('Warning: FLAG_UserProvidedRefSeq is SET but no reference sequence is provided. Computing ref sequence from time point 1. ')
    FLAG_UserProvidedRefSeq = false;
elseif(FLAG_UserProvidedRefSeq == true && nargin == 4)
    userRefSequence = varargin{1};
end

FLAG_SaveFile = true;
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

runAnalysisCode = true;
if(exist([dirNameAnalysis 'Analysis_Misc' chosenSlash], 'dir') == 0)
    disp('Can not find Analysis_Misc folder. Run preprocessing Step 1.')
    disp(' Exiting without running Analysis code.')
    runAnalysisCode = false;
end

if(runAnalysisCode)
    disp('Running analysis Step 1 - ')
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
        end

    end

    %----------------------------------------------------------------------
    % 2. Calculate frequencies from MSA
    %----------------------------------------------------------------------
    
    uniqueTimePoints = unique(timePointVecFromFileName);
    numUniqueTimePoints = length(uniqueTimePoints);

    timeStep = uniqueTimePoints(2:end) - uniqueTimePoints(1:end-1);
    
    fprintf('Calculating one and two point probabilities...')
    for tp = 1:numUniqueTimePoints
        thisTP = uniqueTimePoints(tp);
        numbSampleThisTP = sum(timePointVecFromFileName == thisTP);
        
        % load MSA, calculate 1,2,3,4 pt frequencies from it
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

            %------------------------------------------------------------------
            % check if header is consistent with MPL format and extract
            % info from header

            % header consistency check
            tempHeader = Header{1};
            if(contains(tempHeader, 'freq: ') == 1 )    
                %disp('pass.') 
            else
                fprintf(['Header format of ' thisFileNameHUFiles])
                disp(' not in proper format. Cannot run MPL analysis code.') 
                disp(' ')
                disp('Make sure the headers in FASTA files are in proper format for MPL code. See readme.txt')
                break
            end    

            %------------------------------------------------------------------
            % get information from header
            %[thisSeqTimeVec, thisSeqFreqVec, thisSeqNumReadsVec] = getDayFreqNumReedsFromHeader(Header);
            [thisSeqTimeVec, thisSeqFreqVecTemp, thisSeqNumReadsVec] = getInfoFromHeader(Header);
            
            timePointVec = unique(thisSeqTimeVec);
            numTimePoints = length(timePointVec);

            %------------------------------------------------------------------
            % use numReeds to find frequency of each seq in the MSA
            %if(sum(thisSeqNumReadsVec == -1) > 0)
            if(FLAG_useFreqEntry == true)
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


            % make compressed MSA 
            msaNT = repmat(' ', size(seqNT_All, 2), length(seqNT_All{1}));
            for k = 1:size(seqNT_All, 2)
               msaNT(k,:) = seqNT_All{k};
            end


            msaNTInt = nt2int(msaNT);
            numSitesNT = size(msaNT,2);
            numUniqueSeqs = size(msaNT,1);
        
            if(tp == 1)
                freqCountNT = ntFreqCountCompMSA(msaNTInt, thisSeqFreqVec);
                [~, consensusSeqTP1] = max(freqCountNT);
            end
            
            if(tp == 1 && bsamp == 1)
                q11_all = zeros(numSitesNT, numSitesNT ,numUniqueTimePoints);
            end
            %------------------------------------------------------------------
            % 2.1 make binary MSA 

            %  possible states at each site: [0 1], consensus founder:0, mut:1
            %  each site represented by 1 
            %disp(' ')
            %fprintf('Converting NT MSA to binary...')
            
            numSitesNTExteded = numSitesNT;%numSitesNT*numNTAtEachSite;
            msaNT_Bin = -1*ones(numUniqueSeqs, numSitesNTExteded);
            
            if(FLAG_UserProvidedRefSeq == true)
                referenceSequenceForConversion = nt2int(userRefSequence);
            elseif(FLAG_UserProvidedRefSeq == false)
                referenceSequenceForConversion = consensusSeqTP1;
            end
            for j = 1:numSitesNT
                msaNT_Bin(:,j) = referenceSequenceForConversion(j) ~= msaNTInt(:,j);
            end
            
            %------------------------------------------------------------------
            % 2.2 Calculate 1,2,3,4 point probabilities, mutational flux
            
            % Calculate 1-point freqs
            q = thisSeqFreqVec*msaNT_Bin;
            q(q > 1) = 1;
            q(q < 0) = 0;
            indOfNonZeroQs = find(q);
            indOfNonZeroQsLen = length(indOfNonZeroQs);
            if(bsamp == 1)
                q_sp = sparse(q);
            else
                q_sp = q_sp + sparse(q);
            end        
            % needed for filtering based on 2,3,4 consecutive poly timepoints
            if(tp == 1 && bsamp == 1)
               q_All = zeros(numUniqueTimePoints, numSitesNT); 
               indOfNonZeroQs_All = [];
            end
            q_All(tp,:) = q_All(tp,:) + q;

            indOfNonZeroQs_All = unique([indOfNonZeroQs_All indOfNonZeroQs]);    


            %------------------------------------------------------------------
            % 2.2.1 Calculate 2-point freqs
            % q11_vec contains all non-zero 2 points in vector format
            % q11_sp is the sparse matrix representation of q11

            q11_vec = zeros(1, indOfNonZeroQsLen*indOfNonZeroQsLen);
    
            for kk = 1:indOfNonZeroQsLen
                thisInd = indOfNonZeroQs(kk);
                thisTemp = (thisSeqFreqVec.*msaNT_Bin(:,thisInd)')*msaNT_Bin(:,indOfNonZeroQs);
                thisTemp(thisTemp < 0) = 0;
                thisTemp(thisTemp > 1) = 1;

                if(sum(thisTemp < 0) > 0)
                    disp('Calculated q11 is less than 0...')
                    thisTemp(thisTemp < 0)
                    pause
                end

                if(sum(thisTemp > 1) > 0)
                    disp('Calculated q11 is greater than 1...')
                    thisTemp(thisTemp > 1) - 1
                    pause
                end

                q11_vec(1, (kk-1)*indOfNonZeroQsLen + 1:kk*indOfNonZeroQsLen) = thisTemp;
            end

            firstInd = repmat(indOfNonZeroQs, 1, indOfNonZeroQsLen);
            secondInd = reshape(repmat(indOfNonZeroQs, indOfNonZeroQsLen,1),1, indOfNonZeroQsLen*indOfNonZeroQsLen);
            if(bsamp == 1)
                q11_sp = sparse(firstInd, secondInd, q11_vec, numSitesNTExteded, numSitesNTExteded);
            else
                q11_sp = q11_sp + sparse(firstInd, secondInd, q11_vec, numSitesNTExteded, numSitesNTExteded);
            end
            q11_all(:,:,tp) = full(q11_sp);
            %-----------------------for Epistasis only---------------------
            % 2.2.2 Calculate 3 and 4-point freqs
            if(FLAG_Epi == true)
    
                % 2.2.2.1 reshape q11 in vector form to concatenate with q (needed
                % for inferring epistasis) 

                % first and second rows of twoPointSiteIDs contain the site IDs of
                % 2-point freq entries.

                twoPointSiteIDs = [secondInd(firstInd > secondInd);
                         firstInd(firstInd > secondInd)];
                numTwoPointSiteIDs = size(twoPointSiteIDs, 2);
                numTotalTwoPointEntries = numSitesNTExteded*(numSitesNTExteded-1)/2;

                ww = reshape(repmat((numSitesNT - (1:numSitesNT-1))', 1, 1)', 1, (numSitesNT-1));
                ww(1) = 0;
                www = cumsum(ww);

                if(~isempty(twoPointSiteIDs))
                    twoPointSiteIDs = [twoPointSiteIDs; www(twoPointSiteIDs(1,:))];
                    indOfNonZeroTwoPointInB_sp = sum(twoPointSiteIDs(2:3,:)) - 1;
                else
                    indOfNonZeroTwoPointInB_sp = [];
                end



                % reshape the upper triangle of q11_sp in vetor form to make the
                % vetor of single and double mutant freqs
                
                b_vec = zeros(1, numSitesNTExteded*(numSitesNTExteded-1)/2);

                % replace this with
                startInd = 1;
                for i = 1:numSitesNTExteded-1
                    stopInd = startInd + numSitesNTExteded-i - 1;
                    b_vec(startInd:stopInd) = q11_sp(i,i+1:end);
                    startInd = stopInd + 1;
                end
                
                if(bsamp == 1)
                    b_sp = sparse(b_vec);
                else
                    b_sp = b_sp + sparse(b_vec);
                end

                % 3.3 calculate 3-point freqs
                q111_vec = zeros(1, indOfNonZeroQsLen*numTwoPointSiteIDs);
                firstInd_q111 = zeros(1, indOfNonZeroQsLen*numTwoPointSiteIDs);
                secondInd_q111 = zeros(1, indOfNonZeroQsLen*numTwoPointSiteIDs);
                counter111 = 0;
                for kk = 1:indOfNonZeroQsLen
                    firstSiteInd = indOfNonZeroQs(kk);
                    for j = 1:numTwoPointSiteIDs
                        secondSiteInd = twoPointSiteIDs(1,j);
                        thirdSiteInd = twoPointSiteIDs(2,j);
                        counter111 = counter111 + 1;
                        q111_vec(counter111) = (thisSeqFreqVec.*msaNT_Bin(:,firstSiteInd)')*(msaNT_Bin(:,secondSiteInd).*msaNT_Bin(:,thirdSiteInd));

                        firstInd_q111(counter111) = firstSiteInd;
                        secondInd_q111(counter111) = indOfNonZeroTwoPointInB_sp(j);
                    end
                end
                
                if(bsamp == 1)
                    q111_sp = sparse(firstInd_q111(1:counter111), secondInd_q111(1:counter111), q111_vec(1:counter111), numSitesNTExteded, ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded);
                else
                    q111_sp = q111_sp + sparse(firstInd_q111(1:counter111), secondInd_q111(1:counter111), q111_vec(1:counter111), numSitesNTExteded, ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded);
                end
                numTwoPointSiteIDs;

                % 3.4 calculate 4-point freqs        
                q1111_vec = zeros(1, numTwoPointSiteIDs*numTwoPointSiteIDs);
                firstInd_q1111 = zeros(1, numTwoPointSiteIDs*numTwoPointSiteIDs);
                secondInd_q1111 = zeros(1, numTwoPointSiteIDs*numTwoPointSiteIDs);
                counter1111 = 0;
                for kk = 1:numTwoPointSiteIDs
                    firstSiteInd = twoPointSiteIDs(1,kk);
                    secondSiteInd = twoPointSiteIDs(2,kk);
                    for j = 1:numTwoPointSiteIDs
                        thirdSiteInd = twoPointSiteIDs(1,j);
                        fourthSiteInd = twoPointSiteIDs(2,j);

                        counter1111 = counter1111 + 1;
                        q1111_vec(counter1111) = (thisSeqFreqVec.*msaNT_Bin(:,firstSiteInd)')*(msaNT_Bin(:,secondSiteInd).*msaNT_Bin(:,thirdSiteInd).*msaNT_Bin(:,fourthSiteInd));

                        firstInd_q1111(counter1111) = indOfNonZeroTwoPointInB_sp(kk);
                        secondInd_q1111(counter1111) = indOfNonZeroTwoPointInB_sp(j);
                    end
                end
                        
                if(bsamp == 1)
                    q1111_sp = sparse(firstInd_q1111(1:counter1111), secondInd_q1111(1:counter1111), q1111_vec(1:counter1111), ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded, ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded);
                else
                    q1111_sp = q1111_sp + sparse(firstInd_q1111(1:counter1111), secondInd_q1111(1:counter1111), q1111_vec(1:counter1111), ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded, ((numSitesNTExteded+1)*numSitesNTExteded/2) - numSitesNTExteded);
                end
            end
        
        end
        
        % distance term for for recomb
        mtxDim = numSitesNT;
        vecLen = mtxDim*(mtxDim-1)/2;
        distVec = zeros(vecLen,1);

        counter = 0;
        for i = 1:(mtxDim-1)
            counter = counter + 1;
            distVec(counter:(counter+mtxDim-(i + 1))) = 1:(mtxDim-i);
            counter = counter + (mtxDim) - (i+1);
        end
        %-----------------------end for Epistasis only---------------------

        %----------------------------------------------------------------------
        % normalize to account for multiple bootstrap samples (freq values in q
        % and q11 should be in the range 0-1, but these will be in range
        % 0-numbSampleThisTP because of summing over all bootstrap samples.
        q_All(tp,:) = q_All(tp,:)/numbSampleThisTP;
        q_sp = q_sp/numbSampleThisTP;
        q11_sp = q11_sp/numbSampleThisTP;
        if(FLAG_Epi == true)
            q111_sp = q111_sp/numbSampleThisTP;
            q1111_sp = q1111_sp/numbSampleThisTP;
            b_sp = b_sp/numbSampleThisTP;
        end
        
        if((sum(q_sp > 1) > 0))
            q_sp(q_sp > 1) - 1
        
            pause
        end
        if(sum(q_sp < 0) > 0)
            q_sp(q_sp < 0)
        
            pause
        end
        
        
        %-----------------------------------------------------------------
        % 3. Calculate terms required for MPL estimate
        %-----------------------------------------------------------------
        %fprintf('Calculating mutation flux term...')
        xTemp1 = repmat(1:numSitesNT, 1, numSitesNT);
        xTemp2 = reshape(repmat(1:numSitesNT, numSitesNT, 1),1, numSitesNT*numSitesNT);
        xTemp = [xTemp2(xTemp1 > xTemp2);
                 xTemp1(xTemp1 > xTemp2)];

        if(FLAG_stratonovich == false)
            if(FLAG_linearInt == false)
                if(tp == 1)
                    % MPL
                    tempintCovMtx = (q11_sp - q_sp'*q_sp);
                    if(abs(tempintCovMtx) > 1)
                        disp('abs(tempintCovMtx) > 1')
                        tempintCovMtx(abs(tempintCovMtx) > 1)
                        pause
                    end
                    tempintCovMtx(abs(tempintCovMtx) < 1e-10) = 0;
                    tempintCovMtx(abs(tempintCovMtx) > 1) = 1;
                    intCovMtx = tempintCovMtx*timeStep(1);
        %                 whos q_sp
        %                 whos intCovMtx
        %                 pause
                    vEst1_noTS = sparse(1 - q_sp');
                    vEst2_noTS = - q_sp';
                    vEst1 = vEst1_noTS*timeStep(1);
                    vEst2 = vEst2_noTS*timeStep(1);
                    q_0 = q_sp';
                    % MPL-Epi
                    if(FLAG_Epi == true)
                        %intCovMtxEpi = [(q11_sp - q_sp'*q_sp) covE_sp; covE_sp' covF_sp]*timeStep(1);
                        covE_sp = q111_sp - q_sp'*b_sp;
                        covF_sp = q1111_sp - b_sp'*b_sp;
                        intCovMtxEpi = [tempintCovMtx covE_sp; covE_sp' covF_sp]*timeStep(1);
                        tempSpVec3 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(1,:)), numTotalTwoPointEntries, 1);
                        tempSpVec4 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(2,:)), numTotalTwoPointEntries, 1);

                        vEst3_Epi = (tempSpVec3 - b_sp')*timeStep(1); % [(tempSpVec1 + tempSpVec2 - 2*b_sp)']*timeStep(1);
                        vEst4_Epi = (tempSpVec4 - b_sp')*timeStep(1); %- [q_sp'; 2*b_sp']*timeStep(1);
                        vEst5_Epi = - b_sp'*timeStep(1);
                        
                        %cov4RecombTerm_sp_old = cov4RecombTerm_sp;
                        
                        qiqjVec = sparse(matUpperTriuToVec(q_sp'*q_sp));
                        cov4RecombTerm_sp = distVec.*(b_sp' - qiqjVec');
                        cov4RecombTerm_Epi = cov4RecombTerm_sp*timeStep(1);
                        qBar_0 = [q_sp'; b_sp'];
                    end
                elseif(tp == numUniqueTimePoints)
                    % MPL
                    q_T = q_sp';
                    % MPL-Epi
                    if(FLAG_Epi == true)
                        qBar_T = [q_sp';b_sp'];
                    end
                else
                    % MPL
                    tempintCovMtx = (q11_sp - q_sp'*q_sp);
                    if(abs(tempintCovMtx) > 1)
                        disp('abs(tempintCovMtx) > 1')
                        tempintCovMtx(abs(tempintCovMtx) > 1)
                        pause
                    end
                    tempintCovMtx(abs(tempintCovMtx) < 1e-10) = 0;
                    tempintCovMtx(abs(tempintCovMtx) > 1) = 1;
                    
                    intCovMtx = intCovMtx +  tempintCovMtx*timeStep(tp);
                    vEst1_noTS = sparse(1 - q_sp');
                    vEst2_noTS = - q_sp';
                    vEst1 = vEst1 + vEst1_noTS*timeStep(tp);
                    vEst2 = vEst2 + vEst2_noTS*timeStep(tp);
                    
                    % MPL-Epi
                    if(FLAG_Epi == true)
                        %intCovMtxEpi = intCovMtxEpi +  [(q11_sp - q_sp'*q_sp) covE_sp; covE_sp' covF_sp]*timeStep(t);
                        covE_sp = q111_sp - q_sp'*b_sp;
                        covF_sp = q1111_sp - b_sp'*b_sp;
                        intCovMtxEpi = intCovMtxEpi +  [tempintCovMtx covE_sp; covE_sp' covF_sp]*timeStep(tp);
                        tempSpVec3 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(1,:)), numTotalTwoPointEntries, 1);
                        tempSpVec4 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(2,:)), numTotalTwoPointEntries, 1);
                        vEst3_Epi = vEst3_Epi + (tempSpVec3 - b_sp')*timeStep(tp); % [(tempSpVec1 + tempSpVec2 - 2*b_sp)']*timeStep(1);
                        vEst4_Epi = vEst4_Epi + (tempSpVec4 - b_sp')*timeStep(tp); %- [q_sp'; 2*b_sp']*timeStep(1);
                        vEst5_Epi = vEst5_Epi - b_sp'*timeStep(tp);

                        %cov4RecombTerm_sp_old = cov4RecombTerm_sp;
                        
                        qiqjVec = sparse(matUpperTriuToVec(q_sp'*q_sp));
                        cov4RecombTerm_sp = distVec.*(b_sp' - qiqjVec');
                        
                        cov4RecombTerm_Epi = cov4RecombTerm_Epi + cov4RecombTerm_sp*timeStep(tp);
                    end
                end
            elseif(FLAG_linearInt == true) % linear interpolation
                guessOfNonZeroEntriesOfIntCovMtx = ceil((numSitesNTExteded*numSitesNTExteded)*0.07);
                if(FLAG_Epi == true)
                    guessOfNonZeroEntriesOfIntCovMtxEpi = ceil((numSitesNTExteded*(numSitesNTExteded+1)/2)*(numSitesNTExteded*(numSitesNTExteded+1)/2)*0.001);
                    guessOfNonZeroEntriesOfvEst345Epi = ceil(numSitesNTExteded*(numSitesNTExteded-1)/2*0.75);
                end
                if(tp == 1)
                    % MPL
                    last_q11_sp = q11_sp;
                    last_q_sp = q_sp';
                    q_0 = q_sp';                

                    vEst1 = sparse(1, 1, 0, numSitesNTExteded, 1, numSitesNTExteded);
                    vEst2 = sparse(1, 1, 0, numSitesNTExteded, 1, numSitesNTExteded);
                    intCovMtx = spalloc(numSitesNTExteded, numSitesNTExteded,guessOfNonZeroEntriesOfIntCovMtx);
                    cov4RecombTerm_Epi = sparse(numSitesNTExteded*(numSitesNTExteded-1)/2, 1, ceil(numSitesNTExteded*(numSitesNTExteded-1)/2*0.1));
                    
                    % MPL-Epi
                    if(FLAG_Epi == true)
                        last_q111_sp = q111_sp;
                        last_q1111_sp = q1111_sp;
                        last_b_sp = b_sp';
                        qBar_0 = [q_sp'; b_sp'];
                        
                        last_x11_sp = [last_q11_sp last_q111_sp; last_q111_sp' last_q1111_sp];
                        last_x_sp = [last_q_sp; last_b_sp];
                        
                        intCovMtxEpi = spalloc(numSitesNTExteded*(numSitesNTExteded+1)/2, numSitesNTExteded*(numSitesNTExteded+1)/2, guessOfNonZeroEntriesOfIntCovMtxEpi);
                        vEst3_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                        vEst4_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                        vEst5_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                    end
                else
                    % MPL
                    this_q11_sp = q11_sp;
                    this_q_sp = q_sp';
                    q_T = q_sp'; % so that finally q_T is set for t == numFileNameToAnalyzeCell

                    q_sp_mid_point = (last_q_sp + this_q_sp)/2;

                    thisCovMtx = (last_q11_sp + this_q11_sp)/2 - (last_q_sp*last_q_sp'/3 + this_q_sp*this_q_sp'/3 + last_q_sp*this_q_sp'/6 + this_q_sp*last_q_sp'/6);
                    thisCovMtx(abs(thisCovMtx) < 1e-10) = 0;
                    if(sum(sum(abs(thisCovMtx) > 1)) > 0)
                        disp('Warning: some entries of abs(lastCovMtx) > 1')
                        thisCovMtx(abs(thisCovMtx) > 1)
                        pause
                        thisCovMtx(abs(thisCovMtx) > 1) = 1;
                    end

                    intCovMtx = intCovMtx + thisCovMtx*timeStep(tp-1);

                    vEst1 = vEst1 + (sparse(1 - q_sp_mid_point))*timeStep(tp-1);
                    vEst2 = vEst2 + (-q_sp_mid_point)*timeStep(tp-1);

                    last_q11_sp = this_q11_sp;
                    last_q_sp = this_q_sp;
                    
                    % MPL-Epi
                    if(FLAG_Epi == true)
                        this_q111_sp = q111_sp;
                        this_q1111_sp = q1111_sp;
                        this_b_sp = b_sp';
                        qBar_T = [this_q_sp; this_b_sp]; % so that finally q_T is set for t == numFileNameToAnalyzeCell
                        
                        this_x11_sp = [this_q11_sp this_q111_sp; this_q111_sp' this_q1111_sp];
                        this_x_sp = [this_q_sp; this_b_sp];
                                                
                        thisCovMtxEpi = (last_x11_sp + this_x11_sp)/2 - (last_x_sp*last_x_sp'/3 + this_x_sp*this_x_sp'/3 + last_x_sp*this_x_sp'/6 + this_x_sp*last_x_sp'/6);
                                                
                        intCovMtxEpi = intCovMtxEpi + thisCovMtxEpi*timeStep(tp-1);
                        
                        tempSpVec3 = sparse(1:numTotalTwoPointEntries, 1, (last_q_sp(xTemp(1,:)) + this_q_sp(xTemp(1,:)))/2, numTotalTwoPointEntries, 1);
                        tempSpVec4 = sparse(1:numTotalTwoPointEntries, 1, (last_q_sp(xTemp(2,:)) + this_q_sp(xTemp(2,:)))/2, numTotalTwoPointEntries, 1);
                        b_sp_mid_point = (last_b_sp + this_b_sp)/2;
                        
                        vEst3_Epi = vEst3_Epi + (tempSpVec3 - b_sp_mid_point)*timeStep(tp-1);
                        vEst4_Epi = vEst4_Epi + (tempSpVec4 - b_sp_mid_point)*timeStep(tp-1);
                        vEst5_Epi = - b_sp_mid_point*timeStep(tp-1);
                        
                        recTerm2_Epi = sparse(matUpperTriuToVec(last_q_sp*last_q_sp'/3));
                        recTerm3_Epi = sparse(matUpperTriuToVec(this_q_sp*this_q_sp'/3));
                        recTerm4_Epi = sparse(matUpperTriuToVec(this_q_sp*last_q_sp'/6));
                        recTerm5_Epi = sparse(matUpperTriuToVec(last_q_sp*this_q_sp'/6));
                        
                        cov4RecombTerm_sp = distVec.*(b_sp_mid_point - (recTerm2_Epi' + recTerm3_Epi' + recTerm4_Epi' + recTerm5_Epi'));
                        cov4RecombTerm_Epi = cov4RecombTerm_Epi + cov4RecombTerm_sp*timeStep(tp-1);
                        
                        last_b_sp = this_b_sp;
                        last_q111_sp = this_q111_sp;
                        last_q1111_sp = this_q1111_sp;
                        last_x_sp = this_x_sp;
                        last_x11_sp = this_x11_sp;
                    end
                end    
            end
        elseif(FLAG_stratonovich == true)
            guessOfNonZeroEntriesOfIntCovMtx = ceil((numSitesNTExteded*numSitesNTExteded)*0.07);
            if(FLAG_Epi == true)
                guessOfNonZeroEntriesOfIntCovMtxEpi = ceil((numSitesNTExteded*(numSitesNTExteded+1)/2)*(numSitesNTExteded*(numSitesNTExteded+1)/2)*0.001);
                guessOfNonZeroEntriesOfvEst345Epi = ceil(numSitesNTExteded*(numSitesNTExteded-1)/2*0.75);
            end
            if(tp == 1) % tell Saqib bhai
                % MPL
                lastCovMtx = (q11_sp - q_sp'*q_sp);
                lastCovMtx(abs(lastCovMtx) < 1e-10) = 0;
                if(sum(sum(abs(lastCovMtx) > 1)) > 0)
                    disp('Warning: some entries of abs(lastCovMtx) > 1')
                    pause
                    lastCovMtx(abs(lastCovMtx) > 1) = 1;
                end
                
                q_0 = q_sp';
                last_vEst1_noTS = sparse(1 - q_sp');
    
                last_vEst2_noTS = - q_sp';
                vEst1 = sparse(1, 1, 0, numSitesNTExteded, 1, numSitesNTExteded);
                vEst2 = sparse(1, 1, 0, numSitesNTExteded, 1, numSitesNTExteded);
                intCovMtx = spalloc(numSitesNTExteded, numSitesNTExteded,guessOfNonZeroEntriesOfIntCovMtx);
                cov4RecombTerm_Epi = sparse(numSitesNTExteded*(numSitesNTExteded-1)/2, 1, ceil(numSitesNTExteded*(numSitesNTExteded-1)/2*0.1));
                
                % MPL-Epi
                if(FLAG_Epi == true)
                    %lastCovMtxEpi = [(q11_sp - q_sp'*q_sp) covE_sp; covE_sp' covF_sp];
                    covE_sp = q111_sp - q_sp'*b_sp;
                    covF_sp = q1111_sp - b_sp'*b_sp;
                    lastCovMtxEpi = [lastCovMtx covE_sp; covE_sp' covF_sp];
                    qBar_0 = [q_sp'; b_sp'];
                    tempSpVec3 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(1,:)), numTotalTwoPointEntries, 1);
                    tempSpVec4 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(2,:)), numTotalTwoPointEntries, 1);
                    last_vEst3_noTS_Epi = (tempSpVec3 - b_sp');
                    last_vEst4_noTS_Epi = (tempSpVec4 - b_sp');
                    last_vEst5_noTS_Epi = - b_sp';
                    
                    if(last_vEst3_noTS_Epi(abs(last_vEst3_noTS_Epi) < 1e-6))
                        disp('Warning: last_vEst3_noTS_Epi(abs(last_vEst3_noTS_Epi) < 1e-6)')
                        asdf
                    end
                    if(last_vEst4_noTS_Epi(abs(last_vEst4_noTS_Epi) < 1e-6))
                        disp('Warning: last_vEst4_noTS_Epi(abs(last_vEst4_noTS_Epi) < 1e-6)')
                        asdf
                    end
                    
                    intCovMtxEpi = spalloc(numSitesNTExteded*(numSitesNTExteded+1)/2, numSitesNTExteded*(numSitesNTExteded+1)/2, guessOfNonZeroEntriesOfIntCovMtxEpi);
                    vEst3_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                    vEst4_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                    vEst5_Epi = sparse(1, 1, 0, numSitesNTExteded*(numSitesNTExteded-1)/2, 1, guessOfNonZeroEntriesOfvEst345Epi);
                    
                    qiqjVec = sparse(matUpperTriuToVec(q_sp'*q_sp));
                    cov4RecombTerm_sp = distVec.*(b_sp' - qiqjVec');
                    last_cov4RecombTerm_sp_noTS_Epi = cov4RecombTerm_sp;
                end
            else
                % MPL
                thisCovMtx = (q11_sp - q_sp'*q_sp);
                thisCovMtx(abs(thisCovMtx) < 1e-10) = 0;
                if(sum(sum(abs(lastCovMtx) > 1)) > 0)
                    disp('Warning: some entries of abs(lastCovMtx) > 1')
                    pause
                    lastCovMtx(abs(lastCovMtx) > 1) = 1;
                end
                %thisCovMtx = (q11_sp - q_sp'*q_sp);
                thisCovMtxMidPoint = 0.5*thisCovMtx + 0.5*lastCovMtx;
                intCovMtx = intCovMtx +  thisCovMtxMidPoint*timeStep(tp-1);



                sumOfAbsIntCovEntriesDiag1 = diag(intCovMtx);
                if(sum(sumOfAbsIntCovEntriesDiag1 < 0) > 0)
                    temp = find(sumOfAbsIntCovEntriesDiag1 < 0);
                    disp('Error: sumOfAbsIntCovEntriesDiag1 < 0')
                    disp(['      at ' num2str(temp)])

                    q_All(:,temp)
                    intCovMtx(temp, temp)
                    thisCovMtx(temp, temp)
                    lastCovMtx(temp, temp)
                    temp
                    q11_sp(temp)
                    q_sp(temp) - 1
                    pause
                end

                lastCovMtx = thisCovMtx;
                q_T = q_sp'; % so that finally q_T is set for t == numFileNameToAnalyzeCell

                this_vEst1_noTS = sparse(1 - q_sp');
                this_vEst1_noTS_MidPoint = 0.5*this_vEst1_noTS + 0.5*last_vEst1_noTS;
                vEst1 = vEst1 + this_vEst1_noTS_MidPoint*timeStep(tp-1);

                this_vEst2_noTS = - q_sp';
                this_vEst2_noTS_MidPoint = 0.5*this_vEst2_noTS + 0.5*last_vEst2_noTS;
                vEst2 = vEst2 + this_vEst2_noTS_MidPoint*timeStep(tp-1);
                last_vEst1_noTS = this_vEst1_noTS;
                last_vEst2_noTS = this_vEst2_noTS;
                
                % MPL-Epi
                if(FLAG_Epi == true)
                    %thisCovMtxEpi = [(q11_sp - q_sp'*q_sp) covE_sp; covE_sp' covF_sp];
                    covE_sp = q111_sp - q_sp'*b_sp;
                    covF_sp = q1111_sp - b_sp'*b_sp;
                    thisCovMtxEpi = [thisCovMtx covE_sp; covE_sp' covF_sp];
                    thisCovMtxMidPointEpi = 0.5*thisCovMtxEpi + 0.5*lastCovMtxEpi;
                    intCovMtxEpi = intCovMtxEpi +  thisCovMtxMidPointEpi*timeStep(tp-1);

                    lastCovMtxEpi = thisCovMtxEpi;
                    qBar_T = [q_sp';b_sp']; % so that finally q_T is set for t == numDataFiles

                    tempSpVec3 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(1,:)), numTotalTwoPointEntries, 1);
                    this_vEst3_noTS_Epi = (tempSpVec3 - b_sp');
                    this_vEst3_noTS_MidPointEpi = 0.5*this_vEst3_noTS_Epi + 0.5*last_vEst3_noTS_Epi;
                    vEst3_Epi = vEst3_Epi + this_vEst3_noTS_MidPointEpi*timeStep(tp-1);

                    if(this_vEst3_noTS_Epi(abs(this_vEst3_noTS_Epi) < 1e-6))
                        disp('Warning: this_vEst3_noTS_Epi(abs(this_vEst3_noTS_Epi) < 1e-6)')
                        asdf
                    end
                    if(this_vEst3_noTS_MidPointEpi(abs(this_vEst3_noTS_MidPointEpi) < 1e-6))
                        disp('Warning: this_vEst3_noTS_MidPointEpi(abs(this_vEst3_noTS_MidPointEpi) < 1e-6)')
                        asdf
                    end
                    
                    tempSpVec4 = sparse(1:numTotalTwoPointEntries, 1, q_sp(xTemp(2,:)), numTotalTwoPointEntries, 1);
                    this_vEst4_noTS_Epi = (tempSpVec4 - b_sp');
                    this_vEst4_noTS_MidPointEpi = 0.5*this_vEst4_noTS_Epi + 0.5*last_vEst4_noTS_Epi;
                    vEst4_Epi = vEst4_Epi + this_vEst4_noTS_MidPointEpi*timeStep(tp-1);

                    if(this_vEst4_noTS_Epi(abs(this_vEst4_noTS_Epi) < 1e-6))
                        disp('Warning: this_vEst4_noTS_Epi(abs(this_vEst4_noTS_Epi) < 1e-6)')
                        asdf
                    end
                    
                    this_vEst5_noTS_Epi = - b_sp';
                    this_vEst5_noTS_MidPointEpi = 0.5*this_vEst5_noTS_Epi + 0.5*last_vEst5_noTS_Epi;
                    vEst5_Epi = vEst5_Epi + this_vEst5_noTS_MidPointEpi*timeStep(tp-1);
                    
                    last_vEst3_noTS_Epi = this_vEst3_noTS_Epi;
                    last_vEst4_noTS_Epi = this_vEst4_noTS_Epi;
                    last_vEst5_noTS_Epi = this_vEst5_noTS_Epi;

                    qiqjVec = sparse(matUpperTriuToVec(q_sp'*q_sp));
                    cov4RecombTerm_sp = distVec.*(b_sp' - qiqjVec');
                    
                    this_cov4RecombTerm_sp_noTS_Epi = cov4RecombTerm_sp;
                    this_cov4RecombTerm_sp_noTS_MidPointEpi = 0.5*this_cov4RecombTerm_sp_noTS_Epi + 0.5*last_cov4RecombTerm_sp_noTS_Epi;
                    cov4RecombTerm_Epi = cov4RecombTerm_Epi + this_cov4RecombTerm_sp_noTS_MidPointEpi*timeStep(tp-1);
                    
                    last_cov4RecombTerm_sp_noTS_Epi = this_cov4RecombTerm_sp_noTS_Epi;
                end
            end
        end
    end
    disp('done.')


    %----------------------------------------------------------------------
    % 4. Set file names to save analysis results
    %----------------------------------------------------------------------
    
    % check accessibility of estimates and save to file
    if(FLAG_stratonovich == true)
        convention = 'Stratonovich';
    elseif(FLAG_stratonovich == false && FLAG_linearInt == false)
        convention = 'Ito';
    elseif(FLAG_stratonovich == false && FLAG_linearInt == true)
        convention = 'LinearInter';
    end
    
    thisFileNameHUFiles = fileNameHUFilesCell{1};
    indOfDash = strfind(thisFileNameHUFiles, '_');
    
    fileNameSelEst = ['SelEst_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.txt'];
    fileNameSelEstSL = ['SelEstSL_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.txt'];
    if(FLAG_Epi == true)
        fileNameSelEstEpi = ['SelEstEpi_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '_' num2str(priorConstEpi) '.txt'];
    end
    fileNameIntCovMtx = ['IntCovMtx_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.txt'];
    fileNameAllTrajs = ['AllTrajs_' thisFileNameHUFiles(1:indOfDash(2)-1) '.txt'];
    fileNameAllTrajsR = ['AllTrajsWithTimeInfo_' thisFileNameHUFiles(1:indOfDash(2)-1) '.txt'];
    fileNameAccessEst = ['AccessibilityMPL_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '.txt'];
    fileNameAccessEstEpi = ['AccessibilityMPLEpi_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '.txt'];
    
    if(FLAG_Epi == false)
        fileNameTestWorkspace = ['TestWorkspace_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.mat'];
    elseif(FLAG_Epi == true)
        fileNameTestWorkspace = ['TestWorkspace_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '_' num2str(priorConstEpi) '.mat'];
    end
        
    if(FLAG_Epi == false)
        if(exist([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEst], 'file') == 2)
           disp('Accessibility file already exists. Not computing accessibility...') 
           FLAG_MarkAccessibility = false;
        end
    elseif(FLAG_Epi == true)
        cond1_temp = exist([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEst], 'file') == 2;
        cond2_temp = exist([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEstEpi], 'file') == 2;
        if(cond1_temp && cond2_temp)
           disp('Accessibility file already exists. Not computing accessibility...') 
           FLAG_MarkAccessibility = false;
        end
    end
    
    %----------------------------------------------------------------------
    % 4.1 Calculate accessibitily
    % right now, can not handle too large matricies
    if(FLAG_MarkAccessibility && (size(intCovMtx, 1) < 5000))
       disp(' ')
       disp('Calculate accessibility...NOTE: Accessibility code works for additive model with < 5000 poly sites') 
       
       [selcSites, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(intCovMtx);
       
       wellCondSitesLen = length(wellCondSites);
       illCondSitesLen = length(illCondSites);
       numClusters = length(cluster);
       
       if(isempty(wellCondSitesLen))
           wellCondSitesLen = 0;
       end
       
       if(isempty(illCondSitesLen))
           illCondSitesLen = 0;
       end
       
       if(~isempty(numClusters))
           clusterLengths = zeros(1, numClusters);
           for cc = 1:numClusters
               clusterLengths(cc) = length(cluster{cc});
           end
           numColAccDataSave = 2 + numClusters;
       else
           numClusters = 0;
           clusterLengths = 0;
           numColAccDataSave = 2;
       end
       
       maxEntriesTemp = max([wellCondSitesLen, illCondSitesLen, clusterLengths]);
           
       accDataMtxToWrite = zeros(maxEntriesTemp, numColAccDataSave);
       accDataMtxToWrite(:,1) = [wellCondSites zeros(1, maxEntriesTemp-wellCondSitesLen)]';
       accDataMtxToWrite(:,2) = [illCondSites zeros(1, maxEntriesTemp-illCondSitesLen)]';
       for cc = 1:numClusters
           accDataMtxToWrite(:,2+cc) = [cluster{cc} zeros(1, maxEntriesTemp-clusterLengths(cc))]';
       end
       
       if(exist([dirNameAnalysis 'Estimates' chosenSlash], 'dir') == 0)
            mkdir([dirNameAnalysis 'Estimates' chosenSlash])
        end
       dlmwrite([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEst], accDataMtxToWrite);
       
       % for Epi
       if(FLAG_Epi == true)
           if(size(intCovMtxEpi, 1) < 50)
               [selcSitesEpi, wellCondSitesEpi, illCondSitesEpi, clusterEpi] = findIndColsOfHmat(intCovMtxEpi);

               wellCondSitesEpiLen = length(wellCondSitesEpi);
               illCondSitesEpiLen = length(illCondSitesEpi);
               numClustersEpi = length(clusterEpi);

               if(isempty(wellCondSitesEpiLen))
                   wellCondSitesEpiLen = 0;
               end

               if(isempty(illCondSitesEpiLen))
                   illCondSitesEpiLen = 0;
               end

               if(~isempty(numClustersEpi))
                   clusterEpiLengths = zeros(1, numClustersEpi);
                   for cc = 1:numClustersEpi
                       clusterEpiLengths(cc) = length(clusterEpi{cc});
                   end
                   numColAccDataSaveEpi = 2 + numClustersEpi;
               else
                   numClustersEpi = 0;
                   clusterEpiLengths = 0;
                   numColAccDataSaveEpi = 2;
               end

               maxEntriesTempEpi = max([wellCondSitesEpiLen, illCondSitesEpiLen, clusterEpiLengths]);

               accDataMtxToWriteEpi = zeros(maxEntriesTempEpi, numColAccDataSaveEpi);
               accDataMtxToWriteEpi(:,1) = [wellCondSitesEpi zeros(1, maxEntriesTempEpi-wellCondSitesEpiLen)]';
               accDataMtxToWriteEpi(:,2) = [illCondSitesEpi zeros(1, maxEntriesTempEpi-illCondSitesEpiLen)]';
               for cc = 1:numClustersEpi
                   accDataMtxToWriteEpi(:,2+cc) = [clusterEpi{cc} zeros(1, maxEntriesTempEpi-clusterEpiLengths(cc))]';
               end

               if(exist([dirNameAnalysis 'Estimates' chosenSlash], 'dir') == 0)
                    mkdir([dirNameAnalysis 'Estimates' chosenSlash])
               end
               dlmwrite([dirNameAnalysis 'Estimates' chosenSlash  fileNameAccessEstEpi], accDataMtxToWriteEpi);
           else
               disp('IntCovMtEpi size > 50x50...case not handled, Accessibility not computed.')
           end
       end
    end
    
    %----------------------------------------------------------------------
    % 5. Calculating MPL estimates
    %----------------------------------------------------------------------
    
    % load mutation vectors
    fileNameMutVec = [thisFileNameHUFiles(1:indOfDash(2)-1) '_mutVecs.txt'];
    
    mutationVectors = dlmread([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVec]);
    mutVecWT2Mut = mutationVectors(1,:)';
    mutVecMut2WT = mutationVectors(2,:)';

    if(FLAG_Epi == true)
        fileNameMutVecEpi = [thisFileNameHUFiles(1:indOfDash(2)-1) '_mutVecsEpi.txt'];
        mutationVectorsEpi = dlmread([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameMutVecEpi]);
        mutVecWT2Mut_j = mutationVectorsEpi(1,:)';
        mutVecWT2Mut_i = mutationVectorsEpi(2,:)';
        mutVecMut2WT_iPlusj = mutationVectorsEpi(3,:)';
    end

    % MPL Estimate in sparse matrix format
    vEstMu = vEst1.*mutVecWT2Mut + vEst2.*mutVecMut2WT;
    regMtx = sparse(1:numSitesNTExteded, 1:numSitesNTExteded, priorConstSC*ones(numSitesNTExteded,1), numSitesNTExteded, numSitesNTExteded); 
    numer = (q_T - q_0 - vEstMu);
    denom = (intCovMtx + regMtx);
    selEst = denom\numer;
    denomSL = diag(diag(intCovMtx + regMtx));
    selEstSL = denomSL\numer;
    
    if(FLAG_troubleShoot == 1)
        selEstNoMu = denom\(q_T - q_0);
        selEstSLNoMu = denomSL\(q_T - q_0);
    end

    if(FLAG_Epi == true)
        vEstMuEpiTemp = vEst3_Epi.*mutVecWT2Mut_j + vEst4_Epi.*mutVecWT2Mut_i + vEst5_Epi.*mutVecMut2WT_iPlusj;
        vEstMuEpi = [vEstMu; vEstMuEpiTemp];
        etaEstREpi = [sparse(zeros(numSitesNT, 1)); cov4RecombTerm_Epi]*recombProb;

        regVec = [priorConstSC*ones(numSitesNTExteded,1);  priorConstEpi*ones((numSitesNTExteded*(numSitesNTExteded-1)/2),1)];
    
        regMtxEpi = sparse(1:(numSitesNTExteded*(numSitesNTExteded+1)/2), 1:(numSitesNTExteded*(numSitesNTExteded+1)/2), regVec, (numSitesNTExteded*(numSitesNTExteded+1)/2), (numSitesNTExteded*(numSitesNTExteded+1)/2));
        numerEpi = (qBar_T - qBar_0 - vEstMuEpi - etaEstREpi);
        denomEpi = (intCovMtxEpi + regMtxEpi);
        selEstEpi = denomEpi\numerEpi;
    
%        numerEpi1 = (qBar_T - qBar_0);
%        numerEpi2 = (- vEstMuEpi);
%        numerEpi3 = (- etaEstREpi);

%        sumOfIntCovEntriesEpi = (sum(intCovMtxEpi(1:numSitesNT, (numSitesNT+1):end), 2));
%        sumOfAbsIntCovEntriesEpi = (sum(abs(intCovMtxEpi(1:numSitesNT, (numSitesNT+1):end)), 2));
%        sumOfAbsIntCovEntriesLink = (sum(abs(intCovMtxEpi(1:numSitesNT, 2:numSitesNT)), 2));
%        sumOfAbsIntCovEntriesDiag = diag(intCovMtxEpi(1:numSitesNT,1:numSitesNT));
%        if(sum(sumOfAbsIntCovEntriesDiag < 0) > 0)
%            temp = find(sumOfAbsIntCovEntriesDiag < 0);
%            disp('Error: sumOfAbsIntCovEntriesDiag < 0')
%            disp(['      at ' num2str(temp)])
% 
%            intCovMtxEpi(temp, temp)
%            q_All(:,temp)
%            intCovMtx(temp, temp);
%            (q11_sp - q_sp'*q_sp)
%            pause
%        end
    
%        if(FLAG_troubleShoot == true)
%            fprintf('Calculating the inverse of the int. cov. matrix...')
%            invIntCovMtxEpi = inv(denomEpi);
%            disp('done.')
%        end
 
    end
    
    
    %-----------------------------------------------------------------
    % 6. Save estimates, other data to file
    %-----------------------------------------------------------------
    
    if(FLAG_SaveFile)
        fprintf('Saving selection coefficient estimates to file...')
        
        if(exist([dirNameAnalysis 'Estimates' chosenSlash], 'dir') == 0)
            mkdir([dirNameAnalysis 'Estimates' chosenSlash])
        end
        
        
        if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEst], 'file') == 2)
           delete([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEst]) 
        end
        if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSL], 'file') == 2)
           delete([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSL]) 
        end
        dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEst], full(selEst))
        dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSL], full(selEstSL))
        
        if(FLAG_SaveIntCovMtx == 1)
            if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameIntCovMtx], 'file') == 2)
               delete([dirNameAnalysis 'Estimates' chosenSlash fileNameIntCovMtx]) 
            end
            dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameIntCovMtx], full(intCovMtx))
        end
        q_All = round(q_All*10000)/10000;
        
        if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajs], 'file') == 2)
           delete([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajs]) 
        end
        dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajs], q_All);
        
        if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajsR], 'file') == 2)
           delete([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajsR]) 
        end
        dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajsR], [uniqueTimePoints' q_All]);
        
        if(FLAG_Epi == true)
            dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstEpi], full(selEstEpi))
        end
    
        if(FLAG_troubleShoot == true)
            if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameTestWorkspace], 'file') == 2)
               delete([dirNameAnalysis 'Estimates' chosenSlash fileNameTestWorkspace]) 
            end
            save([dirNameAnalysis 'Estimates' chosenSlash fileNameTestWorkspace])
            
            fileNameSelEstNoMu = ['SelEstNoMu_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.txt'];
            fileNameSelEstSLNoMu = ['SelEstSLNoMu_' thisFileNameHUFiles(1:indOfDash(2)-1) '_' convention '_gamma' num2str(priorConstSC) '.txt'];
            if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstNoMu], 'file') == 2)
               delete([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstNoMu]) 
            end
            
            dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstNoMu], full(selEstNoMu))
            if(exist([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSLNoMu], 'file') == 2)
               delete([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSLNoMu]) 
            end
            dlmwrite([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSLNoMu], full(selEstSLNoMu))
        end
        disp('done.')
    end
end