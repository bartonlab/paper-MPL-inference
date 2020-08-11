clc
clear all
close all
warning off


%========================== INITIALIZATION ================================

% -------------------------- User specified -------------------------------
setAll = 9860001%[8550001 9550001 9650001 9750001 9850001 9860001 9950001];
strainsAll = 5%[1 5 5 5 5 5 5];
%Tused = 300;
tjSamplingSchemeStr = 'scheme1';
ngSamplingSchemeStr = 'schemeD';
dTSamplingSchemeStr = 'scheme33';
% -- ngSamplingSchemeStr --
% schemeA: pBino = 0.0075
% schemeB: pBino = 0.0087
% schemeC: pBino = 0.0095
% schemeD: pBino = 0.0139 <---Data mean
% schemeE: pBino = 0.015
% schemeF: pBino = 0.0187
% schemeG: pBino = 0.02
% -- dTSamplingSchemeStr --
% scheme2:  weight1 = 0.92, const = 118 (small mean dTVec)
% scheme3:  weight1 = 0.87, const = 120 (small mean dTVec)
% scheme33:  weight1 = 0.87, const = 120 <------------------ Data mean
% scheme55:  weight1 = 0.80, const = 124 (large mean dTVec)

allConventions = 3%[1 2 3];
allConvLen = length(allConventions);

for ss = 1:length(setAll)
    thisSet = setAll(ss)

    numStrainsInInitialPop = strainsAll(ss);

    for c = 1:length(allConventions)
        % chose convention 1: Ito, 2: Stratonovich, 3: Linear interpolation
        % Stratonovich has increased robustness to sampling effects than Ito
        % Linear interpolation has most increased robustness to sampling effects
        setConvention = allConventions(c);

        priorConstSC = 5; % this is the strength of the SC regularization term

        thisGenomicSegStartInd = 1; % this is the starting location of the first NT of the protein in the whole genome
        thisGenomicSegStopInd = 50; % this is the ending location of the last NT of the protein in the whole genome

        % for GAG use these 
        % thisGenomicSegStartInd = 790;
        % thisGenomicSegStopInd = 2289;

        % this file contains the NT-to-NT mutation probability. It must be located
        % in the folder .../MPL Pipeline/Data_Misc/MutationProbabilities/
        fileNameContainingMutProb = 'MutProb_HIV_SyntheticData_1eminus4.txt';
        recombProb = 0;%1e-4; % recombination probability

        FLAG_UserProvidedRefSeq = true; % SET: user provides reference sequence in ACGT- form
        referenceSequence = repmat('A', 1, 50);
        FLAG_binaryApprox = true; % SET: use binary approximation (only binary approximation wroks currently)
        numNT = 5; % specify the number of NT that 'can' occur in the provided fasta files. Default value is 5 (ACGT-)

        FLAG_MarkAccessibility = false; % KEEP this FALSE for the time being, Accessibility code needs to be checked

        FLAG_SaveIntCovMtx = false; % SET: will save Integrated Covariance matrix (for debugging only)
        FLAG_useFreqEntry = true;
        FLAG_troubleShoot = true; % SET: saves SelEstNoMu and SelEstSLNoMu
        FLAG_Epi = false; % SET: use MPL with epistasis, UNSET: MPL with epistasis not used

        textCell{1} = ['dirNamesSet' num2str(thisSet) '_ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_'];
        %--------------------------------------------------------------------------


        % ------------------------- AUTO INITIALIZATION ---------------------------
        % NO USER INPUT REQUIRED
        if(setConvention == 1)
            FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = false;
        elseif(setConvention == 2)
            FLAG_stratonovich = true;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = false;
        elseif(setConvention == 3)
            FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = true;
        end

        % this file will contain the names of .fasta files to analyze suing the
        % AnalysisMPL_shortRead code. This file will be generated autotomatically
        % in preprocessingStep1. Here we just need to specify the name. 
        fileNameFastaFilesWithHU = 'fastaFilesHU.txt'; 
        meFastaFileToAnalyze = 'fastaFilesToAnalyze.txt'; % right now, these fasta files need to be generated on laptop

        FLAG_firstSeqIsRef = true; % set: 1st sequence of every fasta file is reference sequence
        mainDir = pwd;
        if(ispc)
            chosenSlash = '\';
        elseif(isunix)
            chosenSlash = '/';
        else
            disp('Error: system is not unix and not PC...')
            pause
        end
        dirNameTemp123 = 'dirNameFiles';
        dirNameStr1Files = [mainDir chosenSlash 'Data_Misc' chosenSlash dirNameTemp123 chosenSlash];


        fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, textCell);
        numPat = length(fileNamesListThisDir);
        %--------------------------------------------------------------------------
        if(numPat == 0)
            disp('NumPat = 0. Check initialization settings and run again.')
        end


        % ========================== BEGIN PROCESSING =============================

        for pat = 1:numPat

            fileNameContainingDirPath = [dirNameStr1Files fileNamesListThisDir{pat}];
            indOfDash = strfind(fileNameContainingDirPath, '_');
            indOfDot = strfind(fileNameContainingDirPath, '.');
            patID = fileNameContainingDirPath(indOfDash(end-1)+1:indOfDash(end)-1);
            thisProt = fileNameContainingDirPath(indOfDash(end)+1:indOfDot(end)-1);

            disp('-----------------------------------------------------------------')
            disp(' ')
            disp(['Patient: ' patID])
            disp(['Protein: ' thisProt])
            FLAG_Skip = false;
            if(strcmp(patID, 'p3') && strcmp(thisProt, 'p6'))
                FLAG_Skip = true;
            end
            if(FLAG_Skip == false)
                %analysisStep1_v2(fileNameContainingDirPath, priorConstSC, FLAG_stratonovich, FLAG_MarkAccessibility, FLAG_UserProvidedRefSeq, FLAG_SaveIntCovMtx, FLAG_useFreqEntry, FLAG_troubleShoot, FLAG_linearInt);
                priorConst = priorConstSC;
                FLAG_vector = [FLAG_stratonovich;
                               FLAG_MarkAccessibility;
                               FLAG_UserProvidedRefSeq;
                               FLAG_SaveIntCovMtx;
                               FLAG_useFreqEntry;
                               FLAG_troubleShoot;
                               FLAG_linearInt];
                % no need to specify  FLAG_Epi for MPL without epistasis analyss
                analysisStep1_v2(fileNameContainingDirPath, priorConst, FLAG_vector, referenceSequence);
            end
        end
    end
end

%end