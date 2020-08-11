clc
clear all
close all
warning off

%perSiteSelctionSelc is loaded from a heterogeneous samplied data set


%========================== INITIALIZATION ================================

% -------------------------- User specified -------------------------------
% dataSet
setAll = 9860001%[8550001 9550001 9650001 9750001 9850001 9860001 9950001];
strainsAll = 5%[1 5 5 5 5 5 5];
%Tused = 300;
tjSamplingSchemeStr = 'scheme1';
allConventions = 3%[1 2 3];
allConvLen = length(allConventions);

%--------------- 
% selcSites based on:
setAll_unEven = setAll;
strainsAll_unEven = strainsAll;
%Tused_unEven = Tused;
ngSamplingSchemeStr_unEven = 'schemeD';
dTSamplingSchemeStr_unEven = 'scheme33';
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
[dirNameDataTemp, dirNameAnalysisTemp, dirNameResultsTemp] = setDirNamesMPLPipeline('Set_dirNames_MPL_SimData1.txt');

ResultFolder = [dirNameResultsTemp 'Set' num2str(setAll_unEven) '/ng_' ngSamplingSchemeStr_unEven '_dT_' dTSamplingSchemeStr_unEven '_Tused_' tjSamplingSchemeStr '_initStr' num2str(strainsAll_unEven) '/'];
thisSetConvSelcSitesFile = ['SelcSites.txt'];
%-------------------

for ss = 1:length(setAll)
    thisSet = setAll(ss)

    numStrainsInInitialPop = strainsAll(ss);
    
    selcSitesAll = dlmread([ResultFolder thisSetConvSelcSitesFile]);

    for c = 1:length(allConventions)

        % chose convention 1: Ito, 2: Stratonovich, 3: Linear interpolation
        % Stratonovich has increased robustness to sampling effects than Ito
        % Linear interpolation has most increased robustness to sampling effects
        setConvention = allConventions(c);

        priorConstSC = 5; % this is the strength of the regularization term

        dTStep = 1;
        ng = 1000;


        textCell{1} = ['dirNamesSet' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dTStep) '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_' ];
        ResultFolder = [dirNameResultsTemp 'Set' num2str(thisSet) '/ng' num2str(ng) '_dT' num2str(dTStep) '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];

        FLAG_SaveIntovMtx = false; % SET: will save Integrated Covariance matrix (for debugging only)

        FLAG_firstSeqIsRef = true; % set: 1st sequence of every fasta file is reference sequence
        FLAG_Epi = false; % SET: use MPL with epistasis, UNSET: MPL with epistasis not used

        % ---------- NO user input required for followng initializations ----------

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

        % Use stratonovich for increased robustness to sampling effects
        if(FLAG_stratonovich == true)
            convention = 'Stratonovich';
        elseif(FLAG_stratonovich == false && FLAG_linearInt == false)
            convention = 'Ito';
        elseif(FLAG_stratonovich == false && FLAG_linearInt == true)
            convention = 'LinearInter';
        end

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

        synNonSynThisProt = [];
        selEstMPLThisProt = [];
        selEstSLThisProt = [];

        selEstMPLSynAll = [];
        selEstMPLNonSynAll = [];
        selEstSLSynAll = [];
        selEstSLNonSynAll = [];

        aucLinkItr = zeros(numPat,2);
        aucLinkNoMuItr = zeros(numPat,2);
        aucUnLinkItr = zeros(numPat,2);
        aucUnLinkNoMuItr = zeros(numPat,2);
        aucLinkLinUniItr = zeros(numPat,2);
        %====================== END INITIALIZATION ================================

        for pat = 1:numPat

            fileNameContainingDirPath = fileNamesListThisDir{pat};
            indOfDash = strfind(fileNameContainingDirPath, '_');
            indOfDot = strfind(fileNameContainingDirPath, '.');
            patID = fileNameContainingDirPath(indOfDash(end-1)+1:indOfDash(end)-1);
            thisProt = fileNameContainingDirPath(indOfDash(end)+1:indOfDot(end)-1);
            fileNameContainingDirPath = [dirNameStr1Files fileNamesListThisDir{pat}];

            disp('-----------------------------------------------------------------')
            disp(' ')
            disp(['Patient: ' patID])
            disp(['Protein: ' thisProt])
            FLAG_Skip = false;
            if(FLAG_Skip == false)
                [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);

                fileNameSynNonSyn = [patID '_' thisProt '_SynNonSyn.txt'];
                synNonSynAll = dlmread([dirNameAnalysis 'Analysis_Misc' chosenSlash fileNameSynNonSyn]);


                fileNameSelEst = ['SelEst_' patID '_' thisProt '_' convention '_gamma' num2str(priorConstSC) '.txt'];
                selEstMPLAll = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEst]);

                fileNameSelEstNoMu = ['SelEstNoMu_' patID '_' thisProt '_' convention '_gamma' num2str(priorConstSC) '.txt'];
                selEstMPLNoMuAll = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstNoMu]);

                fileNameSelEstSL = ['SelEstSL_' patID '_' thisProt '_' convention '_gamma' num2str(priorConstSC) '.txt'];
                selEstSLAll = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSL]);

                fileNameSelEstSLNoMu = ['SelEstSLNoMu_' patID '_' thisProt '_' convention '_gamma' num2str(priorConstSC) '.txt'];
                selEstSLNoMuAll = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameSelEstSLNoMu]);

                fileNameAllTrajs = ['AllTrajs_' patID '_' thisProt '.txt'];
                allTrajs = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajs]);
                polySitesIndicator = (sum(allTrajs > 0 & allTrajs < 1) > 0);


                fileNameIntCovMtx = ['IntCovMtx_' patID '_' thisProt '_' convention '_gamma' num2str(priorConstSC) '.txt'];
                %intCovMtx = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameIntCovMtx]);


                %---------- this only use syn and non syn ---------------------
                synNonSyn = synNonSynAll((synNonSynAll == 0 | synNonSynAll == 1) & polySitesIndicator == 1);
                selEstMPL = selEstMPLAll((synNonSynAll == 0 | synNonSynAll == 1) & polySitesIndicator == 1);
                selEstSL = selEstSLAll((synNonSynAll == 0 | synNonSynAll == 1) & polySitesIndicator == 1);

                selEstMPLSyn = selEstMPLAll((synNonSynAll == 0) & polySitesIndicator == 1);
                selEstMPLNonSyn = selEstMPLAll((synNonSynAll == 1) & polySitesIndicator == 1);
                selEstSLSyn = selEstSLAll((synNonSynAll == 0) & polySitesIndicator == 1);
                selEstSLNonSyn = selEstSLAll((synNonSynAll == 1) & polySitesIndicator == 1);

                selEstMPLSynAll = [selEstMPLSynAll; selEstMPLSyn];
                selEstMPLNonSynAll = [selEstMPLNonSynAll; selEstMPLNonSyn];
                selEstSLSynAll = [selEstSLSynAll; selEstSLSyn];
                selEstSLNonSynAll = [selEstSLNonSynAll; selEstSLNonSyn];


        %--------------------------------------------------------------------------
              
                newDir = [dirNameDataTemp 'Set' num2str(thisSet) '/ng' num2str(ng) '_dT' num2str(dTStep) '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/MC' num2str(pat) '/synth/perSiteSelc/'];
               
                perSiteSelction = dlmread([newDir 'perSiteSelction.txt']);

        %--------------------------------------------------------------------------

                posOnlyItrTemp = perSiteSelction > 0 ;%== max(perSiteSelction);
                negOnlyItrTemp = perSiteSelction < 0; %== min(perSiteSelction);
                
                %selcSites = polySitesIndicator;
                selcSites = selcSitesAll(pat,:);
                selcSites = logical(selcSites);
                
                perSiteSelctionSelc = perSiteSelction(selcSites);
                posOnlySelcItrTemp = perSiteSelctionSelc > 0;% == max(perSiteSelctionSelc);
                negOnlySelcItrTemp = perSiteSelctionSelc < 0;% == min(perSiteSelctionSelc);

                sigmaEstOutLink = selEstMPLAll;
                sigmaEstOutLinkNoMu = selEstMPLNoMuAll;%zeros(size(selEstMPLAll, 1), size(selEstMPLAll, 2));
                sigmaEstOutUnLink = selEstSLAll;
                sigmaEstOutUnLinkNoMu = selEstSLNoMuAll;%zeros(size(selEstMPLAll, 1), size(selEstMPLAll, 2));
                
                sigmaEstOutLinkSelc = sigmaEstOutLink(selcSites);
                sigmaEstOutLinkNoMuSelc = sigmaEstOutLinkNoMu(selcSites);
                sigmaEstOutUnLinkSelc = sigmaEstOutUnLink(selcSites);
                sigmaEstOutUnLinkNoMuSelc = sigmaEstOutUnLinkNoMu(selcSites);

%                 [~,~,~, aucLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink, 1);
%                 [~,~,~, aucLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink, 1);
% 
%                 [~,~,~, aucLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLinkNoMu, 1);
%                 [~,~,~, aucLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLinkNoMu, 1);
% 
%                 [~,~,~, aucUnLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLink, 1);
%                 [~,~,~, aucUnLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLink, 1);
% 
%                 [~,~,~, aucUnLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);
%                 [~,~,~, aucUnLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);

                [aucLinkItrTemp, aucLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLink, sigmaEstOutLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
                [aucUnLinkItrTemp, aucUnLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutUnLink, sigmaEstOutUnLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
                [aucLinkNoMuItrTemp, aucLinkNoMuSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLinkNoMu, sigmaEstOutLinkNoMuSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
                [aucUnLinkNoMuItrTemp, aucUnLinkNoMuSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutUnLinkNoMu, sigmaEstOutUnLinkNoMuSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);

                aucLinkItr(pat,:) = aucLinkItrTemp;
                aucLinkNoMuItr(pat,:) = aucLinkNoMuItrTemp;
                aucUnLinkItr(pat,:) = aucUnLinkItrTemp;
                aucUnLinkNoMuItr(pat,:) = aucUnLinkNoMuItrTemp;
                
                aucLinkSelcItr(pat,:) = aucLinkSelcItrTemp;
                aucLinkNoMuSelcItr(pat,:) = aucLinkNoMuSelcItrTemp;
                aucUnLinkSelcItr(pat,:) = aucUnLinkSelcItrTemp;
                aucUnLinkNoMuSelcItr(pat,:) = aucUnLinkNoMuSelcItrTemp;

            else
                disp('....Skipping this patient protein combination...')
            end
        end
        
        disp('-------------------------------')
        disp(' mean AUROC ')
        allAuc = [mean(aucLinkItr(aucLinkItr(:,1) ~= -1, 1)) mean(aucLinkItr(aucLinkItr(:,2) ~= -1, 2));
                  mean(aucLinkNoMuItr(aucLinkNoMuItr(:,1) ~= -1, 1)) mean(aucLinkNoMuItr(aucLinkNoMuItr(:,2) ~= -1, 2));
                  mean(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1)) mean(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuItr(aucUnLinkNoMuItr(:,1) ~= -1, 1)) mean(aucUnLinkNoMuItr(aucUnLinkNoMuItr(:,2) ~= -1, 2))]
        disp('-------------------------------')
        
        allAucSelc = [mean(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1, 1)) mean(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1, 2));
                  mean(aucLinkNoMuSelcItr(aucLinkNoMuSelcItr(:,1) ~= -1, 1)) mean(aucLinkNoMuSelcItr(aucLinkNoMuSelcItr(:,2) ~= -1, 2));
                  mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1, 1)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuSelcItr(aucUnLinkNoMuSelcItr(:,1) ~= -1, 1)) mean(aucUnLinkNoMuSelcItr(aucUnLinkNoMuSelcItr(:,2) ~= -1, 2))]

        % thisResultFile is same as before, diffSelcSites only effects
        % selcSites variable
        thisResultFile = ['AUROC_results_priorConstSC' num2str(priorConstSC) '_' convention '.txt'];
        thisResultSelcFile = ['AUROC_Selc_results_priorConstSC' num2str(priorConstSC) '_' convention '_selcBasedOn_ngScheme_' ngSamplingSchemeStr_unEven '_dtScheme_' dTSamplingSchemeStr_unEven '.txt'];

        if(exist(ResultFolder, 'dir') == 0)
            mkdir(ResultFolder)
        end

        if(exist([ResultFolder thisResultFile], 'file') == 2)
            delete([ResultFolder thisResultFile])
        end
        if(exist([ResultFolder thisResultSelcFile], 'file') == 2)
            delete([ResultFolder thisResultSelcFile])
        end
        
        dlmwrite([ResultFolder thisResultFile], [aucLinkItr aucLinkNoMuItr aucUnLinkItr aucUnLinkNoMuItr])
        dlmwrite([ResultFolder thisResultSelcFile], [aucLinkSelcItr aucLinkNoMuSelcItr aucUnLinkSelcItr aucUnLinkNoMuSelcItr])

    end
end



