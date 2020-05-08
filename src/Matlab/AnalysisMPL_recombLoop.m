%% Run Maximum Path Likelihood (MPLin) method on GT data


% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates one and two point frequencies
%      3. Find estimates using MPL

% Last updated 31-Oct 2017
% Author: M Saqib Sohail

clc

clear all
close all

saveFile = 1;
%L = 50; 
allSets = [755 75500001 7550001 755001];
%thisSet = 755001;%35500001;%755001;%355001;%561101%55%5617011%315;%67;%52;%67;%5611%561701%5991%315%5617%7%5613%315%56%5611001%5601%5991%53900%6510;%6500%1990;%1991%1034;%7%58%7;
numItr = 100;%90;
allTused = 1000;%[300 500 700 1000 1300];
for kkw = 1:length(allSets)
    
thisSet = allSets(kkw);

useSameT = 0; % flag that controls if 
              %     1: usable T will be loaded from data
              %     0: supplied by user in variable Tused
% Tused = 1000%2500;%2500;%10000;
% Tstart = 1%500;%1500%1001;
loadDataOption = 1;
% 1 : load data from full trajectories information
% 2 : load data from sampled trajectories information (*.dat file)

priorConst = 1;% 1/10;%1/10;%1/N^2;

runLinearInt = 0; % 1: run linear interpolation code, 0: skip it
calcLikelihood = 0;
fileNameContainingDirPath = 'dirNames_recomb.txt';
getSysParam_long;

%Tstart = 31;
dT = 10;
ng = 100;
Tused = 1000;%allTused(kkw);
Tstart = 1;


if(useSameT == 1) % use the T thats in the data file, i.e., use all generations
else 
    %Tend = Tused + Tstart - 1;
    Tend = Tused + Tstart;
end

plotFigs = 0;

timeLoadDataMPL = zeros(1, numItr);
timeProcDataMPL = zeros(1, numItr);
timeRunAlgoMPL = zeros(1, numItr);
timeWholeCodeMPL = zeros(1, numItr);
%
for thisItr = 1:numItr%1:numItr
    tic

%%  1. load data (in .mat or .dat format)
%--------------------------------------------------------------------------

    thisItr
    mainDir = pwd;    
    if(ispc)
        chosenSlash = '\';
    elseif(isunix)
        chosenSlash = '/';
    else
        disp('Error: system si not unix and not PC...')
        pause
    end
    if(classesOfSites == 2)
        selTypeName = 'PosSel';
    elseif(classesOfSites == 3)
        selTypeName = 'PosDelSel';
    end

    [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    
    
    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end
    
    disp('Loading file...')
    if(loadDataOption == 1)
        % load data from full trajectories
        if(recombination == 0)
            fileName = ['WFsim_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '.mat'];
        elseif(recombination == 1)
            fileName = ['WFsim_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '.mat'];
        end
        load([dirNameData fileName], 'N', 'D', 'masterStrainList', 'masterTimeFreq', 'perSiteSelction', 'sitesUnderSelection', 'muVal');

        masterStrainList = masterStrainList(:,1:Lin);
    elseif(loadDataOption == 2)
        % load data from sampled trajectories
        if(Tstart == 1)
            fileName = ['WFsim_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
        else
            fileName = ['WFsim_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
        end
        dataIn = dlmread([dirNameData fileName]);
    end
    
    timeLoadDataMPL(thisItr) = toc;
    warning off
%% 2. calculates one and two point frequencies
%--------------------------------------------------------------------------
    if(loadDataOption == 1)
        % Reconstruct MSA from raw data files, perform finite sampling by 
        % selecting ng individuals from the MSA and find sampled 1 and 2 point
        % frequencies
        %numGen = (Tend - Tstart + 1)/dT; % number of generation in the whole msa
        numGen = (Tend - Tstart)/dT; % number of generation in the whole msa

        % later can make a switch to chose less than N samples too to simulate what
        % happens in FLU samples 
        numSamplesPerGen = N; 
        numSamplesPerGenSelected = ng;

        samplingTimePoints = Tstart:dT:Tend;
        numSamplingPoints = length(samplingTimePoints);
        numGen = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);
        q11 = -1*ones(Lin, Lin, numSamplingPoints);

        %randSelctIndAll = zeros(numSamplingPoints, numSamplesPerGenSelected);
        % this section calculates q_t for each site, for the synthetic protein, 
        % 2 allele per site case this is just the sum of occurances of ones
        fprintf('Calculate 1 and 2 point probabilites...')
        for t = 1:numSamplingPoints

            thisSamplingTimePoint = samplingTimePoints(t);
            strainsThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,1);
            freqStrainThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,2);
            thisMSATemp = -1*ones(numSamplesPerGen, Lin);
            count1 = 0;
            for k = 1:length(strainsThisTimePoint)
                count1 = count1 + 1;
                thisMSATemp(count1:count1+freqStrainThisTimePoint(k)-1,:) = repmat(masterStrainList(strainsThisTimePoint(k),:), freqStrainThisTimePoint(k), 1);
                count1 = count1 + freqStrainThisTimePoint(k) - 1;
            end
            fileNameRandPerm = ['RandPerm_N' num2str(N) '_t' num2str(t)];
            load([dirNameRandPermFiles fileNameRandPerm],'randPermN')
            temp4051 = randPermN;
            randSelctInd = temp4051(1:numSamplesPerGenSelected);
            %randSelctIndAll(t, :) = randSelctInd;

            thisMSA = thisMSATemp(randSelctInd,:);
            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
            thisMSALogical = logical(thisMSA);
            for l = 1:Lin
                % multiply the lth column of MSA with all MSA to find q11 for (l,:)
                % pairs of cov matrix
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);

                % sum for the lth row of the covariance matrix
                q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;
            end
        end

        disp('done')

        clear qijAtTimeTk;
        clear masterStrainFreq;
    elseif(loadDataOption == 2)
        AllMSATimeVec = dataIn(:,1) + 1;
        AllMSAFreqIn = dataIn(:,2);
        AllMSAIn = dataIn(:,3:end);
        samplingTimePoints = unique(dataIn(:,1) + 1)';
        numSamplingPoints = length(samplingTimePoints);
        numSamplesPerGen = Nin;
        numSamplesPerGenSelected = ng;
        numGen = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);
        q11 = -1*ones(Lin, Lin, numSamplingPoints);
        
        % this section calculates q_t for each site, for the synthetic protein, 
        % 2 allele per site case this is just the sum of occurances of ones
        fprintf('Calculate 1 and 2 point probabilites...')
        for t = 1:numSamplingPoints

            thisSamplingTimePoint = samplingTimePoints(t);
            thisTimePointSelcRows = AllMSATimeVec == thisSamplingTimePoint;
            strainsThisTimePoint = AllMSAIn(thisTimePointSelcRows,:);
            freqStrainThisTimePoint = AllMSAFreqIn(thisTimePointSelcRows);
            thisMSA = -1*ones(numSamplesPerGenSelected, Lin);
            count1 = 0;
            for k = 1:size(strainsThisTimePoint,1)
                count1 = count1 + 1;
                thisMSA(count1:count1+freqStrainThisTimePoint(k)-1,:) = repmat(strainsThisTimePoint(k,:), freqStrainThisTimePoint(k), 1);
                count1 = count1 + freqStrainThisTimePoint(k) - 1;
            end

            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
            thisMSALogical = logical(thisMSA);
            for l = 1:Lin
                % multiply the lth column of MSA with all MSA to find q11 for (l,:)
                % pairs of cov matrix
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);
                % sum for the lth row of the covariance matrix
                q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;
            end
        end
        disp('done')
    end
    timeProcDataMPL(thisItr) = toc;
    % plot tranjectory at each site
    if(plotFigs)
        for l = 1:Lin
           figure(l)
           subplot(2,1,1)
           plot(q(:,l), 'b.-')
           xlabel('time points')
           ylabel('Frequency Count')
           title(['Site number: ' num2str(l)])
        end
    end
    
    

%% 3. Find estimates using MPL
%     vectors with mu

    
    fprintf('Calculating selection coefficients estimates...')
    model = 'linkDiff';
    
    

    startTime = 1;
    stopTime = numGen;
    
    sumAijWithLink = zeros(Lin,Lin);
    vecCovEntries = zeros((Lin*(Lin-1)/2), numSamplingPoints);
    for t = 1:numSamplingPoints%startTime:stopTime-1
        tempDiagEntries = q(t,:).*(1 - q(t,:));
        temp200 = zeros(Lin,Lin);
        for thisSite = 1:Lin
            % remember q11 is just the 2 point, i.e., q_ij, we still
            % need ot subtract q_i*q_j from it
            % this makes q_i*q_i q_i*q_j ... q_i*q_k ( will remove
            % diagonal later)
            temp100 = q(t,:).*q(t,thisSite);

            % except for (i,i) entry, rest are q_ij - q_i*q_j
            temp200(thisSite,:) = q11(thisSite, :, t) - temp100;
        end
        covMtxThisTime = ~(eye(Lin)).*temp200 + diag(tempDiagEntries);
        sumAijWithLink = sumAijWithLink + covMtxThisTime;
        polySitesThisTP = (q(t,:) > 0 & q(t,:) < 1);
        [avgLD, allPairInVec] = getAvgCov(covMtxThisTime(polySitesThisTP,polySitesThisTP));
        avgLDPerTime(t) = avgLD;
        temp1 = triu(covMtxThisTime, 1);
        temp2 = temp1(temp1 ~= 0);
        req_len = Lin*(Lin-1)/2;
        temp2_len = length(temp2);
        numZeroInsrt = req_len - temp2_len;
        if(numZeroInsrt > 0)
            temp2 = [temp2; zeros(numZeroInsrt, 1)];
        end
        vecCovEntries(:,t) = temp2;
    end
    
    tempTimeSumCovProcessig = toc;

    % 6.1.2.1 calc estimates of selection and errors
    sumAijUnLink = diag(diag(sumAijWithLink));


    v_EstOutLink = dT*muVal*(sum((ones(Lin,stopTime-startTime+1-1) - 2*q(startTime:stopTime-1,:)'),2));
    
    sigmaEstOutLink = (sumAijWithLink*dT + priorConst*eye(Lin))\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
    sigmaEstOutLinkNoMu = (sumAijWithLink*dT + priorConst*eye(Lin))\(q(stopTime,:) - q(startTime,:))';
    sigmaEstOutUnLink = (sumAijUnLink*dT + priorConst*eye(Lin))\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
    sigmaEstOutUnLinkNoMu = (sumAijUnLink*dT + priorConst*eye(Lin))\(q(stopTime,:) - q(startTime,:))';
    sigmaEstOutLink_LinUni = zeros(Lin,1);
    
    if(calcLikelihood == 1)
        likelihoodLink = 0;
        likelihoodUnLink = 0;
        likelihoodLinkZero = 0;
        likelihoodUnLinkZero = 0;
        likelihoodLink50 = 0;
        likelihoodUnLink50 = 0;
        likelihoodLink100 = 0;
        likelihoodUnLink100 = 0;
        likelihoodLinkRand = 0;
        likelihoodLink50Rand = 0;
        likelihoodLink100Rand = 0;
        
        temp9090 = randperm(Lin);
        sigmaEstOutLinkRand = sigmaEstOutLink(temp9090);
        selcSites50 = ~(abs(max(q) - min(q)) <= 0.05 & (max(q) == 1 | min(q) == 0))'; % ~ rejection criterion
        selcSites100 = ~(abs(max(q) - min(q)) <= 0.1 & (max(q) == 1 | min(q) == 0))'; % ~ rejection criterion
        for t = 1:numSamplingPoints%startTime:stopTime-1
            tempDiagEntries = q(t,:).*(1 - q(t,:));
            temp200 = zeros(Lin,Lin);
            for thisSite = 1:Lin
                % remember q11 is just the 2 point, i.e., q_ij, we still
                % need ot subtract q_i*q_j from it
                % this makes q_i*q_i q_i*q_j ... q_i*q_k ( will remove
                % diagonal later)
                temp100 = q(t,:).*q(t,thisSite);
                % except for (i,i) entry, rest are q_ij - q_i*q_j
                temp200(thisSite,:) = q11(thisSite, :, t) - temp100;
            end
            covMtxThisTime = ~(eye(Lin)).*temp200 + diag(tempDiagEntries);
           %sitesUsedForLLCalc = q(t,:)' > 10/Nin & q(t,:)' < 1 - 10/Nin;
            sitesUsedForLLCalc = q(t,:)' > 0 & q(t,:)' < 1;
            sitesUsedForLLCalc50 = sitesUsedForLLCalc & selcSites50;
            sitesUsedForLLCalc100 = sitesUsedForLLCalc & selcSites100;
            
            temp10 = (q(t+1,sitesUsedForLLCalc)' - q(t,sitesUsedForLLCalc)' - dT*(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)*sigmaEstOutLink(sitesUsedForLLCalc) - muVal*(1 - 2*q(t,sitesUsedForLLCalc)')));
            temp20 = covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc) + 0.000001*eye(sum(sitesUsedForLLCalc));
            temp30 = temp10'*(temp20\temp10)/dT; % no need for dT as it is constant and can be neglected
            likelihoodLink = likelihoodLink + sum(temp30) + log(1/sqrt(det(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)+ 0.000001*eye(sum(sitesUsedForLLCalc)))));
            temp10z = (q(t+1,sitesUsedForLLCalc)' - q(t,sitesUsedForLLCalc)' - dT*(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)*zeros(sum(sitesUsedForLLCalc),1) - muVal*(1 - 2*q(t,sitesUsedForLLCalc)')));
            temp20z = covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc) + 0.000001*eye(sum(sitesUsedForLLCalc));
            temp30z = temp10z'*(temp20z\temp10z)/dT; % no need for dT as it is constant and can be neglected
            likelihoodLinkZero = likelihoodLinkZero + sum(temp30z) + log(1/sqrt(det(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)+ 0.000001*eye(sum(sitesUsedForLLCalc)))));
             
            temp11 = (q(t+1,sitesUsedForLLCalc)' - q(t,sitesUsedForLLCalc)' - dT*(diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)))*sigmaEstOutUnLink(sitesUsedForLLCalc) - muVal*(1 - 2*q(t,sitesUsedForLLCalc)')));
            temp21 = diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc) + 0.000001*eye(sum(sitesUsedForLLCalc))));
            temp31 = temp11'*(temp21\temp11)/dT; % no need for dT as it is constant and can be neglected
            likelihoodUnLink = likelihoodUnLink + sum(temp31) + log(1/sqrt(det(diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)+ 0.000001*eye(sum(sitesUsedForLLCalc)))))));
            temp11z = (q(t+1,sitesUsedForLLCalc)' - q(t,sitesUsedForLLCalc)' - dT*(diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)))*zeros(sum(sitesUsedForLLCalc),1) - muVal*(1 - 2*q(t,sitesUsedForLLCalc)')));
            temp21z = diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc) + 0.000001*eye(sum(sitesUsedForLLCalc))));
            temp31z = temp11z'*(temp21z\temp11z)/dT; % no need for dT as it is constant and can be neglected
            likelihoodUnLinkZero = likelihoodUnLinkZero + sum(temp31z) + log(1/sqrt(det(diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)+ 0.000001*eye(sum(sitesUsedForLLCalc)))))));
             
            
%             
%             temp2020 = q(t,sitesUsedForLLCalc)' + covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)*sigmaEstOutLink(sitesUsedForLLCalc) + muVal*(1 - 2*q(t,sitesUsedForLLCalc)');
%             likelihoodLink = likelihoodLink + sum(N*q(t,sitesUsedForLLCalc)'.*log(temp2020));
%             temp2020z = q(t,sitesUsedForLLCalc)' + covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)*zeros(sum(sitesUsedForLLCalc),1) + muVal*(1 - 2*q(t,sitesUsedForLLCalc)');
%             likelihoodLinkZero = likelihoodLinkZero + sum(N*q(t,sitesUsedForLLCalc)'.*log(temp2020z));

%             temp2020_50 = q(t,sitesUsedForLLCalc50)' + covMtxThisTime(sitesUsedForLLCalc50,sitesUsedForLLCalc50)*sigmaEstOutLink(sitesUsedForLLCalc50) + muVal*(1 - 2*q(t,sitesUsedForLLCalc50)');
%             likelihoodLink50 = likelihoodLink50 + sum(N*q(t,sitesUsedForLLCalc50)'.*log(temp2020_50));
%             temp2020_100 = q(t,sitesUsedForLLCalc100)' + covMtxThisTime(sitesUsedForLLCalc100,sitesUsedForLLCalc100)*sigmaEstOutLink(sitesUsedForLLCalc100) + muVal*(1 - 2*q(t,sitesUsedForLLCalc100)');
%             likelihoodLink100 = likelihoodLink100 + sum(N*q(t,sitesUsedForLLCalc100)'.*log(temp2020_100));
%             
%             temp2021 = q(t,sitesUsedForLLCalc)' + diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)))*sigmaEstOutUnLink(sitesUsedForLLCalc) + muVal*(1 - 2*q(t,sitesUsedForLLCalc)');
%             likelihoodUnLink = likelihoodUnLink + sum(N*q(t,sitesUsedForLLCalc)'.*log(temp2021));
%             temp2021z = q(t,sitesUsedForLLCalc)' + diag(diag(covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)))*zeros(sum(sitesUsedForLLCalc),1) + muVal*(1 - 2*q(t,sitesUsedForLLCalc)');
%             likelihoodUnLinkZero = likelihoodUnLinkZero + sum(N*q(t,sitesUsedForLLCalc)'.*log(temp2021z));
%             temp2021_50 = q(t,sitesUsedForLLCalc50)' + diag(diag(covMtxThisTime(sitesUsedForLLCalc50,sitesUsedForLLCalc50)))*sigmaEstOutUnLink(sitesUsedForLLCalc50) + muVal*(1 - 2*q(t,sitesUsedForLLCalc50)');
%             likelihoodUnLink50 = likelihoodUnLink50 + sum(N*q(t,sitesUsedForLLCalc50)'.*log(temp2021_50));
%             temp2021_100 = q(t,sitesUsedForLLCalc100)' + diag(diag(covMtxThisTime(sitesUsedForLLCalc100,sitesUsedForLLCalc100)))*sigmaEstOutUnLink(sitesUsedForLLCalc100) + muVal*(1 - 2*q(t,sitesUsedForLLCalc100)');
%             likelihoodUnLink100 = likelihoodUnLink100 + sum(N*q(t,sitesUsedForLLCalc100)'.*log(temp2021_100));
%             
%             temp2000 = q(t,sitesUsedForLLCalc)' + covMtxThisTime(sitesUsedForLLCalc,sitesUsedForLLCalc)*sigmaEstOutLinkRand(sitesUsedForLLCalc) + muVal*(1 - 2*q(t,sitesUsedForLLCalc)');
%             likelihoodLinkRand = likelihoodLinkRand + sum(N*q(t,sitesUsedForLLCalc)'.*log(temp2000));
%             temp2000_50 = q(t,sitesUsedForLLCalc50)' + covMtxThisTime(sitesUsedForLLCalc50,sitesUsedForLLCalc50)*sigmaEstOutLinkRand(sitesUsedForLLCalc50) + muVal*(1 - 2*q(t,sitesUsedForLLCalc50)');
%             likelihoodLink50Rand = likelihoodLink50Rand + sum(N*q(t,sitesUsedForLLCalc50)'.*log(temp2000_50));
%             temp2000_100 = q(t,sitesUsedForLLCalc100)' + covMtxThisTime(sitesUsedForLLCalc100,sitesUsedForLLCalc100)*sigmaEstOutLinkRand(sitesUsedForLLCalc100) + muVal*(1 - 2*q(t,sitesUsedForLLCalc100)');
%             likelihoodLink100Rand = likelihoodLink100Rand + sum(N*q(t,sitesUsedForLLCalc100)'.*log(temp2000_100));
%             
        end
        % no need for - sign as we put alt hypthesis (Link) first
        deltaLogLikeli = 2*(likelihoodLink - likelihoodUnLink); 
        
        deltaLogLikeliMPL = 2*(likelihoodLink - likelihoodLinkZero); 
        deltaLogLikeliSL = 2*(likelihoodUnLink - likelihoodUnLinkZero); 
        
        
        deltaLogLikeli50 = 2*(likelihoodLink50 - likelihoodUnLink50);
        deltaLogLikeli100 = 2*(likelihoodLink100 - likelihoodUnLink100);
        
        deltaLogLikeliRand = 2*(likelihoodLink - likelihoodLinkRand); 
        deltaLogLikeli50Rand = 2*(likelihoodLink50 - likelihoodLink50Rand);
        deltaLogLikeli100Rand = 2*(likelihoodLink100 - likelihoodLink100Rand);
        
        
    end
    
    q11Temp = q11;
    if(runLinearInt == 1)
        % LinearInter_Nsites_Dec2016;
        LinearInter_Nsites_Jan2017_1;

        startTime_LinUni = 1;
        stopTime_LinUni = size(q_LinUni,1);%dTByDTNew*stopTime;

        sumAijWithLink_LinUni = zeros(Lin,Lin);
        
        for t_LinUni = startTime_LinUni:stopTime_LinUni-1
            tempDiagEntries_LinUni = q_LinUni(t_LinUni,:).*(1 - q_LinUni(t_LinUni,:));
            temp200_LinUni = zeros(Lin,Lin);
            for thisSite = 1:Lin
                % remember q11 is just the 2 point, i.e., q_ij, we still
                % need ot subtract q_i*q_j from it
                % this makes q_i*q_i q_i*q_j ... q_i*q_k ( will remove
                % diagonal later)
                temp100_LinUni = q_LinUni(t_LinUni,:).*q_LinUni(t_LinUni,thisSite);

                % except for (i,i) entry, rest are q_ij - q_i*q_j
                temp200_LinUni(thisSite,:) = q11_LinUni(thisSite, :, t_LinUni) - temp100_LinUni;
            end
            covMtxThisTime_LinUni = ~(eye(Lin)).*temp200_LinUni + diag(tempDiagEntries_LinUni);
            sumAijWithLink_LinUni = sumAijWithLink_LinUni + covMtxThisTime_LinUni;
        end
        tempTimeSumCovProcessig_LinUni = toc;


        v_EstOutLink_LinUni = muVal*(sum((ones(Lin,stopTime_LinUni-startTime_LinUni+1-1) - 2*q_LinUni(startTime_LinUni:stopTime_LinUni-1,:)'),2));
        sigmaEstOutLink_LinUni = (sumAijWithLink_LinUni*dT_LinUni + priorConst*eye(Lin))\(q_LinUni(stopTime_LinUni,:) - q_LinUni(startTime_LinUni,:) - v_EstOutLink_LinUni')';

    end    
    disp('done.')
    timeRunAlgoMPL(thisItr) = toc;

%%
%     thisRun = 1;
    if(loadDataOption == 1)        
        if(Tstart == 1)
            fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        else
             fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        end
    elseif(loadDataOption == 2)
        % in this case, only change the file extension
        fileNameSave = [fileName(1:end-4) '.mat'];
    end

%     fileExist = exist([dirNameAnalysis fileNameSave], 'file');
%     while(fileExist == 2)
%         thisRun = thisRun + 1;
%         fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model  '_filter_' filterStr '_' num2str(thisRun) '_M3' fileName(end-3:end)];
%         fileExist = exist([dirNameAnalysis fileNameSave], 'file');
%     end

    timeWholeCodeMPL(thisItr) = toc;
    timeWholeCodeMPL(thisItr)
    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysis fileNameSave])
    end
    

% deltaLogLikeli
% deltaLogLikeli50
% deltaLogLikeli100
% pause(0.5)
% [likelihoodLink likelihoodUnLink likelihoodLinkRand]
%[deltaLogLikeli deltaLogLikeli50 deltaLogLikeli100]
%[deltaLogLikeliRand deltaLogLikeli50Rand deltaLogLikeli100Rand]
%  pause(5)
%[likelihoodLink likelihoodLinkZero likelihoodLink-likelihoodLinkZero  likelihoodUnLink likelihoodUnLinkZero likelihoodUnLink-likelihoodUnLinkZero]
 %[likelihoodLink likelihoodUnLink likelihoodLink-likelihoodUnLink]
 
% pause(1)
end
end


Save_AnalysisData_forPlot;