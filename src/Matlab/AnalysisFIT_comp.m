%% Analysis of GT data using Feder's FIT method
%
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates 1-point frequencies
%      3. Classify sites using Feder's method and estimate ML estimate and save
% Last updated 26-Nov 2017
% Author: M Saqib Sohail



% 
clc
clear all
close all

saveFile = 1%1;
thisSet = 'medium_simple';%'medium_complex';
numItr = 1%90;
useSameT = 0; % flag that controls if 
              %     1: usable T will be loaded from data
              %     0: supplied by user in variable Tused
loadDataOption = 2;
% 1 : load data from full trajectories information
% 2 : load data from sampled trajectories information (*.dat file)
fileNameContainingDirPath = 'dirNames.txt';
getSysParam_4;

if(useSameT == 1) % use the T thats in the data file, i.e., use all generations
else 
    T = Tused + Tstart;
end

timeWholeCodeFeder = zeros(1, numItr);

%%
for thisItr = 1:numItr
    tic

%%  1. Load data   
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
        
    [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData thisSet chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];
    dirNameAnalysisFeder = [dirNameAnalysis 'Feder' chosenSlash];

    if(exist(dirNameAnalysisFeder, 'dir') == 0)
        mkdir(dirNameAnalysisFeder)        
    end

    disp('Loading file...')
    if(loadDataOption == 1)
        % load data from full trajectories
    elseif(loadDataOption == 2)
        % load data from sampled trajectories
        fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];
        dataIn = dlmread([dirNameData fileName]);
    end

    warning off
%% 2. calculates 1-point frequencies
%--------------------------------------------------------------------------
    if(loadDataOption == 1)
        % Reconstruct MSA from raw data files, perform finite sampling by 
        % selecting ng individuals from the MSA and find sampled 1 and 2 point
        % frequencies
        numGen = (T - Tstart)/dT; % number of generation in the whole msa

        % later can make a switch to chose less than N samples too to simulate what
        % happens in FLU samples 
        numSamplesPerGen = Nin;
        numSamplesPerGenSelected = ng;

        samplingTimePoints = Tstart:dT:T;
        numSamplingPoints = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);

        randSelctIndAll = zeros(numSamplingPoints, numSamplesPerGenSelected);
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
            fileNameRandPerm = ['RandPerm_N' num2str(Nin) '_t' num2str(t)];
            load([dirNameRandPermFiles fileNameRandPerm],'randPermN')
            temp4051 = randPermN;
            randSelctInd = temp4051(1:numSamplesPerGenSelected);
            randSelctIndAll(t, :) = randSelctInd;

            thisMSA = thisMSATemp(randSelctInd,:);
            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
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
        end
        disp('done')
    end
    

%% 3. Classify sites using Feder's method and estimate ML estimate and save
    
    sigmaEstFeder = -1*ones(Lin,1);
    
    meanYi = zeros(1, Lin);
    Ssq = zeros(1, Lin);
    degFree = zeros(1, Lin);
    for l = 1:Lin
        yNumer = q(2:end,l) - q(1:end-1,l);
        yDenomTemp = sqrt(2*dT*q(:,l).*(1 - q(:,l)));
        yDenom = yDenomTemp(1:end-1);
        YiThisSite = yNumer./yDenom;
        validInd = ~isnan(YiThisSite) & YiThisSite ~= Inf & YiThisSite ~=-Inf;
        meanYi(l) = sum(YiThisSite(validInd))/sum(validInd);
        Ssq(l) = sum((YiThisSite(validInd) - meanYi(l)).^2)/(sum(validInd)-1);
        degFree(l) = sum(validInd) - 1;
    end
    
    t_FI = meanYi./sqrt(Ssq/numSamplingPoints);
    
    % 95% confidence interval
    CIup = 0.0975;
    CIlow = 0.025;
    tDistUpLimitSigLev = tinv(CIup,degFree);
    tDistlowLimitSigLev = tinv(CIlow,degFree);
    
    
%     t_FI > tDistUpLimitSigLev;
%     t_FI < tDistlowLimitSigLev;
    
    disp('done.')
    
    model = 'linkDiff';
    
    if(loadDataOption == 1)        
        if(Tstart == 1)
            fileNameSave = ['Feder_' fileName(1:end-4) '_Tend' num2str(T) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        else
             fileNameSave = ['Feder_' fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(T) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        end
    elseif(loadDataOption == 2)
        % in this case, only change the file extension
        fileNameSave = ['Feder_' fileName(1:end-4) '_output.mat'];
    end
    
    timeWholeCodeFeder(thisItr) = toc;
    timeWholeCodeFeder(thisItr)
    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysisFeder fileNameSave])
    end
end
