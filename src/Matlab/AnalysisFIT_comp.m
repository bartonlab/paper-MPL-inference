%% Analysis of GT data using Feder's FIT method
%
%%
% Run the script in MATLAB (press F5). The code will prompt the user for
% any required action. The output of this code is estimates of selection
% coefficients stored in the file:
% FIT_medium_simple_old_collected_extended_Tend1000.csv 

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates 1-point frequencies
%      3. Classify sites using Feder's method and estimate ML estimate and save
% Written 20-Nov 2017
% Author: M Saqib Sohail

% Last updated: 16-July 2019
%               -automated assignment of various variable
%               -estimates saved to .csv file

%% 
clc
clear all
close all

saveFile = 1%1;

repeatInput = 1;
disp('-------------------------------------------------------------------')
disp(' ')
disp(' This code will run the FIT method on ground truth data. ')
disp(' ')
disp('Dataset to analyze: (1) Medium simple, (2) Medium complex');
prompt = ' ';
while(repeatInput == 1)
    str = input(prompt,'s');
    if(length(str) == 1 && str == '1')
        repeatInput = 0;
        thisSet = 'medium_simple';
    elseif(length(str) == 1 && str == '2')
        repeatInput = 0;
        thisSet = 'medium_complex';
    else
        disp('Unexpected input, please pren ''1'' or ''2''. Dataset to analyze: (1) Medium simple, (2) Medium complex');
        prompt = ' ';
    end
end


numItr = 100;
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
        
    makeDirNames();
    [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
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
    % check if output dat file already exists
    if(thisItr == 1)
        dirNameSaveCSVFile = [];
        if(strcmp(thisSet, 'medium_simple'))
            fileNamePaperDataCompFIT = ['FIT_medium_simple_old_collected_extended_Tend' num2str(T) '.csv'];
        elseif(strcmp(thisSet, 'medium_complex'))
            fileNamePaperDataCompFIT = ['FIT_medium_complex_old_collected_extended_Tend' num2str(T) '.csv'];
        end

        if(exist([dirNameSaveCSVFile fileNamePaperDataCompFIT], 'file') == 2)
            disp(['The following file already exists ' dirNameSaveCSVFile fileNamePaperDataCompFIT])
            prompt = 'Overwrite? (y) yes, (n) no. Your choice: ';
            repeatInput2 = 1;
            while(repeatInput2 == 1)
                str = input(prompt,'s');
                if(length(str) == 1 && str == 'y')
                    repeatInput2 = 0;

                elseif(length(str) == 1 && str == 'n')
                    repeatInput2 = 0;
                    disp('FIT not run. Rename existing file and run this code again.')
                    break
                else
                    disp('Unexpected input, please pren ''y'' or ''n''. Overwrite? (y) yes, (n) no. Your choice: ');
                    prompt = ' ';
                end
            end        
        end
    end
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
    
    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysisFeder fileNameSave])
    end
    
    sigmaEstOutFeder = t_FI;
    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);
    [~,~,~, aucFederItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFeder, 1);
    [~,~,~, aucFederItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFeder, 1);
    aucFederItr(thisItr,:) = aucFederItrTemp;

    if(thisItr == 1)
        str = ' ';
        str = [str ',trajectory,method,t0,T,ns,deltat,runtime'];
        for l = 1:Lin
            str = [str ',s' num2str(l-1)];
        end
        str = [str ',AUROC_ben,AUROC_del'];
        for l = 1:Lin
            str = [str ',ds' num2str(l-1)];
        end
        str = [str '\n'];
        fid = fopen(fileNamePaperDataCompFIT,'wt');
        fprintf(fid, str);
        fclose(fid);
    end

    inpDataFIT = [num2str(thisItr-1) ',' num2str(thisItr-1) ',FIT,' num2str(Tstart) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeWholeCodeFeder(thisItr))];
    for l = 1:Lin
        inpDataFIT = [inpDataFIT ',' num2str(sigmaEstOutFeder(l))];
    end
    inpDataFIT = [inpDataFIT ',' num2str(aucFederItrTemp(1)) ',' num2str(aucFederItrTemp(2))];
    for l = 1:Lin
        inpDataFIT = [inpDataFIT ',' num2str(perSiteSelction(l) - sigmaEstOutFeder(l))];
    end
    inpDataFIT = [inpDataFIT '\n'];
    
    fid = fopen(fileNamePaperDataCompFIT,'a');
    fprintf(fid, inpDataFIT);
    fclose(fid);
end
disp(['Output saved to ...' chosenSlash 'Matlab' chosenSlash  fileNamePaperDataCompFIT])