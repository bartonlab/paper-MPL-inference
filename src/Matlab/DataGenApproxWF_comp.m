%% Converts GT test data to format readable by Ferrer's ApproxWF code
%
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates 1-point frequencies
%      3. Filter and convert data to Ferrer's ApproxWF format
% Last updated 29-Nov 2017
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


plotFigs = 0;
L = Lin;

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
    
    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
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
%% 2. extract and and two point frequencies
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
    
    %% 3. Filter and convert data to Ferrer's ApproxWF format
    fprintf('Converting data to Ferrers format...')
    
    qOrig = q;
    noiseThreshold = 0.05; % the frequency cutoff for conidering a trajectory to be not sampling noise
    disp(['Filtering trajectories with a sampling noise threshold of ' num2str(noiseThreshold) '...'])
    FiltTrajFoll;
    
    for l = 1:L
        % for the lth site, 
        %    i. select trajs that go from 0 to 0 only if
        %             highest is > 0.05
        %    ii. select trajs that go from 0 to 1 
        %    iii. select trajs that go from 1 to 1 only if
        %             lowest is < 0.95
        %    iv. select trajs that go from 1 to 0 
        % 
        thisSiteAllTrajs = allSitesAllTrajs(allSitesAllTrajs(:,1) == l,:);
        minVal = min(thisSiteAllTrajs(:,4:end)')';
        maxVal = max(thisSiteAllTrajs(:,4:end)')';
        thisSiteAllTrajs = [ thisSiteAllTrajs thisSiteAllTrajs(:,1) == l & (abs(maxVal - minVal) > noiseThreshold)];
        
        for tt = 1:size(thisSiteAllTrajs,1)
            if(min(thisSiteAllTrajs(tt,4:5)) == 0 && thisSiteAllTrajs(tt,end) == 0)
                q(thisSiteAllTrajs(tt,2):thisSiteAllTrajs(tt,3),l) = 0;
            elseif(min(thisSiteAllTrajs(tt,4:5)) == 1 && thisSiteAllTrajs(tt,end) == 0)
                q(thisSiteAllTrajs(tt,2):thisSiteAllTrajs(tt,3),l) = 1;
            end
        end
        
        if(sum(q(:,l)) == 0)
            q(:,l) = qOrig(:,l);
        end
    end
    disp('done.')

    % save data
    dirNameDataFerrer = [dirNameData 'Ferrer' chosenSlash];
    if(exist(dirNameDataFerrer, 'dir') == 0)
        mkdir(dirNameDataFerrer)
    end
        
    if(loadDataOption == 1)        
    elseif(loadDataOption == 2)
        % in this case, only change the file extension
        fileNameSave = fileName(1:end-4);
        fileNameFerrer = ['Ferrer_' fileNameSave '.txt'];
    end

    if(saveFile == 1)
        disp('Saving file in format of Ferrer...') 
        fid = fopen([dirNameDataFerrer fileNameFerrer],'wt');
        
        %fid = fopen('thisTestFile.txt', 'wt');
        
        % save 1st line
        for ll = 1:(L+1)
            if(ll == 1)
                strToSave = 'time';
            else
                strToSave = ['L' num2str(ll-1)];
            end
            
            if(ll ~= L+1)
                formatSpec = '%s\t';
            else
                formatSpec = '%s\n';
            end
            
            fprintf(fid, formatSpec, strToSave);
        end
        
        formatSpec = '%s\n';
        qN = round(q*ng);
        ngStr = num2str(ng);
        for tt = 1:numSamplingPoints
            thisTimeStr = [];
            for ll = 1:(L+1)
                if(ll == 1)
                    thisTimeStr = [num2str(samplingTimePoints(tt))];
                else
                    thisTimeStr = [thisTimeStr ' ' num2str(qN(tt,ll-1)) '/' ngStr];
                end
            end
            fprintf(fid, formatSpec, thisTimeStr);
        end
        fclose(fid)
    end
    toc
end
