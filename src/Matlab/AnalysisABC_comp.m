%% Analysis of GT data using Foll's WFABC_1 aad WFABC_2 programs
%
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. run WFABC_1 and WFABC_2 
%      3. calculate the ML estimate and save
% Last updated 26-Nov 2017
% Author: M Saqib Sohail


clc
clear all
close all


thisSet = 'medium_simple';%'medium_complex';
numItr = 1%100;
useSameT = 0; % flag that controls if 
              %     1: usable T will be loaded from data
              %     0: supplied by user in variable Tuse
fileNameContainingDirPath = 'dirNames.txt';
getSysParam_4;

actualT = T/1000;


if(useSameT == 1) % use the T thats in the data file, i.e., use all generations
else 
    T = Tused + Tstart;
end

timeToStartWFABC1 = zeros(1, numItr);
timeWFABC1 = zeros(1, numItr);
timeWFABC2 = zeros(1, numItr);
timeWholeCodeFoll = zeros(1, numItr);

for thisItr = 1:numItr
    tic
%%  1. Load data (in .mat or .dat format)
    
    thisItr
    dirNameScriptFile = pwd;    
    if(ispc)
        chosenSlash = '\';
    elseif(isunix)
        chosenSlash = '/';
    else
        display('Error: system si not unix and not PC...')
        pause
    end
    
    [dirNameData, dirNameAnalysis, dirNameABC] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData thisSet chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];

    
    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end

    dirNameDataFoll = [dirNameData 'Foll' chosenSlash];
    dirNameAnalysisFoll = [dirNameAnalysis 'Foll' chosenSlash];
    if(exist(dirNameAnalysisFoll, 'dir') == 0)
        mkdir(dirNameAnalysisFoll)        
    end
    
    display('Loading file...') 
    
    fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];
    
    fileNameFoll = ['Foll_' fileName(1:end-4) '.txt'];
    fileNameFollFlipVec = [fileNameFoll(1:end-4) '_flip.dat'];
    flipVec = dlmread([dirNameDataFoll fileNameFollFlipVec]);
    dirNameDataFollLinux = checkDirName(dirNameDataFoll);
    dirNameAnalysisFollLinux = checkDirName(dirNameAnalysisFoll);
    
    mainDir = pwd;
    
    % copy the WFABC source files into the data directory. From the data
    % directory, we will move the analysis files to the analysis directory.
    if(thisItr == 1)
        sourceDirName = dirNameABC;%'/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/WFABC_v1.1/WFABC_v1.1/binaries/Linux/';
        sourceDirNameLinux = checkDirName(sourceDirName);

        destDirNameLinux = checkDirName(dirNameDataFoll);
        wfabcFile1 = 'wfabc_1';
        wfabcFile2 = 'wfabc_2';
        commandCopyFiles1 = ['cp ' sourceDirNameLinux wfabcFile1 ' ' destDirNameLinux];
        [statusCopyFiles1, cmdoutCopyFiles1] = system(commandCopyFiles1);
        disp(cmdoutCopyFiles1)

        commandCopyFiles2 = ['cp ' sourceDirNameLinux chosenSlash wfabcFile2 ' ' destDirNameLinux];
        [statusCopyFiles2, cmdoutCopyFiles2] = system(commandCopyFiles2);
        disp(cmdoutCopyFiles2)
    end
    
    cd(dirNameDataFoll)

    %% 2. run WFABC_1 and WFABC_2 
    
    % run WFABC_1 command of Foll
    command1 = ['./wfabc_1 ' fileNameFoll];
    timeToStartWFABC1(thisItr) = toc;
    [status1, cmdout1] = system(command1);
    timeWFABC1(thisItr) = toc;
    disp(cmdout1)
    
    
    % run WFABC_2 command of Foll
    command2 = ['./wfabc_2 -fixed_N ' num2str(Nin) ' -ploidy 1 -min_s -1 -max_s 0.5 ' fileNameFoll];
    
    [status2, cmdout2] = system(command2);
    timeWFABC2(thisItr) = toc;
    disp(cmdout2)
  
    [fileNameFoll(1:end-4) '_posterior_s.txt'];
    
    % copy all 4 analysis files to the analysis folder
    if(exist([fileNameFoll(1:end-4) '_Ne_bootstrap.txt'], 'file') == 2)
        commandCopyAnalysisFile1 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_Ne_bootstrap.txt ' dirNameAnalysisFollLinux];
        [statusCopyAnalysisFile1, cmdoutCopyAnalysisFile1] = system(commandCopyAnalysisFile1);
        disp(cmdoutCopyAnalysisFile1)
    end
    
    if(exist([fileNameFoll(1:end-4) '_obs_stats.txt'], 'file') == 2)
        commandCopyAnalysisFile2 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_obs_stats.txt ' dirNameAnalysisFollLinux];
        [statusCopyAnalysisFile2, cmdoutCopyAnalysisFile2] = system(commandCopyAnalysisFile2);
        disp(cmdoutCopyAnalysisFile2)
    end
    
    if(exist([fileNameFoll(1:end-4) '_posterior_s.txt'], 'file') == 2)
        commandCopyAnalysisFile3 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_posterior_s.txt ' dirNameAnalysisFollLinux];
        [statusCopyAnalysisFile3, cmdoutCopyAnalysisFile3] = system(commandCopyAnalysisFile3);
        disp(cmdoutCopyAnalysisFile3)
    end
    
    if(exist([fileNameFoll(1:end-4) '_posterior_N.txt'], 'file') == 2)
        commandCopyAnalysisFile4 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_posterior_N.txt ' dirNameAnalysisFollLinux];
        [statusCopyAnalysisFile4, cmdoutCopyAnalysisFile4] = system(commandCopyAnalysisFile4);
        disp(cmdoutCopyAnalysisFile4)
    end
    
    cd(mainDir)

    posterior_s = dlmread([dirNameAnalysisFoll fileNameFoll(1:end-4) '_posterior_s.txt']);
    temp2000 = repmat(flipVec', 1, size(posterior_s, 2));
    posterior_s = posterior_s.*temp2000;
    sigmaEstOutFoll_mean = mean(posterior_s,2)';
    sigmaEstOutFoll_median = median(posterior_s,2)';
    
    %% 3. calculate the ML estimate and save
    % find ML of the numerical posterior distribution
    leftEdge = -1;
    rightEdge = 0.5;
    numBins = 100;
    step = (rightEdge-leftEdge)/numBins;
    edges = leftEdge:step:rightEdge;
    
    sigmaEstOutFoll_ML = zeros(Lin,1);
    
    for ll = 1:Lin
        histValues = histcounts(posterior_s(ll,:), edges);
        [valllMax, indllMax] = max(histValues);
        sigmaEstOutFoll_ML(ll) = edges(indllMax) + step/2;
    end

    timeWholeCodeFoll(thisItr) = toc;
    timeWholeCodeFoll(thisItr)
        save([dirNameAnalysisFoll fileNameFoll(1:end-4) '_posterior_s.mat'])%, ...
            %'posterior_s', 'sigmaEstOutFoll_mean', 'sigmaEstOutFoll_median', 'sigmaEstOutFoll_ML', 'timeWholeCodeFoll')
    if(thisItr == numItr)
        commandRemoveFile1 = ['rm ' destDirNameLinux wfabcFile1];
        [statusRemoveFile1, cmdoutRemoveFile1] = system(commandRemoveFile1);
        commandRemoveFile2 = ['rm ' destDirNameLinux wfabcFile2];
        [statusRemoveFile2, cmdoutRemoveFile2] = system(commandRemoveFile2);
    end
end
%%
disp('Time required to run the code')
[timeToStartWFABC1;
timeWFABC1;
timeWFABC2;
timeWholeCodeFoll]
