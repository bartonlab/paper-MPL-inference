%% Analysis of GT data using Ferrer's ApproxWF method
%
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. run ApproxWF method  and save estimate
% Last updated 26-Nov 2017
% Author: M Saqib Sohail

clc
clear all
close all
% set 1990, set 1991 31:100,

thisSet = 'medium_simple';%'medium_complex';
numItr = 1;
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

timeWholeCodeFerrer = zeros(1, numItr);

for thisItr = 1:numItr
    tic
%%  1. Load data   
    
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
    
    [dirNameData, dirNameAnalysis, ~, dirNameApproxWF] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData thisSet chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];

    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end
    

    dirNameDataFerrer = [dirNameData 'Ferrer' chosenSlash];
    dirNameAnalysisFerrer = [dirNameAnalysis 'Ferrer' chosenSlash];
    if(exist(dirNameAnalysisFerrer, 'dir') == 0)
        mkdir(dirNameAnalysisFerrer)        
    end
    
    display('Loading file...') 
    
    fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];

    %load([dirNameAnalysis fileNameForMu], 'muVal');

    fileNameFerrer = ['Ferrer_' fileName(1:end-4) '.txt'];

    dirNameDataFerrerLinux = checkDirName(dirNameDataFerrer);
    dirNameAnalysisFerrerLinux = checkDirName(dirNameAnalysisFerrer);
    
    mainDir = pwd;

    %% 2. run ApproxWF method and save estimate
    % copy the WFABC source files into the data directory
    sourceDirName = dirNameApproxWF;%'/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Ferrer/Approxwf/';
    sourceDirNameLinux = checkDirName(sourceDirName);
    dataDirNameLinux = checkDirName(dirNameDataFerrer);
    
    
    destDirNameLinux = dirNameAnalysisFerrerLinux;
    ApproxWFFile1 = 'ApproxWF';
    if(thisItr == 1)
        commandCopyFiles1 = ['cp ' sourceDirNameLinux ApproxWFFile1 ' ' destDirNameLinux];
        [statusCopyFiles1, cmdoutCopyFiles1] = system(commandCopyFiles1);
        disp(cmdoutCopyFiles1)
    end

    disp('Running Ferrer''s method...')
    cd(dirNameAnalysisFerrer)
    commandRunApproxWF = ['./ApproxWF task=estimate h=0.5 N=' num2str(Nin) ' mutRate=' num2str(muVal) ' MCMClength=10000' ' loci=' dataDirNameLinux fileNameFerrer ' outName=' fileNameFerrer(1:end-4) '_output.txt verbose'];
    timeToStartFerrer(thisItr) = toc;
    [statusRunApproxWF, cmdoutRunApproxWF] = system([commandRunApproxWF]);
    timeFerrer(thisItr) = toc;
    disp(cmdoutRunApproxWF)
        
    
 %%
    cd(mainDir)
    
    FerrerDataTable = readtable([dirNameAnalysisFerrer fileNameFerrer(1:end-4) '_output.txt']);
    [maxLL, maxLLInd] = max(FerrerDataTable{:,2});
    sigmaEstOutFerrer = FerrerDataTable{maxLLInd,3:end}';

    timeWholeCodeFerrer(thisItr) = toc;
    timeWholeCodeFerrer(thisItr)
%     timeWholeCode(thisItr) = toc;
%     timeWholeCode(thisItr)
    save([dirNameAnalysisFerrer fileNameFerrer(1:end-4) '_output.mat'])%, ...
        %'sigmaEstOutFerrer', 'timeWholeCodeFerrer')

    if(thisItr == numItr)
        commandRemoveFile1 = ['rm ' destDirNameLinux ApproxWFFile1];
        [statusRemoveFile1, cmdoutRemoveFile1] = system(commandRemoveFile1);
    end
end
