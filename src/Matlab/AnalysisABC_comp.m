%% Analysis of GT data using Foll's WFABC_1 aad WFABC_2 programs
%
%%
% Run the script in MATLAB (press F5). The code will prompt the user for
% any required action. The output of this code is estimates of selection
% coefficients stored in the file:
% ABC_medium_simple_old_collected_extended_Tend1000.csv 

% this code :
%      1. loads data in .dat format
%      2. run WFABC_1 and WFABC_2 
%      3. calculate the ML estimate and save
% Last updated 26-Nov 2017
% Author: M Saqib Sohail
%
% Last updated: 16-July 2019
%               -automated assignment of various variable
%               -estimates saved to .csv file
%%
clc
clear all
close all


repeatInput = 1;
disp('-------------------------------------------------------------------')
disp(' ')
disp(' This code will run WFABCv1.1 script on ground truth data. ')
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
errorExit = 1;
for thisItr = 1:numItr
    tic
%%  1. Load data in .mat or .dat format)
    
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
    
    makeDirNames();
    [dirNameData, dirNameAnalysis, dirNameABC] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];

    
    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end

    dirNameDataFoll = [dirNameData 'Foll' chosenSlash];
    dirNameAnalysisFoll = [dirNameAnalysis 'Foll' chosenSlash];
    if(exist(dirNameAnalysisFoll, 'dir') == 0)
        mkdir(dirNameAnalysisFoll)        
    end
    
    % check if output .csv file already exists
    if(thisItr == 1)
        dirNameSaveCSVFile = [];
        if(strcmp(thisSet, 'medium_simple'))
            fileNamePaperDataCompABC = ['ABC_medium_simple_old_collected_extended_Tend' num2str(T) '.csv'];
        elseif(strcmp(thisSet, 'medium_complex'))
            fileNamePaperDataCompABC = ['ABC_medium_complex_old_collected_extended_Tend' num2str(T) '.csv'];
        end

        if(exist([dirNameSaveCSVFile fileNamePaperDataCompABC], 'file') == 2)
            disp(['The following file already exists ' dirNameSaveCSVFile fileNamePaperDataCompABC])
            prompt = 'Overwrite? (y) yes, (n) no. Your choince: ';
            repeatInput2 = 1;
            while(repeatInput2 == 1)
                str = input(prompt,'s');
                if(length(str) == 1 && str == 'y')
                    repeatInput2 = 0;

                elseif(length(str) == 1 && str == 'n')
                    repeatInput2 = 0;
                    disp('WFABC not run. Rename existing file and run this code again.')
                    break
                else
                    disp('Unexpected input, please pren ''y'' or ''n''. Overwrite? (y) yes, (n) no. Your choince: ');
                    prompt = ' ';
                end
            end        
        end
    end
    
    % load data to process
    disp('Loading file...') 
    
    fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];
    
    fileNameFoll = ['Foll_' fileName(1:end-4) '.txt'];
    fileNameFollFlipVec = [fileNameFoll(1:end-4) '_flip.dat'];
    if(exist([dirNameDataFoll fileNameFollFlipVec], 'file') ~= 2)
        disp('Error: data file not found.')
        disp(' ')
        disp('Possible cause: Run DataGenWFABC_comp.m prior to running this Analysis script.')
        break
    end
    flipVec = dlmread([dirNameDataFoll fileNameFollFlipVec]);
    dirNameDataFollLinux = checkDirName(dirNameDataFoll);
    dirNameAnalysisFollLinux = checkDirName(dirNameAnalysisFoll);
    
    mainDir = pwd;
    if(exist(dirNameABC, 'dir') ~= 7)
        disp('Download WFABC from http://jjensenlab.org/software and instal it in directory')
        disp('.../paper-MPL-inference-master/src/external/')
        disp('Verify the installation directory name matches the one given in file dirNames.txt line 12')
        break
    end
    
    % copy the WFABC source files into the data directory. From the data
    % directory, we will move the analysis files to the analysis directory.
    if(thisItr == 1)
        sourceDirName = dirNameABC;
        sourceDirNameLinux = checkDirName(sourceDirName);
        destDirNameLinux = checkDirName(dirNameDataFoll);
        
        sourceDirNameWindows_hasSpace = checkDirNameWindows(sourceDirName);
        destDirNameWindows_hasSpace = checkDirNameWindows(dirNameDataFoll);
        
        if(ispc)
            if(sourceDirNameWindows_hasSpace == true || destDirNameWindows_hasSpace == true)
                disp('Error: Directory names can not have space in them. Check:')
                disp(sourceDirName)
                disp(dirNameDataFoll)
                break
            end
        end
        
        wfabcFile1 = 'wfabc_1';
        wfabcFile2 = 'wfabc_2';
        if(ispc)
            commandCopyFiles1 = ['copy ' sourceDirName wfabcFile1 '.exe ' dirNameDataFoll(1:end-1)];
        else
            commandCopyFiles1 = ['cp ' sourceDirNameLinux wfabcFile1 ' ' destDirNameLinux];
        end
        
        [statusCopyFiles1, cmdoutCopyFiles1] = system(commandCopyFiles1);
        disp(cmdoutCopyFiles1)

        if(ispc)
            commandCopyFiles2 = ['copy ' sourceDirName wfabcFile2 '.exe ' dirNameDataFoll];
        else
            commandCopyFiles2 = ['cp ' sourceDirNameLinux chosenSlash wfabcFile2 ' ' destDirNameLinux];
        end
        [statusCopyFiles2, cmdoutCopyFiles2] = system(commandCopyFiles2);
        disp(cmdoutCopyFiles2)
    end
    
    cd(dirNameDataFoll)

    %% 2. run WFABC_1 and WFABC_2 
    
    % run WFABC_1 command of Foll
    if(ispc)
        command1 = ['wfabc_1 ' fileNameFoll];
    else
        command1 = ['./wfabc_1 ' fileNameFoll];
    end
    
    timeToStartWFABC1(thisItr) = toc;
    [status1, cmdout1] = system(command1);
    timeWFABC1(thisItr) = toc;
    disp(cmdout1)
    
    
    % run WFABC_2 command of Foll
    if(ispc)
        command2 = ['wfabc_2 -fixed_N ' num2str(Nin) ' -ploidy 1 -min_s -1 -max_s 0.5 ' fileNameFoll];
    else
        command2 = ['./wfabc_2 -fixed_N ' num2str(Nin) ' -ploidy 1 -min_s -1 -max_s 0.5 ' fileNameFoll];
    end
    
    
    [status2, cmdout2] = system(command2);
    timeWFABC2(thisItr) = toc;
    disp(cmdout2)
  
    [fileNameFoll(1:end-4) '_posterior_s.txt'];
    
    % copy all 4 analysis files to the analysis folder
    if(exist([fileNameFoll(1:end-4) '_Ne_bootstrap.txt'], 'file') == 2)
        if(ispc)
            commandCopyAnalysisFile1 = ['move ' dirNameDataFoll fileNameFoll(1:end-4) '_Ne_bootstrap.txt ' dirNameAnalysisFoll(1:end-1)];
        else
            commandCopyAnalysisFile1 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_Ne_bootstrap.txt ' dirNameAnalysisFollLinux];
        end
        [statusCopyAnalysisFile1, cmdoutCopyAnalysisFile1] = system(commandCopyAnalysisFile1);
        disp(cmdoutCopyAnalysisFile1)
    end
    
    if(exist([fileNameFoll(1:end-4) '_obs_stats.txt'], 'file') == 2)
        if(ispc)
            commandCopyAnalysisFile2 = ['move ' dirNameDataFoll fileNameFoll(1:end-4) '_obs_stats.txt ' dirNameAnalysisFoll(1:end-1)];
        else
            commandCopyAnalysisFile2 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_obs_stats.txt ' dirNameAnalysisFollLinux];
        end
        [statusCopyAnalysisFile2, cmdoutCopyAnalysisFile2] = system(commandCopyAnalysisFile2);
        disp(cmdoutCopyAnalysisFile2)
    end
    
    if(exist([fileNameFoll(1:end-4) '_posterior_s.txt'], 'file') == 2)
        if(ispc)
            commandCopyAnalysisFile3 = ['move ' dirNameDataFoll fileNameFoll(1:end-4) '_posterior_s.txt ' dirNameAnalysisFoll(1:end-1)];
        else
            commandCopyAnalysisFile3 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_posterior_s.txt ' dirNameAnalysisFollLinux];
        end
        [statusCopyAnalysisFile3, cmdoutCopyAnalysisFile3] = system(commandCopyAnalysisFile3);
        disp(cmdoutCopyAnalysisFile3)
    end
    
    if(exist([fileNameFoll(1:end-4) '_posterior_N.txt'], 'file') == 2)
        if(ispc)
            commandCopyAnalysisFile4 = ['move ' dirNameDataFoll fileNameFoll(1:end-4) '_posterior_N.txt ' dirNameAnalysisFoll(1:end-1)];
        else
            commandCopyAnalysisFile4 = ['mv ' dirNameDataFollLinux fileNameFoll(1:end-4) '_posterior_N.txt ' dirNameAnalysisFollLinux];
        end
        [statusCopyAnalysisFile4, cmdoutCopyAnalysisFile4] = system(commandCopyAnalysisFile4);
        disp(cmdoutCopyAnalysisFile4)
    end
    
    cd(mainDir)

    % check if some sites were flipped, adjust sign to accomodate
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
    disp([' Time required to analyze this iteration: ' num2str(timeWholeCodeFoll(thisItr))])
        save([dirNameAnalysisFoll fileNameFoll(1:end-4) '_posterior_s.mat'])%, ...
            %'posterior_s', 'sigmaEstOutFoll_mean', 'sigmaEstOutFoll_median', 'sigmaEstOutFoll_ML', 'timeWholeCodeFoll')
    if(thisItr == numItr)
        if(ispc)
            commandRemoveFile1 = ['delete ' destDirNameLinux wfabcFile1 '.exe'];
        else
            commandRemoveFile1 = ['rm ' destDirNameLinux wfabcFile1];
        end
        
        [statusRemoveFile1, cmdoutRemoveFile1] = system(commandRemoveFile1);
        if(ispc)
            commandRemoveFile2 = ['delete ' destDirNameLinux wfabcFile2 '.exe'];
        else
            commandRemoveFile2 = ['rm ' destDirNameLinux wfabcFile2];
        end
        
        [statusRemoveFile2, cmdoutRemoveFile2] = system(commandRemoveFile2);
    end
    
    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);
    [~,~,~, aucFoll_MLItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFoll_ML, 1);
    [~,~,~, aucFoll_MLItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFoll_ML, 1);

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
        fid = fopen(fileNamePaperDataCompABC,'wt');
        fprintf(fid, str);
        fclose(fid);
    end

    inpDataABC = [num2str(thisItr-1) ',' num2str(thisItr-1) ',ABC,' num2str(Tstart) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeWholeCodeFoll(thisItr))];
    for l = 1:Lin
        inpDataABC = [inpDataABC ',' num2str(sigmaEstOutFoll_ML(l))];
    end
    inpDataABC = [inpDataABC ',' num2str(aucFoll_MLItrTemp(1)) ',' num2str(aucFoll_MLItrTemp(2))];
    for l = 1:Lin
        inpDataABC = [inpDataABC ',' num2str(perSiteSelction(l) - sigmaEstOutFoll_ML(l))];
    end
    inpDataABC = [inpDataABC '\n'];
    
    
    fid = fopen([dirNameSaveCSVFile fileNamePaperDataCompABC],'a');
    fprintf(fid, inpDataABC);
    fclose(fid);
    errorExit = 0;
end
if(errorExit == 1)
else
    disp(['Output saved to ...' chosenSlash 'Matlab' chosenSlash  fileNamePaperDataCompABC])
end
