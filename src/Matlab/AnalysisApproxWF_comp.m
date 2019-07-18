%% Analysis of GT data using Ferrer's ApproxWF method
%
%%
% Run the script in MATLAB (press F5). The code will prompt the user for
% any required action. The output of this code is estimates of selection
% coefficients stored in the file:
% ApproxWF_medium_simple_old_collected_extended_Tend1000.csv 

% this code :
%      1. loads data (in .mat or .dat format)
%      2. run ApproxWF method  and save estimate
% Written: 26-Nov 2017
% Author: M Saqib Sohail

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
disp(' This code will run ApproxWF script on ground truth data. ')
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

timeWholeCodeFerrer = zeros(1, numItr);
errorExit = 1;
for thisItr = 1:numItr
    tic
%%  1. Load data   
    
    thisItr
    dirNameScriptFile = pwd;    
    if(ispc)
        chosenSlash = '\';
        disp('Error: This code is not compatible with a Windows based system. Try running on a Linux based system.')
        errorExit = 1;
        break
    elseif(isunix)
        chosenSlash = '/';
    else
        display('Error: system si not unix and not PC...')
        pause
    end
    
    makeDirNames();
    [dirNameData, dirNameAnalysis, ~, dirNameApproxWF] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];

    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end
    

    dirNameDataFerrer = [dirNameData 'Ferrer' chosenSlash];
    dirNameAnalysisFerrer = [dirNameAnalysis 'Ferrer' chosenSlash];
    if(exist(dirNameAnalysisFerrer, 'dir') == 0)
        mkdir(dirNameAnalysisFerrer)        
    end
    
    % check if output dat file already exists
    if(thisItr == 1)
        dirNameSaveCSVFile = [];
        if(strcmp(thisSet, 'medium_simple'))
            fileNamePaperDataCompApproxWF = ['ApproxWF_medium_simple_old_collected_extended_Tend' num2str(T) '.csv'];
        elseif(strcmp(thisSet, 'medium_complex'))
            fileNamePaperDataCompApproxWF = ['ApproxWF_medium_complex_old_collected_extended_Tend' num2str(T) '.csv'];
        end

        if(exist([dirNameSaveCSVFile fileNamePaperDataCompApproxWF], 'file') == 2)
            disp(['The following file already exists ' dirNameSaveCSVFile fileNamePaperDataCompApproxWF])
            prompt = 'Overwrite? (y) yes, (n) no. Your choince: ';
            repeatInput2 = 1;
            while(repeatInput2 == 1)
                str = input(prompt,'s');
                if(length(str) == 1 && str == 'y')
                    repeatInput2 = 0;

                elseif(length(str) == 1 && str == 'n')
                    repeatInput2 = 0;
                    disp('WFApproxWF not run. Rename existing file and run this code again.')
                    break
                else
                    disp('Unexpected input, please pren ''y'' or ''n''. Overwrite? (y) yes, (n) no. Your choince: ');
                    prompt = ' ';
                end
            end        
        end
    end
    
    display('Loading file...') 
    
    fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];

    %load([dirNameAnalysis fileNameForMu], 'muVal');

    fileNameFerrer = ['Ferrer_' fileName(1:end-4) '.txt'];

    dirNameDataFerrerLinux = checkDirName(dirNameDataFerrer);
    dirNameAnalysisFerrerLinux = checkDirName(dirNameAnalysisFerrer);
    
    mainDir = pwd;

    %% 2. run ApproxWF method and save estimate
    % copy the WFApproxWF source files into the data directory
    sourceDirName = dirNameApproxWF;
    sourceDirNameLinux = checkDirName(sourceDirName);
    dataDirNameLinux = checkDirName(dirNameDataFerrer);
    
    if(exist(dirNameDataFerrer, 'dir') ~= 7)
        disp('Error: data file not in Format compatible with ApproxWF. Run DataGenApproxWF_comp.m before running this code.')
        break
    end
    
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
    
    
    
    if(exist([dirNameAnalysisFerrer fileNameFerrer(1:end-4) '_output.txt'], 'file') ~= 2)
        disp('Error: data file not found.')
        disp(' ')
        disp('       Possible cause:')
        disp('       1. Verify the ApproxWF package is properly installed and functioning (see help files')
        disp('          provided with the ApproxWF package')
        disp('       2. Note that the ApproxWF code needs to be compiled before first use by running the')
        disp('          following command as explained in the installation guide of ApproxWF package')
        disp('                   g++ -std=c++11 -O3 -o ApproxWF *.cpp -fopenmp -lz ')
        disp('       3. Verify names of data and ApproxWF directories are correctly specified in dirNames.txt.')
        disp('          The Approxwf package should be in directory .../paper-MPL-inference-master/src/external/Approxwf/')
        disp('       4. Run DataGenApproxWF_comp.m prior to running this Analysis script.')
        break
    end
    
    FerrerDataTable = readtable([dirNameAnalysisFerrer fileNameFerrer(1:end-4) '_output.txt']);
    [maxLL, maxLLInd] = max(FerrerDataTable{:,2});
    sigmaEstOutFerrer = FerrerDataTable{maxLLInd,3:end}';

    timeWholeCodeFerrer(thisItr) = toc;
    
    save([dirNameAnalysisFerrer fileNameFerrer(1:end-4) '_output.mat'])%, ...

    if(thisItr == numItr)
        commandRemoveFile1 = ['rm ' destDirNameLinux ApproxWFFile1];
        [statusRemoveFile1, cmdoutRemoveFile1] = system(commandRemoveFile1);
    end
    
    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);
    [~,~,~, aucFerrerItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFerrer, 1);
    [~,~,~, aucFerrerItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFerrer, 1);
    aucFerrerItr(thisItr,:) = aucFerrerItrTemp;

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
        fid = fopen(fileNamePaperDataCompApproxWF,'wt');
        fprintf(fid, str);
        fclose(fid);
    end

    inpDataApproxWF = [num2str(thisItr-1) ',' num2str(thisItr-1) ',ApproxWF,' num2str(Tstart) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeWholeCodeFerrer(thisItr))];

    for l = 1:Lin
        inpDataApproxWF = [inpDataApproxWF ',' num2str(sigmaEstOutFerrer(l))];
    end
    inpDataApproxWF = [inpDataApproxWF ',' num2str(aucFerrerItrTemp(1)) ',' num2str(aucFerrerItrTemp(2))];
    for l = 1:Lin
        inpDataApproxWF = [inpDataApproxWF ',' num2str(perSiteSelction(l) - sigmaEstOutFerrer(l))];
    end
    inpDataApproxWF = [inpDataApproxWF '\n'];
    fid = fopen(fileNamePaperDataCompApproxWF,'a');
    fprintf(fid, inpDataApproxWF);
    fclose(fid);
    errorExit = 0;
end
if(errorExit == 1)
else
    disp(['Output saved to ...' chosenSlash 'Matlab' chosenSlash  fileNamePaperDataCompApproxWF])
end