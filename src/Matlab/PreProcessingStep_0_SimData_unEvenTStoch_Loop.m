
% prepares  dirNames_X_Y.txt files needed in PreProcessingStep1

% usage: data folder should have the following structure
%        Data dir:  /.../dir1/[Patient_ID]/[Protein_name]/[FileNames].fasta
%        Analysis Dir should be different from data directiry


clc
clear all
close all

%-------------------- USER CONTROLLED INITIALIZATION ----------------------
FLAG_SaveFile = 1;
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



for ss = 1:length(setAll)
thisSet = setAll(ss);%9950001;%55;%560;%755;

%ngAll = 100%[10 20 50 80 100];
%dTStepAll = 10%[5 10 20 50 100 200];

%for nn = 1:length(ngAll)
%for tt = 1:length(dTStepAll)

%dTStep = dTStepAll(tt);
%ng = ngAll(nn);
numStrainsInInitialPop = strainsAll(ss);%5;%3; % number of strains in the initial population

[dirNameDataTemp, dirNameAnalysisTemp, dirNameResultsTemp] = setDirNamesMPLPipeline('Set_dirNames_MPL_SimData1.txt');

% These variables need to be set for given dataset                   'ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused' num2str(Tused)
dataDirNameMain = [dirNameDataTemp 'Set' num2str(thisSet) '/ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];
%dataDirNameMain = ['/local1/staff/ee/mssohail/Matlab Codes 2/MPL Pipeline/Sim_Data/dT' num2str(dTStep) '/Set' num2str(thisSet) '/'];

analysisDirNameMain = [dirNameAnalysisTemp 'Set' num2str(thisSet) '/ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];
%analysisDirNameMain = ['/local1/staff/ee/mssohail/Matlab Codes 2/MPL Pipeline/Sim_Analysis/dT' num2str(dTStep) '/Set' num2str(thisSet) '/'];


str1 = ['dirNamesSet' num2str(thisSet) '_ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_'];
%str1 = ['dirNamesSet' num2str(thisSet) 'dT' num2str(dTStep) '_'];
%--------------------------------------------------------------------------

% ------------------------- AUTO INITIALIZATION ---------------------------
% NO USER INPUT REQUIRED
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
if(exist(dirNameStr1Files, 'dir') == 0)
    mkdir(dirNameStr1Files)
end

allProts{1} = 'synth';

numProts = length(allProts);

for prot = 1:numProts
    thisProt = allProts{prot};

    dirNamePatList = getFolderContent(dataDirNameMain, 'dir');
    
    numPat = length(dirNamePatList);
    if(numPat == 0)
        disp('NumPat = 0. Check initialization settings and run again.')
    end

    %%
    for pat = 1:numPat
        disp(['Patient: ' dirNamePatList{pat}])
        dirNameDataThisPat = [dataDirNameMain dirNamePatList{pat} chosenSlash thisProt chosenSlash];
        dirNameAnalysisThisPat = [analysisDirNameMain dirNamePatList{pat} chosenSlash thisProt chosenSlash];

        if(exist(dirNameAnalysisThisPat, 'dir') == 0)
            mkdir(dirNameAnalysisThisPat)
        end
        fileNameContainingDirPath = [str1 dirNamePatList{pat} '_' thisProt '.txt'];

        if(FLAG_SaveFile == 1)
            fprintf('Saving .txt file containing names of data and analysis folders...')

            if(exist([dirNameStr1Files fileNameContainingDirPath], 'file') == 2)
                delete([dirNameStr1Files fileNameContainingDirPath])
            end
            fileID = fopen([dirNameStr1Files fileNameContainingDirPath],'w');
            for f = 1:2
               if(f == 1)
                  fprintf(fileID,'%s\n', '% This file specifies directory paths to Data, Analysis');
                  fprintf(fileID,'%s\n', '% the syntax is /dir1/dir2/dir3/ ');
                  fprintf(fileID,'%s\n', '% for windows based system, the code will automatically reformat the path');
                  fprintf(fileID,'%s\n', '%');
                  fprintf(fileID,'%s\n', '% -------------------------------------------------------------------------');
               end
               if(f == 1)
                   fprintf(fileID,'%s\n',['dirNameData=' dirNameDataThisPat]);
               elseif(f == 2)
                   fprintf(fileID,'%s\n',['dirAnalysisName=' dirNameAnalysisThisPat]);
               end
            end
            fclose(fileID);
            disp('done.')
        elseif(FLAG_SaveFile == 0)
            disp('Warning: .txt data file not saved as FLAG_SaveFile flag not set.')
        else
            disp('Error: case undefined')
            pause
        end
    end
end
end
%end
%end