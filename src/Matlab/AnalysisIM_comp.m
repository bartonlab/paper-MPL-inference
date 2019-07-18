%% Implements the linked-diffusion method with linear interpolation for 
% 3 class of sites (deleterious, beneficial, neutral) sims 
%%
% Run the script in MATLAB (press F5). The code will prompt the user for
% any required action. The output of this code is estimates of selection
% coefficients stored in the file:
% Det_medium_simple_old_collected_extended_Tend1000.csv 

%  1. Load data   
%  2. Reconstruct MSA from raw data files, perform finite sampling by 
%     selecting x_g individuals from the MSA and find sampled 1 and 2 point
%     frequencies
%  6. Estimation...SW diff with linkage for Linear interpolated freq 
%     vectors with mu

% Written: 10-Feb 2017
% Author: M Saqib Sohail

% Last updated: 16-July 2019
%               -automated assignment of various variable
%               -estimates saved to .csv file
%%

clc
clear all
close all

debuggingDisplay = 0;
saveFile = 1;

repeatInput = 1;
disp('-------------------------------------------------------------------')
disp(' ')
disp(' This code will analyze ground truth data using Illingworth''s method. ')
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

% optimizeOption = 2; % 1: maxSA runs, 2: tolVal = 0.01*1/StallIterLimIn
initCondTrue = 0; % 3: init with 10% fluct around true, 2: init all zero, 
                  % 1: init true, 0: init random
plotFigs = 0;
SARuns = 100000;
StallIterLimIn = 10000;
numOptRuns = 5; % number of optimization runs
TolFunIn = 0.0001/StallIterLimIn;
dispStr = 'diagnose';%'diagnose';%'off';
fileNameContainingDirPath = 'dirNames.txt';
getSysParam_4;


if(useSameT == 1) % use the T thats in the data file, i.e., use all generations
else 
    Tend = Tused + Tstart;
end
if(classesOfSites == 2)
    selTypeName = 'PosSel';
elseif(classesOfSites == 3)
    selTypeName = 'PosDelSel';
end

timeLoadDataIlling = zeros(1, numItr);
timeProcDataIlling = zeros(1, numItr);
timeRunAlgoIlling = zeros(1, numItr);
timeWholeCodeIlling = zeros(1, numItr);
timeEachSAItrIlling = zeros(numOptRuns, numItr);

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
    [dirNameData, dirNameAnalysis, dirNameDet] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];
    

    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end

    dirNameAnalysisIlling = [dirNameAnalysis 'Illing' chosenSlash];
    if(exist(dirNameAnalysisIlling, 'dir') == 0)
        mkdir(dirNameAnalysisIlling)
    end

    % check if output dat file already exists
    if(thisItr == 1)
        dirNameSaveCSVFile = [];
        if(strcmp(thisSet, 'medium_simple'))
            fileNamePaperDataCompDet = ['Det_medium_simple_old_collected_extended_Tend' num2str(T) '.csv'];
        elseif(strcmp(thisSet, 'medium_complex'))
            fileNamePaperDataCompDet = ['Det_medium_complex_old_collected_extended_Tend' num2str(T) '.csv'];
        end

        if(exist([dirNameSaveCSVFile fileNamePaperDataCompDet], 'file') == 2)
            disp(['The following file already exists ' dirNameSaveCSVFile fileNamePaperDataCompDet])
            prompt = 'Overwrite? (y) yes, (n) no. Your choice: ';
            repeatInput2 = 1;
            while(repeatInput2 == 1)
                str = input(prompt,'s');
                if(length(str) == 1 && str == 'y')
                    repeatInput2 = 0;

                elseif(length(str) == 1 && str == 'n')
                    repeatInput2 = 0;
                    disp('IM not run. Rename existing file and run this code again.')
                    break
                else
                    disp('Unexpected input, please pren ''y'' or ''n''. Overwrite? (y) yes, (n) no. Your choice: ');
                    prompt = ' ';
                end
            end        
        end
    end
    
    
    disp('Loading file...')
    if(loadDataOption == 1)
        % load data from full trajectories
    elseif(loadDataOption == 2)
        % load data from sampled trajectories
        fileName = ['wfsim_' thisSet '_' num2str(thisItr-1) '_T' num2str(Tend) '_ns' num2str(ng) '_dt' num2str(dT) '.dat'];
        dataIn = dlmread([dirNameData fileName]);
    end

    timeLoadDataIlling(thisItr) = toc;
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
        numSamplesPerGen = N; 
        numSamplesPerGenSelected = ng;
        
        samplingTimePoints = Tstart:dT:T;
        numSamplingPoints = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);
        q01 = -1*ones(Lin, Lin, numSamplingPoints);
        q10 = -1*ones(Lin, Lin, numSamplingPoints);
        q11 = -1*ones(Lin, Lin, numSamplingPoints);

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
            fileNameRandPerm = ['RandPerm_N' num2str(N) '_t' num2str(t)];
            load([dirNameRandPermFiles fileNameRandPerm],'randPermN')
            temp4051 = randPermN;
            randSelctInd = temp4051(1:numSamplesPerGenSelected);
            randSelctIndAll(t, :) = randSelctInd;

            thisMSA = thisMSATemp(randSelctInd,:);
            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
            thisMSALogical = logical(thisMSA);
            for l = 1:Lin
                % multiply the lth column of MSA with all MSA to find q11 for (l,:)
                % pairs of cov matrix
                tempq01 = (repmat(~thisMSALogical(:,l),1, Lin).* thisMSALogical);
                tempq10 = (repmat(thisMSALogical(:,l),1, Lin).* ~thisMSALogical);
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);

                % sum for the lth row of the covariance matrix
                q01(l,:,t) = sum(tempq01)./numSamplesPerGenSelected;
                q10(l,:,t) = sum(tempq10)./numSamplesPerGenSelected;
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
        q01 = -1*ones(Lin, Lin, numSamplingPoints);
        q10 = -1*ones(Lin, Lin, numSamplingPoints);
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
                tempq01 = (repmat(~thisMSALogical(:,l),1, Lin).* thisMSALogical);
                tempq10 = (repmat(thisMSALogical(:,l),1, Lin).* ~thisMSALogical);
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);

                % sum for the lth row of the covariance matrix
                q01(l,:,t) = sum(tempq01)./numSamplesPerGenSelected;
                q10(l,:,t) = sum(tempq10)./numSamplesPerGenSelected;
                q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;
            end
        end
        disp('done')
    end
    
    
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
    
    

    
    qOrig = q;
    N = Nin;
    fprintf('Filtering trajectories...')
    FiltTrajIlling;
    disp('done.')
    model = 'Illing';
    q = qTemp3;
    
    timePointsToDrop = 0;
%%


    % =====  C. correct sampling problems and again find q, qij  =====
    %--------------------------------------------------------------------------
    % L : length of sequence
    % T : total number of generations

    % this code makes sampled MSA from the output of the WF simulation, outer
    % loop same as that used in step A, this one also calculates qij


    qijAtTimeTk = cell(Lin,Lin);
    qijAtTimeTk(:) = {zeros(2)};

    qijAll = cell(1, numSamplingPoints);
    qijAll(:) ={qijAtTimeTk};


    % this section calculates q_t for each site, for the synthetic protein, 
    % 2 allele per site case this is just the sum of occurances of ones


    for t = 1:(numSamplingPoints -timePointsToDrop)
    %    t
        thisSamplingTimePoint = samplingTimePoints(t);
       % calculate qij for bialellic (0,1) system
        qijAtTimeTk = cell(Lin,Lin);
        qijAtTimeTk(:) = {zeros(2)};
        % calculate qij
        for i = 1:Lin
            for j = 1:Lin
                if(~(i==j))
                    q11Temp333 = q11(i,j,t);
                    q10Temp333 = q10(i,j,t);
                    q01Temp333 = q01(i,j,t);
                    q00Temp333 = 1 - (q01Temp333 + q10Temp333 + q11Temp333);
                    % early check

                    if((validSelCoeff(t,j) == 1) && (validQ(t,i) == 1))
                         if(q00Temp333 == 0 && q01Temp333 == 0)
                            if(debuggingDisplay ==1)
                               disp(' Error 0: q00 and q01 are both zero! (early warning)')
                               zTemp2 = [q00Temp333 q01Temp333; q10Temp333 q11Temp333]
                               disp('i    j    t')
                               [i j t]
                               disp('changing to')
                            end
                            if(q00Temp333 == 0)
                                q00Temp333 = min(1 - qTemp3(t,i), 1 - qTemp3(t,j));
                            end
                            if(q01Temp333 == 0)
                                q01Temp333 = min(1 - qTemp3(t,i), 1 - qTemp3(t,j));
                            end
                            zTemp2 = [q00Temp333 q01Temp333; q10Temp333 q11Temp333];
    %                         if(sum(zTemp2(:)) > 1)
                               zTemp = [q00Temp333; q01Temp333; q10Temp333; q11Temp333];%zTemp2(:);
                               [zTempMaxVal, zTempMaxInd] = max(zTemp);
                               zTemp(zTempMaxInd) = 1 - sum(zTemp(setdiff(1:4,zTempMaxInd)));%zTempMaxVal - sum(zTemp(1:2));
                               q00Temp333 = zTemp(1);
                               q01Temp333 = zTemp(2);
                               q10Temp333 = zTemp(3);
                               q11Temp333 = zTemp(4);

                               if(q00Temp333 + q01Temp333 + q10Temp333 + q11Temp333 ~= 1)
                                   q00Temp333
                                   q01Temp333
                                   q10Temp333
                                   q11Temp333
                                   pause
                               end
    %                        end
                            if(debuggingDisplay ==1)
                                [q00Temp333 q01Temp333; q10Temp333 q11Temp333]
                                pause
                            end
                        end
                        if(q11Temp333 == 0 && q10Temp333 == 0)
                            if(debuggingDisplay ==1)
                               disp(' Error 1: q11 and q10 are both zero! (early warning)')
                               zTemp2 = [q00Temp333 q01Temp333; q10Temp333 q11Temp333]
                               disp('i    j    t')
                               [i j t]
                               disp('changing to')
                            end
                            if(q11Temp333 == 0)
                                q11Temp333 = min(qTemp3(t,i), qTemp3(t,j));
                            end
                            if(q10Temp333 == 0)
                                q10Temp333 = min(qTemp3(t,i), 1 - qTemp3(t,j));
                            end
                            zTemp2 = [q00Temp333 q01Temp333; q10Temp333 q11Temp333];
    %                         if(sum(zTemp2(:)) > 1)
                               zTemp = [q00Temp333; q01Temp333; q10Temp333; q11Temp333];%zTemp = zTemp2(:);
                               [zTempMaxVal, zTempMaxInd] = max(zTemp);
                               zTemp(zTempMaxInd) = 1 - sum(zTemp(setdiff(1:4,zTempMaxInd)));%zTempMaxVal - sum(zTemp([2 4]));
                               q00Temp333 = zTemp(1);
                               q01Temp333 = zTemp(2);
                               q10Temp333 = zTemp(3);
                               q11Temp333 = zTemp(4);
    %                        end
                               if(q00Temp333 + q01Temp333 + q10Temp333 + q11Temp333 ~= 1)
                                   q00Temp333
                                   q01Temp333
                                   q10Temp333
                                   q11Temp333
                                   pause
                               end
                            if(debuggingDisplay ==1)
                               [q00Temp333 q01Temp333; q10Temp333 q11Temp333]
                               pause
                            end
                        end


                    qijAtTimeTk_temp = [q00Temp333 q01Temp333; q10Temp333 q11Temp333];

                    qijAtTimeTk{i,j} = qijAtTimeTk_temp;
    %                 qijAtTimeTk{j,i} = qijAtTimeTk_temp';
                    end
                end

            end 
        end
        qijAll{t} = qijAtTimeTk;

    end
    clear qijAtTimeTk;
    disp('done')
    toc

    

    % make the gMat matrix for faster calculation in likelihood function
    gMat = zeros(Lin,Lin, (numSamplingPoints - timePointsToDrop));

    for t = 1:(numSamplingPoints - timePointsToDrop)
    %    t
         qijAtTimeTk = qijAll{t};
         for i = 1:Lin
             % if qi < threshold, sigmaEst_eff = 0, pp. 992
             if(validQ(t,i) == 1) 
            %disp('y')
                 tempSum = 0;
                 for j = 1:Lin

                    % qj < threshold, ignore its linkage effect pp. 992
                    if(~(i==j) && (validSelCoeff(t,j) == 1)) 
                        thisQ = qijAtTimeTk{i,j};

                        q00Temp444 = thisQ(1,1);
                        q01Temp444 = thisQ(1,2);
                        q10Temp444 = thisQ(2,1);
                        q11Temp444 = thisQ(2,2);


                        if(q00Temp444 == 0 && q01Temp444 == 0)
                            disp(' Error 0: q00 and q01 are both zero!')
                            [q00Temp444 q01Temp444; q10Temp444 q11Temp444]
                            disp('i   j   t')
                            [i j t]
                            pause
                        end
                        if(q11Temp444 == 0 && q10Temp444 == 0)
                            disp(' Error 1: q11 and q10 are both zero!')
                            [q00Temp444 q01Temp444; q10Temp444 q11Temp444]
                            disp('i   j   t')
                            [i j t]
                            pause
                        end

                        gMat(i,j, t) = (q11Temp444/(q11Temp444 + q10Temp444) - q01Temp444/(q01Temp444 + q00Temp444));
                    end
                 end
             end
         end
    end

    timeProcDataIlling(thisItr) = toc;

    %% optmization

    lastValidTimePoint = numGen - timePointsToDrop;


    sigmaEstOutLinkIllingAll = zeros(Lin, numOptRuns);
    fvalAll = zeros(1, numOptRuns);
    



    for k = 1:numOptRuns
        overAllLogLikelihood = @(sigmaEstIn) func1_v3_mex(sigmaEstIn, numSamplingPoints, Lin, timePointsToDrop, lastValidTimePoint, dT, ng, trajPerSite, validQ, validSelCoeff, q, siteTrajAll, gMat);
    
        lowerBound = -1*10^0*ones(1, Lin);
        upperBound = 1*10^0*ones(1, Lin);
        initialUpLim = max(perSiteSelction);
        initialLowLim = min(perSiteSelction);
        if(initCondTrue == 1)
            initialSigmaEstIn = perSiteSelction;
            initCond = 'True';
        elseif(initCondTrue == 0)
            initialSigmaEstIn = (initialUpLim - initialLowLim)*rand(1,Lin) + initialLowLim;
            initCond = 'Rand';
        elseif(initCondTrue == 2)
            initialSigmaEstIn = zeros(1,Lin);
            initCond = 'Zero';
        elseif(initCondTrue == 3)
            fluctuation = 0.1*( 2*max(perSiteSelction)*randn(1,Lin) - max(perSiteSelction)*ones(1,Lin) );
            initialSigmaEstIn = perSiteSelction + fluctuation;
            initCond = 'TrueFluc10';
        elseif(initCondTrue == 4)
            fluctuation = 0.25*( 2*max(perSiteSelction)*randn(1,Lin) - max(perSiteSelction)*ones(1,Lin) );
            initialSigmaEstIn = perSiteSelction + fluctuation;
            initCond = 'TrueFluc25';
        elseif(initCondTrue == 5)
            fluctuation = 0.5*( 2*max(perSiteSelction)*randn(1,Lin) - max(perSiteSelction)*ones(1,Lin) );
            initialSigmaEstIn = perSiteSelction + fluctuation;
            initCond = 'TrueFluc50';
        elseif(initCondTrue == 6)
            modelMPL = 'linkDiff';
            temp4356 = strfind(fileName, '_'); 
            fileNameSS = [fileName(1:temp4356(end-2)) modelMPL '_Set' num2str(thisSet) '_itr' num2str(thisItr) '.mat'];
           load([dirNameAnalysis fileNameSS], 'sigmaEstOutUnLink');
            initCond = 'SS';
        end


        % simulated annealing
        %==================
        optAlgo = 'SA';
        optionsSA = saoptimset('simulannealbnd');
    
        optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'TolFun',TolFunIn,'StallIterLim',StallIterLimIn,'Display',dispStr);
    
        [sigmaEstOutLinkIlling, fval] = simulannealbnd(overAllLogLikelihood, initialSigmaEstIn, lowerBound, upperBound, optionsSA);

    
        sigmaEstOutLinkIllingAll(:,k) = sigmaEstOutLinkIlling;
        fvalAll(k) = fval;
        timeEachSAItrIlling(k,thisItr) = toc;
        disp(['Total time required by run ' num2str(k) ' of iter ' num2str(thisItr) ': ' num2str(timeEachSAItrIlling(k,thisItr))])
    end
    if(numOptRuns > 1)
        sigmaEstOutLinkIlling = mean(sigmaEstOutLinkIllingAll,2)';
        sigmaEstOutLinkIlling = sigmaEstOutLinkIlling';
    else
        sigmaEstOutLinkIlling = sigmaEstOutLinkIllingAll;
    end
    timeRunAlgoIlling(thisItr) = toc;
    model = [model optAlgo initCond];
%%
%     thisRun = 1;
    if(loadDataOption == 1)        
        if(Tstart == 1)
            fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        else
             fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        end
    elseif(loadDataOption == 2)
        
        temp4356 = strfind(fileName, '_');
        fileNameSave = [fileName(1:temp4356(end-2)) model '_Set' num2str(thisSet) '_itr' num2str(thisItr) '.mat'];
    end

    timeWholeCodeIlling(thisItr) = toc;
    
    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysisIlling fileNameSave])
    end
    
    
    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);
    [~,~,~, aucLinkIllingItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLinkIlling, 1);
    [~,~,~, aucLinkIllingItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLinkIlling, 1);

    aucLinkIllingItr(thisItr,:) = aucLinkIllingItrTemp;

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

        fid = fopen(fileNamePaperDataCompDet,'wt');
        fprintf(fid, str);
        fclose(fid);
    end


    inpDataDet = [num2str(thisItr-1) ',' num2str(thisItr-1) ',Det,' num2str(Tstart) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeWholeCodeIlling(thisItr))];
    for l = 1:Lin
        inpDataDet = [inpDataDet ',' num2str(sigmaEstOutLinkIlling(l))];
    end
    inpDataDet = [inpDataDet ',' num2str(aucLinkIllingItrTemp(1)) ',' num2str(aucLinkIllingItrTemp(2))];
    for l = 1:Lin
        inpDataDet = [inpDataDet ',' num2str(perSiteSelction(l) - sigmaEstOutLinkIlling(l))];
    end
    inpDataDet = [inpDataDet '\n'];
    fid = fopen(fileNamePaperDataCompDet,'a');
    fprintf(fid, inpDataDet);
    fclose(fid);
end
disp(['Output saved to ...' chosenSlash 'Matlab' chosenSlash  fileNamePaperDataCompDet])