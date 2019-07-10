%% Implements the linked-diffusion method with linear interpolation for 
% 3 class of sites (deleterious, beneficial, neutral) sims 

%  1. Load data   
%  2. Reconstruct MSA from raw data files, perform finite sampling by 
%     selecting x_g individuals from the MSA and find sampled 1 and 2 point
%     frequencies
%  6. Estimation...SW diff with linkage for Linear interpolated freq 
%     vectors with mu
%  7. Calculate all y12 values and LD (both measures of linkage)

% this was called untitled53.m previously

% Last updated Feb 2017
% Author: M Saqib Sohail




clc
clear all
close all

debuggingDisplay = 0;
saveFile = 0%1;

thisSet = 'medium_complex';
numItr = 1%90;
useSameT = 0; % flag that controls if 
              %     1: usable T will be loaded from data
              %     0: supplied by user in variable Tused

%Tused = 1000;%%2500;%300;%10000;
%Tstart = 1%1500%1001;
loadDataOption = 2%1;
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

    
    [dirNameData, dirNameAnalysis, dirNameABC] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData thisSet chosenSlash];
    dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];
    

    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end

    dirNameAnalysisIlling = [dirNameAnalysis 'Illing' chosenSlash];
    if(exist(dirNameAnalysisIlling, 'dir') == 0)
        mkdir(dirNameAnalysisIlling)
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


%fprintf('Calculate cleaned 1 and 2 point probabilities...')
%  we drop last time points at this stage and not before
%parfor t = 1:(numSamplingPoints -timePointsToDrop)
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

% % this is input for the likelihood calculation function
% % find indices of trajectory start ans stop from validQ
% siteTrajAll = zeros(maxTrajPerSite, 2 ,Lin);
% trajPerSite = zeros(1,Lin);
% for l = 1:L
%     thisSiteAllValidTraj = validQ(:,l) == 1;
%     
%     temp12 = thisSiteAllValidTraj(2:end) - thisSiteAllValidTraj(1:end-1);
%     alltrajStartThisSite = find(temp12 == 1) + 1;
%     alltrajEndThisSite = find(temp12 == -1);
%     enteriesThisSite = length(alltrajStartThisSite);
%     siteTrajAll(1:enteriesThisSite,:,l) = [alltrajStartThisSite alltrajEndThisSite];
%     trajPerSite(l) = enteriesThisSite;
% end

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
%parfor k = 1:numOptRuns



for k = 1:numOptRuns
    

    
    overAllLogLikelihood = @(sigmaEstIn) func1_v3_mex(sigmaEstIn, numSamplingPoints, Lin, timePointsToDrop, lastValidTimePoint, dT, ng, trajPerSite, validQ, validSelCoeff, q, siteTrajAll, gMat);
    %overAllLogLikelihood = @(sigmaEstIn) func1_v4(sigmaEstIn, numSamplingPoints, Lin, timePointsToDrop, dT, ng, exp700, expm700, trajPerSite, validQ, q, siteTrajAll, gMat);
    %overAllLogLikelihood = @(sigmaEstIn) func1(sigmaEstIn ,numSamplingPoints, Lin, timePointsToDrop, lastValidTimePoint, qijAll, validQ, validSelCoeff, dT, q, ng);


    
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
    
    % lowerBound = [-1*10^-2*ones(1, 5) -1*10^-3*ones(1,15)];
    % upperBound = [1*10^-2*ones(1, 5) 1*10^-3*ones(1,15)];

    % lowerBound = zeros(1, Lin);
    % upperBound = zeros(1, Lin);
    % lbVal1 = -1*10^-2;
    % lbVal2 = -1*10^-3;
    % ubVal1 = 1*10^-2;
    % ubVal2 = 1*10^-3;
    
    % % find which set of sites fixate atleast once, allow theam higher
    % % coeeficents
    % setNonFixate = [];
    % for l = 1:L
    %     
    %     numOnes_temp = sum(q(1:lastValidTimePoint,l) == 1);
    %     numZeros_temp = sum(q(1:lastValidTimePoint,l) == 0);
    %     if(numOnes_temp == 0 || numZeros_temp == 0) % did not fixate
    %         setNonFixate = [setNonFixate l];
    %     end
    % end
    % setFixate = setdiff(1:Lin,setNonFixate);
    % setNonFixate;
    % 
    % 
    % 
    % lowerBound(setFixate) = lbVal1;
    % lowerBound(setNonFixate) = lbVal2;
    % upperBound(setFixate) = ubVal1;
    % upperBound(setNonFixate) = ubVal2;


    % simulated annealing
    %==================
    optAlgo = 'SA';
    optionsSA = saoptimset('simulannealbnd');
%    optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'Display',dispStr);
    optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'TolFun',TolFunIn,'StallIterLim',StallIterLimIn,'Display',dispStr);
%     if(optimizeOption == 1)
%         optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'Display',dispStr);
%     elseif(optimizeOption == 2)
%         optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'TolFun',TolFunIn,'StallIterLim',StallIterLimIn,'Display',dispStr);
%     end
    %optionsSA = saoptimset(optionsSA,'TolFun',1e-3, 'InitialTemperature',100,'Display', 'diagnose');
    %optionsSA = saoptimset(optionsSA,'InitialTemperature',100,'Display', 'diagnose');
    % optionsSA.TolFun = 1e-1;
    [sigmaEstOutLinkIlling, fval] = simulannealbnd(overAllLogLikelihood, initialSigmaEstIn, lowerBound, upperBound, optionsSA);

    % % GA
    % %==================
    % optAlgo = 'GA';
    % numberOfVariables = L;
    % optionsGA = gaoptimset('CrossoverFraction',0.9,'PopulationSize',300,...
    %             'Generations',350,'EliteCount',1,'InitialPopulation',[],...
    %             'CrossoverFcn',@crossoverheuristic,'MutationFcn',{@mutationuniform,0.15},... 
    %             'TolFun',1e-12,'SelectionFcn',@selectionuniform,...
    %             'PopInitRange',[lowerBound;upperBound]);
    % [sigmaEstOutLinkIlling_GA, fval_GA] = ga(overAllLogLikelihood,numberOfVariables,[],[],[],[],lowerBound,upperBound,[],optionsGA);

    % % CLPSO
    % %==================
    % optAlgo = 'CLPSO';
    % [g_best_clpso,cost_gb_clpso] = optimizeWithCLPSO(Lin, lowerBound, upperBound, overAllLogLikelihood);

    % sigmaEstOutLinkIlling = g_best_clpso(end,:);
    % fval = cost_gb_clpso(end);
    % % HPSO
    % %==================
%     optAlgo = 'HPSO';
%     [g_best_hpso,cost_gb_hpso] = optimizeWithHPSO(Lin, lowerBound, upperBound, overAllLogLikelihood);
% 
%     sigmaEstOutLinkIlling = g_best_hpso(end,:);
%     fval = cost_gb_hpso(end);
%      
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

%     fileExist = exist([dirNameAnalysis fileNameSave], 'file');
%     while(fileExist == 2)
%         thisRun = thisRun + 1;
%         fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model  '_filter_' filterStr '_' num2str(thisRun) '_M3' fileName(end-3:end)];
%         fileExist = exist([dirNameAnalysis fileNameSave], 'file');
%     end

    timeWholeCodeIlling(thisItr) = toc;
    timeWholeCodeIlling(thisItr)
    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysisIlling fileNameSave])
    end
    
end
