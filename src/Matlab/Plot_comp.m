%% Method comparison 
% Comparison between MPLin, SS, ABC, FIT, ApproxWF, Deterministic also
% aggregate's results in .CSV file 
%%

% this code :
%      1. load data
%      2. calculate NRMSE, AUROC per iteration
%      3. aggregate results in a .CSV file
%      4. calculate and plot AUROC, Bubble plot, violin plot

% Last updated 31-Oct 2017
% Author: M Saqib Sohail


clc
close all
clear all

init_ColorNames_Settings;
% set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'DefaultTextFontName','CMU Serif Roman')
% set(0,'DefaultAxesFontSize',15)
% set(0,'DefaultTextFontSize',15)

% Setting the default renderer to "Painters" instead of "OpenGL"
set(0, 'DefaultFigureRenderer', 'painters');
 
% Setting font and size
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)
saveFigs = 0;
fileNameContainingDirPath = 'dirNames.txt';

thisSet = 'medium_simple';%'medium_complex';
modelIlling = 'IllingSASS';%Rand';%'IllingSAZero';% 'IllingSATrue'; 'IllingSARand';
getSysParam_4;

numItr = 3%90%90%250;
T = Tused + Tstart - 1;

actualT = T/1000;
posSelc = selVal/2/Nin;
negSelc = delSelVal/2/Nin;
lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};
numTps = Tused/dT;

dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    display('Error: system si not unix and not PC...')
    pause
end

[~, dirNameAnalysis, dirNameABC] = loadDirNames(fileNameContainingDirPath);
dirNameAnalysis = [dirNameAnalysis thisSet chosenSlash];

trueClassROC_AllItr = [];
sigmaEstOutLink_AllItr = [];
sigmaEstOutUnLink_AllItr = [];
sigmaEstOutUnLinkNoMu_AllItr = [];
sigmaEstOutLink_LinUni_AllItr = [];
sigmaEstOutLinkNoMu_AllItr = [];
sigmaEstOutFoll_mean_AllItr = [];
sigmaEstOutFoll_median_AllItr = [];
sigmaEstOutFoll_ML_AllItr = [];
sigmaEstOutFeder_AllItr = [];
sigmaEstOutFerrer_AllItr = [];
sigmaEstOutLinkIlling_AllItr = [];

perSiteSelction_AllItr  = [];
siteToExculeCuzOfLimitedPolyTps_AllItr = [];
aucLinkItr = zeros(numItr,2);
aucLinkNoMuItr = zeros(numItr,2);
aucUnLinkItr = zeros(numItr,2);
aucUnLinkNoMuItr = zeros(numItr,2);
aucLinkLinUniItr = zeros(numItr,2);
aucFoll_meanItr = zeros(numItr,2);
aucFoll_medianItr = zeros(numItr,2);
aucFoll_MLItr = zeros(numItr,2);
aucFederItr = zeros(numItr,2);
aucFerrerItr = zeros(numItr,2);
aucLinkIllingItr = zeros(numItr,2);
if(strcmp(thisSet, 'small_simple'))
    fileNamePaperDataCompABC = 'ABC_small_simple_collected_extended.csv';
    fileNamePaperDataCompFIT = 'FIT_small_simple_collected_extended.csv';
    fileNamePaperDataCompApproxWF = 'ApproxWF_small_simple_collected_extended.csv';
    fileNamePaperDataCompDet = 'Det_small_simple_collected_extended.csv';
elseif(strcmp(thisSet, 'small_complex'))
    fileNamePaperDataCompABC = 'ABC_small_complex_collected_extended.csv';
    fileNamePaperDataCompFIT = 'FIT_small_complex_collected_extended.csv';
    fileNamePaperDataCompApproxWF = 'ApproxWF_small_complex_collected_extended.csv';
    fileNamePaperDataCompDet = 'Det_small_complex_collected_extended.csv';
elseif(strcmp(thisSet, 'medium_simple'))
    fileNamePaperDataCompABC = 'ABC_medium_simple_collected_extended.csv';
    fileNamePaperDataCompFIT = 'FIT_medium_simple_collected_extended.csv';
    fileNamePaperDataCompApproxWF = 'ApproxWF_medium_simple_collected_extended.csv';
    fileNamePaperDataCompDet = 'Det_medium_simple_collected_extended.csv';
elseif(strcmp(thisSet, 'medium_complex'))
    fileNamePaperDataCompABC = 'ABC_medium_comple_collected_extended.csv';
    fileNamePaperDataCompFIT = 'FIT_medium_comple_collected_extended.csv';
    fileNamePaperDataCompApproxWF = 'ApproxWF_medium_comple_collected_extended.csv';
    fileNamePaperDataCompDet = 'Det_medium_comple_collected_extended.csv';
end

for itr = 1:numItr
    %% 1. load data
    fileName = ['wfsim_' thisSet '_' num2str(itr-1) '_T' num2str(T) '_ns' num2str(ng) '_dt' num2str(dT) '_output.mat'];
    %load([dirNameAnalysis fileName], 'L', 'perSiteSelction', ...
    %         'sigmaEstOutLink', 'sigmaEstOutLink_LinUni', 'sigmaEstOutLinkNoMu',...
    %         'sigmaEstOutUnLink', 'sigmaEstOutUnLinkNoMu', 'q', 'priorConst', 'timeWholeCodeMPL');
    
    sigmaEstOutLink = zeros(Lin,1);
    sigmaEstOutLink_LinUni = zeros(Lin,1);
    sigmaEstOutLinkNoMu = zeros(Lin,1);
    sigmaEstOutUnLink = zeros(Lin,1);
    sigmaEstOutUnLinkNoMu = zeros(Lin,1);
    timeMPLItr(itr) = 0%timeWholeCodeMPL(itr);
    % identify neutral, pos/neg selection
    trueClassROC = 3*ones(1, Lin);
    trueClassROC(perSiteSelction == posSelc) = 4;
    trueClassROC(perSiteSelction == negSelc) = 5;
    
    trueClassROC_AllItr = [trueClassROC_AllItr trueClassROC];
    sigmaEstOutLink_AllItr = [sigmaEstOutLink_AllItr sigmaEstOutLink'];
    sigmaEstOutUnLink_AllItr = [sigmaEstOutUnLink_AllItr sigmaEstOutUnLink'];
    sigmaEstOutUnLinkNoMu_AllItr = [sigmaEstOutUnLinkNoMu_AllItr sigmaEstOutUnLinkNoMu'];
    sigmaEstOutLink_LinUni_AllItr = [sigmaEstOutLink_LinUni_AllItr sigmaEstOutLink_LinUni'];
    sigmaEstOutLinkNoMu_AllItr = [sigmaEstOutLinkNoMu_AllItr sigmaEstOutLinkNoMu'];
    perSiteSelction_AllItr = [perSiteSelction_AllItr perSiteSelction];
    
    
    fileNameFoll = ['Foll_' fileName(1:end-10) 'posterior_s.mat'];
    load([dirNameAnalysis 'Foll' chosenSlash fileNameFoll], 'sigmaEstOutFoll_mean', ...
        'sigmaEstOutFoll_median', 'sigmaEstOutFoll_ML', 'timeWholeCodeFoll');%, 'posterior_s');

    timeFollItr(itr) = timeWholeCodeFoll(itr);
    sigmaEstOutFoll_mean = sigmaEstOutFoll_mean';
    sigmaEstOutFoll_median = sigmaEstOutFoll_median';
    
    sigmaEstOutFoll_mean_AllItr = [sigmaEstOutFoll_mean_AllItr sigmaEstOutFoll_mean'];
    sigmaEstOutFoll_median_AllItr = [sigmaEstOutFoll_median_AllItr sigmaEstOutFoll_median'];
    sigmaEstOutFoll_ML_AllItr = [sigmaEstOutFoll_ML_AllItr sigmaEstOutFoll_ML'];
    
    
    fileNameFeder = ['Feder_' fileName(1:end-10) 'output.mat'];
    load([dirNameAnalysis 'Feder' chosenSlash fileNameFeder], 'sigmaEstFeder', 't_FI', 'timeWholeCodeFeder');
    
    timeFederItr(itr) = timeWholeCodeFeder(itr);
    sigmaEstOutFeder = sigmaEstFeder';
    sigmaEstOutFeder = t_FI';
    sigmaEstOutFeder_AllItr = [sigmaEstOutFeder_AllItr sigmaEstOutFeder'];
    
    fileNameFerrer = ['Ferrer_' fileName(1:end-10) 'output.mat'];
    
    load([dirNameAnalysis 'Ferrer' chosenSlash fileNameFerrer], 'sigmaEstOutFerrer', 'timeWholeCodeFerrer');
    timeFerrerItr(itr) = timeWholeCodeFerrer(itr);
%sigmaEstOutFerrer = zeros(Lin,1);
%timeFerrerItr(itr) = 0;%526.34;%466.435; % 1990: 526.34 (5), 1991: 466.435 (30)
    sigmaEstOutFerrer_AllItr = [sigmaEstOutFerrer_AllItr sigmaEstOutFerrer'];
    

    fileNameIlling = ['Illing_' fileName(1:end-10) 'output.mat'];
    %load([dirNameAnalysis 'Illing' chosenSlash fileNameIlling], 'sigmaEstOutLinkIlling', 'timeWholeCodeIlling');
    %timeIllingItr(itr) = timeWholeCodeIlling(itr);
sigmaEstOutLinkIlling = zeros(Lin,1);
timeIllingItr(itr) = 0;
    sigmaEstOutLinkIlling_AllItr = [sigmaEstOutLinkIlling_AllItr sigmaEstOutLinkIlling'];

    
    iMax = 0;%numTps*0.02;
    siteToExculeCuzOfLimitedPolyTps = zeros(1, Lin);
    for i = 1:iMax
        temp = sum(q > 0 & q < 1) == i;
        siteToExculeCuzOfLimitedPolyTps = siteToExculeCuzOfLimitedPolyTps | temp;
    end
    
    siteToExculeCuzOfLimitedPolyTps_AllItr = [siteToExculeCuzOfLimitedPolyTps_AllItr siteToExculeCuzOfLimitedPolyTps];
    
    %% 2. calculate NRMSE, AUROC
    nrmse_3class_Link_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutLink').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_LinkNoMu_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutLinkNoMu').^2)/sum(abs(perSiteSelction).^2));
    nrmse_LinUni_3class_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutLink_LinUni').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_UnLink_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutUnLink').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_UnLinkNoMu_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutUnLinkNoMu').^2)/sum(abs(perSiteSelction).^2));

    nrmse_3class_Foll_mean_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutFoll_mean').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_Foll_median_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutFoll_median').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_Foll_ML_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutFoll_ML').^2)/sum(abs(perSiteSelction).^2));

    nrmse_3class_Feder_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutFeder').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_Ferrer_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutFerrer').^2)/sum(abs(perSiteSelction).^2));
    nrmse_3class_LinkIlling_itr(itr) = sqrt(sum(abs(perSiteSelction - sigmaEstOutLinkIlling').^2)/sum(abs(perSiteSelction).^2));

    nrmse_3class_Link_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLink(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_LinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLinkNoMu(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_LinUni_3class_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLink_LinUni(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_UnLink_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutUnLink(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_UnLinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutUnLinkNoMu(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));

    nrmse_3class_Foll_mean_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutFoll_mean(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_Foll_median_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutFoll_median(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_Foll_ML_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutFoll_ML(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));

    nrmse_3class_Feder_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutFeder(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_Ferrer_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutFerrer(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
    nrmse_3class_LinkIlling_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLinkIlling(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));

    nrmse_3class_Link_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLink(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_LinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLinkNoMu(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_LinUni_3class_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLink_LinUni(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_UnLink_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutUnLink(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_UnLinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutUnLinkNoMu(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));

    nrmse_3class_Foll_mean_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutFoll_mean(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_Foll_median_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutFoll_median(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_Foll_ML_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutFoll_ML(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));

    nrmse_3class_Feder_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutFeder(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_Ferrer_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutFerrer(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
    nrmse_3class_LinkIlling_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLinkIlling(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));

    rmse_3class_Link_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLink(trueClassROC == 3)').^2));
    rmse_3class_LinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLinkNoMu(trueClassROC == 3)').^2));
    rmse_LinUni_3class_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLink_LinUni(trueClassROC == 3)').^2));
    rmse_3class_UnLink_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutUnLink(trueClassROC == 3)').^2));
    rmse_3class_UnLinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutUnLinkNoMu(trueClassROC == 3)').^2));

    rmse_3class_Foll_mean_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutFoll_mean(trueClassROC == 3)').^2));
    rmse_3class_Foll_median_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutFoll_median(trueClassROC == 3)').^2));
    rmse_3class_Foll_ML_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutFoll_ML(trueClassROC == 3)').^2));

    rmse_3class_Feder_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutFeder(trueClassROC == 3)').^2));
    rmse_3class_Ferrer_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutFerrer(trueClassROC == 3)').^2));
    rmse_3class_LinkIlling_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLinkIlling(trueClassROC == 3)').^2));

    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);
    [~,~,~, aucLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink, 1);
    [~,~,~, aucLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink, 1);

    [~,~,~, aucLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLinkNoMu, 1);
    [~,~,~, aucLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLinkNoMu, 1);

    [~,~,~, aucUnLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLink, 1);
    [~,~,~, aucUnLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLink, 1);

    [~,~,~, aucUnLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);
    [~,~,~, aucUnLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);

    [~,~,~, aucLinkLinUniItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink_LinUni, 1);
    [~,~,~, aucLinkLinUniItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink_LinUni, 1);
    
    [~,~,~, aucFoll_meanItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFoll_mean, 1);
    [~,~,~, aucFoll_meanItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFoll_mean, 1);

    [~,~,~, aucFoll_medianItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFoll_median, 1);
    [~,~,~, aucFoll_medianItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFoll_median, 1);

    [~,~,~, aucFoll_MLItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFoll_ML, 1);
    [~,~,~, aucFoll_MLItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFoll_ML, 1);

    [~,~,~, aucFederItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFeder, 1);
    [~,~,~, aucFederItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFeder, 1);

    [~,~,~, aucFerrerItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutFerrer, 1);
    [~,~,~, aucFerrerItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutFerrer, 1);
    
    [~,~,~, aucLinkIllingItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLinkIlling, 1);
    [~,~,~, aucLinkIllingItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLinkIlling, 1);

    aucLinkItr(itr,:) = aucLinkItrTemp;
    aucLinkNoMuItr(itr,:) = aucLinkNoMuItrTemp;
    aucUnLinkItr(itr,:) = aucUnLinkItrTemp;
    aucUnLinkNoMuItr(itr,:) = aucUnLinkNoMuItrTemp;
    aucLinkLinUniItr(itr,:) = aucLinkLinUniItrTemp;
    
    aucFoll_meanItr(itr,:) = aucFoll_meanItrTemp;
    aucFoll_medianItr(itr,:) = aucFoll_medianItrTemp;
    aucFoll_MLItr(itr,:) = aucFoll_MLItrTemp;
    
    aucFederItr(itr,:) = aucFederItrTemp;
    aucFerrerItr(itr,:) = aucFerrerItrTemp;
    aucLinkIllingItr(itr,:) = aucLinkIllingItrTemp;
    
    %% 3. aggregate results in a .CSV file
    if(itr == 1)
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
        fid = fopen(fileNamePaperDataCompFIT,'wt');
        fprintf(fid, str);
        fclose(fid);
        fid = fopen(fileNamePaperDataCompApproxWF,'wt');
        fprintf(fid, str);
        fclose(fid);
        fid = fopen(fileNamePaperDataCompDet,'wt');
        fprintf(fid, str);
        fclose(fid);
    end
    inpDataABC = [num2str(itr-1) ',' num2str(itr-1) ',ABC,' num2str(Tstart-1) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeFollItr(itr))];
    inpDataFIT = [num2str(itr-1) ',' num2str(itr-1) ',FIT,' num2str(Tstart-1) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeFederItr(itr))];
    inpDataApproxWF = [num2str(itr-1) ',' num2str(itr-1) ',ApproxWF,' num2str(Tstart-1) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeFerrerItr(itr))];
    inpDataDet = [num2str(itr-1) ',' num2str(itr-1) ',Det,' num2str(Tstart-1) ',' num2str(T) ',' num2str(ng) ',' num2str(dT) ',' num2str(timeFerrerItr(itr))];
    for l = 1:Lin
        inpDataABC = [inpDataABC ',' num2str(sigmaEstOutFoll_ML(l))];
        inpDataFIT = [inpDataFIT ',' num2str(sigmaEstOutFeder(l))];
        inpDataApproxWF = [inpDataApproxWF ',' num2str(sigmaEstOutFerrer(l))];
        inpDataDet = [inpDataDet ',' num2str(sigmaEstOutLinkIlling(l))];
    end
    inpDataABC = [inpDataABC ',' num2str(aucFoll_MLItrTemp(1)) ',' num2str(aucFoll_MLItrTemp(2))];
    inpDataFIT = [inpDataFIT ',' num2str(aucFederItrTemp(1)) ',' num2str(aucFederItrTemp(2))];
    inpDataApproxWF = [inpDataApproxWF ',' num2str(aucFerrerItrTemp(1)) ',' num2str(aucFerrerItrTemp(2))];
    inpDataDet = [inpDataDet ',' num2str(aucLinkIllingItrTemp(1)) ',' num2str(aucLinkIllingItrTemp(2))];
    for l = 1:Lin
        inpDataABC = [inpDataABC ',' num2str(perSiteSelction(l) - sigmaEstOutFoll_ML(l))];
        inpDataFIT = [inpDataFIT ',' num2str(perSiteSelction(l) - sigmaEstOutFeder(l))];
        inpDataApproxWF = [inpDataApproxWF ',' num2str(perSiteSelction(l) - sigmaEstOutFerrer(l))];
        inpDataDet = [inpDataDet ',' num2str(perSiteSelction(l) - sigmaEstOutLinkIlling(l))];
    end
    inpDataABC = [inpDataABC '\n'];
    inpDataFIT = [inpDataFIT '\n'];
    inpDataApproxWF = [inpDataApproxWF '\n'];
    inpDataDet = [inpDataDet '\n'];
    fid = fopen(fileNamePaperDataCompABC,'a');
    fprintf(fid, inpDataABC);
    fclose(fid);
    fid = fopen(fileNamePaperDataCompFIT,'a');
    fprintf(fid, inpDataFIT);
    fclose(fid);
    fid = fopen(fileNamePaperDataCompApproxWF,'a');
    fprintf(fid, inpDataApproxWF);
    fclose(fid);
    fid = fopen(fileNamePaperDataCompDet,'a');
    fprintf(fid, inpDataDet);
    fclose(fid);

end
%% 4. calculate and plot AUROC, Bubble plot, violin plot
nrmse_3class_Link_itr_mean = mean(nrmse_3class_Link_itr);
nrmse_3class_LinkNoMu_itr_mean = mean(nrmse_3class_LinkNoMu_itr);
nrmse_LinUni_3class_itr_mean = mean(nrmse_LinUni_3class_itr);
nrmse_3class_UnLink_itr_mean = mean(nrmse_3class_UnLink_itr);
nrmse_3class_UnLinkNoMu_itr_mean = mean(nrmse_3class_UnLinkNoMu_itr);

nrmse_3class_Foll_mean_itr_mean = mean(nrmse_3class_Foll_mean_itr);
nrmse_3class_Foll_median_itr_mean = mean(nrmse_3class_Foll_median_itr);
nrmse_3class_Foll_ML_itr_mean = mean(nrmse_3class_Foll_ML_itr);

nrmse_3class_Feder_itr_mean = mean(nrmse_3class_Feder_itr);
nrmse_3class_Ferrer_itr_mean = mean(nrmse_3class_Ferrer_itr);
nrmse_3class_LinkIlling_itr_mean = mean(nrmse_3class_LinkIlling_itr);


nrmse_3class_Link_itrPos_mean = mean(nrmse_3class_Link_itrPos);
nrmse_3class_LinkNoMu_itrPos_mean = mean(nrmse_3class_LinkNoMu_itrPos);
nrmse_LinUni_3class_itrPos_mean = mean(nrmse_LinUni_3class_itrPos);
nrmse_3class_UnLink_itrPos_mean = mean(nrmse_3class_UnLink_itrPos);
nrmse_3class_UnLinkNoMu_itrPos_mean = mean(nrmse_3class_UnLinkNoMu_itrPos);

nrmse_3class_Foll_mean_itrPos_mean = mean(nrmse_3class_Foll_mean_itrPos);
nrmse_3class_Foll_median_itrPos_mean = mean(nrmse_3class_Foll_median_itrPos);
nrmse_3class_Foll_ML_itrPos_mean = mean(nrmse_3class_Foll_ML_itrPos);

nrmse_3class_Feder_itrPos_mean = mean(nrmse_3class_Feder_itrPos);
nrmse_3class_Ferrer_itrPos_mean = mean(nrmse_3class_Ferrer_itrPos);
nrmse_3class_LinkIlling_itrPos_mean = mean(nrmse_3class_LinkIlling_itrPos);

nrmse_3class_Link_itrNeg_mean = mean(nrmse_3class_Link_itrNeg);
nrmse_3class_LinkNoMu_itrNeg_mean = mean(nrmse_3class_LinkNoMu_itrNeg);
nrmse_LinUni_3class_itrNeg_mean = mean(nrmse_LinUni_3class_itrNeg);
nrmse_3class_UnLink_itrNeg_mean = mean(nrmse_3class_UnLink_itrNeg);
nrmse_3class_UnLinkNoMu_itrNeg_mean = mean(nrmse_3class_UnLinkNoMu_itrNeg);

nrmse_3class_Foll_mean_itrNeg_mean = mean(nrmse_3class_Foll_mean_itrNeg);
nrmse_3class_Foll_median_itrNeg_mean = mean(nrmse_3class_Foll_median_itrNeg);
nrmse_3class_Foll_ML_itrNeg_mean = mean(nrmse_3class_Foll_ML_itrNeg);

nrmse_3class_Feder_itrNeg_mean = mean(nrmse_3class_Feder_itrNeg);
nrmse_3class_Ferrer_itrNeg_mean = mean(nrmse_3class_Ferrer_itrNeg);
nrmse_3class_LinkIlling_itrNeg_mean = mean(nrmse_3class_LinkIlling_itrNeg);

rmse_3class_Link_itrNeu_mean = mean(rmse_3class_Link_itrNeu);
rmse_3class_LinkNoMu_itrNeu_mean = mean(rmse_3class_LinkNoMu_itrNeu);
rmse_LinUni_3class_itrNeu_mean = mean(rmse_LinUni_3class_itrNeu);
rmse_3class_UnLink_itrNeu_mean = mean(rmse_3class_UnLink_itrNeu);
rmse_3class_UnLinkNoMu_itrNeu_mean = mean(rmse_3class_UnLinkNoMu_itrNeu);

rmse_3class_Foll_mean_itrNeu_mean = mean(rmse_3class_Foll_mean_itrNeu);
rmse_3class_Foll_median_itrNeu_mean = mean(rmse_3class_Foll_median_itrNeu);
rmse_3class_Foll_ML_itrNeu_mean = mean(rmse_3class_Foll_ML_itrNeu);

rmse_3class_Feder_itrNeu_mean = mean(rmse_3class_Feder_itrNeu);
rmse_3class_Ferrer_itrNeu_mean = mean(rmse_3class_Ferrer_itrNeu);
rmse_3class_LinkIlling_itrNeu_mean = mean(rmse_3class_LinkIlling_itrNeu);

nrmse_3class_mean = [nrmse_3class_Link_itr_mean nrmse_3class_LinkNoMu_itr_mean ...
                     nrmse_3class_UnLink_itr_mean nrmse_3class_UnLinkNoMu_itr_mean ...
                     nrmse_3class_Foll_ML_itr_mean nrmse_3class_Feder_itr_mean ...
                     nrmse_3class_Ferrer_itr_mean nrmse_3class_LinkIlling_itr_mean];

nrmse_3class_Pos_mean = [nrmse_3class_Link_itrPos_mean nrmse_3class_LinkNoMu_itrPos_mean ...
                         nrmse_3class_UnLink_itrPos_mean nrmse_3class_UnLinkNoMu_itrPos_mean ...
                         nrmse_3class_Foll_ML_itrPos_mean nrmse_3class_Feder_itrPos_mean ...
                         nrmse_3class_Ferrer_itrPos_mean nrmse_3class_LinkIlling_itrPos_mean];

nrmse_3class_Neg_mean = [nrmse_3class_Link_itrNeg_mean nrmse_3class_LinkNoMu_itrNeg_mean ...
                         nrmse_3class_UnLink_itrNeg_mean nrmse_3class_UnLinkNoMu_itrNeg_mean ...
                         nrmse_3class_Foll_ML_itrNeg_mean nrmse_3class_Feder_itrNeg_mean ...
                         nrmse_3class_Ferrer_itrNeg_mean nrmse_3class_LinkIlling_itrNeg_mean];

rmse_3class_Neu_mean = [rmse_3class_Link_itrNeu_mean rmse_3class_LinkNoMu_itrNeu_mean ...
                        rmse_3class_UnLink_itrNeu_mean rmse_3class_UnLinkNoMu_itrNeu_mean ...
                        rmse_3class_Foll_ML_itrNeu_mean rmse_3class_Feder_itrNeu_mean ...
                        rmse_3class_Ferrer_itrNeu_mean rmse_3class_LinkIlling_itrNeu_mean];

sigmaEstOutLink_AllItr_FULL = sigmaEstOutLink_AllItr;
sigmaEstOutUnLink_AllItrr_FULL = sigmaEstOutUnLink_AllItr;
sigmaEstOutUnLinkNoMu_AllItrr_FULL = sigmaEstOutUnLinkNoMu_AllItr;
sigmaEstOutLink_LinUni_AllItrr_FULL = sigmaEstOutLink_LinUni_AllItr;
sigmaEstOutLinkNoMu_AllItrr_FULL = sigmaEstOutLinkNoMu_AllItr;
perSiteSelction_AllItrr_FULL = perSiteSelction_AllItr;

sigmaEstOutFoll_mean_AllItr_FULL = sigmaEstOutFoll_mean_AllItr;
sigmaEstOutFoll_median_AllItr_FULL = sigmaEstOutFoll_median_AllItr;
sigmaEstOutFoll_ML_AllItr_FULL = sigmaEstOutFoll_ML_AllItr;
sigmaEstOutFeder_AllItr_FULL = sigmaEstOutFeder_AllItr;
sigmaEstOutFerrer_AllItr_FULL = sigmaEstOutFerrer_AllItr;
sigmaEstOutLinkIlling_AllItr_FULL = sigmaEstOutLinkIlling_AllItr;

% siteToExculeCuzOfLimitedPolyTps_AllItr_ind = find(siteToExculeCuzOfLimitedPolyTps_AllItr);
% indToExclude = find(sigmaEstOutLink_AllItr == -1/0.0001/dT);
% indToInclude = setdiff(1:length(trueClassROC_AllItr), [indToExclude siteToExculeCuzOfLimitedPolyTps_AllItr_ind]);
% trueClassROC_AllItr = trueClassROC_AllItr(indToInclude);
% sigmaEstOutLink_AllItr = sigmaEstOutLink_AllItr(indToInclude);
% sigmaEstOutUnLink_AllItr = sigmaEstOutUnLink_AllItr(indToInclude);
% sigmaEstOutUnLinkNoMu_AllItr = sigmaEstOutUnLinkNoMu_AllItr(indToInclude);
% sigmaEstOutLink_LinUni_AllItr = sigmaEstOutLink_LinUni_AllItr(indToInclude);
% sigmaEstOutLinkNoMu_AllItr = sigmaEstOutLinkNoMu_AllItr(indToInclude);
% perSiteSelction_AllItr = perSiteSelction_AllItr(indToInclude);


% estimation performance

% nrmse_3class_Link = sqrt(sum(abs(perSiteSelction_AllItr - sigmaEstOutLink_AllItr).^2)/sum(abs(perSiteSelction_AllItr).^2));
% nrmse_3class_LinkNoMu = sqrt(sum(abs(perSiteSelction_AllItr - sigmaEstOutLinkNoMu_AllItr).^2)/sum(abs(perSiteSelction_AllItr).^2));
% nrmse_LinUni_3class = sqrt(sum(abs(perSiteSelction_AllItr - sigmaEstOutLink_LinUni_AllItr).^2)/sum(abs(perSiteSelction_AllItr).^2));
% nrmse_3class_UnLink = sqrt(sum(abs(perSiteSelction_AllItr - sigmaEstOutUnLink_AllItr).^2)/sum(abs(perSiteSelction_AllItr).^2));
% nrmse_3class_UnLinkNoMu = sqrt(sum(abs(perSiteSelction_AllItr - sigmaEstOutUnLinkNoMu_AllItr).^2)/sum(abs(perSiteSelction_AllItr).^2));
% 
% 
% % only positive sites
posOnly = perSiteSelction_AllItr == max(perSiteSelction_AllItr);
% nrmse_3class_LinkPosOnly = sqrt(sum(abs(perSiteSelction_AllItr(posOnly) - sigmaEstOutLink_AllItr(posOnly)).^2)/sum(abs(perSiteSelction_AllItr(posOnly)).^2));
% nrmse_3class_LinkNoMuPosOnly = sqrt(sum(abs(perSiteSelction_AllItr(posOnly) - sigmaEstOutLinkNoMu_AllItr(posOnly)).^2)/sum(abs(perSiteSelction_AllItr(posOnly)).^2));
% nrmse_LinUni_3classPosOnly = sqrt(sum(abs(perSiteSelction_AllItr(posOnly) - sigmaEstOutLink_LinUni_AllItr(posOnly)).^2)/sum(abs(perSiteSelction_AllItr(posOnly)).^2));
% nrmse_3class_UnLinkPosOnly = sqrt(sum(abs(perSiteSelction_AllItr(posOnly) - sigmaEstOutUnLink_AllItr(posOnly)).^2)/sum(abs(perSiteSelction_AllItr(posOnly)).^2));
% nrmse_3class_UnLinkNoMuPosOnly = sqrt(sum(abs(perSiteSelction_AllItr(posOnly) - sigmaEstOutUnLinkNoMu_AllItr(posOnly)).^2)/sum(abs(perSiteSelction_AllItr(posOnly)).^2));
% 
% % only negative sites
negOnly = perSiteSelction_AllItr == min(perSiteSelction_AllItr);
% nrmse_3class_LinkNegOnly = sqrt(sum(abs(perSiteSelction_AllItr(negOnly) - sigmaEstOutLink_AllItr(negOnly)).^2)/sum(abs(perSiteSelction_AllItr(negOnly)).^2));
% nrmse_3class_LinkNoMuNegOnly = sqrt(sum(abs(perSiteSelction_AllItr(negOnly) - sigmaEstOutLinkNoMu_AllItr(negOnly)).^2)/sum(abs(perSiteSelction_AllItr(negOnly)).^2));
% nrmse_LinUni_3classNegOnly = sqrt(sum(abs(perSiteSelction_AllItr(negOnly) - sigmaEstOutLink_LinUni_AllItr(negOnly)).^2)/sum(abs(perSiteSelction_AllItr(negOnly)).^2));
% nrmse_3class_UnLinkNegOnly = sqrt(sum(abs(perSiteSelction_AllItr(negOnly) - sigmaEstOutUnLink_AllItr(negOnly)).^2)/sum(abs(perSiteSelction_AllItr(negOnly)).^2));
% nrmse_3class_UnLinkNoMuNegOnly = sqrt(sum(abs(perSiteSelction_AllItr(negOnly) - sigmaEstOutUnLinkNoMu_AllItr(negOnly)).^2)/sum(abs(perSiteSelction_AllItr(negOnly)).^2));
% 
% 
% % only neutral sites
% neutOnly = perSiteSelction_AllItr == 0;
% rmse_3class_LinkNeutOnly = sqrt(sum(abs(perSiteSelction_AllItr(neutOnly) - sigmaEstOutLink_AllItr(neutOnly)).^2)/sum(neutOnly));
% rmse_3class_LinkNoMuNeutOnly = sqrt(sum(abs(perSiteSelction_AllItr(neutOnly) - sigmaEstOutLinkNoMu_AllItr(neutOnly)).^2)/sum(neutOnly));
% rmse_LinUni_3classNeutOnly = sqrt(sum(abs(perSiteSelction_AllItr(neutOnly) - sigmaEstOutLink_LinUni_AllItr(neutOnly)).^2)/sum(neutOnly));
% rmse_3class_UnLinkNeutOnly = sqrt(sum(abs(perSiteSelction_AllItr(neutOnly) - sigmaEstOutUnLink_AllItr(neutOnly)).^2)/sum(neutOnly));
% rmse_3class_UnLinkNoMuNeutOnly = sqrt(sum(abs(perSiteSelction_AllItr(neutOnly) - sigmaEstOutUnLinkNoMu_AllItr(neutOnly)).^2)/sum(neutOnly));


%

[xLink{1}, yLink{1}, zLink{1}, aucLink(1)] = perfcurve(double(posOnly), sigmaEstOutLink_AllItr, 1);
[xLink{2}, yLink{2}, zLink{2}, aucLink(2)] = perfcurve(double(~negOnly), sigmaEstOutLink_AllItr, 1);

[xLinkNoMu{1}, yLinkNoMu{1}, zLinkNoMu{1}, aucLinkNoMu(1)] = perfcurve(double(posOnly), sigmaEstOutLinkNoMu_AllItr, 1);
[xLinkNoMu{2}, yLinkNoMu{2}, zLinkNoMu{2}, aucLinkNoMu(2)] = perfcurve(double(~negOnly), sigmaEstOutLinkNoMu_AllItr, 1);

[xUnLink{1}, yUnLink{1}, zUnLink{1}, aucUnLink(1)] = perfcurve(double(posOnly), sigmaEstOutUnLink_AllItr, 1);
[xUnLink{2}, yUnLink{2}, zUnLink{2}, aucUnLink(2)] = perfcurve(double(~negOnly), sigmaEstOutUnLink_AllItr, 1);

[xUnLinkNoMu{1}, yUnLinkNoMu{1}, zUnLinkNoMu{1}, aucUnLinkNoMu(1)] = perfcurve(double(posOnly), sigmaEstOutUnLinkNoMu_AllItr, 1);
[xUnLinkNoMu{2}, yUnLinkNoMu{2}, zUnLinkNoMu{2}, aucUnLinkNoMu(2)] = perfcurve(double(~negOnly), sigmaEstOutUnLinkNoMu_AllItr, 1);

[xLinkLinUni{1}, yLinkLinUni{1}, zLinkLinUni{1}, aucLinkLinUni(1)] = perfcurve(double(posOnly), sigmaEstOutLink_LinUni_AllItr, 1);
[xLinkLinUni{2}, yLinkLinUni{2}, zLinkLinUni{2}, aucLinkLinUni(2)] = perfcurve(double(~negOnly), sigmaEstOutLink_LinUni_AllItr, 1);

[xFoll_mean{1}, yFoll_mean{1}, zFoll_mean{1}, aucFoll_mean(1)] = perfcurve(double(posOnly), sigmaEstOutFoll_mean_AllItr, 1);
[xFoll_mean{2}, yFoll_mean{2}, zFoll_mean{2}, aucFoll_mean(2)] = perfcurve(double(~negOnly), sigmaEstOutFoll_mean_AllItr, 1);

[xFoll_median{1}, yFoll_median{1}, zFoll_median{1}, aucFoll_median(1)] = perfcurve(double(posOnly), sigmaEstOutFoll_median_AllItr, 1);
[xFoll_median{2}, yFoll_median{2}, zFoll_median{2}, aucFoll_median(2)] = perfcurve(double(~negOnly), sigmaEstOutFoll_median_AllItr, 1);

[xFoll_ML{1}, yFoll_ML{1}, zFoll_ML{1}, aucFoll_ML(1)] = perfcurve(double(posOnly), sigmaEstOutFoll_ML_AllItr, 1);
[xFoll_ML{2}, yFoll_ML{2}, zFoll_ML{2}, aucFoll_ML(2)] = perfcurve(double(~negOnly), sigmaEstOutFoll_ML_AllItr, 1);

[xFeder{1}, yFeder{1}, zFeder{1}, aucFeder(1)] = perfcurve(double(posOnly), sigmaEstOutFeder_AllItr, 1);
[xFeder{2}, yFeder{2}, zFeder{2}, aucFeder(2)] = perfcurve(double(~negOnly), sigmaEstOutFeder_AllItr, 1);

[xFerrer{1}, yFerrer{1}, zFerrer{1}, aucFerrer(1)] = perfcurve(double(posOnly), sigmaEstOutFerrer_AllItr, 1);
[xFerrer{2}, yFerrer{2}, zFerrer{2}, aucFerrer(2)] = perfcurve(double(~negOnly), sigmaEstOutFerrer_AllItr, 1);

[xLinkIlling{1}, yLinkIlling{1}, zLinkIlling{1}, aucLinkIlling(1)] = perfcurve(double(posOnly), sigmaEstOutLinkIlling_AllItr, 1);
[xLinkIlling{2}, yLinkIlling{2}, zLinkIlling{2}, aucLinkIlling(2)] = perfcurve(double(~negOnly), sigmaEstOutLinkIlling_AllItr, 1);
figure
plot(xLink{1},yLink{1},lineCol{1}, 'color', cornFlowerBlue, 'LineWidth', 1, 'MarkerFaceColor', cornFlowerBlue)
hold on
plot(xLink{2},yLink{2},lineCol{2}, 'color', rosyBrown, 'LineWidth', 1, 'MarkerFaceColor', rosyBrown)
plot(xUnLink{1},yUnLink{1},'--', 'color', cornFlowerBlue, 'LineWidth', 1, 'MarkerFaceColor', cornFlowerBlue)
hold on
plot(xUnLink{2},yUnLink{2},'--', 'color', rosyBrown, 'LineWidth', 1, 'MarkerFaceColor', rosyBrown)
plot(xUnLinkNoMu{1},yUnLinkNoMu{1},':', 'color', cornFlowerBlue, 'LineWidth', 1, 'MarkerFaceColor', cornFlowerBlue)
hold on
plot(xUnLinkNoMu{2},yUnLinkNoMu{2},':', 'color', rosyBrown, 'LineWidth', 1, 'MarkerFaceColor', rosyBrown)

xlabel('False positive rate FPR')
ylabel('True positive rate TPR')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...%'YTick'       , 0:0.2:1, ...
  'LineWidth'   , 1        );
aucLink
aucLinkNoMu
aucUnLink
aucUnLinkNoMu
aucLinkLinUni

aucLinkAve = mean(aucLinkItr)
aucLinkNoMuAve = mean(aucLinkNoMuItr)
aucUnLinkAve = mean(aucUnLinkItr)
aucUnLinkNoMuAve = mean(aucUnLinkNoMuItr)
aucFoll_meanAve = mean(aucFoll_meanItr)
aucFoll_medianAve = mean(aucFoll_medianItr)
aucFoll_MLAve = mean(aucFoll_MLItr)
aucFederAve = mean(aucFederItr)
aucFerrerAve = mean(aucFerrerItr)
aucLinkIllingAve = mean(aucLinkIllingItr)

aucLinkSTD = std(aucLinkItr)
aucLinkNoMuSTD = std(aucLinkNoMuItr)
aucUnLinkSTD = std(aucUnLinkItr)
aucUnLinkNoMuSTD = std(aucUnLinkNoMuItr)
aucFoll_meanSTD = std(aucFoll_meanItr)
aucFoll_medianSTD = std(aucFoll_medianItr)
aucFoll_MLSTD = std(aucFoll_MLItr)
aucFederSTD = std(aucFederItr)
aucFerrerSTD = std(aucFerrerItr)
aucLinkIllingSTD = std(aucLinkIllingItr)





%%
sigmaEstOutLink_AllItr_average = mean(reshape(sigmaEstOutLink_AllItr, Lin, numItr)')
% sigmaEstOutUnLink_AllItr = [];
% sigmaEstOutUnLinkNoMu_AllItr = [];
% sigmaEstOutLink_LinUni_AllItr = [];
% sigmaEstOutLinkNoMu_AllItr = [];
% perSiteSelction_AllItr  = [];

% sigmaEstOutLink_AllItr_AllReplicates = (reshape(sigmaEstOutLink_AllItr, Lin, numItr)');
% sigmaEstOutLink_AllItr_average3Rep = mean(sigmaEstOutLink_AllItr_AllReplicates([11 19 29],:))
% 
% sqrt(sum(abs(perSiteSelction(1:10) - sigmaEstOutLink_AllItr_average(1:10)).^2)./sum(abs(perSiteSelction(1:10)).^2))
% sqrt(sum(abs(perSiteSelction(1:10) - sigmaEstOutLink_AllItr_average3Rep(1:10)).^2)./sum(abs(perSiteSelction(1:10)).^2))





%Labels = {'MPL', 'MPL no \mu', 'SS', 'SS no \mu', 'ABC ML', 'ABC mean', 'ABC median'}

comparisonPosAve = [aucLinkAve(1) aucLinkNoMuAve(1) aucUnLinkAve(1) aucUnLinkNoMuAve(1) aucFoll_MLAve(1) aucFederAve(1) aucFerrerAve(1) aucLinkIllingAve(1)];
comparisonNegAve = [aucLinkAve(2) aucLinkNoMuAve(2) aucUnLinkAve(2) aucUnLinkNoMuAve(2) aucFoll_MLAve(2) aucFederAve(2) aucFerrerAve(2) aucLinkIllingAve(2)];

% comparisonPosSTD = [aucLinkSTD(1) aucLinkNoMuSTD(1) aucUnLinkSTD(1) aucUnLink+NoMuSTD(1) aucFoll_MLSTD(1) aucFederSTD(1) aucFerrerSTD(1) aucLinkIllingSTD(1)];
% comparisonNegSTD = [aucLinkSTD(2) aucLinkNoMuSTD(2) aucUnLinkSTD(2) aucUnLinkNoMuSTD(2) aucFoll_MLSTD(2) aucFederSTD(2) aucFerrerSTD(2) aucLinkIllingSTD(2)];

Labels = {'MPL', 'MPL no \mu', 'SS', 'SS no \mu', 'ABC', 'Feder', 'Ferrer', 'Illing'}

figure
bar([comparisonPosAve])%, [cornFlowerBlue;cornFlowerBlue;cornFlowerBlue;cornFlowerBlue])
alpha 0.6
map6Paired = brewermap(6,'Paired');
map12Paired = brewermap(12,'Paired');
thisColorMap = map12Paired;
% colormap(thisColorMap)% map6Paired)


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTickLabel', Labels, ...'YLim', [0.5 1.001], ...  'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0.4:0.1:1, ...
  'LineWidth', 1)
ylabel('AUROC')
axis([0.5 8.5 0.50 1])
grid on
title(['Set : ' num2str(thisSet)])
title('Detecting beneficial sites')
% hold on
% errorbar(1:length(comparisonPosAve), comparisonPosAve, comparisonPosSTD, 'ko')
if(Lin == 10)
    siteNumAlphab = 'Ten';
elseif(Lin == 50)
    siteNumAlphab = 'Fifty';
elseif(Lin == 500)
    siteNumAlphab = 'FiveHund';
end
if(saveFigs)
   figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3class_Set' num2str(thisSet) '_AUROC_ben'];
%   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8.6 6.45])% ,[0 0 8 6])
%   set(gcf, 'renderer', 'painters');
   print(figname, '-dpng','-r400')
%    print(figname, '-depsc')
end

figure
bar([comparisonNegAve])%, [cornFlowerBlue;cornFlowerBlue;cornFlowerBlue;cornFlowerBlue])
alpha 0.6
map6Paired = brewermap(6,'Paired');
map12Paired = brewermap(12,'Paired');
map2Paired = brewermap(2,'Paired');
thisColorMap = map12Paired;
% colormap(thisColorMap)% map6Paired)


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTickLabel', Labels, ...'YLim', [0.5 1.001], ...  'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0.4:0.1:1, ...
  'LineWidth', 1)
ylabel('AUROC')
axis([0.5 8.5 0.5 1])
title(['Set : ' num2str(thisSet)])
grid on
title('Detecting deleterious sites')
xlabel('Generations')
% hold on
% errorbar(1:length(comparisonNegAve), comparisonNegAve, comparisonNegSTD, 'ko')

if(saveFigs)
   figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3class_Set' num2str(thisSet) '_AUROC_del'];
%   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8.6 6.45])% ,[0 0 8 6])
%   set(gcf, 'renderer', 'painters');
   print(figname, '-dpng','-r400')
%    print(figname, '-depsc')
end

drawBubblePlot;
thisColorMap = map2Paired;
mainDir = pwd;
cd('/local/staff/ee/mssohail/Matlab Codes/gramm-master')
drawViolinPlotTrajWork;
cd(mainDir)
%
if(saveFigs)
   figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3class_Set' num2str(thisSet) '_violin'];
%   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8.6 6.45])% ,[0 0 8 6])
%   set(gcf, 'renderer', 'painters');
   print(figname, '-dpng','-r400')
%    print(figname, '-depsc')
end

% ----------------------------------------------------------------------- %



