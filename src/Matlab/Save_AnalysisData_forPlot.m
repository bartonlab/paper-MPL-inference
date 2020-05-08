numItr = 100%90%90%250;

set(0,'DefaultAxesFontName','Arial')% Helvetica
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',7)
set(0,'DefaultTextFontSize',7)

init_ColorNames_Settings;
npg1 = color_scheme_npg(1,:);
npg4 = color_scheme_npg(4,:);
npg5 = color_scheme_npg(5,:);
npg7 = color_scheme_npg(7,:);
npg8 = color_scheme_npg(8,:);
set1_2 = color_scheme_set1(2,:);
% set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'DefaultTextFontName','CMU Serif Roman')
% set(0,'DefaultAxesFontSize',15)
% set(0,'DefaultTextFontSize',15)

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

saveFigs = 0;
saveFile = 0;
%Lin = 10;
%allSets = [55 25500002 2550002 255002];
%allSets = [55 35500001 3550001 355001];


allSets = [755 75500001 7550001 755001];
allTused = 1000;%[300 500 700 1000 1300];
numAllSets = length(allSets);
numAllTused = length(allTused);

meanAucLinkPosMtx = zeros(numAllSets, numAllTused);
meanAucUnLinkPosMtx = zeros(numAllSets, numAllTused);
meanAucLinkNegMtx = zeros(numAllSets, numAllTused);
meanAucUnLinkNegMtx = zeros(numAllSets, numAllTused);
meanSelcSitesMtx1 = zeros(numAllSets, numAllTused);
meanSelcSitesMtx5 = zeros(numAllSets, numAllTused);
meanSelcSitesMtx10 = zeros(numAllSets, numAllTused);
meanSelcSitesMtx25 = zeros(numAllSets, numAllTused);
meanSumOfAbsOfAllOffDiagTermsICM = zeros(numAllSets, numAllTused);
meanSumOfAbsOfAllOffDiagTermsICorrM = zeros(numAllSets, numAllTused);
StdOffDiagVecICorrM_IterWise = zeros(numAllSets, numItr);

stdOffDiagVecICorrM = zeros(numAllSets, 1);
recValVec = zeros(numAllSets, 1);
noiseThresh = 0.05;

for ss = 1:numAllSets
thisSet = allSets(ss);

fileNameContainingDirPath = 'dirNames_recomb.txt';
recombination = 0;
recVal = 0;
getSysParam_long;
%Tstart = 31;
Tstart = 1;

for tt = 1:numAllTused
Tused = allTused(tt);%1300;
dT = 10;
ng = 100;
step = 0.002;

edges = [-1 -0.1:step:0.1] + step/2;
if(thisSet == 330)
    step = 0.01;
    edges = [-1 -0.5:step:0.5] + step/2;
end


Tend = Tused + Tstart;

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

if(classesOfSites == 2)
    selTypeName = 'PosSel';
elseif(classesOfSites == 3)
    selTypeName = 'PosDelSel';
end
[~, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);

dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];



trueClassROC_AllItr = [];
sigmaEstOutLink_AllItr = [];
sigmaEstOutUnLink_AllItr = [];
sigmaEstOutUnLinkNoMu_AllItr = [];
sigmaEstOutLink_LinUni_AllItr = [];
sigmaEstOutLinkNoMu_AllItr = [];
perSiteSelction_AllItr  = [];
siteToExculeCuzOfLimitedPolyTps_AllItr = [];
aucLinkItr = zeros(numItr,2);
aucLinkNoMuItr = zeros(numItr,2);
aucUnLinkItr = zeros(numItr,2);
aucUnLinkNoMuItr = zeros(numItr,2);
aucLinkLinUniItr = zeros(numItr,2);
xLinkPos = zeros(numItr, 51);
yLinkPos = zeros(numItr, 51);
xLinkNeg = zeros(numItr, 51);
yLinkNeg = zeros(numItr, 51);
selcSitesItr1 = zeros(numItr,1);
selcSitesItr5 = zeros(numItr,1);
selcSitesItr10 = zeros(numItr,1);
selcSitesItr25 = zeros(numItr,1);

offDiagVecICM = zeros(numItr,1225);
offDiagVecICorrM = zeros(numItr,1225);


sumOfAbsOfAllOffDiagTermsICM = zeros(numItr,1);
sumOfDiagTermsICM = zeros(numItr,1);

sumOfAbsOfAllOffDiagTermsICorrM = zeros(numItr,1);
sumOfDiagTermsICorrM = zeros(numItr,1);
vecCovEntriesAll = zeros(1225*101*100, 1);
for itr = 1:numItr
    selcSites1 = [];
    selcSites5 = [];
    selcSites10 = [];
    selcSites25 = [];
    if(recombination == 0)
        if(Tstart == 1)
            fileName = ['WFsim_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        else        
            fileName = ['WFsim_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        end
    elseif(recombination == 1)
        if(Tstart == 1)
            fileName = ['WFsim_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        else        
            fileName = ['WFsim_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        end
    end
    load([dirNameAnalysis fileName], 'L', 'perSiteSelction', 'perSiteSelctionSelc', ...
             'sigmaEstOutLink', 'sigmaEstOutLink_LinUni', 'sigmaEstOutLinkNoMu', 'priorConst', ...
             'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc', 'sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'sigmaEstOutUnLink', 'sigmaEstOutUnLinkNoMu', 'q', 'priorConst', 'timeWholeCodeMPL', 'q_LinUni', ...
             'sumAijWithLink', 'dT', 'avgLDPerTime', 'vecCovEntries');

    % identify neutral, pos/neg selection
    trueClassROC = 3*ones(1, Lin);
    trueClassROC(perSiteSelction == posSelc) = 4;
    trueClassROC(perSiteSelction == negSelc) = 5;
    
    selcSites1 = ~(abs(max(q) - min(q)) <= 0.01 & (max(q) == 1 | min(q) == 0)); % ~ rejection criterion
    selcSites5 = ~(abs(max(q) - min(q)) <= 0.05 & (max(q) == 1 | min(q) == 0)); % ~ rejection criterion
    selcSites10 = ~(abs(max(q) - min(q)) <= 0.10 & (max(q) == 1 | min(q) == 0)); % ~ rejection criterion
    selcSites25 = ~(abs(max(q) - min(q)) <= 0.25 & (max(q) == 1 | min(q) == 0)); % ~ rejection criterion
    
    trueClassROC_AllItr = [trueClassROC_AllItr trueClassROC];
    sigmaEstOutLink_AllItr = [sigmaEstOutLink_AllItr sigmaEstOutLink'];
    sigmaEstOutUnLink_AllItr = [sigmaEstOutUnLink_AllItr sigmaEstOutUnLink'];
    sigmaEstOutUnLinkNoMu_AllItr = [sigmaEstOutUnLinkNoMu_AllItr sigmaEstOutUnLinkNoMu'];
    sigmaEstOutLink_LinUni_AllItr = [sigmaEstOutLink_LinUni_AllItr sigmaEstOutLink_LinUni'];
    sigmaEstOutLinkNoMu_AllItr = [sigmaEstOutLinkNoMu_AllItr sigmaEstOutLinkNoMu'];
    perSiteSelction_AllItr = [perSiteSelction_AllItr perSiteSelction];
        
    
    nrmse_3class_Link_itrPos(itr) = sqrt(sum(abs(perSiteSelction(1:DAll) - sigmaEstOutLink(1:DAll)').^2)/sum(abs(perSiteSelction(1:DAll)).^2));
    nrmse_3class_LinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(1:DAll) - sigmaEstOutLinkNoMu(1:DAll)').^2)/sum(abs(perSiteSelction(1:DAll)).^2));
    nrmse_LinUni_3class_itrPos(itr) = sqrt(sum(abs(perSiteSelction(1:DAll) - sigmaEstOutLink_LinUni(1:DAll)').^2)/sum(abs(perSiteSelction(1:DAll)).^2));
    nrmse_3class_UnLink_itrPos(itr) = sqrt(sum(abs(perSiteSelction(1:DAll) - sigmaEstOutUnLink(1:DAll)').^2)/sum(abs(perSiteSelction(1:DAll)).^2));
    nrmse_3class_UnLinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(1:DAll) - sigmaEstOutUnLinkNoMu(1:DAll)').^2)/sum(abs(perSiteSelction(1:DAll)).^2));

    nrmse_3class_Link_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(DAll+1:DeleAll+DAll) - sigmaEstOutLink(DAll+1:DeleAll+DAll)').^2)/sum(abs(perSiteSelction(DAll+1:DeleAll+DAll)).^2));
    nrmse_3class_LinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(DAll+1:DeleAll+DAll) - sigmaEstOutLinkNoMu(DAll+1:DeleAll+DAll)').^2)/sum(abs(perSiteSelction(DAll+1:DeleAll+DAll)).^2));
    nrmse_LinUni_3class_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(DAll+1:DeleAll+DAll) - sigmaEstOutLink_LinUni(DAll+1:DeleAll+DAll)').^2)/sum(abs(perSiteSelction(DAll+1:DeleAll+DAll)).^2));
    nrmse_3class_UnLink_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(DAll+1:DeleAll+DAll) - sigmaEstOutUnLink(DAll+1:DeleAll+DAll)').^2)/sum(abs(perSiteSelction(DAll+1:DeleAll+DAll)).^2));
    nrmse_3class_UnLinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(DAll+1:DeleAll+DAll) - sigmaEstOutUnLinkNoMu(DAll+1:DeleAll+DAll)').^2)/sum(abs(perSiteSelction(DAll+1:DeleAll+DAll)).^2));

    nrmse_3class_Link_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(DeleAll+DAll+1:end) - sigmaEstOutLink(DeleAll+DAll+1:end)').^2));
    nrmse_3class_LinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(DeleAll+DAll+1:end) - sigmaEstOutLinkNoMu(DeleAll+DAll+1:end)').^2));
    nrmse_LinUni_3class_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(DeleAll+DAll+1:end) - sigmaEstOutLink_LinUni(DeleAll+DAll+1:end)').^2));
    nrmse_3class_UnLink_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(DeleAll+DAll+1:end) - sigmaEstOutUnLink(DeleAll+DAll+1:end)').^2));
    nrmse_3class_UnLinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(DeleAll+DAll+1:end) - sigmaEstOutUnLinkNoMu(DeleAll+DAll+1:end)').^2));

    posOnlyItrTemp = perSiteSelction == max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction == min(perSiteSelction);

    [xLinkPosTemp1010,yLinkPosTemp1010,~, aucLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink, 1);
    [xLinkNegTemp1010,yLinkNegTemp1010,~, aucLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink, 1);


    xLinkPos(itr,:) = [zeros(1, Lin + 1 - length(xLinkPosTemp1010)) xLinkPosTemp1010'];
    yLinkPos(itr,:) = [zeros(1, Lin + 1 - length(yLinkPosTemp1010)) yLinkPosTemp1010'];
    xLinkNeg(itr,:) = [zeros(1, Lin + 1 - length(xLinkNegTemp1010)) xLinkNegTemp1010'];
    yLinkNeg(itr,:) = [zeros(1, Lin + 1 - length(yLinkNegTemp1010)) yLinkNegTemp1010'];
    %------------------------------------------------------------------------------

    [~,~,~, aucLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLinkNoMu, 1);
    [~,~,~, aucLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLinkNoMu, 1);

    [~,~,~, aucUnLinkItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLink, 1);
    [~,~,~, aucUnLinkItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLink, 1);

    [~,~,~, aucUnLinkNoMuItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);
    [~,~,~, aucUnLinkNoMuItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutUnLinkNoMu, 1);

    [~,~,~, aucLinkLinUniItrTemp(1)] = perfcurve(double(posOnlyItrTemp), sigmaEstOutLink_LinUni, 1);
    [~,~,~, aucLinkLinUniItrTemp(2)] = perfcurve(double(~negOnlyItrTemp), sigmaEstOutLink_LinUni, 1);

    aucLinkItr(itr,:) = aucLinkItrTemp;
    aucLinkNoMuItr(itr,:) = aucLinkNoMuItrTemp;
    aucUnLinkItr(itr,:) = aucUnLinkItrTemp;
    aucUnLinkNoMuItr(itr,:) = aucUnLinkNoMuItrTemp;
    aucLinkLinUniItr(itr,:) = aucLinkLinUniItrTemp;
    selcSitesItr1(itr) = sum(selcSites1);
    selcSitesItr5(itr) = sum(selcSites5);
    selcSitesItr10(itr) = sum(selcSites10);
    selcSitesItr25(itr) = sum(selcSites25);
    
    
    sumAijWithLinkdT = sumAijWithLink*dT;
    sumOfAbsOfAllOffDiagTermsICM(itr) = sum(sum(abs((sumAijWithLinkdT - diag(diag(sumAijWithLinkdT))))));
    sumOfDiagTermsICM(itr) = sum(diag(sumAijWithLinkdT));
    
    Btemp = triu(sumAijWithLinkdT,1);
    Ctemp = Btemp(Btemp ~= 0);
    numZeroPadding = 1225 - length(Ctemp);
    if(numZeroPadding > 0)
        Ctemp = [Ctemp; zeros(numZeroPadding,1)];
    end
    offDiagVecICM(itr,:) = Ctemp;
    
    sumAijWithLinkdTCorr = corrcov(sumAijWithLinkdT);
    sumOfAbsOfAllOffDiagTermsICorrM(itr) = sum(sum(abs((sumAijWithLinkdTCorr - diag(diag(sumAijWithLinkdTCorr))))));
    sumOfDiagTermsICorrM(itr) = sum(diag(sumAijWithLinkdTCorr));
    
    Btemp = triu(sumAijWithLinkdTCorr,1);
    Ctemp = Btemp(Btemp ~= 0);
    numZeroPadding = 1225 - length(Ctemp);
    if(numZeroPadding > 0)
        Ctemp = [Ctemp; zeros(numZeroPadding,1)];
    end
    offDiagVecICorrM(itr,:) = Ctemp;
    
    StdOffDiagVecICorrM_IterWise(ss,itr) = std(offDiagVecICorrM(itr,:));
    
    
    avgLDPerTimeItr(itr,:) = avgLDPerTime;
    
    oneTirElements = 1225*101;
    vecCovEntriesAll(oneTirElements*(itr-1)+1:oneTirElements*itr) = vecCovEntries(:);
end
avgLDPerTimeItrSet{ss} = avgLDPerTimeItr;
if(ss == 1)
    useItr = 15;%10;
elseif(ss == 2)
    useItr = 15;%20;
elseif(ss == 3)
    useItr = 40;
elseif(ss == 4)
    useItr = 50;
end
% figure
% h = histogram(abs(offDiagVecICorrM(useItr,:)), 0:0.05:1)
% axis([0 1 0 500])


figure
h = histogram((offDiagVecICM), -165:10:165)
maxVal(ss) = max(max(abs(offDiagVecICM)))
minVal(ss) = min(min((offDiagVecICM)))
set(gca, ...
    'YScale', 'log')
%h = histogram((offDiagVecICorrM), -1:0.01:1);
histValues(ss,:) = h.Values;
histEdges(ss,:) = h.BinEdges(2:end);

% figure
% hSepTP = histogram((vecCovEntriesAll), -0.25:0.005:0.25)
% maxValSepTP(ss) = max(max((vecCovEntriesAll)))
% minValSepTP(ss) = min(min((vecCovEntriesAll)))
% 
% vecCovEntriesAllCell{ss} = vecCovEntriesAll;
% set(gca, ...
%     'YScale', 'log')
% %h = histogram((offDiagVecICorrM), -1:0.01:1);
% histValuesSepTP(ss,:) = hSepTP.Values;
% histEdgesSepTP(ss,:) = hSepTP.BinEdges(2:end);





% stdOffDiagVecICorrM(ss) = std(offDiagVecICorrM(:));
% meanOffDiagVecICorrM(ss) = mean(abs(offDiagVecICorrM(:)));
% medianOffDiagVecICorrM(ss) = median(abs(offDiagVecICorrM(:)));

% set(gca,...
%     'YScale', 'log')
% axis([-1 1 0 1e5])
meanAucLinkTemp = mean(aucLinkItr);
meanAucUnLinkTemp = mean(aucUnLinkItr);
meanAucLinkPosMtx(ss,tt) = meanAucLinkTemp(1);
meanAucUnLinkPosMtx(ss,tt) = meanAucUnLinkTemp(1);
meanAucLinkNegMtx(ss,tt) = meanAucLinkTemp(2);
meanAucUnLinkNegMtx(ss,tt) = meanAucUnLinkTemp(2);
meanSelcSitesMtx1(ss,tt) = mean(selcSitesItr1);
meanSelcSitesMtx5(ss,tt) = mean(selcSitesItr5);
meanSelcSitesMtx10(ss,tt) = mean(selcSitesItr10);
meanSelcSitesMtx25(ss,tt) = mean(selcSitesItr25);
meanSumOfAbsOfAllOffDiagTermsICM(ss,tt) = mean(sumOfAbsOfAllOffDiagTermsICM);
meanSumOfAbsOfAllOffDiagTermsICorrM(ss,tt) = mean(sumOfAbsOfAllOffDiagTermsICorrM);
% meanAucLinkSets(ss,:) = mean(aucLinkItr);
% meanAucUnLinkSets(ss,:) = mean(aucUnLinkItr);

thisGroup = 1;
Std_only_si_Ben(ss,tt) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)));
thisGroup = 2;
Std_only_si_Del(ss,tt) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)));

thisGroup = 1;
SE_only_si_Ben(ss,tt) = Std_only_si_Ben(ss,tt)/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));
thisGroup = 2;
SE_only_si_Del(ss,tt) = Std_only_si_Del(ss,tt)/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));

end
recValVec(ss) = recVal;

fileNameRecombPlot = ['Set' num2str(thisSet) '_AUC_COVMtx.csv'];

if(exist(fileNameRecombPlot, 'file') == 2)
    delete(fileNameRecombPlot)
end
csvwrite(fileNameRecombPlot, [aucLinkItr offDiagVecICM])
end