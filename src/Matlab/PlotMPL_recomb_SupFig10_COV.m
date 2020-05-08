%%


clc
close all
clear all

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


allSets = [755 75500001 7550001 755001];

numAllSets = length(allSets);

meanAucLinkPosMtx = zeros(numAllSets, 1);
meanAucUnLinkPosMtx = zeros(numAllSets, 1);
meanAucLinkNegMtx = zeros(numAllSets, 1);
meanAucUnLinkNegMtx = zeros(numAllSets, 1);
meanSelcSitesMtx1 = zeros(numAllSets, 1);
meanSelcSitesMtx5 = zeros(numAllSets, 1);
meanSelcSitesMtx10 = zeros(numAllSets, 1);
meanSelcSitesMtx25 = zeros(numAllSets, 1);
meanSumOfAbsOfAllOffDiagTermsICM = zeros(numAllSets, 1);
meanSumOfAbsOfAllOffDiagTermsICorrM = zeros(numAllSets, 1);
StdOffDiagVecICorrM_IterWise = zeros(numAllSets, numItr);

stdOffDiagVecICorrM = zeros(numAllSets, 1);
recValVec = zeros(numAllSets, 1);
noiseThresh = 0.05;

for ss = 1:numAllSets
thisSet = allSets(ss);
fileNameRecombPlot = ['Set' num2str(thisSet) '_AUC_COVMtx.csv'];

recombination = 0;
recVal = 0;
getSysParam_long;

Tstart = 1;


Tused = 1000;
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

% if(classesOfSites == 2)
%     selTypeName = 'PosSel';
% elseif(classesOfSites == 3)
%     selTypeName = 'PosDelSel';
% end
% [~, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
% 
%dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];



         
temp = csvread(fileNameRecombPlot);
pause(1)
aucLinkItr = temp(:,[1 2]);
offDiagVecICM  = temp(:,3:end);
clear temp

figure
h = histogram((offDiagVecICM), -165:10:165)
maxVal(ss) = max(max(abs(offDiagVecICM)))
minVal(ss) = min(min((offDiagVecICM)))
set(gca, ...
    'YScale', 'log')
%h = histogram((offDiagVecICorrM), -1:0.01:1);
histValues(ss,:) = h.Values;
histEdges(ss,:) = h.BinEdges(2:end);


meanAucLinkTemp = mean(aucLinkItr);
%meanAucUnLinkTemp = mean(aucUnLinkItr);
meanAucLinkPosMtx(ss) = meanAucLinkTemp(1);
%meanAucUnLinkPosMtx(ss) = meanAucUnLinkTemp(1);
meanAucLinkNegMtx(ss) = meanAucLinkTemp(2);
%meanAucUnLinkNegMtx(ss) = meanAucUnLinkTemp(2);
%meanSelcSitesMtx1(ss) = mean(selcSitesItr1);
%meanSelcSitesMtx5(ss) = mean(selcSitesItr5);
%meanSelcSitesMtx10(ss) = mean(selcSitesItr10);
%meanSelcSitesMtx25(ss) = mean(selcSitesItr25);
%meanSumOfAbsOfAllOffDiagTermsICM(ss) = mean(sumOfAbsOfAllOffDiagTermsICM);
%meanSumOfAbsOfAllOffDiagTermsICorrM(ss) = mean(sumOfAbsOfAllOffDiagTermsICorrM);
% meanAucLinkSets(ss,:) = mean(aucLinkItr);
% meanAucUnLinkSets(ss,:) = mean(aucUnLinkItr);

thisGroup = 1;
Std_only_si_Ben(ss) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)));
thisGroup = 2;
Std_only_si_Del(ss) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)));

thisGroup = 1;
SE_only_si_Ben(ss) = Std_only_si_Ben(ss)/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));
thisGroup = 2;
SE_only_si_Del(ss) = Std_only_si_Del(ss)/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));

recValVec(ss) = recVal;


end








%%
close all
xDim = 17.4;
yDim = 7-1.83;


fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.5;
rightMargin = 0.25;
bottomMargin = 1;
topMargin = 0.5;
hgap1 = 0.25;
hgap2 = 2;
height1 = (yDim - bottomMargin - topMargin);
%height2 = (yDim - bottomMargin - topMargin);
width1 = (xDim - leftMargin - rightMargin - 1*hgap1 - hgap2)/6;
width2 = width1*4;


ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1+width1+hgap2 bottomMargin width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
%%
color_scheme1 = brewermap(100,'Blues');
color_scheme3 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;

for k = 1:2
    if(k == 1)
        barData = meanAucLinkPosMtx;
        errBarData = Std_only_si_Ben;
        yLabelStr = 'Mean AUROC';% (beneficial)';
        titleStr = 'Beneficial';
    elseif(k == 2)
        barData = meanAucLinkNegMtx;
        errBarData = Std_only_si_Del;
        yLabelStr = ' ';%'Mean AUROC (deleterious)';
        titleStr = 'Deleterious'
    end
    axes(ha(k))
    bb = bar(barData) % this is only so axes is set by bar plot and not errorbar plot
    hold on
    xCord = [1 2 3 4];
    yCord = barData;
    errorbar(xCord, yCord, errBarData, 'k', 'LineStyle', 'none', 'CapSize', 0)
    bb(1).FaceColor = color_scheme31(50,:);
    bb(1).BarWidth = 0.6;

    %colormap(myColorMap3(2,:));
    %xTickLabelTemp = flip(numStrainsInInitialPopAll);
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
      'XTickLabel'  , {'0' '10^{-5}' '10^{-4}' '10^{-3}'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
      'YTick'  , 0.5:0.1:1, ...
      'LineWidth', 0.5, ...
      'FontSize', 6)
    if(k == 2)
       set(gca, ...
           'YTickLabel', ' ')
        xlab = xlabel('Recombination probability, r')
        set(xlab,'Units','normalized');
        set(xlab,'position',get(xlab,'position') + [-0.6 0 0]);
    end
    axis([0.5 4.5 0.5 1])
    ylabel(yLabelStr)
    
    title(titleStr)
end



%save RecombMPL

colorList{1} = color_scheme11(100,:);
colorList{2} = color_scheme11(80,:);
colorList{3} = color_scheme11(60,:);
colorList{4} = color_scheme11(40,:);

axes(ha(3))
histValues(histValues == 0) = 0.1;
for ss = 4:-1:1
    semilogy(histEdges, (histValues(ss,:)), 'color', colorList{ss}, 'LineWidth', 2)
    if(ss == 4)
        hold on
        axis([-160 160 0.5 1e6])
    end
end
xlabel('Magnitude of absolute off-diagonal entries of integrated covariance matrix')
ylabel('Count of pairs')
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...'XTickLabel'  , Nin*recValVec, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
      'YTick'  , [1 1e2 1e4 1e6], ...
      'LineWidth', 0.5, ...
      'FontSize', 6)
pause(1)

dimDummy = [0.1 0.1 0.1 0.1];
line1 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color',colorList{1})
line1.Units = 'centimeter';
line1.X = [8.45 9.15]+6;
line1.Y = (yDim-topMargin-0.35)*ones(1,2);
line1.LineWidth = 2;

line2 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color',colorList{2})
line2.Units = 'centimeter';
line2.X = [8.45 9.15]+6;
line2.Y = (yDim-topMargin-0.65-0.01)*ones(1,2);
line2.LineWidth = 2;

line3 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color',colorList{3})
line3.Units = 'centimeter';
line3.X = [8.45 9.15]+6;
line3.Y = (yDim-topMargin-0.95-0.02)*ones(1,2);
line3.LineWidth = 2;

line4 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color',colorList{4})
line4.Units = 'centimeter';
line4.X = [8.45 9.15]+6;
line4.Y = (yDim-topMargin-1.25-0.03)*ones(1,2);
line4.LineWidth = 2;

textLeg1 = annotation('textbox',dimDummy,'String','r = 0','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [9.2+6 yDim-topMargin-0.6 5 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontSize = 6;

textLeg2 = annotation('textbox',dimDummy,'String','r = 10^{-5}','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [9.2+6 yDim-topMargin-0.9 5 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontSize = 6;

textLeg3 = annotation('textbox',dimDummy,'String','r = 10^{-4}','FitBoxToText','on')
textLeg3.Units = 'centimeter';
textLeg3.Position = [9.2+6 yDim-topMargin-1.2 5 0.5];
textLeg3.LineStyle = 'none';
textLeg3.FontSize = 6;

textLeg4 = annotation('textbox',dimDummy,'String','r = 10^{-3}','FitBoxToText','on')
textLeg4.Units = 'centimeter';
textLeg4.Position = [9.2+6 yDim-topMargin-1.5 5 0.5];
textLeg4.LineStyle = 'none';
textLeg4.FontSize = 6;

pause(2)

textA = annotation('textbox',dimDummy,'String','a','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25 yDim-topMargin 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;

textB = annotation('textbox',dimDummy,'String','b','FitBoxToText','on')
textB.Units = 'centimeter';
textB.Position = [0.15+2*width1+hgap1+hgap2 yDim-topMargin 5 0.5];
textB.LineStyle = 'none';
textB.FontWeight = 'bold';
textB.FontSize = 9;




%%

pause(2)
saveFigs = 1;
if(saveFigs == 1)
   figname = ['figs-sim-recombination_COV'];
   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])% ,[0 0 8 6])
   set(gcf, 'renderer', 'painters');
   print(figname, '-dpng','-r400')
   set(gcf, 'PaperSize', [xDim yDim])
   print(figname, '-dpdf', '-fillpage')
   %print(figname, '-depsc', '-r400')
end

