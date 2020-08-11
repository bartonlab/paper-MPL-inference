clc
clear all
%close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

saveFigs = 1;
itrs = 1000000;
maxBinEdge = 100;
xBinomial = binornd(1000,[0.0139*ones(1, itrs)]);
fig0 = figure;
h1 = histogram(xBinomial,0:1:maxBinEdge);
valsBinomial = h1.Values;
binWidthBinomial = h1.BinWidth;
binEdgesBinomial = h1.BinEdges;
close(fig0);
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;


xDim = 17.4;
yDim = 2*3.67+1.5 + 1.5 +1.6;


fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.5;
rightMargin = 0.5;
bottomMargin = 1+1.6;
topMargin = 0.5;
hgap1 = 1.75;%1.5;
hgap2 = 1.5;
hgap2 = 2;
vgap1 = 1.5;
height1 = (yDim - bottomMargin - topMargin - vgap1)/2;
width1 = (xDim - leftMargin - rightMargin - 2*hgap1)/3;
width2 = (xDim - leftMargin - rightMargin -hgap2)/2;

ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+height1+vgap1 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+height1+vgap1 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1) bottomMargin+height1+vgap1 width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');

            
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+width2+hgap2 bottomMargin width2 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');



            
valsNewTemp = valsBinomial/sum(valsBinomial*binWidthBinomial);
binEdgesNew = (binEdgesBinomial(1:end-1) + binEdgesBinomial(2:end))/2;
valsNew = valsNewTemp/sum(binEdgesNew.*valsNewTemp);
xSpec = [binEdgesNew(1) binEdgesNew binEdgesNew(end)];
ySpec = [0 valsNew 0];
axes(ha(1))
plot(binEdgesNew, valsNew, 'color', color_scheme11(20,:))
fill(xSpec, ySpec, color_scheme11(80,:))
axis([0 40 0 0.008])
ylabel('Density')
xlabel('Number of sequences per time point, n_s')
%title('Binomial with n = 1000, p = 0.0139')

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'YTick', 1e-3*(0:2:8), ...
  'YTickLabel'  , {'0' '0.002' '0.004' '0.006' '0.008'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...'YTick'  , 0.5:0.1:1, ...'YTickLabel'  , ' ', ...'YTick'  , 0:0.2:1, ...
  'LineWidth', 0.5, ...
  'FontSize', 6)

%-----------------------------------------------------------
% delta t samp dist
itrs = 1000000;
maxBinEdge_dT = 200;
k = 3.5;
theta = 8.4;
x1 = gamrnd(k*ones(1, 0.87*itrs),theta);
k = 3;
theta = 2;
const = 120;
gamma_mean = (k*theta) + const;
xTemp = gamrnd(k*ones(1, 0.13*itrs),theta) + const;
xTempUp = xTemp(xTemp > gamma_mean);
xTempLow = xTemp(xTemp <= gamma_mean);

x2Low = xTempUp - 2*(xTempUp - gamma_mean);
x2Up = xTempLow + 2*(gamma_mean - xTempLow);

x2 = [x2Up x2Low];
x = [x1 x2];
fig2 = figure
h2 = histogram(x,0:1:maxBinEdge_dT);



vals_dT = h2.Values
binWidth_dT = h2.BinWidth
binEdges_dT = h2.BinEdges;
close(fig2);
color_scheme1 = brewermap(100,'Blues');
color_scheme3 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;

valsNew_dT = vals_dT/sum(vals_dT*binWidth_dT);
binEdgesNew_dT = (binEdges_dT(1:maxBinEdge_dT) + binEdges_dT(2:(maxBinEdge_dT+1)))/2;
axes(ha(2))
plot(binEdgesNew_dT, valsNew_dT, 'color', color_scheme11(20,:))
xSpec_dT = [binEdgesNew_dT(1) binEdgesNew_dT binEdgesNew_dT(end)];
ySpec_dT = [0 valsNew_dT 0];
fill(xSpec_dT, ySpec_dT, color_scheme11(80,:))

ylabel('Density')
xlabel('Time sampling step, \Deltat_k')
%title(['Gamma with k = ' num2str(k) ', \theta = ' num2str(theta)])
%title(['Mixture of two Gamma distribution'])
axis([0 150 0 0.03])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...'XTickLabel'  , {'0' '10^{-5}' '10^{-4}' '10^{-3}'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...'YTick'  , 0.5:0.1:1, ...'YTickLabel'  , ' ', ...'YTick'  , 0:0.2:1, ...
  'LineWidth', 0.5, ...
  'FontSize', 6)






%------------------------------------------------------------------------
% plot Trajectory length
itrs = 1000000*5;
itrs_A = round(3/14*itrs);
itrs_B = round(11/14*itrs);
maxBinEdgeTused = 350;
k3 = 5.5%2.8;
theta3 = 7.2%8;
const3_1 = 5;
x3 = gamrnd(k3*ones(1, itrs_A),theta3) + const3_1;

k3 = 15;
theta3 = 8;
const3 = 143;
gamma_mean3 = (k3*theta3) + const3;
xTemp = gamrnd(k3*ones(1, itrs_B),theta3) + const3;
xTempUp = xTemp(xTemp > gamma_mean3);
xTempLow = xTemp(xTemp <= gamma_mean3);
% 
x4Low = xTempUp - 2*(xTempUp - gamma_mean3);
x4Up = xTempLow + 2*(gamma_mean3 - xTempLow);
% 
x4 = [x4Up x4Low];

x = round([x3 x4]);
fig3 = figure
h3 = histogram(x,0:1:maxBinEdgeTused);



valsTused = h3.Values
binWidthTused = h3.BinWidth
binEdgesTused = h3.BinEdges;
close(fig3);
color_scheme1 = brewermap(100,'Blues');
color_scheme3 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;

valsNewTused = valsTused/sum(valsTused*binWidthTused);
binEdgesNewTused = (binEdgesTused(1:maxBinEdgeTused) + binEdgesTused(2:(maxBinEdgeTused+1)))/2;

% plot sampling dist T
axes(ha(3))
plot(binEdgesNewTused, valsNewTused, 'color', color_scheme11(20,:))
xSpecTused = [binEdgesNewTused(1) binEdgesNewTused binEdgesNewTused(end)];
ySpecTused = [0 valsNewTused 0];
fill(xSpecTused, ySpecTused, color_scheme11(80,:))

ylabel('Density')
xlabel('Generations used for inference, T')
%title(['Gamma with k = ' num2str(k) ', \theta = ' num2str(theta)])
%title(['Mixture of two Gamma distribution'])
%axis([0 100 0 0.015])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
  'XTick', 0:100:400, ...'XTickLabel'  , {'0' '10^{-5}' '10^{-4}' '10^{-3}'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...'YTick'  , 0.5:0.1:1, ...'YTickLabel'  , ' ', ...'YTick'  , 0:0.2:1, ...
  'LineWidth', 0.5, ...
  'FontSize', 6)

%---------------------------------------------------
%%
pat = 6%9%6%50;

thisSet = 9860001%[8550001 9550001 9650001 9750001 9850001 9860001 9950001];
numStrainsInInitialPop = 5%[1 5 5 5 5 5 5];
%Tused = 300;
tjSamplingSchemeStr = 'scheme1';
ngSamplingSchemeStr = 'schemeD';
dTSamplingSchemeStr = 'scheme33';

setConvention = 3;

priorConstSC = 5; % this is the strength of the SC regularization term

textCell{1} = ['dirNamesSet' num2str(thisSet) '_ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_'];
[dirNameDataTemp, dirNameAnalysisTemp, dirNameResultsTemp] = setDirNamesMPLPipeline('Set_dirNames_MPL_SimData1.txt');

ResultFolder = [dirNameResultsTemp 'Set' num2str(thisSet) '/ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];
if(setConvention == 1)
    FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
    FLAG_linearInt = false;
elseif(setConvention == 2)
    FLAG_stratonovich = true;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
    FLAG_linearInt = false;
elseif(setConvention == 3)
    FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
    FLAG_linearInt = true;
end

% Use stratonovich for increased robustness to sampling effects
if(FLAG_stratonovich == true)
    convention = 'Stratonovich';
elseif(FLAG_stratonovich == false && FLAG_linearInt == false)
    convention = 'Ito';
elseif(FLAG_stratonovich == false && FLAG_linearInt == true)
    convention = 'LinearInter';
end

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
fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, textCell);

fileNameContainingDirPath = fileNamesListThisDir{pat};
indOfDash = strfind(fileNameContainingDirPath, '_');
indOfDot = strfind(fileNameContainingDirPath, '.');
patID = fileNameContainingDirPath(indOfDash(end-1)+1:indOfDash(end)-1);
thisProt = fileNameContainingDirPath(indOfDash(end)+1:indOfDot(end)-1);
fileNameContainingDirPath = [dirNameStr1Files fileNamesListThisDir{pat}];
[dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);




                     
fileNameAllTrajs = ['AllTrajsWithTimeInfo_' patID '_' thisProt '.txt'];
dataInTemp = dlmread([dirNameAnalysis 'Estimates' chosenSlash fileNameAllTrajs]);
samplingTimePoints = dataInTemp(:,1);
q = dataInTemp(:,2:end);

axes(ha(4))


% %------ solid lines coloured based on type, dashed lines for sampling time -----
plot(samplingTimePoints, q(:,21:50), '-', 'color', color_scheme31(50,:), 'LineWidth', 1)
hold on
plot(samplingTimePoints, q(:,1:10), '-', 'color', color_scheme21(99,:), 'LineWidth', 1)
plot(samplingTimePoints, q(:,11:20), '-', 'color', color_scheme11(99,:), 'LineWidth', 1)
for i = 2:length(samplingTimePoints)
    plot(samplingTimePoints(i)*ones(1, 11), 0:0.1:1, '--', 'color', color_scheme31(40,:), 'LineWidth', 0.5)
end

xlabel('Generation')
ylabel('Frequency')
axis([0 samplingTimePoints(end) 0 1])
%title('Heterogeneous sampled trajectories')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...'XTickLabel'  , {'0' '10^{-5}' '10^{-4}' '10^{-3}'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...'YTick'  , 0.5:0.1:1, ...'YTickLabel'  , ' ', ...
  'YTick'  , 0:0.2:1, ...
  'LineWidth', 0.5, ...
  'FontSize', 6)

dimDummy = [0.1 0.1 0.1 0.1];
line1 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color', color_scheme21(99,:))
line1.Units = 'centimeter';
line1.X = leftMargin + [0.5 1.5] ;
line1.Y = (1.5)*ones(1,2);
line1.LineWidth = 2;

line2 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color', color_scheme11(99,:))
line2.Units = 'centimeter';
line2.X = leftMargin + [0.5 1.5] ;
line2.Y = (1.1)*ones(1,2);
line2.LineWidth = 2;

line3 = annotation('line',[0.1 0.1], [0.1 0.1], 'Color', color_scheme31(50,:))
line3.Units = 'centimeter';
line3.X = leftMargin + [0.5 1.5] ;
line3.Y = (0.7)*ones(1,2);
line3.LineWidth = 2;

textLeg01 = annotation('textbox',dimDummy,'String','Beneficial','FitBoxToText','on')
textLeg01.Units = 'centimeter';
textLeg01.Position = [leftMargin+1.5 1.3 4 0.5];
textLeg01.LineStyle = 'none';
textLeg01.FontSize = 8;

textLeg02 = annotation('textbox',dimDummy,'String','Deleterious','FitBoxToText','on')
textLeg02.Units = 'centimeter';
textLeg02.Position = [leftMargin+1.5 0.9 4 0.5];
textLeg02.LineStyle = 'none';
textLeg02.FontSize = 8;

textLeg03 = annotation('textbox',dimDummy,'String','Neutral','FitBoxToText','on')
textLeg03.Units = 'centimeter';
textLeg03.Position = [leftMargin+1.5 0.5 4 0.5];
textLeg03.LineStyle = 'none';
textLeg03.FontSize = 8;



%-------------------------------------------------------------------



%========================== INITIALIZATION ================================

% -------------------------- User specified -------------------------------
% dataSet
FLAG_loadDiffSelcSites = true; % if set, it will pick selcSites of 'perfect case' based on heterogenous sampling
setAll = [9860001];
strainsAll = [5];

all_TusedCode = [1 1];
%unEvenSamplingAll = [false true]; 
ngSamplingSchemeStr = 'schemeD';
dTSamplingSchemeStr = 'scheme33';
% -- ngSamplingSchemeStr --
% schemeA: pBino = 0.008
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
numAllSets = length(setAll);
numAllTused = length(all_TusedCode);

dTStep = 1; % for good sampling case
ng = 1000; % for good sampling case

meanAucBen_si_MPL = -1*ones(numAllSets, numAllTused);
meanAucDel_si_MPL = -1*ones(numAllSets, numAllTused);
meanAucBen_si_SL = -1*ones(numAllSets, numAllTused);
meanAucDel_si_SL = -1*ones(numAllSets, numAllTused);
meanAucBen_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
meanAucDel_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
meanAucBen_si_SL_unEven = -1*ones(numAllSets, numAllTused);
meanAucDel_si_SL_unEven = -1*ones(numAllSets, numAllTused);

meanAucBenSelc_si_MPL = -1*ones(numAllSets, numAllTused);
meanAucDelSelc_si_MPL = -1*ones(numAllSets, numAllTused);
meanAucBenSelc_si_SL = -1*ones(numAllSets, numAllTused);
meanAucDelSelc_si_SL = -1*ones(numAllSets, numAllTused);
meanAucBenSelc_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
meanAucDelSelc_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
meanAucBenSelc_si_SL_unEven = -1*ones(numAllSets, numAllTused);
meanAucDelSelc_si_SL_unEven = -1*ones(numAllSets, numAllTused);

stdAucBen_si_MPL = -1*ones(numAllSets, numAllTused);
stdAucDel_si_MPL = -1*ones(numAllSets, numAllTused);
stdAucBen_si_SL = -1*ones(numAllSets, numAllTused);
stdAucDel_si_SL = -1*ones(numAllSets, numAllTused);
stdAucBen_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
stdAucDel_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
stdAucBen_si_SL_unEven = -1*ones(numAllSets, numAllTused);
stdAucDel_si_SL_unEven = -1*ones(numAllSets, numAllTused);

stdAucBenSelc_si_MPL = -1*ones(numAllSets, numAllTused);
stdAucDelSelc_si_MPL = -1*ones(numAllSets, numAllTused);
stdAucBenSelc_si_SL = -1*ones(numAllSets, numAllTused);
stdAucDelSelc_si_SL = -1*ones(numAllSets, numAllTused);
stdAucBenSelc_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
stdAucDelSelc_si_MPL_unEven = -1*ones(numAllSets, numAllTused);
stdAucBenSelc_si_SL_unEven = -1*ones(numAllSets, numAllTused);
stdAucDelSelc_si_SL_unEven = -1*ones(numAllSets, numAllTused);

dataLinkPosCell = cell(numAllSets, numAllTused);
dataLinkNegCell = cell(numAllSets, numAllTused);
dataUnLinkPosCell = cell(numAllSets, numAllTused);
dataUnLinkNegCell = cell(numAllSets, numAllTused);

dataLinkSelcPosCell = cell(numAllSets, numAllTused);
dataLinkSelcNegCell = cell(numAllSets, numAllTused);
dataUnLinkSelcPosCell = cell(numAllSets, numAllTused);
dataUnLinkSelcNegCell = cell(numAllSets, numAllTused);


dataLinkPosCell_unEven = cell(numAllSets, numAllTused);
dataLinkNegCell_unEven = cell(numAllSets, numAllTused);
dataUnLinkPosCell_unEven = cell(numAllSets, numAllTused);
dataUnLinkNegCell_unEven = cell(numAllSets, numAllTused);

dataLinkSelcPosCell_unEven = cell(numAllSets, numAllTused);
dataLinkSelcNegCell_unEven = cell(numAllSets, numAllTused);
dataUnLinkSelcPosCell_unEven = cell(numAllSets, numAllTused);
dataUnLinkSelcNegCell_unEven = cell(numAllSets, numAllTused);


ngVec_mean = zeros(numAllSets, numAllTused);
dTVec_mean = zeros(numAllSets, numAllTused);
ngVec_unEven_mean = zeros(numAllSets, numAllTused);
dTVec_unEven_mean = zeros(numAllSets, numAllTused);
trajLen_unEven_mean = zeros(numAllSets, numAllTused);
for ss = 1:numAllSets
    thisSet = setAll(ss)

    %unEvenSampling = unEvenSamplingAll(ss);
    numStrainsInInitialPop = strainsAll(ss);    
    
    % chose convention 1: Ito, 2: Stratonovich, 3: Linear interpolation
    % Stratonovich has increased robustness to sampling effects than Ito
    % Linear interpolation has most increased robustness to sampling effects
    setConvention = 3;
    
    priorConstSC = 5; % this is the strength of the SC regularization term
    
    for tt = 1:numAllTused
        if(all_TusedCode(tt) == 1)
            tjSamplingSchemeStr = 'scheme1';
        elseif(all_TusedCode(tt) == 2)
            tjSamplingSchemeStr = 'scheme2';
        end
        %Tused = TusedAll(tt);
                
        textCell{1} = ['dirNamesSet' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dTStep) '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_' ];
        textCell_unEven{1} = ['dirNamesSet' num2str(thisSet) '_ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '_'];
        
        FLAG_SaveIntovMtx = false; % SET: will save Integrated Covariance matrix (for debugging only)

        FLAG_firstSeqIsRef = true; % set: 1st sequence of every fasta file is reference sequence
        FLAG_Epi = false; % SET: use MPL with epistasis, UNSET: MPL with epistasis not used

        [dirNameDataTemp, dirNameAnalysisTemp, dirNameResultsTemp] = setDirNamesMPLPipeline('Set_dirNames_MPL_SimData1.txt');

        ResultFolder = [dirNameResultsTemp 'Set' num2str(thisSet) '/ng' num2str(ng) '_dT' num2str(dTStep) '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];
        ResultFolder_unEven = [dirNameResultsTemp 'Set' num2str(thisSet) '/ng_' ngSamplingSchemeStr '_dT_' dTSamplingSchemeStr '_Tused_' tjSamplingSchemeStr '_initStr' num2str(numStrainsInInitialPop) '/'];
        % ---------- NO user input required for followng initializations ----------

        if(setConvention == 1)
            FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = false;
        elseif(setConvention == 2)
            FLAG_stratonovich = true;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = false;
        elseif(setConvention == 3)
            FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
            FLAG_linearInt = true;
        end

        % Use stratonovich for increased robustness to sampling effects
        if(FLAG_stratonovich == true)
            convention = 'Stratonovich';
        elseif(FLAG_stratonovich == false && FLAG_linearInt == false)
            convention = 'Ito';
        elseif(FLAG_stratonovich == false && FLAG_linearInt == true)
            convention = 'LinearInter';
        end

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

        fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, textCell);
        fileNamesListThisDir_unEven = findFileNamesWithGivenText(dirNameStr1Files, textCell_unEven);
        numPat = length(fileNamesListThisDir);
        numPat_unEven = length(fileNamesListThisDir_unEven);
        
        if(numPat ~= numPat_unEven)
           disp('Warning: Number of MC runs is not the same for good and sparse sampling cases.')
           disp(['         Computing AUROC over ' num2str(numPat) ' MC runs'])
           pause(1)
        end

        %====================== END INITIALIZATION ================================
        
        thisResultFile = ['AUROC_results_priorConstSC' num2str(priorConstSC) '_' convention '.txt'];
        
        dataIn = dlmread([ResultFolder thisResultFile]);
        dataIn_unEven = dlmread([ResultFolder_unEven thisResultFile]);
        
        aucLinkItr = dataIn(:,[1 2]);
        aucLinkNoMuItr = dataIn(:,[3 4]);
        aucUnLinkItr = dataIn(:,[5 6]);
        aucUnLinkNoMuItr = dataIn(:,[7 8]);
        
        aucLinkItr_unEven = dataIn_unEven(:,[1 2]);
        aucLinkNoMuItr_unEven = dataIn_unEven(:,[3 4]);
        aucUnLinkItr_unEven = dataIn_unEven(:,[5 6]);
        aucUnLinkNoMuItr_unEven = dataIn_unEven(:,[7 8]);
        
        thisResultSelcFile = ['AUROC_Selc_results_priorConstSC' num2str(priorConstSC) '_' convention '.txt'];
        
        
        if(FLAG_loadDiffSelcSites == false)
            dataInSelc = dlmread([ResultFolder thisResultSelcFile]);
            dataInSelc_unEven = dlmread([ResultFolder_unEven thisResultSelcFile]);
        else
            dataInSelc = dlmread([ResultFolder ['AUROC_Selc_results_priorConstSC' num2str(priorConstSC) '_' convention '_selcBasedOn_ngScheme_' ngSamplingSchemeStr '_dtScheme_' dTSamplingSchemeStr '.txt']]);
            dataInSelc_unEven = dlmread([ResultFolder_unEven thisResultSelcFile]);
        end
        
        aucLinkSelcItr = dataInSelc(:,[1 2]);
        aucLinkNoMuSelcItr = dataInSelc(:,[3 4]);
        aucUnLinkSelcItr = dataInSelc(:,[5 6]);
        aucUnLinkNoMuSelcItr = dataInSelc(:,[7 8]);
        
        aucLinkSelcItr_unEven = dataInSelc_unEven(:,[1 2]);
        aucLinkNoMuSelcItr_unEven = dataInSelc_unEven(:,[3 4]);
        aucUnLinkSelcItr_unEven = dataInSelc_unEven(:,[5 6]);
        aucUnLinkNoMuSelcItr_unEven = dataInSelc_unEven(:,[7 8]);
        


        allAuc = [mean(aucLinkItr(aucLinkItr(:,1) ~= -1, 1)) mean(aucLinkItr(aucLinkItr(:,2) ~= -1, 2));
                  mean(aucLinkNoMuItr(aucLinkNoMuItr(:,1) ~= -1, 1)) mean(aucLinkNoMuItr(aucLinkNoMuItr(:,2) ~= -1, 2));
                  mean(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1)) mean(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuItr(aucUnLinkNoMuItr(:,1) ~= -1, 1)) mean(aucUnLinkNoMuItr(aucUnLinkNoMuItr(:,2) ~= -1, 2))];
              
        allAuc_unEven = [mean(aucLinkItr_unEven(aucLinkItr_unEven(:,1) ~= -1, 1)) mean(aucLinkItr_unEven(aucLinkItr_unEven(:,2) ~= -1, 2));
                  mean(aucLinkNoMuItr_unEven(aucLinkNoMuItr_unEven(:,1) ~= -1, 1)) mean(aucLinkNoMuItr_unEven(aucLinkNoMuItr_unEven(:,2) ~= -1, 2));
                  mean(aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,1) ~= -1, 1)) mean(aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuItr_unEven(aucUnLinkNoMuItr_unEven(:,1) ~= -1, 1)) mean(aucUnLinkNoMuItr_unEven(aucUnLinkNoMuItr_unEven(:,2) ~= -1, 2))];      
    
        meanAucBen_si_MPL(ss,tt) = allAuc(1, 1);
        meanAucDel_si_MPL(ss,tt) = allAuc(1, 2);
        meanAucBen_si_SL(ss,tt) = allAuc(3, 1);
        meanAucDel_si_SL(ss,tt) = allAuc(3, 2);
        
        meanAucBen_si_MPL_unEven(ss,tt) = allAuc_unEven(1, 1);
        meanAucDel_si_MPL_unEven(ss,tt) = allAuc_unEven(1, 2);
        meanAucBen_si_SL_unEven(ss,tt) = allAuc_unEven(3, 1);
        meanAucDel_si_SL_unEven(ss,tt) = allAuc_unEven(3, 2);
        
        stdAucBen_si_MPL(ss,tt) = std(aucLinkItr(aucLinkItr(:,1) ~= -1, 1));
        stdAucDel_si_MPL(ss,tt) = std(aucLinkItr(aucLinkItr(:,2) ~= -1, 2));
        stdAucBen_si_SL(ss,tt) = std(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1));
        stdAucDel_si_SL(ss,tt) = std(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2));
        
        stdAucBen_si_MPL_unEven(ss,tt) = std(aucLinkItr_unEven(aucLinkItr_unEven(:,1) ~= -1, 1));
        stdAucDel_si_MPL_unEven(ss,tt) = std(aucLinkItr_unEven(aucLinkItr_unEven(:,2) ~= -1, 2));
        stdAucBen_si_SL_unEven(ss,tt) = std(aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,1) ~= -1, 1));
        stdAucDel_si_SL_unEven(ss,tt) = std(aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,2) ~= -1, 2));
        
        allAucSelc = [mean(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1, 1)) mean(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1, 2));
                  mean(aucLinkNoMuSelcItr(aucLinkNoMuSelcItr(:,1) ~= -1, 1)) mean(aucLinkNoMuSelcItr(aucLinkNoMuSelcItr(:,2) ~= -1, 2));
                  mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1, 1)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuSelcItr(aucUnLinkNoMuSelcItr(:,1) ~= -1, 1)) mean(aucUnLinkNoMuSelcItr(aucUnLinkNoMuSelcItr(:,2) ~= -1, 2))];
              
        allAucSelc_unEven = [mean(aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,1) ~= -1, 1)) mean(aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,2) ~= -1, 2));
                  mean(aucLinkNoMuSelcItr_unEven(aucLinkNoMuSelcItr_unEven(:,1) ~= -1, 1)) mean(aucLinkNoMuSelcItr_unEven(aucLinkNoMuSelcItr_unEven(:,2) ~= -1, 2));
                  mean(aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,1) ~= -1, 1)) mean(aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,2) ~= -1, 2));
                  mean(aucUnLinkNoMuSelcItr_unEven(aucUnLinkNoMuSelcItr_unEven(:,1) ~= -1, 1)) mean(aucUnLinkNoMuSelcItr_unEven(aucUnLinkNoMuSelcItr_unEven(:,2) ~= -1, 2))];      
                    
        meanAucBenSelc_si_MPL(ss,tt) = allAucSelc(1, 1);
        meanAucDelSelc_si_MPL(ss,tt) = allAucSelc(1, 2);
        meanAucBenSelc_si_SL(ss,tt) = allAucSelc(3, 1);
        meanAucDelSelc_si_SL(ss,tt) = allAucSelc(3, 2);
        
        meanAucBenSelc_si_MPL_unEven(ss,tt) = allAucSelc_unEven(1, 1);
        meanAucDelSelc_si_MPL_unEven(ss,tt) = allAucSelc_unEven(1, 2);
        meanAucBenSelc_si_SL_unEven(ss,tt) = allAucSelc_unEven(3, 1);
        meanAucDelSelc_si_SL_unEven(ss,tt) = allAucSelc_unEven(3, 2);
        
        stdAucBenSelc_si_MPL(ss,tt) = std(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1, 1));
        stdAucDelSelc_si_MPL(ss,tt) = std(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1, 2));
        stdAucBenSelc_si_SL(ss,tt) = std(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1, 1));
        stdAucDelSelc_si_SL(ss,tt) = std(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1, 2));
        
        stdAucBenSelc_si_MPL_unEven(ss,tt) = std(aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,1) ~= -1, 1));
        stdAucDelSelc_si_MPL_unEven(ss,tt) = std(aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,2) ~= -1, 2));
        stdAucBenSelc_si_SL_unEven(ss,tt) = std(aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,1) ~= -1, 1));
        stdAucDelSelc_si_SL_unEven(ss,tt) = std(aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,2) ~= -1, 2));
        
        dataLinkPosCell{ss,tt} = aucLinkItr(aucLinkItr(:,1) ~= -1, 1);
        dataLinkNegCell{ss,tt} = aucLinkItr(aucLinkItr(:,2) ~= -1, 2);
        dataUnLinkPosCell{ss,tt} = aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1);
        dataUnLinkNegCell{ss,tt} = aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2);

        dataLinkSelcPosCell{ss,tt} = aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1, 1);
        dataLinkSelcNegCell{ss,tt} = aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1, 2);
        dataUnLinkSelcPosCell{ss,tt} = aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1, 1);
        dataUnLinkSelcNegCell{ss,tt} = aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1, 2);
        
        dataLinkPosCell_unEven{ss,tt} = aucLinkItr_unEven(aucLinkItr_unEven(:,1) ~= -1, 1);
        dataLinkNegCell_unEven{ss,tt} = aucLinkItr_unEven(aucLinkItr_unEven(:,2) ~= -1, 2);
        dataUnLinkPosCell_unEven{ss,tt} = aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,1) ~= -1, 1);
        dataUnLinkNegCell_unEven{ss,tt} = aucUnLinkItr_unEven(aucUnLinkItr_unEven(:,2) ~= -1, 2);

        dataLinkSelcPosCell_unEven{ss,tt} = aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,1) ~= -1, 1);
        dataLinkSelcNegCell_unEven{ss,tt} = aucLinkSelcItr_unEven(aucLinkSelcItr_unEven(:,2) ~= -1, 2);
        dataUnLinkSelcPosCell_unEven{ss,tt} = aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,1) ~= -1, 1);
        dataUnLinkSelcNegCell_unEven{ss,tt} = aucUnLinkSelcItr_unEven(aucUnLinkSelcItr_unEven(:,2) ~= -1, 2);
        
        
 
    end
end



% all sites
dataLinkPosCellPlot{1} = dataLinkPosCell{1,1};
dataLinkPosCellPlot{2} = dataLinkPosCell_unEven{1,1};
dataLinkPosCellPlot{3} = [0 0 0 0 0.1 -0.1];
dataLinkPosCellPlot{4} = dataLinkPosCell{1,2};
dataLinkPosCellPlot{5} = dataLinkPosCell_unEven{1,2};

dataLinkNegCellPlot{1} = dataLinkNegCell{1,1};
dataLinkNegCellPlot{2} = dataLinkNegCell_unEven{1,1};
dataLinkNegCellPlot{3} = [0 0 0 0 0.1 -0.1];
dataLinkNegCellPlot{4} = dataLinkNegCell{1,2};
dataLinkNegCellPlot{5} = dataLinkNegCell_unEven{1,2};


% Only Poly sites
dataLinkSelcPosCellPlot{1} = dataLinkSelcPosCell{1,1};
dataLinkSelcPosCellPlot{2} = dataLinkSelcPosCell_unEven{1,1};
dataLinkSelcPosCellPlot{3} = [0 0 0 0 0.1 -0.1];
dataLinkSelcPosCellPlot{4} = dataLinkSelcPosCell{1,2};
dataLinkSelcPosCellPlot{5} = dataLinkSelcPosCell_unEven{1,2};

dataLinkSelcNegCellPlot{1} = dataLinkSelcNegCell{1,1};
dataLinkSelcNegCellPlot{2} = dataLinkSelcNegCell_unEven{1,1};
dataLinkSelcNegCellPlot{3} = [0 0 0 0 0.1 -0.1];
dataLinkSelcNegCellPlot{4} = dataLinkSelcNegCell{1,2};
dataLinkSelcNegCellPlot{5} = dataLinkSelcNegCell_unEven{1,2};


dataLinkSelcCellPlot{1} = dataLinkSelcPosCell{1,1};
dataLinkSelcCellPlot{2} = dataLinkSelcPosCell_unEven{1,1};
dataLinkSelcCellPlot{3} = [0 0 0 0 0.1 -0.1];
dataLinkSelcCellPlot{4} = dataLinkSelcNegCell{1,1};
dataLinkSelcCellPlot{5} = dataLinkSelcNegCell_unEven{1,1};

axes(ha(5))

bar(1, mean(dataLinkSelcCellPlot{1}), 'FaceColor', color_scheme31(30,:))
hold on
bar(2, mean(dataLinkSelcCellPlot{2}), 'FaceColor', color_scheme11(80,:))
bar(4, mean(dataLinkSelcCellPlot{4}), 'FaceColor', color_scheme31(30,:))
bar(5, mean(dataLinkSelcCellPlot{5}), 'FaceColor', color_scheme11(80,:))
xCord = [1 2 4 5];
yCord = [mean(dataLinkSelcCellPlot{1}) mean(dataLinkSelcCellPlot{2}) mean(dataLinkSelcCellPlot{4}) mean(dataLinkSelcCellPlot{5})];
errBarData = [std(dataLinkSelcCellPlot{1}) std(dataLinkSelcCellPlot{2}) std(dataLinkSelcCellPlot{4}) std(dataLinkSelcCellPlot{5})];
errorbar(xCord, yCord, errBarData, 'k', 'LineStyle', 'none', 'CapSize', 0)

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
      'XTickLabel'  , {'Beneficial' 'Deleterious'}, ...'XTickLabel'  , {' ' '5' ' ' '10' ' ' '20'}, ...
      'XTick'       , [1.5 4.5], ...
      'YTick'  , 0.5:0.1:1, ...
      'LineWidth', 0.5, ...
      'FontSize', 6)
axis([0 6 0.5 1])
%title()
ylabel('Mean AUROC')

dimDummy = [0 0 0 0];
box1 = annotation('rectangle',dimDummy,'FaceColor',color_scheme31(30,:))
box1.Units = 'centimeter';
box1.Position = [10 1.55 0.485 0.28];

box2 = annotation('rectangle',dimDummy,'FaceColor',color_scheme11(80,:))
box2.Units = 'centimeter';
box2.Position = [10 1.05 0.485 0.28];

textLeg1 = annotation('textbox',dimDummy,'String','Perfect sampling','FitBoxToText','on')
textLeg1.Units = 'centimeter';
textLeg1.Position = [10.6 1.5 4 0.5];
textLeg1.LineStyle = 'none';
textLeg1.FontSize = 8;

textLeg2 = annotation('textbox',dimDummy,'String','Heterogeneous sampling','FitBoxToText','on')
textLeg2.Units = 'centimeter';
textLeg2.Position = [10.6 1.0 5 0.5];
textLeg2.LineStyle = 'none';
textLeg2.FontSize = 8;


%---------------------------------------------------------------------

textA = annotation('textbox',dimDummy,'String','a','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25 yDim-topMargin 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;

textA = annotation('textbox',dimDummy,'String','b','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25+width1+hgap1-0.1 yDim-topMargin 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;

textA = annotation('textbox',dimDummy,'String','c','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25+2*(width1+hgap1)-0.1 yDim-topMargin 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;

textA = annotation('textbox',dimDummy,'String','d','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25 yDim-topMargin-height1-vgap1 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;

textA = annotation('textbox',dimDummy,'String','e','FitBoxToText','on')
textA.Units = 'centimeter';
textA.Position = [0.25+width2+hgap2-0.1 yDim-topMargin-height1-vgap1 5 0.5];
textA.LineStyle = 'none';
textA.FontWeight = 'bold';
textA.FontSize = 9;


if(saveFigs == 1)
   figname = ['fig-supp12'];
   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])% ,[0 0 8 6])
   set(gcf, 'renderer', 'painters');
   %print(figname, '-dpng','-r400')
   set(gcf, 'PaperSize', [xDim yDim])
   print([dirNameResultsTemp figname], '-dpdf', '-fillpage')   
end

