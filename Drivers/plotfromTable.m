%% Read and combine data from tables
clc;
clear;
savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\20200126_benchmarking';
saveName = '20200126_benchmarking';
fNames1P = '1pFstim*.mat';
fNames2P = '2pFstim*.mat';

% Load results table 1P and 2P
fList1P = dir(fullfile(savepath,fNames1P));
lib1P = struct;
T_well_1P = table;
for i = 1:length(fList1P)
    lib1P(i).T_well_all = load(fullfile(fList1P(i).folder,fList1P(i).name),'T_well_all').T_well_all;
    T_well_1P = [T_well_1P;lib1P(i).T_well_all]; 
end

fList2P = dir(fullfile(savepath,fNames2P));
lib2P = struct;
T_well_2P = table;
for i = 1:length(fList2P)
    lib2P(i).T_well_all = load(fullfile(fList2P(i).folder,fList2P(i).name),'T_well').T_well;
    T_well_2P = [T_well_2P;lib2P(i).T_well_all]; 
end

%% Standardized the names
pruneName = @(x) regexprep(x,'-cyOFP.*','');
stname = @(x) libplot.fzmatch(x);
prunedTName = rowfun(pruneName,T_well_1P,'InputVariables','ConstructName');
T_well_1P = [rowfun(stname,prunedTName),T_well_1P];
T_well_1P.Properties.VariableNames{'Var1'} = 'Legal Name';

prunedTName = rowfun(pruneName,T_well_2P,'InputVariables','ConstructName');
T_well_2P = [rowfun(stname,prunedTName),T_well_2P];
T_well_2P.Properties.VariableNames{'Var1'} = 'Legal Name';

%% Add kinetics and group rows by name 
selpp1P = {'Legal Name','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','TgtBrightnessFinalMean','TgtBleachingRatioMean','AveTraceGRratioMean',...
'dFF0_shortStim','dFF0_longStim','Bleachingratio_total','TgtAUC_sig_NormaltoRandT',...
'Areableach','TgtBrightnessMeanBleach','RefBrightnessMeanBleach','GRratio_pixelMeanBleach',...
'GRratio_meanMeanBleach','Bleachingratio_totalBleach','TgtAUCNormaltoRandTBleach','Kinetics'};
selpp2P = {'Legal Name','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','TgtBrightnessFinalMean','TgtBleachingRatioMean','AveTraceGRratioMean',...
'dFF0_shortStim','dFF0_longStim','Bleachingratio_total','Bleachingratio_2_1s','TgtAUC_sig_NormaltoRandT','Kinetics'};
T_well_1P.Kinetics = T_well_1P.dFF0_shortStim./T_well_1P.dFF0_longStim;
T_well_2P.Kinetics = T_well_2P.dFF0_shortStim./T_well_2P.dFF0_longStim;
OnePhoton = libplot.groupWells(T_well_1P,'groupOn','Legal Name','selectProperties',selpp1P);
TwoPhoton = libplot.groupWells(T_well_2P,'groupOn','Legal Name','selectProperties',selpp2P);
Benchmarking_nangroup = outerjoin(OnePhoton,TwoPhoton,'LeftKeys','Legal Name','RightKeys','Legal Name','MergeKeys',true);
Benchmarking_nangroup = sortrows(Benchmarking_nangroup,'Legal Name','ascend');
Benchmarking_nangroup.Properties.RowNames = Benchmarking_nangroup.('Legal Name');

% Filling missing values in dFF0 column
T_well_1Pfilled = T_well_1P;
T_well_2Pfilled = T_well_2P;
T_well_1Pfilled.Kinetics = T_well_1Pfilled.dFF0_shortStim./T_well_1Pfilled.dFF0_longStim;
T_well_1Pfilled.Kinetics = T_well_1Pfilled.dFF0_shortStim./T_well_1Pfilled.dFF0_longStim;
ppTBF = {'dFF0_shortStim','dFF0_longStim','Kinetics'};
T_well_1Pfilled = fillmissing(T_well_1Pfilled,'constant',0,'DataVariables',ppTBF);
T_well_2Pfilled = fillmissing(T_well_2Pfilled,'constant',0,'DataVariables',ppTBF);
OnePhoton = libplot.groupWells(T_well_1Pfilled,'groupOn','Legal Name','selectProperties',selpp1P);
TwoPhoton = libplot.groupWells(T_well_2Pfilled,'groupOn','Legal Name','selectProperties',selpp2P);
Benchmarking_group = outerjoin(OnePhoton,TwoPhoton,'LeftKeys','Legal Name','RightKeys','Legal Name','MergeKeys',true);
Benchmarking_group = sortrows(Benchmarking_group,'Legal Name','ascend');
Benchmarking_group.Properties.RowNames = Benchmarking_group.('Legal Name');

save(fullfile(savepath,strcat(saveName,'_T_well_nanstats_grouplvl.mat')),'Benchmarking_nangroup');
save(fullfile(savepath,strcat(saveName,'_T_well_stats_grouplvl.mat')),'Benchmarking_group');

%% Bar Plot
saveformat = '.fig';
pltTable = Benchmarking_nangroup;
dt = {'Area_Mean_OnePhoton','dFF0_shortStim_Mean_OnePhoton', 'dFF0_longStim_Mean_OnePhoton', ...
    'GRratioMean_Mean_OnePhoton', 'Bleachingratio_total_Mean_OnePhoton',...
    'AveTraceGRratioMean_Mean_TwoPhoton'};

close all
% Figure 1: 1P dFF0 
SortOn = 'dFF0_longStim_Mean_OnePhoton'; % Input Sorting Criteria
cmap = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980 ];
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'dFF01P',saveformat));
[b,T_well_stats_sorted] = libplot.formattedBar(pltTable, 'dFF0_longStim_Mean_OnePhoton', ...
    'dFF0_shortStim_Mean_OnePhoton', 'pp1std','dFF0_longStim_STD_OnePhoton',...
    'pp2std', 'dFF0_shortStim_STD_OnePhoton', 'sort_on', SortOn, 'Name','Legal Name','ylabel',...
    'Response Amplitude \DeltaF/F0','colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

% Figure 2: GR ratio 1P
SortOn = 'GRratioMean_Mean_OnePhoton';
cmap = [0.4660, 0.6740, 0.1880];
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'GRratio1P',saveformat));
libplot.formattedBar(pltTable,'GRratioMean_Mean_OnePhoton','pp1std','GRratioMean_STD_OnePhoton', ...
    'sort_on', SortOn, 'sort_direction','descend','Name','Legal Name','ylabel','Normalized GR ratio ',...
    'norm2ctrl','ASAP1','colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

% figure 3 1P Photostability 
SortOn = 'AveTraceGRratioMean_Mean_OnePhoton';
cmap = [60,179,113;152,251,152]./256;
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'PS',saveformat));
libplot.formattedBar(pltTable,'Bleachingratio_total_Mean_OnePhoton','Bleachingratio_total_Mean_TwoPhoton',...
    'pp1std','Bleachingratio_total_STD_OnePhoton', 'pp2std', 'Bleachingratio_total_STD_TwoPhoton', ...
    'sort_on', SortOn, 'Name','Legal Name','ylabel','Bleaching ratio','colormap', ...
    cmap, 'datatip',dt, 'savePath',savefigpath);

% figure4 1P Kinetics bar 
SortOn = 'Kinetics_Mean_OnePhoton';
cmap = [30,144,255; 135,206,250]/256;
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'kinetics',saveformat));
libplot.formattedBar(pltTable,'Kinetics_Mean_OnePhoton','Kinetics_Mean_TwoPhoton', ...
    'pp1std','Kinetics_STD_OnePhoton', 'pp2std','Kinetics_STD_TwoPhoton', ...
     'sort_on', SortOn,'Name','Legal Name', 'sort_direction','descend','ylabel','kinetics (short/long)',...
     'colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

 % figure 5 2P dFF0 
SortOn = 'dFF0_longStim_Mean_TwoPhoton'; % Input Sorting Criteria
cmap = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980 ];
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'dFF02P',saveformat));
[b,T_well_stats_sorted] = libplot.formattedBar(pltTable, 'dFF0_longStim_Mean_TwoPhoton', ...
    'dFF0_shortStim_Mean_TwoPhoton', 'pp1std','dFF0_longStim_STD_TwoPhoton',...
    'pp2std', 'dFF0_shortStim_STD_TwoPhoton', 'sort_on', SortOn, 'Name','Legal Name','ylabel',...
    'Response Amplitude \DeltaF/F0','colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

%% Scatter plot



%% Load Ephys
dirEphys = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\Ephys';
NameEphys = 'ThreeTest_*_wASAP1.mat';
fListEphys = dir(fullfile(dirEphys,NameEphys));
for i = 1:length(fListEphys)
    load(fullfile(fListEphys(i).folder,fListEphys(i).name));
end

% Standardize and sort by the names
Train_AP = [rowfun(stname,Train_AP(:,{'Name'})),Train_AP];
Train_AP.Properties.VariableNames{'Var1'} = 'Legal Name';
Train_AP.Properties.RowNames = Train_AP.('Legal Name');

Transition_curve = [rowfun(stname,Transition_curve(:,{'Name'})),Transition_curve];
Transition_curve.Properties.VariableNames{'Var1'} = 'Legal Name';
Transition_curve.Properties.RowNames = Transition_curve.('Legal Name');

Single_AP = [rowfun(stname,Single_AP(:,{'Name'})),Single_AP];
Single_AP.Properties.VariableNames{'Var1'} = 'Legal Name';
Single_AP.Properties.RowNames = Single_AP.('Legal Name');

%% Compare with Ephys
ShortStimMean = 'dFF0_shortStim_Mean_OnePhoton';
ShortStimStd = 'dFF0_shortStim_STD_OnePhoton';
LongStimMean = 'dFF0_longStim_Mean_OnePhoton';
LongStimStd = 'dFF0_longStim_STD_OnePhoton';
CtrlFP = 'cyOFP (One-photon)';
sel = Single_AP.('Legal Name')([1:3,5:7]);
saveext = 'fig';

figure(1) % Compare Ephys with Fstim
clf
hold on

EphysdFF0Mean = [];
EphysdFF0Std = [];
FstimLongdFF0Mean = [];
FstimLongdFF0Std = [];
for i = 1:length(sel)
    EphysdFF0Mean(i) = Transition_curve{sel(i),'V30Mean'};
    EphysdFF0Std(i) = Transition_curve{sel(i),'V30Std'};
    FstimLongdFF0Mean(i) = Benchmarking_group{sel(i),LongStimMean};
    FstimLongdFF0Std(i) = Benchmarking_group{sel(i),LongStimStd};
    errorbar(-EphysdFF0Mean(i),-FstimLongdFF0Mean(i),FstimLongdFF0Std(i),...
        FstimLongdFF0Std(i),EphysdFF0Std(i),EphysdFF0Std(i),'.','MarkerSize',30)
end
x = -EphysdFF0Mean;
y = -FstimLongdFF0Mean;
[k, R2] = libplot.lnrFitting(x,y,false,true,1);
xlabel({'Ephys VClamp = +30mV','-\DeltaF/F0'})
ylabel({'Fstim = 100ms@100Hz','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis([0 inf 0 inf])
axis square
legend(Transition_curve.Name(sel),'location','best')
set(gcf, 'Position',  [50, 10, 500, 500])
libplot.saveEXT(fullfile(savepath,strcat(saveName,'_VClamp_',CtrlFP,'.',saveext)));

figure(2)
clf
hold on
EphysAPdFF0Mean = [];
EphysAPdFF0Std = [];
FstimShortdFF0Mean = [];
FstimShortdFF0Std = [];
for i = 1:length(sel) %1:size(Single_AP,1)
    EphysAPdFF0Mean(i) = 1-Single_AP{sel(i),'SpikeMean'};
    EphysAPdFF0Std(i) = Single_AP{sel(i),'SpikeStd'};
    FstimShortdFF0Mean(i) = Benchmarking_group{sel(i),ShortStimMean};
    FstimShortdFF0Std(i) = Benchmarking_group{sel(i),ShortStimStd};
    errorbar(EphysAPdFF0Mean(i),-FstimShortdFF0Mean(i),FstimShortdFF0Std(i),...
        FstimShortdFF0Std(i),EphysAPdFF0Std(i),EphysAPdFF0Std(i),'.','MarkerSize',30)
end
x = EphysAPdFF0Mean;
y = -FstimShortdFF0Mean;
[k, R2] = libplot.lnrFitting(x,y,false,true,2);
xlabel({'Ephys Single AP','-\DeltaF/F0'})
ylabel({'Fstim = 1ms ','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis([0 inf 0 inf])
axis square
legend(Single_AP.Name(sel),'location','best')
set(gcf, 'Position',  [550, 10, 500, 500])
libplot.saveEXT(fullfile(savepath,strcat(saveName,'_SingleAP_',CtrlFP,'.',saveext)));

figure(3)
clf
hold on
EphysTrainAPdFF0Mean = [];
EphysTrainAPdFF0Std = [];
FstimLongdFF0Mean = [];
FstimLongdFF0Std = [];
for i = 1:length(sel) %size(Train_AP,1)
    EphysTrainAPdFF0Mean(i) =  1-Train_AP{sel(i),'SpikeMean'};
    EphysTrainAPdFF0Std(i) = Train_AP{sel(i),'SpikeStd'};
    FstimLongdFF0Mean(i) = Benchmarking_group{sel(i),LongStimMean};
    FstimLongdFF0Std(i) = Benchmarking_group{sel(i),LongStimStd};
    errorbar(EphysTrainAPdFF0Mean(i),-FstimLongdFF0Mean(i),FstimLongdFF0Std(i),...
        FstimLongdFF0Std(i),EphysTrainAPdFF0Std(i),EphysTrainAPdFF0Std(i),'.','MarkerSize',30)
end
x = EphysTrainAPdFF0Mean;
y = -FstimLongdFF0Mean;
p = fitlm(x,y,'Intercept',false);
[k, R2] = libplot.lnrFitting(x,y,false,true,3);
xlabel({'Ephys train AP','-\DeltaF/F0'})
ylabel({'Fstim = 100ms 100Hz','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis square
axis([0 inf 0 inf])
legend(Train_AP.Name(sel),'location','best')
set(gcf, 'Position',  [1050, 10, 500, 500])
libplot.saveEXT(fullfile(savepath,strcat(saveName,'_TrainAP_',CtrlFP,'.',saveext)));

%% Compare 1P with 2P
% cmap = colormap(colorcube(12));
cmap = [255,223,158;255,194,115;229,105,105;193,85,139;138,73,161;24, 120, 192;...
    48, 168, 120;168,216,234;170,150,218;150,218,170;252,186,211;250,116,79;...
    250,79,85;175,161,100;100,114,175;97,232,49;208,120,240]./256;
pp1 = 'dFF0_shortStim_Mean_OnePhoton';
pp1Std = 'dFF0_shortStim_STD_OnePhoton';
pp2 = 'dFF0_shortStim_Mean_TwoPhoton';
pp2Std = 'dFF0_shortStim_STD_TwoPhoton';
saveext = 'fig';
table1 = Benchmarking_group;
table2 = Benchmarking_group;
sel = table1.('Legal Name')([1:17]);

figure(1) % Compare Ephys with Fstim
clf
hold on
subplot(2,1,1)
hold on
x = [];
xStd = [];
y = [];
yStd = [];
for i = 1:length(sel)
    x(i) = table1{sel(i),pp1};
    xStd(i) = table1{sel(i),pp1Std};
    y(i) = table2{sel(i),pp2};
    yStd(i) = table2{sel(i),pp2Std};
    errorbar(x(i),y(i),yStd(i),yStd(i),xStd(i),xStd(i),'.','MarkerSize',30,...
        'Color',cmap(i,:));
end
[k, R2] = libplot.lnrFitting(x,y,false,true,1);
xlabel({'One-Photon -\DeltaF/F0'})
ylabel({'Two-Photon -\DeltaF/F0'})
title({'One-Photon vs Two-Photon','Short Stimulation'})
axis square
legend(sel,'location','eastoutside')

pp1 = 'dFF0_longStim_Mean_OnePhoton';
pp1Std = 'dFF0_longStim_STD_OnePhoton';
pp2 = 'dFF0_longStim_Mean_TwoPhoton';
pp2Std = 'dFF0_longStim_STD_TwoPhoton';
saveext = 'fig';
table1 = Benchmarking_group;
table2 = Benchmarking_group;

subplot(2,1,2)
hold on
x = [];
xStd = [];
y = [];
yStd = [];
for i = 1:length(sel)
    x(i) = table1{sel(i),pp1};
    xStd(i) = table1{sel(i),pp1Std};
    y(i) = table2{sel(i),pp2};
    yStd(i) = table2{sel(i),pp2Std};
    errorbar(x(i),y(i),yStd(i),yStd(i),xStd(i),xStd(i),'.','MarkerSize',30,...
        'Color',cmap(i,:));
end
[k, R2] = libplot.lnrFitting(x,y,false,true,1);
xlabel({'One-Photon -\DeltaF/F0'})
ylabel({'Two-Photon -\DeltaF/F0'})
title({'One-Photon vs Two-Photon','Long Stimulation'})
axis square
legend(sel,'location','eastoutside')
set(gcf, 'Position',  [50, 10, 900, 900])