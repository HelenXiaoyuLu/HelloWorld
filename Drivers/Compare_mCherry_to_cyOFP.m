%% Compare mCherry and cyOFP benchmarking
% Load and standardize names
clear;
clc
% cyOFP 1P Field stimulation dataset
dir_cyOFP = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200202_benchmarking_plate1\20200202_benchmarking_20X_update20200510';
load(fullfile(dir_cyOFP,'20200202 1PFstim platform characterization 10x vs 20x_T_well_stats v20200512.mat'));
T_group_1PFstim_cyOFP = T_well_stats;
T_group_1PFstim_cyOFP.dFF0_shortStim_Mean = fillmissing(T_group_1PFstim_cyOFP.dFF0_shortStim_Mean,'constant',0);
T_group_1PFstim_cyOFP.dFF0_longStim_Mean = fillmissing(T_group_1PFstim_cyOFP.dFF0_longStim_Mean,'constant',0);

% mCherry 1P Field stimulation dataset
dir_mCherry1PFstim = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_Merge';
load(fullfile(dir_mCherry1PFstim,'20190227 Benchmarking 1P Fstim v20190319'));
T_group_1PFstim_mCherry = T_well_stats;
for i = 1:height(T_group_1PFstim_mCherry)
    T_group_1PFstim_mCherry.Name(i) = libplot.fzmatch(T_group_1PFstim_mCherry.Name(i));
end
T_group_1PFstim_mCherry.Properties.RowNames = T_group_1PFstim_mCherry.Name;
T_group_1PFstim_mCherry.StimShortMean = fillmissing(T_group_1PFstim_mCherry.StimShortMean,'constant',0);
T_group_1PFstim_mCherry.Stim8Mean = fillmissing(T_group_1PFstim_mCherry.Stim8Mean,'constant',0);

% cyOFP 2P Field stimulation dataset
load(fullfile(dir2P,'20200202 2PFstim_output_20200804132138.mat'));
lookup_cyOFP = readtable(fullfile(dir_cyOFP,'All control Green GEVIs (P2).xlsx'),'ReadRowNames',true,'ReadVariableNames',true);
T_well_2P_cyOFP = cyOFPbenchmark2PFstim.T_well;
for i = 1:height(T_well_2P_cyOFP)
    T_well_2P_cyOFP.Name(i) = libplot.fzmatch(lookup_cyOFP.Name(T_well_2P_cyOFP.wellName(i)));
end
T_group_2PFstim_cyOFP = libplot.groupWells(T_well_2P_cyOFP);
T_group_2PFstim_cyOFP.Properties.RowNames = T_group_2PFstim_cyOFP.Name;
save(fullfile(dir_cyOFP,'20200202 cyOFP benchmarking 2PFstim 20200804132138 Trenamed.mat'),'T_well_2P_cyOFP','T_group_2P_cyOFP');

% mCherry 2P Field stimulation dataset
dir_mCherryFstim = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_Merge';
T_fov_2PFstim_mCherry = readtable(fullfile(dir_mCherryFstim,'20190227 Benchmarking 2P Fstim'),'ReadVariableNames',true);
T_well_2PFstim_mCherry = libplot.groupWells(T_fov_2PFstim_mCherry,'groupOn','Construct');
for i = 1:height(T_well_2PFstim_mCherry)
    name = split(T_well_2PFstim_mCherry.Construct{i},'_P');
    T_well_2PFstim_mCherry.Name(i) = libplot.fzmatch(name{1});
end
T_well_2PFstim_mCherry.Properties.RowNames = T_well_2PFstim_mCherry.Construct;
variableNames = {'Construct','UsableRatio_Mean','SNR_Mean','Brightness_Mean',...
    'Background_Mean__Mean','dF_F0_s1__Mean','dF_F0_s2__Mean','Name'};
T_group_2PFstim_mCherry = libplot.groupWells(T_well_2PFstim_mCherry,'selectProperties',variableNames);
T_group_2PFstim_mCherry.Properties.RowNames = T_group_2PFstim_mCherry.Name;
save(fullfile(dir_mCherryFstim,'20190227 Benchmarking 2P Fstim vArchive.mat'),...
    'T_well_2PFstim_mCherry','T_group_2PFstim_mCherry');

% mCherry 2P photostability dataset
dir_mCherryPS = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\2P photobleaching power test';
load(fullfile(dir_mCherryPS,'20190925 2P photobleaching power scan v20200728 T_fovCondition.mat'));
T_fov_2PPs_mCherry_PW30 = T_fovCondition(find(string(T_fovCondition.Settings) =="Power:30 HV:20"),:);;

%% Compare 2P Field stimulation dF/F0
nfig = 1;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_2PFstim_mCherry.Name;
for i = 1:length(nameList)
    scatter(T_group_2PFstim_mCherry.dF_F0_s1__Mean_Mean(nameList(i)),T_group_2PFstim_cyOFP.("dFF0 Short Stim_Mean")(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2PFstim_mCherry.dF_F0_s1__Mean_Mean(nameList);
y = T_group_2PFstim_cyOFP.("dFF0 Short Stim_Mean")(nameList);
xlabel('2P 1ms dF/F0 mCherry')
ylabel('2P 1ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('All benchmarking constructs')

nfig = 2;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_2PFstim_mCherry.dF_F0_s1__Mean_Mean(nameList(i)),T_group_2PFstim_cyOFP.("dFF0 Short Stim_Mean")(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2PFstim_mCherry.dF_F0_s1__Mean_Mean(nameList);
y = T_group_2PFstim_cyOFP.("dFF0 Short Stim_Mean")(nameList);
xlabel('2P 1ms dF/F0 mCherry')
ylabel('2P 1ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,true,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

nfig = 3;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_2PFstim_mCherry.Name;
for i = 1:length(nameList)
    scatter(T_group_2PFstim_mCherry.dF_F0_s2__Mean_Mean(nameList(i)),T_group_2PFstim_cyOFP.("dFF0 Long Stim_Mean")(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2PFstim_mCherry.dF_F0_s2__Mean_Mean(nameList);
y = T_group_2PFstim_cyOFP.("dFF0 Long Stim_Mean")(nameList);
xlabel('2P 10x1.5ms dF/F0 mCherry')
ylabel('2P 10x1.5ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('All benchmarking constructs')

nfig = 4;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_2PFstim_mCherry.dF_F0_s2__Mean_Mean(nameList(i)),T_group_2PFstim_cyOFP.("dFF0 Long Stim_Mean")(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2PFstim_mCherry.dF_F0_s2__Mean_Mean(nameList);
y = T_group_2PFstim_cyOFP.("dFF0 Long Stim_Mean")(nameList);
xlabel('2P 10x1.5ms dF/F0 mCherry')
ylabel('2P 10x1.5ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

%% Compare 2P relative brightness
nfig = 5;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_2Pbrightness_mCherry_plate1.Name;
for i = 1:length(nameList)
    scatter(T_group_2Pbrightness_mCherry_plate1.("GFP2P/RFP1P")(nameList(i)),T_group_2PFstim_cyOFP.GRRatioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2Pbrightness_mCherry_plate1.("GFP2P/RFP1P")(nameList);
y = T_group_2PFstim_cyOFP.GRRatioMean_Mean(nameList);
xlabel('2P GFP/1P RFP mCherry')
ylabel('2P GFP/2P RFP cyOFP')
libplot.lnrFitting(x,y,true,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

nfig = 6;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_2Pbrightness_mCherry_plate1.("GFP2P/RFP1P")(nameList(i)),T_group_2PFstim_cyOFP.GRRatioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_2Pbrightness_mCherry_plate1.("GFP2P/RFP1P")(nameList);
y = T_group_2PFstim_cyOFP.GRRatioMean_Mean(nameList);
xlabel('2P GFP/1P RFP mCherry')
ylabel('2P GFP/2P RFP cyOFP')
libplot.lnrFitting(x,y,true,true,nfig,'polyfit') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

%% Compare 2P photostability

%% Compare 1P Field stimulation
nfig = 11;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_1PFstim_mCherry.Name;
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.StimShortMean(nameList(i)),T_group_1PFstim_cyOFP.dFF0_shortStim_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.StimShortMean(nameList);
y = T_group_1PFstim_cyOFP.dFF0_shortStim_Mean(nameList);
xlabel('1P 1.5ms dF/F0 mCherry')
ylabel('1P 1ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('All benchmarking constructs')

nfig = 12;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.StimShortMean(nameList(i)),T_group_1PFstim_cyOFP.dFF0_shortStim_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.StimShortMean(nameList);
y = T_group_1PFstim_cyOFP.dFF0_shortStim_Mean(nameList);
xlabel('1P 1.5ms dF/F0 mCherry')
ylabel('1P 1ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

nfig = 13;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_1PFstim_mCherry.Name;
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.Stim8Mean(nameList(i)),T_group_1PFstim_cyOFP.dFF0_longStim_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.Stim8Mean(nameList);
y = T_group_1PFstim_cyOFP.dFF0_longStim_Mean(nameList);
xlabel('1P 10x1.5ms dF/F0 mCherry')
ylabel('1P 10x1.5ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('All benchmarking constructs')

nfig = 14;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.Stim8Mean(nameList(i)),T_group_1PFstim_cyOFP.dFF0_longStim_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.Stim8Mean(nameList);
y = T_group_1PFstim_cyOFP.dFF0_longStim_Mean(nameList);
xlabel('1P 10x1.5ms dF/F0 mCherry')
ylabel('1P 10x1.5ms dF/F0 cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

%% Compare 1P relative brightness
nfig = 15;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = T_group_1PFstim_mCherry.Name;
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.GRRatioMean(nameList(i)),T_group_1PFstim_cyOFP.GRratioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.GRRatioMean(nameList);
y = T_group_1PFstim_cyOFP.GRratioMean_Mean(nameList);
xlabel('1P G/R mCherry')
ylabel('1P G/R cyOFP')
libplot.lnrFitting(x,y,false,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('All benchmarking constructs')

nfig = 16;
figure(nfig)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_1PFstim_mCherry.GRRatioMean(nameList(i)),T_group_1PFstim_cyOFP.GRratioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_1PFstim_mCherry.GRRatioMean(nameList);
y = T_group_1PFstim_cyOFP.GRratioMean_Mean(nameList);
xlabel('1P G/R mCherry')
ylabel('1P GRRatioMean cyOFP')
libplot.lnrFitting(x,y,true,true,nfig,'fitlm') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

%% Compare variation in red channel
