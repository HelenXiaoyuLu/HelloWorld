%% Clear 
clearvars -except T_groupFromROI

%% Specify the groups
L1 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\1pFstimLighton_20200202_GEVI_screening_cyOFPset_benchmarking_20X_objective.mat');
L2 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\1pFstimTradition_200126_Fstim_platform_characterization_P1.mat');
L3 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\1pFstimTradition_200126_Fstim_platform_characterization_P2.mat');
savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200202_benchmarking_plate1\20200202_benchmarking_20X_update20200510';
% savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200126_benchmarking';
savepathsplit = split(L2.folderpath,'\');
saveName = savepathsplit{end-1};
% Load archive analysis result from excel
% T_wellArchive = readtable('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\Merge\20190227 1P Fstim 1000Hz archive.xlsx','Sheet','Sheet1');
% concat different data sets
T_well = L1.T_well;
T_well = [L2.T_well;L3.T_well];
% Create standard name column
WellconstructName = cellstr(T_well.ConstructName);
for i = 1:length(WellconstructName)
    Namesplt = split(WellconstructName{i},'.');
    Namesplt2 = split(Namesplt(1),'-cyOFP');
    T_well.Name(i) = string(Namesplt2(1));
end
for i = 1:height(T_well)
    T_well.Name(i) = libplot.fzmatch(T_well.Name{i});
end

%% Well level stats (optional)
T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3','Stim_4'}}, 2); 
T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3'}}, 2); 
T_well.stim_middlemean = mean(T_well{:, {'Stim_4', 'Stim_5' ,'Stim_6'}}, 2); 
T_well.kinetics = T_well.dFF0_shortStim./T_well.dFF0_longStim;
T_well = fillmissing(T_well,'constant',0,'DataVariables',@isnumeric);
    
%% Group Level stats
% replace missing values (nan) with 0
T_well_stats = libplot.groupWells(T_well);
T_well_stats = sortrows(T_well_stats,'Name','ascend');
T_well_stats.Properties.RowNames = T_well_stats.Name;
save(fullfile(savepath,strcat(saveName,'_T_well_stats.mat')),'T_well_stats');

%% Trace (Figure 0)
pltVar = 'TraceNormalizedMean';

f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
t.Padding = 'none';
title(t, 'Well Stats in Construct');
pltVarMean = strcat(pltVar,'_Mean');
pltVarSTD = strcat(pltVar,'_STD');

for i = 1:height(T_well_stats)
    ax = nexttile(t);
    xf = [T_well_stats.(pltVarMean)(i).x; flipud(T_well_stats.(pltVarMean)(i).x)];
    yf = [T_well_stats.(pltVarMean)(i).y + T_well_stats.(pltVarSTD)(i).y; 
          flipud(T_well_stats.(pltVarMean)(i).y - T_well_stats.(pltVarSTD)(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well_stats.(pltVarMean)(i).x, T_well_stats.(pltVarMean)(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well_stats.Name(i), T_well_stats.n_well(i)));
end
saveas(gcf,fullfile(savepath,strcat(saveName,'Well Stats in Construct Trace')));   

%% Group Level Barplot (Figure 1-3, 12)
close all
saveformat = '.fig';
% Figure 1: dFF0
SortOn = 'dFF0_longStim_Mean'; % Input Sorting Criteria
cmap = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980 ];
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'dFF0',saveformat));
dt = {'Area_Mean','dFF0_shortStim_Mean', 'dFF0_longStim_Mean', 'GRratioMean_Mean', ...
    'kinetics_Mean', 'Bleachingratio_total_Mean','Bleachingratio_1s_Mean', };
[b,T_well_stats_sorted] = libplot.formattedBar(T_well_stats, 'dFF0_longStim_Mean', 'dFF0_shortStim_Mean', ...
    'pp1std','dFF0_longStim_STD', 'pp2std', 'dFF0_shortStim_STD', 'sort_on', SortOn, ...
    'ylabel','Response Amplitude \DeltaF/F0','colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

% Figure 2: GR ratio
SortOn = 'GRratioMean_';
cmap = [0.4660, 0.6740, 0.1880];
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'GRratio',saveformat));
libplot.formattedBar(T_well_stats,'GRratioMean_Mean','pp1std','GRratioMean_STD', ...
    'sort_on', SortOn, 'sort_direction','descend', 'ylabel','GR ratio ',...
    'colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

% figure 3 Photostability
SortOn = 'Bleachingratio_total_Mean';
cmap = [60,179,113;152,251,152]./256;
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'PS',saveformat));
libplot.formattedBar(T_well_stats,'Bleachingratio_total_Mean','TgtBleachingRatioMean_Mean',...
    'pp1std','Bleachingratio_total_STD', 'pp2std', 'TgtBleachingRatioMean_STD', ...
    'sort_on', SortOn, 'ylabel','Bleaching ratio','colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

% figure4 Kinetics bar
SortOn = 'kinetics_Mean';
cmap = [30,144,255]/256;
savefigpath = fullfile(savepath,strcat(saveName,'_sortOn_',SortOn,'kinetics',saveformat));
libplot.formattedBar(T_well_stats,'kinetics_Mean', 'pp1std','kinetics_STD', ...
     'sort_on', SortOn, 'sort_direction','descend','ylabel','kinetics (short/long)',...
     'colormap', cmap, 'datatip',dt, 'savePath',savefigpath);

%% Group Level Scatter Plot (Figure 4,13 )
C = colormap(summer);
%C(1:10:160,:)
figure(4)
clf
hold on
Stim5Mean = T_well_stats.Stim5Mean';
Stim5Mean(isnan(Stim5Mean))=0;
StimShortMean = T_well_stats.StimShortMean';
StimShortMean(isnan(StimShortMean))=0;
GRRatio = T_well_stats.GRRatioMean';
AUCMean = T_well_stats.AUCMean;
cmap = colorcube(75);
%AUCMean = normalize(T_well_stats.AUCMean,'range');
% scatter(Stim5Mean',StimShortMean',0.3*AUCMean,normalize(GRRatio,'range'),'filled','MarkerFaceAlpha',0.5)
for i = 1:length(StimShortMean)
    plot(Stim5Mean(i),StimShortMean(i),'.','Color',cmap(i*2-1,:),'MarkerSize',25)%,'MarkerEdgeColor','none')
end
scatter(Stim5Mean',StimShortMean',normalize(GRRatio,'scale')*1000,cmap(1:2:2*length(StimShortMean),:),'filled','MarkerFaceAlpha',1)
errorbar(Stim5Mean', StimShortMean',T_well_stats.StimShortStd',T_well_stats.StimShortStd',T_well_stats.Stim5Std,T_well_stats.Stim5Std,'.k')
xlabel('dF/F0 (100ms 100Hz Stimulations)')
ylabel('dF/F0 (1.5ms Stimulations)')
title({'20190202 1P Fstim','well-level stat'})
legend(T_well_stats.Name,'location','eastoutside')

figure(13)
clf
hold on
StimLongMean = T_well_stats(:, strcat(LongStim,'Mean')).Variables';
StimLongMean(isnan(StimLongMean))=0;
StimShortMean = T_well_stats.StimShortMean';
StimShortMean(isnan(StimShortMean))=0;
GRRatio = T_well_stats.GRRatioMean';
AUCMean = T_well_stats.AUCMean;
cmap = colorcube(75);

for i = 1:length(StimShortMean)
    plot(StimLongMean(i),StimShortMean(i),'.','Color',cmap(i*3-1,:),'MarkerSize',25)%,'MarkerEdgeColor','none')
end
% scatter(StimLongMean',StimShortMean',ones(size(StimShortMean'))*100,cmap(2:3:3*length(StimShortMean)+1,:),'filled','MarkerFaceAlpha',1)
errorbar(StimLongMean', StimShortMean',T_well_stats.StimShortStd',T_well_stats.StimShortStd',T_well_stats(:, strcat(LongStim,'Std')).Variables,T_well_stats(:, strcat(LongStim,'Std')).Variables,'.k')
xlabel('dF/F0 (100ms 100Hz Stimulations)')
ylabel('dF/F0 (1.5ms Stimulations)')
title({'20190202 1P Fstim','well-level stat'})
legend(T_well_stats.Name,'location','eastoutside')

%% Load Archive excel from old pipeline
WellconstructNameArchive = cellstr(T_wellArchive.Name);
for i = 1:length(WellconstructNameArchive)
    Namesplt = split(WellconstructNameArchive{i},' ');
    Namesplt2 = split(Namesplt(1),'.');
    Namesplt3 = split(Namesplt2(1),'_P');
    T_wellArchive.Name(i) = cellstr(Namesplt3(1));
end
T_wellArchive.stim_shortmean = mean(T_wellArchive{:, {'dFF0OfMean_C1S1_', 'dFF0OfMean_C1S2_' ,'dFF0OfMean_C1S3_'}}, 2); 
T_wellArchive.stim_middlemean = mean(T_wellArchive{:, {'dFF0OfMean_C1S4_', 'dFF0OfMean_C1S5_' ,'dFF0OfMean_C1S6_'}}, 2); 

%% Group Level stats from excel
T_well_stats_archive = table();
T_well_stats_archive = libplot.groupWells(T_well_stats_archive);
T_well_stats_archive = sortrows(T_well_stats_archive,'Name','ascend');
save(fullfile(savepath,'T_well_stats_archive.mat'),'T_well_stats_archive')

%% Bar plot from excel (Figure 5-6)
T_well_stats_archive_sorted = T_well_stats_archive(I,:);
LongStim = 'Stim8';

figure(5)
clf
hold on
NameCell = cellstr(T_well_stats_archive_sorted.Name);
X = 1:height(T_well_stats_archive_sorted);
hBar = bar(X,T_well_stats_archive_sorted.GRRatioMean,'BarWidth',0.8,'LineWidth', 1);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_archive_sorted.GRRatioStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel('GR Ratio')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
% title({'20190202 1P GR ratio','well-level stat Archive'})
set(gcf, 'Position',  [110, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_GRRatio_archive.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_GRRatio_archive.fig'))

figure(6)
clf
hold on
NameCell = cellstr(T_well_stats_archive_sorted.Name);
X = 1:numel(T_well_stats_archive_sorted.Name);
hBar = bar(X,[T_well_stats_archive_sorted.StimShortMean,T_well_stats_archive_sorted(:, strcat(LongStim,'Mean')).Variables],'grouped','BarWidth',1,'LineWidth', 1);
Width = hBar(2).XEndPoints - hBar(1).XEndPoints;
errorbar(hBar(1).XEndPoints-0.00*Width,hBar(1).YData,T_well_stats_archive_sorted.StimShortStd,'.k')
errorbar(hBar(2).XEndPoints-0.00*Width,hBar(2).YData,T_well_stats_archive_sorted(:, strcat(LongStim,'Std')).Variables,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Response Amplitude','[\DeltaF/F0]'})
legend('1.0 ms','100 ms','location','best')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.8500, 0.3250, 0.0980];	
hBar(2).FaceColor = [0.9290, 0.6940, 0.1250];	
% title({'20190202 1P Fstim','well-level stat Archive'})
set(gcf, 'Position',  [710, 10, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0_archive.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0_archive.fig'))

%% Pipeline comparison (Figure 7-8)
figure(7)
clf
hold on
scatter(T_well_stats_sorted.StimShortMean,T_well_stats_archive_sorted.StimShortMean,[],[0.8500, 0.3250, 0.0980],'filled','MarkerFaceAlpha',0.9)
scatter(T_well_stats_sorted.Stim8Mean,T_well_stats_archive_sorted.Stim8Mean,[],[0.9290, 0.6940, 0.1250],'filled','MarkerFaceAlpha',0.9)
lnrFitting(T_well_stats_sorted.StimShortMean,T_well_stats_archive_sorted.StimShortMean,false,true,7)
lnrFitting(T_well_stats_sorted.Stim8Mean,T_well_stats_archive_sorted.Stim8Mean,false,true,7)
title('dFF0')
xlabel({'New Pipeline','\DeltaF/F0'})
ylabel({'Old Pipeline','\DeltaF/F0'})
legend('1 ms','100 ms','location','best')
axis square
set(gcf, 'Position',  [10, 10, 400,400])

figure(8)
hold on
% subplot(1,2,2)
h = scatter(T_well_stats_sorted.GRRatioMean,T_well_stats_archive_sorted.GRRatioMean,[],[0.4660, 0.6740, 0.1880],'filled')
lnrFitting(T_well_stats_sorted.GRRatioMean,T_well_stats_archive_sorted.GRRatioMean,false,true,8)
title('GR Ratio')
xlabel({'New Pipeline','Relative Brightness (a.u.)'})
ylabel({'Old Pipeline','Relative Brightness (a.u.)'})
axis square
set(gcf, 'Position',  [10, 10, 400,400])

%% Load excel trace
Trace1P_0227 = readtable('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\SummarizeFigures.xlsx','Sheet','190227 1P photobleaching');
Photostability = Trace1P_0227(end,2:17);
X = table2cell(Photostability)';
X1 = zeros(size(X));
for i = 1:length(X)
    X1(i) = str2double(cell2mat(X(i)));
end
Photostability_0227 = table();
VarNames = Trace1P_0227.Properties.VariableNames(2:17)';
Photostability_0227.Name = VarNames;
Photostability_0227.Remained_Brightess = X1;
figure(3)
subplot(2,1,2)
hold on
NameCell = cellstr(Photostability_0227.Name);
X = 1:16;
hBar = bar(X,Photostability_0227.Remained_Brightess,'BarWidth',0.8,'LineWidth', 1,'FaceAlpha', 0.5);
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel('1 - dF/F0')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
title({'20190227 1P photostabillity','Trace analysis archive'})

%% Compare benchmarking (Figure 9-11)
% Load datasets
% B1 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200202_benchmarking_plate1\20200202_benchmarking_plate1_T_well_stats.mat');
% B2 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_Merge\20190227_benchmarking_Merge_T_well_stats.mat');
% B1_TwellStats = sortrows(B1.T_well_stats,'Name','ascend');
% B2_TwellStats = sortrows(B2.T_well_stats,'Name','ascend');

figure(9)
clf
hold on
B1toB2 = [1:5,7:14,17:19];
scatter(B1_TwellStats(B1toB2,'StimShortMean').Variables,...
    B2_TwellStats.StimShortMean,[],[0.8500, 0.3250, 0.0980],'filled','MarkerFaceAlpha',0.9)
scatter(B1_TwellStats(B1toB2,'Stim5Mean').Variables,...
    B2_TwellStats.Stim8Mean,[],[0.9290, 0.6940, 0.1250],'filled','MarkerFaceAlpha',0.9)
lnrFitting(B1_TwellStats(B1toB2,'StimShortMean').Variables,...
    B2_TwellStats.StimShortMean,false,true,9)
lnrFitting(B1_TwellStats(B1toB2,'Stim5Mean').Variables,...
    B2_TwellStats.Stim8Mean,false,true,9)
title('Compare dFF0 from different constructs')
xlabel({'-GSS-cyOFP','\DeltaF/F0'})
ylabel({'-P2A-mCherry-CAAX','\DeltaF/F0'})
legend('Short Stimulation','Long Stimulation','location','best')
axis square
set(gcf, 'Position',  [10, 10, 400,400])

figure(10)
clf
hold on
% subplot(1,2,2)
scatter(B1_TwellStats(B1toB2,'GRRatioMean').Variables,...
    B2_TwellStats.GRRatioMean,[],[0.4660, 0.6740, 0.1880],'filled')
lnrFitting(B1_TwellStats(B1toB2,'GRRatioMean').Variables,...
    B2_TwellStats.GRRatioMean,false,true,10)
title('Compare GR ratio from different constructs')
xlabel({'-GSS-cyOFP','Relative Brightness (a.u.)'})
ylabel({'-P2A-mCherry-CAAX','Relative Brightness (a.u.)'})
axis square
set(gcf, 'Position',  [10, 10, 400,400])

figure(11)
clf
hold on
% subplot(1,2,2)
cmap = colorcube(75);
%AUCMean = normalize(T_well_stats.AUCMean,'range');
% scatter(Stim5Mean',StimShortMean',0.3*AUCMean,normalize(GRRatio,'range'),'filled','MarkerFaceAlpha',0.5)
for i = 1:min(height(B1_TwellStats),height(B2_TwellStats))
    plot(B1_TwellStats(B1toB2(i),'GRRatioMean').Variables,...
        B2_TwellStats.GRRatioMean(i),'.','Color',cmap(i*2-1,:),'MarkerSize',25)%,'MarkerEdgeColor','none')
end
% scatter(B1_TwellStats(B1toB2,'GRRatioMean').Variables',B2_TwellStats.GRRatioMean',...
%     [],cmap(1:2:2*length(StimShortMean),:),'filled','MarkerFaceAlpha',1)
% 
% scatter(B1_TwellStats(B1toB2,'GRRatioMean').Variables,...
%     B2_TwellStats.GRRatioMean,[],[0.4660, 0.6740, 0.1880],'filled')
lnrFitting(B1_TwellStats(B1toB2,'GRRatioMean').Variables,...
    B2_TwellStats.GRRatioMean,false,true,11)
title('Compare GR ratio from different constructs')
xlabel({'-GSS-cyOFP','Relative Brightness (a.u.)'})
ylabel({'-P2A-mCherry-CAAX','Relative Brightness (a.u.)'})
axis equal
set(gcf, 'Position',  [10, 10, 400,400])
legend(B2_TwellStats.Name,'location','eastoutside')
set(gcf, 'Position',  [10, 10, 400,400])

%% Compare manual roi stats and pipeline well stats
% preprocess
T_groupFromROI = sortrows(T_groupFromROI,'Name','ascend');
T_groupFromROI.StimShortMean = nanmean([T_groupFromROI.Stim_1,T_groupFromROI.Stim_2,T_groupFromROI.Stim_3],2)

%%
figure(1)
clf
hold on
scatter(T_well_stats.StimShortMean,T_groupFromROI.StimShortMean,'filled')
scatter(T_well_stats.Stim8Mean,T_groupFromROI.Stim_8,'filled')
[k1,R21] = libplot.lnrFitting(T_well_stats.StimShortMean,...
    T_groupFromROI.StimShortMean,false,true,1)
[k2,R22] = libplot.lnrFitting(T_well_stats.Stim8Mean,...
    T_groupFromROI.Stim_8,false,true,1)
xlabel('Pipeline analysis dF/F0')
ylabel('Manaual mask dF/F0')
legend('Mean 1ms stim','Mean 100ms stim')
axis square
