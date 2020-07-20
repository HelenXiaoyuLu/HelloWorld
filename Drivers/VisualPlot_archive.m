%% Clear 
clearvars -except T_groupFromROI
%% Specify the groups
L1 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\1pFstimLighton_20200202_GEVI_screening_cyOFPset_benchmarking_20X_objective.mat');
L2 = load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_plate2\20190227 Fstim platform benchmarking 20x plate2 1khzFstim nonTriggered_20200411171118_quant.mat');
savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200202_benchmarking_plate1\20200202_benchmarking_20X_update20200510';
savepathsplit = split(folderpath,'\');
saveName = savepathsplit{end};
% T_wellArchive = readtable('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\Merge\20190227 1P Fstim 1000Hz archive.xlsx','Sheet','Sheet1');
T_well = [L1.T_well;L2.T_well];
% Create standard name column
WellconstructName = cellstr(T_well.ConstructName);
for i = 1:length(WellconstructName)
    Namesplt = split(WellconstructName{i},'.');
    Namesplt2 = split(Namesplt(1),'-cyOFP');
    T_well.Name(i) = string(Namesplt2(1));
end
% Standardize Names
idx = strmatch('JEDI1P',T_well.Name,'exact');
nameCol = strmatch('Name',T_well.Properties.VariableNames,'exact');
T_well(idx,nameCol) = repmat({'JEDI-1P'},length(idx),1)
% replace missing values (nan) with 0
T_well = fillmissing(T_well,'constant',0,'DataVariables',@isnumeric);
%% Well level stats
T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3','Stim_4'}}, 2); 
T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3'}}, 2); 
T_well.stim_middlemean = mean(T_well{:, {'Stim_4', 'Stim_5' ,'Stim_6'}}, 2); 
    
%% Group Level stats
T_well_stats = table();
Wellgroups = findgroups(T_well.Name);
T_well_stats.Name = splitapply(@unique, T_well.Name, Wellgroups);
T_well_stats.n_well = splitapply(@numel, T_well.Name, Wellgroups);
T_well_stats.Stim1Mean = splitapply(@mean, T_well.Stim_1, Wellgroups);
T_well_stats.Stim1Std = splitapply(@std, T_well.Stim_1, Wellgroups);
T_well_stats.Stim2Mean = splitapply(@mean, T_well.Stim_2, Wellgroups);
T_well_stats.Stim2Std = splitapply(@std, T_well.Stim_2, Wellgroups);
T_well_stats.Stim3Mean = splitapply(@mean, T_well.Stim_3, Wellgroups);
T_well_stats.Stim3Std = splitapply(@std, T_well.Stim_3, Wellgroups);
T_well_stats.Stim4Mean = splitapply(@mean, T_well.Stim_4, Wellgroups);
T_well_stats.Stim4Std = splitapply(@std, T_well.Stim_4, Wellgroups);
T_well_stats.Stim5Mean = splitapply(@mean, T_well.Stim_5, Wellgroups);
T_well_stats.Stim5Std = splitapply(@std, T_well.Stim_5, Wellgroups);
%%%%%% Alternatively needed for Xiaoyu-per-2 9 stimulation protocol
T_well_stats.Stim6Mean = splitapply(@mean, T_well.Stim_6, Wellgroups);
T_well_stats.Stim6Std = splitapply(@std, T_well.Stim_6, Wellgroups);
T_well_stats.Stim7Mean = splitapply(@mean, T_well.Stim_7, Wellgroups);
T_well_stats.Stim7Std = splitapply(@std, T_well.Stim_7, Wellgroups);
T_well_stats.Stim8Mean = splitapply(@mean, T_well.Stim_8, Wellgroups);
T_well_stats.Stim8Std = splitapply(@std, T_well.Stim_8, Wellgroups);
T_well_stats.Stim9Mean = splitapply(@mean, T_well.Stim_9, Wellgroups);
T_well_stats.Stim9Std = splitapply(@std, T_well.Stim_9, Wellgroups);
T_well_stats.StimMidMean = splitapply(@mean, T_well.stim_middlemean, Wellgroups);
T_well_stats.StimMidStd = splitapply(@std, T_well.stim_middlemean, Wellgroups);
%%%%%%
T_well_stats.StimShortMean = splitapply(@mean, T_well.stim_shortmean, Wellgroups);
T_well_stats.StimShortStd = splitapply(@std, T_well.stim_shortmean, Wellgroups);
T_well_stats.GRRatioMean = splitapply(@mean, T_well.GRRatioMean, Wellgroups);
T_well_stats.GRRatioStd = splitapply(@std, T_well.GRRatioMean, Wellgroups);
T_well_stats.AUCMean = splitapply(@mean, T_well.AOCMean, Wellgroups);
T_well_stats.AUCStd = splitapply(@std, T_well.AOCMean, Wellgroups);
T_well_stats.KineticsMean = splitapply(@mean, T_well.stim_shortmean./T_well.Stim_8, Wellgroups);
T_well_stats.KineticsStd = splitapply(@std, T_well.stim_shortmean./T_well.Stim_8, Wellgroups);

%%%%% Alternatively, for Eric's 2p library
T_well_stats.TgtBrightnessMean = splitapply(@mean, T_well.TgtBrightnessMean, Wellgroups);
T_well_stats.TgtBrightnessStd = splitapply(@std, T_well.TgtBrightnessMean, Wellgroups);
T_well_stats.RefBrightnessMean = splitapply(@mean, T_well.RefBrightnessMean, Wellgroups);
T_well_stats.RefBrightnessStd = splitapply(@std, T_well.RefBrightnessMean, Wellgroups);
T_well_stats.StimShortMean = splitapply(@mean, T_well.dFF0_shortStim, Wellgroups);
T_well_stats.StimShortStd = splitapply(@std, T_well.dFF0_shortStim, Wellgroups);
T_well_stats.StimLongMean = splitapply(@mean, T_well.dFF0_longStim, Wellgroups);
T_well_stats.StimLongStd = splitapply(@std, T_well.dFF0_longStim, Wellgroups);
T_well_stats.GRRatioMean = splitapply(@mean, T_well.GRratioMean, Wellgroups);
T_well_stats.GRRatioStd = splitapply(@std, T_well.GRratioMean, Wellgroups);
T_well_stats.AUCMean = splitapply(@mean, T_well.TgtAUC_sig_NormaltoRandT, Wellgroups);
T_well_stats.AUCStd = splitapply(@std, T_well.TgtAUC_sig_NormaltoRandT, Wellgroups);
T_well_stats.KineticsMean = splitapply(@mean, T_well.dFF0_shortStim./T_well.dFF0_longStim, Wellgroups);
T_well_stats.KineticsStd = splitapply(@std, T_well.dFF0_shortStim./T_well.dFF0_longStim, Wellgroups);
T_well_stats.Bleachingratio_totalMean = splitapply(@mean, T_well.Bleachingratio_total, Wellgroups);
T_well_stats.Bleachingratio_totalStd = splitapply(@std, T_well.Bleachingratio_total, Wellgroups);
T_well_stats.Bleachingratio_2_1sMean = splitapply(@mean, T_well.Bleachingratio_2_1s, Wellgroups);
T_well_stats.Bleachingratio_2_1sStd = splitapply(@std, T_well.Bleachingratio_2_1s, Wellgroups);
%%%%%%
T_well_stats = sortrows(T_well_stats,'Name','ascend');
save(fullfile(savepath,strcat(saveName,'_T_well_stats.mat')),'T_well_stats')

%% Alternatively, call function
T_well_stats = libplot.groupWells(T_well);

%% Trace (Figure 0)
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
t.Padding = 'none';
title(t, 'Well Stats in Construct');
n_wells = [];
T_well = sortrows(T_well,'Name','ascend');
for i = 1:height(T_well_stats)
    n_wells = [n_wells; ones(T_well_stats.n_well(i),1)*T_well_stats.n_well(i)];
end
T_well_stats.TraceMean = splitapply(@mean, T_well.TraceMean, n_wells, Wellgroups);
% std of FOV traces per construct
T_well_stats.TraceStd = splitapply(@std, T_well.TraceMean, n_wells, Wellgroups);

for i = 1:height(T_well_stats)
    ax = nexttile(t);
    xf = [T_well_stats.TraceMean(i).x; flipud(T_well_stats.TraceMean(i).x)];
    yf = [T_well_stats.TraceMean(i).y + T_well_stats.TraceStd(i).y; 
          flipud(T_well_stats.TraceMean(i).y - T_well_stats.TraceStd(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well_stats.TraceMean(i).x, T_well_stats.TraceMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well_stats.Name(i),T_well_stats.n_well(i)));
end

%% Group Level Barplot (Figure 1-3, 12)
SortOn = 'AUCMean'; % Input Sorting Criteria 
Direction = 'descend';
LongStim = 'StimLong';
[T_well_stats_sorted I]= sortrows(T_well_stats,SortOn,Direction,'MissingPlacement','last');

figure(1) % Relative Brightness
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,T_well_stats_sorted.GRRatioMean,'BarWidth',0.8,'LineWidth', 1);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.GRRatioStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel('GR Ratio')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
% title({'20190202 1P GR ratio','well-level stat New Pipeline'})
set(gcf, 'Position',  [110, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_GRRatio.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_GRRatio.fig'))

figure(2) % Response amplitude
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,[T_well_stats_sorted.StimShortMean,T_well_stats_sorted(:, strcat(LongStim,'Mean')).Variables],'grouped','BarWidth',1,'LineWidth', 1);
Width = hBar(2).XEndPoints - hBar(1).XEndPoints;
errorbar(hBar(1).XEndPoints-0.00*Width,hBar(1).YData,T_well_stats_sorted.StimShortStd,'.k')
errorbar(hBar(2).XEndPoints-0.00*Width,hBar(2).YData,T_well_stats_sorted(:, strcat(LongStim,'Std')).Variables,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Response Amplitude','[\DeltaF/F0]'})
legend('1.0 ms','100 ms','location','best')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.8500, 0.3250, 0.0980];	
hBar(2).FaceColor = [0.9290, 0.6940, 0.1250];	
set(gcf, 'Position',  [710, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0.fig'))

figure(3) % Photostability
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,T_well_stats_sorted.AUCMean,'BarWidth',0.8,'LineWidth', 1,'FaceAlpha', 0.5);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.AUCStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Photostability','[Area-Under-Curve]'})
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
set(gcf, 'Position',  [1310, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_AUC.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_AUC.fig'))

figure(12) %Kinetics bar
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,T_well_stats_sorted.KineticsMean,'BarWidth',0.8,'LineWidth', 1,'FaceAlpha', 0.5);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.KineticsStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Kinetics','dFF0_{1ms}/dFF0_{100ms}'})
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0, 0.4470, 0.7410];	
% title({'20190202 1P photostabillity','well-level stat New Pipeline'})
set(gcf, 'Position',  [1310, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_kinetics_and_dFF0.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_kinetics_and_dFF0.fig'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_kinetics_and_dFF0.png'))

figure(13) % experiment comprehesive figure
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,[T_well_stats_sorted.StimShortMean,T_well_stats_sorted(:, strcat(LongStim,'Mean')).Variables],'grouped','BarWidth',1,'LineWidth', 1);
Width = hBar(2).XEndPoints - hBar(1).XEndPoints;
errorbar(hBar(1).XEndPoints-0.00*Width,hBar(1).YData,T_well_stats_sorted.StimShortStd,'.k')
errorbar(hBar(2).XEndPoints-0.00*Width,hBar(2).YData,T_well_stats_sorted(:, strcat(LongStim,'Std')).Variables,'.k')
errorbar(X,T_well_stats_sorted.KineticsMean,T_well_stats_sorted.KineticsStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Response Amplitude','[\DeltaF/F0]'})
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.8500, 0.3250, 0.0980];	
hBar(2).FaceColor = [0.9290, 0.6940, 0.1250];
scatter(X,T_well_stats_sorted.KineticsMean,[],[0, 0.4470, 0.7410],'filled')
set(gcf, 'Position',  [710, 100, 600, 400])
% title({'20190202 1P Fstim','well-level stat New Pipeline'})
C = get(gca).Children;
legend([C(6),C(5),C(1)],'1.0 ms','100 ms','Kinetics','location','best')
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_1PdFF0.fig'))

figure(14)
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,[T_well_stats_sorted.TgtBrightnessMean,T_well_stats_sorted.RefBrightnessMean],'grouped','BarWidth',1,'LineWidth', 1);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.TgtBrightnessStd,'.k')
errorbar(hBar(2).XEndPoints,hBar(2).YData,T_well_stats_sorted.RefBrightnessStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Brightness (a.u.)'})
legend('Green (Tgt)','Red (Ref)','location','best')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
hBar(2).FaceColor = [0.8500, 0.3250, 0.0980];	
set(gcf, 'Position',  [710, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_AbsBrightness.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_AbsBrightness.fig'))

figure(15) %BleachRatio
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,[T_well_stats_sorted.Bleachingratio_totalMean,T_well_stats_sorted.Bleachingratio_2_1sMean],'grouped','BarWidth',1,'LineWidth', 1,'FaceAlpha', 0.5);
errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.Bleachingratio_totalStd,'.k')
errorbar(hBar(2).XEndPoints,hBar(2).YData,T_well_stats_sorted.Bleachingratio_2_1sStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Brightness (a.u.)'})
legend('Bleachingratio total Mean','Bleachingratio 2.1s Mean','location','best')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0.4660, 0.6740, 0.1880];	
hBar(2).FaceColor = [0, 0.5, 0];	
set(gcf, 'Position',  [710, 100, 600, 400])
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_BleachRatio.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_BleachRatio.fig'))

figure(16)
clf
hold on
NameCell = cellstr(T_well_stats_sorted.Name);
X = 1:height(T_well_stats_sorted);
hBar = bar(X,[T_well_stats_sorted.AUCMean./T_well_stats_sorted.AUCMean(4),...
    T_well_stats_sorted.Bleachingratio_totalMean./T_well_stats_sorted.Bleachingratio_totalMean(4),...
    T_well_stats_sorted.Bleachingratio_2_1sMean./T_well_stats_sorted.Bleachingratio_2_1sMean(4)],...
    'grouped','BarWidth',1,'LineWidth', 1,'FaceAlpha', 0.8);
% errorbar(hBar(1).XEndPoints,hBar(1).YData,T_well_stats_sorted.TgtBrightnessStd,'.k')
% errorbar(hBar(2).XEndPoints,hBar(2).YData,T_well_stats_sorted.Bleachingratio_totalStd,'.k')
% errorbar(hBar(3).XEndPoints,hBar(3).YData,T_well_stats_sorted.Bleachingratio_2_1sStd,'.k')
set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell))
ylabel({'Brightness (a.u.)'})
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
hBar(1).FaceColor = [0, 0.4470, 0.7410];	
hBar(2).FaceColor = [0.8500, 0.3250, 0.0980]	      ;	
hBar(3).FaceColor = [0.9290, 0.6940, 0.1250];
set(gcf, 'Position',  [710, 100, 600, 400])
plot(X,ones(size(X)),'--k')
title('Normalized to ASAP2s')
legend('AUC Mean','Bleachingratio total Mean','Bleachingratio 2.1s Mean','location','best')
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_PhotostabilityNorm.svg'))
saveas(gcf,strcat(savepath,'\',saveName,'_Sorted_on_',SortOn,'_PhotostabilityNorm.fig'))

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
Wellgroups = findgroups(T_wellArchive.Name);
T_well_stats_archive.Name = splitapply(@unique, T_wellArchive.Name, Wellgroups);
T_well_stats_archive.n_well = splitapply(@numel, T_wellArchive.Name, Wellgroups);
T_well_stats_archive.Stim1Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S1_, Wellgroups);
T_well_stats_archive.Stim1Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S1_, Wellgroups);
T_well_stats_archive.Stim2Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S2_, Wellgroups);
T_well_stats_archive.Stim2Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S2_, Wellgroups);
T_well_stats_archive.Stim3Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S3_, Wellgroups);
T_well_stats_archive.Stim3Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S3_, Wellgroups);
T_well_stats_archive.Stim4Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S4_, Wellgroups);
T_well_stats_archive.Stim4Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S4_, Wellgroups);
T_well_stats_archive.Stim5Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S5_, Wellgroups);
T_well_stats_archive.Stim5Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S5_, Wellgroups);
T_well_stats_archive.StimShortMean = splitapply(@mean, T_wellArchive.stim_shortmean, Wellgroups);
T_well_stats_archive.StimShortStd = splitapply(@std, T_wellArchive.stim_shortmean, Wellgroups);
T_well_stats_archive.GRRatioMean = splitapply(@mean, T_wellArchive.RatioOfMeanBrightness_C1_C3_, Wellgroups);
T_well_stats_archive.GRRatioStd = splitapply(@std, T_wellArchive.RatioOfMeanBrightness_C1_C3_, Wellgroups);
%%%%%%%%%
T_well_stats_archive.Stim6Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S6_, Wellgroups);
T_well_stats_archive.Stim6Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S6_, Wellgroups);
T_well_stats_archive.Stim7Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S7_, Wellgroups);
T_well_stats_archive.Stim7Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S7_, Wellgroups);
T_well_stats_archive.Stim8Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S8_, Wellgroups);
T_well_stats_archive.Stim8Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S8_, Wellgroups);
T_well_stats_archive.Stim9Mean = splitapply(@mean, T_wellArchive.dFF0OfMean_C1S9_, Wellgroups);
T_well_stats_archive.Stim9Std = splitapply(@std, T_wellArchive.dFF0OfMean_C1S9_, Wellgroups);
T_well_stats_archive.StimMidMean = splitapply(@mean, T_wellArchive.stim_middlemean, Wellgroups);
T_well_stats_archive.StimMidStd = splitapply(@std, T_wellArchive.stim_middlemean, Wellgroups);
%%%%%%%%%
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
[k1,R21] = lnrFitting(T_well_stats.StimShortMean,...
    T_groupFromROI.StimShortMean,false,true,1)
[k2,R22] = lnrFitting(T_well_stats.Stim8Mean,...
    T_groupFromROI.Stim_8,false,true,1)
xlabel('Pipeline analysis dF/F0')
ylabel('Manaual mask dF/F0')
legend('Mean 1ms stim','Mean 100ms stim')
axis square
%% Functions
function [k, R2] = lnrFitting(x,y,intercept,figureTrue,figNumber)
    p = fitlm(x,y,'Intercept',intercept);   
    k = p.Coefficients(1,1).Variables;
    R2 = p.Rsquared.Adjusted;
    if figureTrue
        figure(figNumber)
        if min(x) <0
            linestart = min(min(x),min(y));
        else
            linestart = 0;
        end
        xplot = linestart:0.01:max(max(x),max(y));        
        plot(xplot,xplot*p.Coefficients(1,1).Variables,'--k')
        text(linestart+0.01,max(y),{strcat('k = ',num2str(k)),strcat('R2 = ',num2str(R2))},'Color',[0 0 0])
    end
end