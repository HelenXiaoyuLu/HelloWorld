%% Load T_plots 
clear
clc
% 2pFstim
% dirp1 = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201014 Evol path';
% dirp2 = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201008 double blind\Benchmarking repeat';
% dirp3 = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200814 benchmarking cyOFP and mCherry';
% dirp4 = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking';
% j1 = load(fullfile(dirp1, '2pFstim 20201014 20201020190736 T_plot T_plot_well'));
% j2 = load(fullfile(dirp2, '2pFstim 20201009 20201009191355 T_plot grouped'));
% j3 = load(fullfile(dirp3, '2pFstim 20200814 20200929223023 T_plot grouped'));
% j4 = load(fullfile(dirp4, 'plate1\2pFstim 20200126 20200814110823 T_plot T_plot_well'));
% j5 = load(fullfile(dirp4, 'plate2\2pFstim 20200126 20200814152346 T_plot T_plot_well'));
% 
% varname = ["ASAP1-EGFP"; "ASAP1"; "ASAP2s"; "ASAP3"; "ASAP2s-H152E"; "JEDI-2P"];
% T = {};
% T{1} = j1.T_plot(varname,:);
% T{2} = j2.T_plot_cyOFP(varname,:);
% T{3} = j3.T_plot_cyOFP(varname,:);
% T{4} = j4.T_plot(varname + "-cyOFP",:);
% T{5} = j5.T_plot(varname + "-cyOFP",:);

%
% 2pFstim
% dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\Compare Benchmakring\corrFilt';
% j1.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20201014 20201029200814 T_out'));
% j2.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20201008 20201029200601 T_out (P1)'));
% j3.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20200814 20201029200356 T_out'));
% j4.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20200126 20201029200154 T_out (P1)'));
% j5.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20200126 20201029201318 T_out (P2)'));

dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\Compare Benchmakring\corrFilt_ver2';
j1.output = load(fullfile(dirp, '2PFieldStimCorrFilt 20201014 20201101004924'));
j2.output = load(fullfile(dirp, '2PFieldStimCorrFilt 20201008 20201101004745 (P1)'));
% j3.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20200814 20201029200356 T_out'));
j4.output = load(fullfile(dirp, '2PFieldStimCorrFilt 20200126 20201101004301 (P1)'));
j5.output = load(fullfile(dirp, '2PFieldStimCorrFilt 20200126 20201101005346 (P2)'));

varname = ["ASAP1-EGFP"; "ASAP1"; "ASAP2s"; "ASAP3"; "ASAP2s-H152E"; "JEDI-2P"];
T = {};
[~, T{1}] = Wellbase_2PFieldStim.Scatter_rand2d(j1);
[~, T{2}] = Wellbase_2PFieldStim.Scatter_rand2d(j2);
% [~, T{3}] = Wellbase_2PFieldStim.Scatter_rand2d(j3);
[~, T{3}] = Wellbase_2PFieldStim.Scatter_rand2d(j4);
[~, T{4}] = Wellbase_2PFieldStim.Scatter_rand2d(j5);
for i = 1:4
    T{i}.Properties.RowNames = T{i}.Name;
    if i == 1
        T{i} = T{i}(varname,:);
    else
        T{i} = T{i}(varname + "-cyOFP",:);
    end   
end

%% Compare experiments
ppmean = "-dF/F0 Short Mean";
ppstd = "dF/F0 Short STD";
y = zeros(numel(varname), numel(T));
for i = 1:numel(T)
    T_plot = T{i};
    pp = T_plot.Properties.VariableNames{contains(T_plot.Properties.VariableNames, ppmean)};
    y(:,i) = T_plot.(pp);
end
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
b = bar(ax,y);
xpos = b(3).XEndPoints;
for i = 1:numel(T)
    T_plot = T{i};
    x = b(i).XEndPoints;
    y = T_plot.(ppmean);
    pp = T_plot.Properties.VariableNames{contains(T_plot.Properties.VariableNames, ppstd)};
    xstd = T_plot.(pp);
    e(i) = errorbar(ax, x, y, xstd, ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
end
l = legend(ax, b, {"20201014", "20201009",...
    "20200814", "20200126 (P1)", ...
    "20200126 (P2)"}, 'Location', 'eastoutside');
l.Title.String = "Mean \pm Standard Deviation";
ax.YAxis.Label.String = ppmean;
ax.Box = false;
ax.XAxis.TickValues = xpos;
ax.XAxis.TickLabels = varname;
ax.XAxis.TickLabelRotation = 90;

%% %% Compare experiments normalized to ASAP1
y = zeros(numel(varname), numel(T));
for i = 1:numel(T)
    T_plot = T{i};
    normidx = contains(T_plot.Properties.RowNames, "ASAP1") & ...
        ~ contains(T_plot.Properties.RowNames, "EGFP"); 
    normval = T_plot.(ppmean)(normidx);
    y(:,i) = T_plot.(ppmean)./ normval;
end
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
b = bar(ax,y);
xpos = b(3).XEndPoints;
for i = 1:numel(T)
    T_plot = T{i};
    normidx = contains(T_plot.Properties.RowNames, "ASAP1") & ...
        ~ contains(T_plot.Properties.RowNames, "EGFP"); 
    x = b(i).XEndPoints;
    normval = T_plot.(ppmean)(normidx);
    y = T_plot.(ppmean)./ normval;
    pp = T_plot.Properties.VariableNames{contains(T_plot.Properties.VariableNames, ppstd)};
    xstd = T_plot.(pp)./ normval;
    e(i) = errorbar(ax, x, y, xstd, ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
end
l = legend(ax, b, {"20201014", "20201009",...
    "20200814", "20200126 (P1)", ...
    "20200126 (P2)"}, 'Location', 'eastoutside');
l.Title.String = "Mean \pm Standard Deviation";
ax.YAxis.Label.String = {ppmean; "normalized to ASAP1"};
ax.Box = false;
ax.XAxis.TickValues = xpos;
ax.XAxis.TickLabels = varname;
ax.XAxis.TickLabelRotation = 90;