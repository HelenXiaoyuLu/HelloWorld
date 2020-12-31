%% load data
clc
clear
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking\plate2';
j1.output = load(fullfile(dirp, "2pFstim 20200126 T_out"));
j2.output = load(fullfile(dirp, "2PFieldStimCorrFilt 20200126 20201101005346 (P2)"));
j3.output = load(fullfile(dirp, "Wellbase_2PFieldStimIlastic 20200126 T_out"));
j4.output = load(fullfile(dirp, "2pFstimSingleCell 20200126 20200902084612"));

%% pull stats from T_fov
[~, T_group1] = Wellbase_2PFieldStim.Scatter_rand2d(j1);
[~, T_group2] = Wellbase_2PFieldStim.Scatter_rand2d(j2);
[~, T_group3] = Wellbase_2PFieldStim.Scatter_rand2d(j3);
T_group4 = j4.output.T_group;
T_group4.Properties.RowNames = T_group4.Group_Name;
T_group4.Properties.VariableNames{11} = 'Brightness Mean';
T_group4.Properties.VariableNames{12} = 'Brightness STD';
% compute dFF0 from group stats
mu = zeros(height(T_group4), 1);
sigma = zeros(height(T_group4), 1);
for i = 1:height(T_group4)
    [mu(i), l]= img.extrema(T_group4.TraceShortMean(i).y, T_group4.TraceShortMean(i).crop(-Inf, 0.99).y, 1);
    sigma(i) = T_group4.TraceShortStd(i).y(l);
end
T_group4.("-dF/F0 Short Mean") = 1 - mu; % -dF/F0
T_group4.("dF/F0 Short STD") = sigma; % -dF/F0  

%% Figure 1: bar plot
cmap = lines;
sel = ["ASAP1-EGFP-cyOFP", "ASAP1-dpOPT-cyOFP", "ASAP1-cyOFP", "ASAP2s-cyOFP", "ASAP3-cyOFP", ...
    "ASAP2s-T207H-cyOFP", "ASAP2s-H152E-cyOFP", "JEDI-1P-cyOFP", "JEDI-2P-cyOFP"];
ppmean = "-dF/F0 Short Mean";
ppstd = "dF/F0 Short STD";

T_group1_sel = sortrows(T_group1(sel,:), ppmean, 'ascend');
T_group2_sel = T_group2(T_group1_sel.Name,:);
T_group3_sel = T_group3(T_group1_sel.Name,:);
T_group4_sel = T_group4(T_group1_sel.Name,:);

y1 = T_group1_sel.(ppmean);
y2 = T_group2_sel.(ppmean);
y3 = T_group3_sel.(ppmean);
y4 = T_group4_sel.(ppmean);
y1std = T_group1_sel.(ppstd);
y2std = T_group2_sel.(ppstd);
y3std = T_group3_sel.(ppstd);
y4std = T_group4_sel.(ppstd);

f = figure('Color', [1, 1, 1]);
f.UserData.singleval = T_group1_sel;   % bind data
f.UserData.corrFilt = T_group2_sel;
f.UserData.Ilastik = T_group3_sel;
f.UserData.Manual = T_group4_sel;
f.UserData.pp = [ppmean, ppstd];

t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
b = bar(ax,[y1, y2, y3, y4]);
b(1).FaceColor = cmap(1,:);
b(1).DisplayName = "singleval";
b(2).FaceColor = cmap(2,:);
b(2).DisplayName = "corrFilt";
b(3).FaceColor = cmap(3,:);
b(3).DisplayName = "Ilastik";
b(4).FaceColor = cmap(5,:);
b(4).DisplayName = "Manual";
y1pos = b(1).XEndPoints;
e(1) = errorbar(ax, y1pos, y1, y1std, ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', 'singleval');
y2pos = b(2).XEndPoints;
e(2) = errorbar(ax, y2pos, y2, y2std,...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', 'corrFilt');
y3pos = b(3).XEndPoints;
e(3) = errorbar(ax, y3pos, y3, y3std,...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', 'Ilastik');
y4pos = b(4).XEndPoints;
e(4) = errorbar(ax, y4pos, y4, y4std,...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', 'Manual');
legend(ax, b, 'Location', 'northwest');
ax.YAxis.Label.String = ppmean;
ax.Box = false;
ax.XAxis.TickLabels = T_group1_sel.Name;
ax.XAxis.TickLabelRotation = 90;
ax.XAxis.Limits = [0.4, height(T_group1_sel) + 0.6];

% Data menu
mData = uimenu('Parent', f, ...
    'Tag', 'Data', ...
    'Text', 'Data');
uimenu('Parent', mData, ...
    'Text', 'Set Reference Variant', ...
    'Tag', 'SetRef', ...
    'MenuSelectedFcn', @setReference);
    




%%
function setReference(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    pp = f.UserData.pp;
    ppmean = pp(1); 
    ppstd = pp(2);
    T_plots = f.UserData;
    T_plot_cell = struct2cell(T_plots);
    T_plot = T_plot_cell{1};
    list = T_plot.Name;
    list = cat(1, "(No Reference)", list);
    [idx, tf] = listdlg('ListString', list, ...
        'Name', 'Reference Variant', ...
        'SelectionMode', 'single', ...
        'ffs', 0, ...
        'ListSize', [250, 160], ...
        'PromptString', 'Please select a reference variant: ');
    if tf == 0  % cancel
        return
    end
    b = findobj(f, '-class', 'matlab.graphics.chart.primitive.Bar');
    e = findobj(f, '-class', 'matlab.graphics.chart.primitive.ErrorBar');
    for i = 1:numel(b)
        mtd =  b(i).DisplayName;
        T_plot = T_plots.(mtd);
        if idx == 1     % original
            normValue = 1;
        else
            normValue = T_plot.(ppmean)(idx - 1);
        end
        b(i).YData = T_plot.(ppmean) ./ normValue;
        e(i).YData = T_plot.(ppmean) ./ normValue;
        e(i).YNegativeDelta = T_plot.(ppstd) ./ normValue;
        e(i).YPositiveDelta = T_plot.(ppstd) ./ normValue;
    end
end