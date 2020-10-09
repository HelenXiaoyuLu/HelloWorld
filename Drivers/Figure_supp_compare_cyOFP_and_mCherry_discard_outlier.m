%% Load jobs 
jobids = [20200831224945, ...  % Reference channel = 1 (green)
          20200929223023, ...  % Reference channel = 2 (red)
          20201002184433];     % 1P Photobleaching (trace only)     
j = db.jobDynamic.pull('ID', jobids(2));

%% Preprocessing
% Well level stats
% [~, T_plot_well] = Wellbase_2PFieldStim.Scatter_rand2dwell(j);
T_plot_well.GEVI = splitapply(@parsename, T_plot_well.Name(:)', 1:height(T_plot_well))';
T_plot_well.refFP = splitapply(@parseref, T_plot_well.Name(:)', 1:height(T_plot_well))';
T_plot_cyOFP_well = T_plot_well(strcmp(T_plot_well.refFP, "cyOFP"),:);
T_plot_mCherry_well = T_plot_well(strcmp(T_plot_well.refFP, "mCherry-CAAX"),:); 
T_plot_mCherry_well(strcmp(T_plot_mCherry_well.GEVI, "ASAP2s-H152E-Q397H"),:) = [];

% Discard outliers
T = T_plot_well;
grouping = findgroups(T.Name);
[T_plot_filtered, T_plot_discard] = discardoutliers(T,grouping, ...
    'method', 'mean',...
    'DataVariables',{'-dF/F0 Short Mean', 'Brightness Mean', 'Photostability Mean'});

% pull group-level stats
weightedMean = @(x, w) sum(x.*w, 'all')./sum(w, 'all');
grouping = findgroups(T_plot_filtered.Name);
T_group = table();
T_group.Name = splitapply(@unique, T_plot_filtered.Name, grouping);
T_group.N = splitapply(@numel, T_plot_filtered.Name, grouping);
T_group.Area = splitapply(@sum, T_plot_filtered.Area, grouping);
T_group.("Brightness Mean") = splitapply(weightedMean, ...
    T_plot_filtered.("Brightness Mean"), T_plot_filtered.Area, grouping);
T_group.("Brightness STD") = splitapply(@std, ...
    T_plot_filtered.("Brightness Mean"), T_plot_filtered.Area, grouping);
T_group.("-dF/F0 Short Mean") = splitapply(weightedMean, ...
    T_plot_filtered.("-dF/F0 Short Mean"), T_plot_filtered.Area, grouping);
T_group.("-dF/F0 Short STD") = splitapply(@std, ...
    T_plot_filtered.("-dF/F0 Short Mean"), T_plot_filtered.Area, grouping);
T_group.("Photostability Mean") = splitapply(weightedMean, ...
    T_plot_filtered.("Photostability Mean"), T_plot_filtered.Area, grouping);
T_group.("Photostability STD") = splitapply(@std, ...
    T_plot_filtered.("Photostability Mean"), T_plot_filtered.Area, grouping);
T_group.("Photostability Mean") = splitapply(weightedMean, ...
    T_plot_filtered.("Photostability Mean"), T_plot_filtered.Area, grouping);
T_group.("Photostability STD") = splitapply(@std, ...
    T_plot_filtered.("Photostability Mean"), T_plot_filtered.Area, grouping);
T_group.("Detectability Index") = splitapply(weightedMean, ...
    T_plot_filtered.("Detectability Index"), T_plot_filtered.Area, grouping);
T_group.("Detectability Index STD") = splitapply(@std, ...
    T_plot_filtered.("Detectability Index"), T_plot_filtered.Area, grouping);
T_group.("AUC unnormalized") = splitapply(weightedMean, ...
    T_plot_filtered.("AUC unnormalized"), T_plot_filtered.Area, grouping);
T_group.("AUC unnormalized STD") = splitapply(@std, ...
    T_plot_filtered.("AUC unnormalized"), T_plot_filtered.Area, grouping);
T_group.GEVI = splitapply(@parsename, T_group.Name(:)', 1:height(T_group))';
T_group.refFP = splitapply(@parseref, T_group.Name(:)', 1:height(T_group))';

% Seperate mCherry and cyOFP
T_plot_cyOFP = T_group(strcmp(T_group.refFP, "cyOFP"),:);
T_plot_mCherry = T_group(strcmp(T_group.refFP, "mCherry-CAAX"),:); 
T_plot_cyOFP.Properties.RowNames = T_plot_cyOFP.GEVI;
T_plot_mCherry.Properties.RowNames = T_plot_mCherry.GEVI;
T_plot_mCherry("ASAP2s-H152E-Q397H",:) = [];

%% Figure 1: compare dF/F0
ppmean = "-dF/F0 Short Mean";
ppstd = "-dF/F0 Short STD";
T_plot_mCherry = sortrows(T_plot_mCherry, ppmean, 'ascend');
T_plot_cyOFP = T_plot_cyOFP(T_plot_mCherry.GEVI,:);
x = T_plot_mCherry.(ppmean);
y = T_plot_cyOFP.(ppmean);
xstd = T_plot_mCherry.(ppstd);
ystd = T_plot_cyOFP.(ppstd);

f = figure('Color', [1, 1, 1]);
f.UserData.PlotData1 = T_plot_cyOFP;   % bind data
f.UserData.PlotData2 = T_plot_mCherry;
% db.watermark(f, j);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
b = bar(ax,[x, y]);
b(1).FaceColor = [0.8500 0.3250 0.0980];
b(1).DisplayName = "mCherry-CAAX";
b(2).FaceColor = [0.9290 0.6940 0.1250];
b(2).DisplayName = "cyOFP";
xpos = b(1).XEndPoints;
e(1) = errorbar(ax, xpos, x, xstd, ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
ypos = b(2).XEndPoints;
e(2) = errorbar(ax, ypos, y, ystd,...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
legend(ax, b, 'Location', 'northwest');
ax.YAxis.Label.String = ppmean;
ax.Box = false;
ax.XAxis.TickLabels = T_plot_mCherry.GEVI;
ax.XAxis.TickLabelRotation = 90;
ax.XAxis.Limits = [0.4, height(T_plot_mCherry) + 0.6];
 
%% Figure 2: scatter compare dF/F0
f = figure('Color', [1, 1, 1]);
f.UserData.PlotData1 = T_plot_cyOFP;   % bind data
f.UserData.PlotData2 = T_plot_mCherry;
% db.watermark(f, j);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
errorbar(x, y, ystd, ystd, xstd, xstd,...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 6, ...
    'Color', 'black', ...
    'DisplayName', '\mu\pm\sigma');
cmap = lines;
for i = 1:numel(x)
    s(i) = scatter(ax, x(i),y(i),[], cmap(i,:),'filled',...
        'DisplayName', T_plot_mCherry.GEVI(i));
end
ax.XAxis.Label.String =  ppmean + " mCherry-CAAX";
ax.YAxis.Label.String =  ppmean + " cyOFP";
% linear fit
[~, ~, ~, l, txt] = libplot.lnrFitting(x, y, true, true, ax,'polyfit');
l.Color = [0.3010 0.7450 0.9330];
l.DisplayName = 'polyfit with outlier';
txt.Color = [0 0.4470 0.7410];

legend(ax, [s,l], 'Location', 'eastoutside');

%% Figure 3: compare photostability
ppmean = "Photostability Mean";
ppstd = "Photostability STD";
T_plot_mCherry = sortrows(T_plot_mCherry, ppmean, 'ascend');
T_plot_cyOFP = T_plot_cyOFP(T_plot_mCherry.GEVI,:);
x = T_plot_mCherry.(ppmean);
y = T_plot_cyOFP.(ppmean);
xstd = T_plot_mCherry.(ppstd);
ystd = T_plot_cyOFP.(ppstd);

f = figure('Color', [1, 1, 1]);
f.UserData.PlotData1 = T_plot_cyOFP;   % bind data
f.UserData.PlotData2 = T_plot_mCherry;
db.watermark(f, j);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
b = bar(ax,[x, y]);
b(1).FaceColor = [0.8500 0.3250 0.0980];
b(1).DisplayName = "mCherry-CAAX";
b(2).FaceColor = [0.9290 0.6940 0.1250];
b(2).DisplayName = "cyOFP";
xpos = b(1).XEndPoints;
e(1) = errorbar(ax, xpos, x, xstd, ...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
ypos = b(2).XEndPoints;
e(2) = errorbar(ax, ypos, y, ystd,...
        'LineStyle', 'none', ...
        'LineWidth', 1, ...
        'CapSize', 6, ...
        'Color', 'black', ...
        'DisplayName', '\mu\pm\sigma');
legend(ax, b, 'Location', 'northwest');
ax.YAxis.Label.String = ppmean;
ax.Box = false;
ax.XAxis.TickLabels = T_plot_mCherry.GEVI;
ax.XAxis.TickLabelRotation = 90;
ax.XAxis.Limits = [0.4, height(T_plot_mCherry) + 0.6];

%% Figure 4: scatter compare photostability
f = figure('Color', [1, 1, 1]);
f.UserData.PlotData1 = T_plot_cyOFP;   % bind data
f.UserData.PlotData2 = T_plot_mCherry;
db.watermark(f, j);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
errorbar(x, y, ystd, ystd, xstd, xstd,...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 6, ...
    'Color', 'black', ...
    'DisplayName', '\mu\pm\sigma');
cmap = lines;
for i = 1:numel(x)
    s(i) = scatter(ax, x(i),y(i),[], cmap(i,:),'filled',...
        'DisplayName', T_plot_mCherry.GEVI(i));
end
ax.XAxis.Label.String =  ppmean + " mCherry-CAAX";
ax.YAxis.Label.String =  ppmean + " cyOFP";
% linear fit
[~, ~, ~, l(1), txt(1)] = libplot.lnrFitting(x, y, true, true, ax,'polyfit');
l(1).Color = [0.3010 0.7450 0.9330];
l(1).DisplayName = 'linear fitting';
txt(1).Color = [0 0.4470 0.7410];
legend(ax, [s,l], 'Location', 'eastoutside');

%% Figure 5: Scatter_rand2dwell
pp1 = "-dF/F0 Short Mean";
pp2 = "Photostability Mean";
cmap = hsv(9);
f = figure('Color', [1, 1, 1]);
f.UserData.PlotData1 = T_plot_cyOFP_well;   % bind data
f.UserData.PlotData2 = T_plot_mCherry_well;
db.watermark(f, j);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
so = gobjects(height(T_plot_mCherry), 1);
sr = gobjects(height(T_plot_mCherry), 1);
for i = 1:height(T_plot_mCherry)
    geviname = T_plot_mCherry.GEVI(i);
    idxo = strcmp(T_plot_cyOFP_well.GEVI, geviname);
    xo = T_plot_cyOFP_well.(pp1)(idxo);
    yo = T_plot_cyOFP_well.(pp2)(idxo);
    idxr = strcmp(T_plot_mCherry_well.GEVI, geviname);
    xr = T_plot_mCherry_well.(pp1)(idxr);
    yr = T_plot_mCherry_well.(pp2)(idxr);
    so(i) = scatter(ax, xo, yo, [], cmap(i,:), 'LineWidth', 2, ...
        'DisplayName', T_plot_mCherry.GEVI(i) + "-cyOFP");
    sr(i) = scatter(ax, xr, yr, [], cmap(i,:), 'filled',...
        'DisplayName', T_plot_mCherry.GEVI(i) + "-mCherry-CAAX");
end
legend(ax, [so, sr], 'Location', 'eastoutside');
ax.XAxis.Label.String =  pp1;
ax.YAxis.Label.String =  pp2;    

% Plot countour
operation = "X.*(Y)";
xvals = T_plot_well.(pp1);
yvals = T_plot_well.(pp2);
[~, ~, h] = libplot.plotcontour(xvals,yvals,'ax',ax,'operation',operation,...
    'note',operation);     

%% Functions
function [gevi, refFP] = parsename(name)
    if contains(name,"-cyOFP") > contains(name,"-mCherry-CAAX")
        gevi = extractBefore(name, "-cyOFP");
    elseif contains(name,"-cyOFP") < contains(name,"-mCherry-CAAX")
        gevi = extractBefore(name, "-mCherry");
    else
        gevi = name;
    end
end

function refFP = parseref(name)
    if contains(name,"-cyOFP") > contains(name,"-mCherry-CAAX")
        refFP = "cyOFP";
    elseif contains(name,"-cyOFP") < contains(name,"-mCherry-CAAX")
        refFP = "mCherry-CAAX";
    else
        refFP = "";
    end
end

function [T_filtered, T_discard] = discardoutliers(T, grouping, p)
    arguments 
        T (:,:) table
        grouping (1,:) double
        p.method (1,:) string = "gesd"
        p.DataVariables (1,:) cell = T.Properties.VariableNames
    end
    mtd = p.method;
    dv = p.DataVariables;
    T_filtered = table();
    T_discard = table();
    for i = 1:length(grouping)
        T_temp = T(find(grouping == i),:);
        [T_temp_filtered, TF] = rmoutliers(T_temp, mtd,  'threshold',1.5, 'DataVariables', dv);
        T_temp_discard = T_temp(TF,:);
        T_filtered = [T_filtered; T_temp_filtered];
        T_discard = [T_discard; T_temp_discard];
    end
end 