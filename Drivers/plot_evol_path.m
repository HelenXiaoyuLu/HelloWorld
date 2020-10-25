%% Plot evol path
% Load jobs
jobid = [20201020190530 , ... % 1p Field stimulation
        20201020190928, ...   % 1p Photobleaching
        20201020190736];      % 2p Field stimulation
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201014 Evol path';
j1.output = load(fullfile(dirp, '1pFstim 20201014 20201020190530')); % 1p Field stimulation
j2.output = load(fullfile(dirp, '1pPS 20201014 20201020190928'));    % 1p Photobleaching
j3.output = load(fullfile(dirp, '2pFstim 20201014 20201020190736')); % 2p Field stimulation
j1.input = '\\stpierrelab7910\E\Image\Xiaoyu\20201014 Evolution path';
j2.input = '\\stpierrelab7910\E\Image\Xiaoyu\20201014 Evolution path';

%% Plot 1P scatter
[~, T_plot_well] = Wellbase_1PFieldStimTrigger.Scatter_rand2dWell(j1);
[~, T_plot_well_PS] = Wellbase_1PPhotobleaching.Bar_AreaUnderCurveWell(j2);
T_plot_well = innerjoin(T_plot_well, T_plot_well_PS,...
        'LeftKeys', 'wellName',...
        'RightKeys', 'wellName', ...
        'RightVariables', {'Area Under the Curve normalized', ...
        'Area Under the Curve unnormalized'}); 
T_plot_well.("G/R x AUC") = T_plot_well.("Relative brightness (GFP/RFP) Mean") .* ...
    T_plot_well.("Area Under the Curve normalized");

%% Scatter plot 
pp1 = "-dFF0 Stim    1ms 60V"; % x value
pp2 = "G/R x AUC"; % y value
grouping = findgroups(T_plot_well.Name);
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
s = gobjects(max(grouping),1);
cmap = jet(max(grouping)); 
for i = 1:max(grouping)
    idx = find(grouping == i);
    s(i) = scatter(T_plot_well.(pp1)(idx), ...
        T_plot_well.(pp2)(idx),'filled',...
        'DisplayName', T_plot_well.Name(idx(1)));
    s(i).CData = repmat(cmap(i,:),numel(idx),1);
end
ax.XAxis.Label.String = pp1;
ax.YAxis.Label.String = pp2;
legend(ax, s, 'location','eastoutside')    

%% Plot 1P scatter group
[~, T_plot] = Wellbase_1PFieldStimTrigger.Scatter_rand2d(j1);
[~, T_plot_PS] = Wellbase_1PPhotobleaching.Bar_AreaUnderCurve(j2);
T_plot = innerjoin(T_plot, T_plot_PS,...
        'LeftKeys', 'Name',...
        'RightKeys', 'Name', ...
        'RightVariables', {'Area Under the Curve normalized', ...
        'Area Under the Curve unnormalized'}); 
T_plot.("G/R x AUC") = T_plot.("Relative brightness (GFP/RFP) Mean") .* ...
    T_plot.("Area Under the Curve normalized");
T_plot.Properties.RowNames = T_plot.Name;

%% Scatter plot 
sel = ["ASAP1", "ASAP2s", "ASAP1-N124I", "ASAP1-N124V", "ASAP1-R406K", ...
    "ASAP1-N124V-R406K", "ASAP2s-H152E", "ASAP2s-N391D", "ASAP2s-Q397H", ...
    "ASAP2s-T207H", "ASAP2s-H152E-T207H-N391D-Q397H", ...
    "JEDI-1P", "JEDI-2P"];
T_plot = T_plot(sel,:);
pp1 = "Detectability Index92.5ms 30V"; % x value
pp2 = "Area Under the Curve normalized"; % y value
f = figure('Color', [1, 1, 1]);
f.UserData.PlotData = T_plot;   % bind data
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
s = gobjects(height(T_plot),1);
cmap = jet(height(T_plot)); 
for i = 1:max(height(T_plot))
    s(i) = scatter(T_plot.(pp1)(i), ...
        T_plot.(pp2)(i),'filled',...
        'DisplayName', T_plot.Name(i),'UserData',i);
    s(i).CData = cmap(i,:);
end
ax.XAxis.Label.String = pp1;
ax.YAxis.Label.String = pp2;
legend(ax, s, 'location','eastoutside')    

% Plot countour
operation = "X.*sqrt(abs(Y))";
xvals = T_plot.(pp1);
yvals = T_plot.(pp2);
[~, ~, h] = libplot.plotcontour(xvals, yvals, 'ax', ax,...
     'operation', operation, 'note', operation);
 
h.UserData = [pp1, pp2, operation];

mData = uimenu('Parent', f, ...
    'Tag', 'Data', ...
    'Text', 'Settings');
uimenu('Parent', mData, ...
    'Text', 'Set Reference Well', ...
    'Tag', 'SetR', ...
    'MenuSelectedFcn', @setReference); 

%% Plot 2P scatter group
[~, T_plot] = Wellbase_2PFieldStim.Scatter_rand2d(j3);
T_plot.Properties.RowNames = T_plot.Name;

%% Scatter 2P group 
sel = ["ASAP1", "ASAP2s", "ASAP1-N124I", "ASAP1-N124V", "ASAP1-R406K", ...
    "ASAP1-N124V-R406K", "ASAP2s-H152E", "ASAP2s-N391D", "ASAP2s-Q397H", ...
    "ASAP2s-T207H", "ASAP2s-H152E-T207H-N391D-Q397H", ...
    "JEDI-1P", "JEDI-2P"];
T_plot = T_plot(sel,:);
pp1 = "Detectability Index Short"; % x value
pp2 = "Photostability Mean"; % y value
f = figure('Color', [1, 1, 1]);
f.UserData.PlotData = T_plot;   % bind data
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
s = gobjects(height(T_plot),1);
cmap = jet(height(T_plot)); 
for i = 1:max(height(T_plot))
    s(i) = scatter(T_plot.(pp1)(i), ...
        T_plot.(pp2)(i),'filled',...
        'DisplayName', T_plot.Name(i),'UserData',i);
    s(i).CData = cmap(i,:);
end
ax.XAxis.Label.String = pp1;
ax.YAxis.Label.String = pp2;
legend(ax, s, 'location','eastoutside')    

% Plot countour
operation = "X.*sqrt(abs(Y))";
xvals = T_plot.(pp1);
yvals = T_plot.(pp2);
[~, ~, h] = libplot.plotcontour(xvals, yvals, 'ax', ax,...
     'operation', operation, 'note', operation);
 
h.UserData = [pp1, pp2, operation];

mData = uimenu('Parent', f, ...
    'Tag', 'Data', ...
    'Text', 'Settings');
uimenu('Parent', mData, ...
    'Text', 'Set Reference Well', ...
    'Tag', 'SetR', ...
    'MenuSelectedFcn', @setReference); 


%%
function setReference(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    ax = f.CurrentAxes;
    T_plot = f.UserData.PlotData;
    list = T_plot.Name;
    list = cat(1, "(No Reference)", list);
    [idx, tf] = listdlg('ListString', list, ...
        'Name', 'Reference Well', ...
        'SelectionMode', 'single', ...
        'ffs', 0, 'ListSize', [250, 160], ...
        'PromptString', 'Please select a reference well: ');
    if tf == 0  % cancel
        return
    end
    s = findobj(f, '-class', 'matlab.graphics.chart.primitive.Scatter');
    h = findobj(f, '-class', 'matlab.graphics.chart.primitive.Contour');
    [~,sidx] = sort(cell2mat({s.UserData}),'ascend');
    s = s(sidx);
    pp1 = h(1).UserData(1);
    pp2 = h(1).UserData(2);
    if idx == 1     % original
        normValueX = 1;
        normValueY = 1;
    else
        normValueX = T_plot.(pp1)(idx - 1);
        normValueY = T_plot.(pp2)(idx - 1);       
        if normValueX ==0 || normValueY ==0
            warning('missing value for selected controls, please select another');
            normValueX = 1;
            normValueY = 1;
        else 
            ax.XAxis.Label.String = strcat('Normalized', pp1);
            ax.YAxis.Label.String = strcat('Normalized', pp2);
        end
    end
    xvals = zeros(1,height(T_plot));
    yvals = zeros(1,height(T_plot));
    for i = 1 : height(T_plot) 
        s(i).XData = T_plot.(pp1)(i)./ normValueX;
        s(i).YData = T_plot.(pp2)(i)./ normValueY;
        xvals(i) = s(i).XData;
        yvals(i) = s(i).YData;
    end
    operation = h(1).UserData(3);
    [~, ~, h1] = libplot.plotcontour(xvals, yvals, 'ax', ax, ...
        'operation', operation, 'note', operation);
    h1.UserData = [pp1, pp2, operation];
    delete(h); % delete the original contour obj
end