%% Plot three-group discrimination
% Load jobs
jobid = [20201018230608 , ... % 1p Field stimulation
        20201009191913, ...   % 1p Photobleaching
        20201009191427];      % 2p Field stimulation
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201008 double blind\Pipeline output';
j1.output = load(fullfile(dirp, '1pFstim 20201008 20201018230823 (P2)')); % 1p Field stimulation
j2.output = load(fullfile(dirp, '1pPS 20201008 20201009191913 (P2)'));    % 1p Photobleaching
j3.output = load(fullfile(dirp, '2pFstim 20201009 20201009191427 (P2)')); % 2p Field stimulation
j1.input = '\\stpierrelab7910\E\Image\Xiaoyu\20201008 Double blind 3 group\Plate 2 double blind 3 group';

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
cmap = lines; 
for i = 1:max(grouping)
    idx = find(grouping == i);
    s(i) = scatter(T_plot_well.(pp1)(idx), ...
        T_plot_well.(pp2)(idx),'filled',...
        'DisplayName', T_plot_well.Name(idx(1)));
    s(i).CData = repmat(cmap(i,:),numel(idx),1)
end
ax.XAxis.Label.String = pp1;
ax.YAxis.Label.String = pp2;
legend(ax, s, 'location','northwest')    