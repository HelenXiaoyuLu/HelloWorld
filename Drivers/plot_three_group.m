clear
clc
j.output = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\Compare Benchmakring\corrFilt_ver2\2PFieldStimCorrFilt 20201008 20201101005147 (P2)');
j.output.input = '\\stpierrelab7910\E\Image\Xiaoyu\20201008 Double blind 3 group\Plate 2 double blind 3 group';
% Wellbase_2PFieldStim.Scatter_rand2dwell(j)

%% 2pFstim 
[~, T_plot_well] = Wellbase_2PFieldStim.Scatter_rand2dwell(j)
pp1 = "-dF/F0 Long Mean"; % x value
pp2 = "Photostability Mean"; % y value
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

%% 1pFstim FOV
T_fov = j.output.T_fov;
xfit = zeros(height(T_fov),4);
rsq = zeros(height(T_fov),1);
for i = 1:height(T_fov)
    tr = T_fov.TraceRaw(i);
    t0 = find(tr.y == max(tr.y));
    xdata = tr.x(t0:end)-tr.x(t0);
    ydata = tr.y(t0:end);
    bsl = math.isbaseline(ydata,'fin');
    bsl(1:100) = 1;
    xdata = (xdata(bsl));
    ydata = (ydata(bsl));
    C = max(ydata);
    [xfit(i,[1,2,4]),rsq(i)] = lsqcurvefit(@(x,xdata) x(1)*exp(-xdata/x(2))+(C-x(1))*exp(-xdata/x(3)),[0.5*C 0.1 30],xdata,ydata);
    xfit(i,3) = C - xfit(i,1);
    rsq(i) = rsq(i)./numel(xdata);
end
T_fov.("Double Exp remove peaks") = xfit;
T_fov.("Fit RESNORM remove peaks") = rsq; %  sum {(FUN(X,XDATA)-YDATA).^2}

% xfit = zeros(height(T_fov),4);
% rsq = zeros(height(T_fov),1);
% for i = 1:height(T_fov)
%     tr = T_fov.TraceRaw(i);
%     t0 = find(tr.y == max(tr.y));
%     xdata = tr.x(t0:end)-tr.x(t0);
%     ydata = tr.y(t0:end);
%     C = max(ydata);
%     [xfit(i,[1,2,4]),rsq(i)] = lsqcurvefit(@(x,xdata) x(1)*exp(-xdata/x(2))+(C-x(1))*exp(-xdata/x(3)),[0.5*C 0.1 30],xdata,ydata);
%     xfit(i,3) = C - xfit(i,1);
%     rsq(i) = rsq(i)./numel(xdata);
% end
% T_fov.("Double Exp raw curve") = xfit;
% T_fov.("Fit RESNORM raw curve") = rsq; %  sum {(FUN(X,XDATA)-YDATA).^2}

% T_fov.("Turning point remove peaks") = log(T_fov.("Double Exp remove peaks")(:,1)./T_fov.("Double Exp remove peaks")(:,3))./(1./T_fov.("Double Exp remove peaks")(:,3)-1./T_fov.("Double Exp remove peaks")(:,4));
% T_fov.("Turning point raw curve") = log(T_fov.("Double Exp raw curve")(:,1)./T_fov.("Double Exp remove peaks")(:,3))./(1./T_fov.("Double Exp remove peaks")(:,3)-1./T_fov.("Double Exp remove peaks")(:,4));

%% Plot sample
% T_fov = sortrows(T_fov,"Fit RESNORM remove peaks",'ascend');
k = 1;
figure()
hold on
tr = T_fov.TraceRaw(k);
t0 = find(tr.y == max(tr.y));
xdata = tr.x(t0:end)-tr.x(t0);
ydata = tr.y(t0:end);
fitwbaseline = T_fov.("Double Exp remove peaks")(k,:);
C = max(ydata);
plot(xdata,ydata)
plot(xdata,fitwbaseline(1)*exp(-xdata/fitwbaseline(2))+fitwbaseline(3)*exp(-xdata/fitwbaseline(4)))
plot(xdata,fitwobaseline(1)*exp(-xdata/fitwobaseline(2))+fitwobaseline(3)*exp(-xdata/fitwobaseline(4)))
legend('raw trace','fit with baseline','fit with raw')
% plot(xdata,-fitwbaseline(1)/fitwbaseline(2)*exp(-xdata/fitwbaseline(2))-fitwbaseline(3)/fitwbaseline(4)*exp(-xdata/fitwbaseline(4)))
% plot(xdata,fitwbaseline(1)/fitwbaseline(2)^2*exp(-xdata/fitwbaseline(2))+fitwbaseline(3)/fitwbaseline(4)^2*exp(-xdata/fitwbaseline(4)))

%% Pull group level stats
T_group = j.output.T_group;
T_well = j.output.T_well;
weightedMean = @(x, w) sum(x.*w, 'all')./sum(w, 'all');
% Pull from FOV to Well
grouping = findgroups(T_fov.Well);
T_fov2Well = table();
T_fov2Well.Well = splitapply(@unique, T_fov.Well, grouping);
T_fov2Well.Area = splitapply(@sum, T_fov.Area, grouping);
% kinectics mean
T_fov2Well.("Double Exp Coeff fast") = splitapply(weightedMean, T_fov.("Double Exp remove peaks")(:,1), T_fov.Area, grouping);
T_fov2Well.("Double Exp Coeff slow") = splitapply(weightedMean, T_fov.("Double Exp remove peaks")(:,3), T_fov.Area, grouping);
T_fov2Well.("Double Exp Tau fast") = splitapply(weightedMean, T_fov.("Double Exp remove peaks")(:,2), T_fov.Area, grouping);
T_fov2Well.("Double Exp Tau slow") = splitapply(weightedMean, T_fov.("Double Exp remove peaks")(:,4), T_fov.Area, grouping);
% kinectics std
T_fov2Well.("Double Exp Coeff fast STD") = splitapply(@std, T_fov.("Double Exp remove peaks")(:,1), T_fov.Area, grouping);
T_fov2Well.("Double Exp Coeff slow STD") = splitapply(@std, T_fov.("Double Exp remove peaks")(:,3), T_fov.Area, grouping);
T_fov2Well.("Double Exp Tau fast STD") = splitapply(@std, T_fov.("Double Exp remove peaks")(:,2), T_fov.Area, grouping);
T_fov2Well.("Double Exp Tau slow STD") = splitapply(@std, T_fov.("Double Exp remove peaks")(:,4), T_fov.Area, grouping);
T_fov2Well.("RMean") = splitapply(weightedMean, T_fov.R, T_fov.Area, grouping);
T_fov2Well.("RStd") = splitapply(@std, T_fov.R, T_fov.Area, grouping);
% Merge
T_well = innerjoin(T_fov2Well, T_well, 'Keys', 'Well');
T_well = innerjoin(T_well, T_group, ...
    'LeftKeys', 'Group', ...
    'RightKeys', 'Group', ...
    'RightVariables', 'Name');


pp1 = "Double Exp Tau fast"; % x value
pp2 = "Double Exp Tau slow"; % y value
pp1std = "Double Exp Tau fast STD"; % x value
pp2std = "Double Exp Tau slow STD"; % y value

grouping = findgroups(T_well.Name);
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
cmap = lines; 

% errorbar(T_well.(pp1),T_well.(pp2),...
%     T_well.(pp2std),T_well.(pp2std),...
%     T_well.(pp1std),T_well.(pp1std),...
%     'LineStyle', 'none',...
%     'DisplayName', '\mu\pm\sigma',...
%     'Color', 'black');

s = gobjects(max(grouping),1);
for i = 1:max(grouping)
    idx = find(grouping == i);
    s(i) = scatter(T_well.(pp1)(idx), ...
        T_well.(pp2)(idx),'filled',...
        'DisplayName', T_well.Name(idx(1)));
    s(i).CData = repmat(cmap(i,:),numel(idx),1)
end
ax.XAxis.Label.String = pp1;
ax.YAxis.Label.String = pp2;
legend(ax, s, 'location','northwest')

%% Plot job output

pp1 = "Double Exp Tau fast"; % x value
pp2 = "Double Exp Tau slow"; % y value
pp1std = "Double Exp Tau fast STD"; % x value
pp2std = "Double Exp Tau slow STD"; % y value
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
cmap = lines;
grouping = findgroups(T_plot.Name);
errorbar(T_plot.(pp1),T_plot.(pp2),...
    T_plot.(pp2std),T_plot.(pp2std),...
    T_plot.(pp1std),T_plot.(pp1std),...
    'LineStyle', 'none',...
    'DisplayName', '\mu\pm\sigma',...
    'Color', 'black');
s = gobjects(max(grouping),1);
for i = 1:max(grouping)
    idx = find(grouping == i);
    s(i) = scatter(T_plot.(pp1)(idx), ...
        T_plot.(pp2)(idx),'filled',...
        'DisplayName', T_plot.Name(idx(1)));
    s(i).CData = repmat(cmap(i,:),numel(idx),1)
end