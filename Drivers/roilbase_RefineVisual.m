% celllevel_RefineVisual
% Goal: generate a trace from each FOV
%% Load Data
if exist('evalout')
    fields = fieldnames(evalout);
    for i = 1:numel(fields)
        assignin('caller', fields{i}, evalout.(fields{i}));
    end
    clear('evalout');
end

%% Fine process of ROI traces
% initiate
n_roi = height(T_roiData);
sevents = stim.segment(0.1);
seventRange = sevents.range();
Coord = cat(1, T_roiInfo.Coord{:});
CoordFOV = cell2mat(T_fovInfo.Coord);
T_roi = T_roiData;
T_roi.ConstructName = T_roiInfo.Name;
T_roi.CoordWell = Coord(:, 2:3);
T_roi.CoordFOV = Coord(:, 4);
% load signal
traceRaw = sig.arbitrary.empty();
for i = 1:n_roi
    s = sig.arbitrary();
    s.setSignal(T_roi{i, 'trace_t'}{1}, T_roiData{i, 'trace'}{1});
    traceRaw(i, 1) = s;
end
T_roi.trace_t = [];
T_roi.trace = [];
T_roi.TraceRaw = traceRaw;

%% Aligmment and Detrend (You need to only run the section according to protocol)
% align using template 1p
traces = sig.arbitrary.empty();
delay = zeros(n_roi, 1);
% pb = ui.progressBar(n_roi);
parfor i = 1:n_roi
    s = T_roi.TraceRaw(i).copy;
    % detrend for alignment
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, 'method', 'normdiv');
    % time alignment
    delay(i) = img.temporalAlignment(Fd, s.x, 'method', 'template', ...
        'template', stim,'delayRange',[0,0.98], 'figshow', false);
    s.lag(delay(i));
    % detrend again to remove response region
    idleIdx = seventRange.t2v(s.x) == 0;
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, ...
        'method', 'div', 'idleIdx', idleIdx);
    % normalize
    Ffinal = math.normbase(Fd, 1:100, 1);
%     if delay(i) == 0
%         Ffinal = nan(size(s.x));
%     end
    s.setSignal(s.x, Ffinal);
    traces(i, 1) = s;
  %  pb.increment();         % increment progress
end
T_roi.Delay = delay;
T_roi.TraceDetrend = traces;

%% Optional Trace Conditioning
% No additional step:
T_roi.TraceFinal = T_roi.TraceDetrend;

%% Resample all FOV Traces to use one common time point
% check if any roi is off too much
for i = 1:n_roi
    a = T_roi.TraceFinal(i);
    if a.x(1)>1
        disp(i)
    end
end
T_roi.TraceFinal(1:end) = T_roi.TraceFinal(1:end).synchronize('dt', 0.001, 'method', 'union');
T_roi.TraceFinal = T_roi.TraceFinal.crop(0.86, 8.0);

%% Pull from ROI to Construct and Wells
% %%%%%%%%% well summary table from roi %%%%%%%%%%%%%%%%%
T_well = table();
% group roi by well
ROIgroups = findgroups(T_roi.CoordWell(:, 1), T_roi.CoordWell(:, 2));
CoordFOV = cell2mat(T_fovInfo.Coord);
FOVgroups = findgroups(CoordFOV(:,2),CoordFOV(:,3));
% name of the well as construct
T_well.Name = splitapply(@unique, T_fovInfo.Name, FOVgroups);
% number of roi per well
T_well.n_roi = splitapply(@numel, T_roi.CoordFOV, ROIgroups);
% number of total pixels per well
T_well.Area = splitapply(@sum, T_roi.area, ROIgroups);
% mean roi trace per well
T_well.TraceMean = splitapply(@mean, T_roi.TraceFinal, T_roi.area, ROIgroups);
% std of roi traces per well
T_well.TraceStd = splitapply(@std, T_roi.TraceFinal, T_roi.area, ROIgroups);

% %%% construct summary table from well %%%%%%%
T_groupFromWell = table();
% group wells by construct
Wellgroups = findgroups(T_well.Name);
% name of the construct
T_groupFromWell.Name = splitapply(@unique, T_well.Name, Wellgroups);
% number of well per construct
T_groupFromWell.n_well = splitapply(@numel, T_well.Name, Wellgroups);
% number of total pixels per construct
T_groupFromWell.Area = splitapply(@sum, T_well.Area, Wellgroups);
% mean well trace per construct
T_groupFromWell.TraceMean = splitapply(@mean, T_well.TraceMean, T_well.Area, Wellgroups);
% std of well traces per construct
T_groupFromWell.TraceStd = splitapply(@std, T_well.TraceMean, T_well.Area, Wellgroups);

% %%% construct summary table from roi %%%
for i = 1:n_roi
    T_roi.ConstructName(i) = T_well.Name(Coord(i,2));
end
T_groupFromROI = table();
% number of ROI per construct
% T_groupFromROI.n_roi = splitapply(@sum, T_well.n_roi, Wellgroups);
% group roi by name of constructs
ROIgroups = findgroups(T_roi.ConstructName);
% name of the construct
T_groupFromROI.Name = splitapply(@unique,T_roi.ConstructName, ROIgroups);
% number of FOV per construct
T_groupFromROI.n_roi = splitapply(@numel, T_roi.ConstructName, ROIgroups);
% number of total pixels per construct
T_groupFromROI.Area = splitapply(@sum, T_roi.area, ROIgroups);
% mean FOV trace per construct
T_groupFromROI.TraceMean = splitapply(@mean, T_roi.TraceFinal, T_roi.area, ROIgroups);
% std of FOV traces per construct
T_groupFromROI.TraceStd = splitapply(@std, T_roi.TraceFinal, T_roi.area, ROIgroups);

% ROIs stats per well -- plot mu and std traces on multi-panel figures
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'ROI Stats in Well');
t.Padding = 'none';
for i = 1:height(T_well)
    ax = nexttile(t);
    xf = [T_well.TraceMean(i).x; flipud(T_well.TraceMean(i).x)];
    yf = [T_well.TraceMean(i).y + T_well.TraceStd(i).y; 
          flipud(T_well.TraceMean(i).y - T_well.TraceStd(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well.TraceMean(i).x, T_well.TraceMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well.Name(i), T_well.n_roi(i)));
end

% well stats per construct -- plot mu and std traces on multi-panel figures
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
t.Padding = 'none';
title(t, 'Well Stats in Construct');
for i = 1:height(T_groupFromWell)
    ax = nexttile(t);
    xf = [T_groupFromWell.TraceMean(i).x; flipud(T_groupFromWell.TraceMean(i).x)];
    yf = [T_groupFromWell.TraceMean(i).y + T_groupFromWell.TraceStd(i).y; 
          flipud(T_groupFromWell.TraceMean(i).y - T_groupFromWell.TraceStd(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_groupFromWell.TraceMean(i).x, T_groupFromWell.TraceMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_groupFromWell.Name(i), T_groupFromWell.n_well(i)));
end

% FOVs stats per construct -- plot mu and std traces on multi-panel figures
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'ROI Stats in Construct');
t.Padding = 'none';
for i = 1:height(T_groupFromROI)
    ax = nexttile(t);
    xf = [T_groupFromROI.TraceMean(i).x; flipud(T_groupFromROI.TraceMean(i).x)];
    yf = [T_groupFromROI.TraceMean(i).y + T_groupFromROI.TraceStd(i).y; 
          flipud(T_groupFromROI.TraceMean(i).y - T_groupFromROI.TraceStd(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_groupFromROI.TraceMean(i).x, T_groupFromROI.TraceMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4, 0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_groupFromROI.Name(i), T_groupFromROI.n_roi(i)));
end

%% Quantification
% response amplitude per construct
respEval = cell(height(T_groupFromROI), 1);
for i = 1:height(T_groupFromROI)
    respEval{i} = sig.evaluator.response.load(T_groupFromROI.TraceMean(i).y, ...
        T_groupFromROI.TraceMean(i).x, ...
        'event', seventRange, 'neighborhood', [ones(1,7)*0.2,ones(1,2)*0.4], 'extrema', true);
end
% plot
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:height(T_groupFromROI)
    ax = nexttile(t);
    respEval{i}.plot(ax);
    title(T_groupFromROI.Name(i));
end
% please check the visualized evaluation quality before continue
% append response amplitude to table
n_stim = seventRange.n_pulse;   % number of stimulation events
for i = 1:n_stim
    dFF0 = zeros(height(T_groupFromROI), 1);
    for j = 1:height(T_groupFromROI)
        dFF0(j) = respEval{j}.event.dFF0(2 * i - 1);
    end
    T_groupFromROI.("Stim_" + i) = dFF0;
end

% response amplitude per well
respEval = cell(height(T_well), 1);
pb = ui.progressBar(height(T_well));
for i = 1:height(T_well)
    respEval{i} = sig.evaluator.response.load(T_well.TraceMean(i).y, ...
        T_well.TraceMean(i).x, ...
        'event', seventRange, 'neighborhood', [ones(1,7)*0.2,ones(1,2)*0.4], 'extrema', true);
    pb.increment();         % increment progress
end
% append response amplitude to table
n_stim = seventRange.n_pulse;   % number of stimulation events
for i = 1:n_stim
    dFF0 = zeros(height(T_well), 1);
    for j = 1:height(T_well)
        dFF0(j) = respEval{j}.event.dFF0(2 * i - 1);
    end
    T_well.("Stim_" + i) = dFF0;
end

% response amplitude per roi
respEval = cell(height(T_roi), 1);
parfor i = 1:height(T_roi)
    respEval{i} = sig.evaluator.response.load(T_roi.TraceFinal(i).y, ...
        T_roi.TraceFinal(i).x, ...
        'event', seventRange, 'neighborhood', [ones(1,7)*0.2,ones(1,2)*0.4], 'extrema', true);
end
for i = 1:n_stim
    dFF0 = zeros(height(T_roi), 1);
    parfor j = 1:height(T_roi)
        dFF0(j) = respEval{j}.event.dFF0(2 * i - 1);
    end
    T_roi.("Stim_" + i) = dFF0;
end


%% Save
save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_plate1\single cell mask results\manualmask_dFF0_quant.mat', 'stim', 'T_roi', 'T_roiInfo','T_roiData', 'T_well', 'T_groupFromROI', 'T_groupFromWell');
