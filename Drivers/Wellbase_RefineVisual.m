% This is the example of sp. 
% 
% Zhuohe Liu, St-Pierre Lab, Apr. 2020
% harry.liu@rice.edu

%% Load Data
% matfile contains the following variables: T_fovData, T_fovInfo, stim
jh = db.job.pull('id','20200410112943');
jh.loadCache();
evalout = jh.result.evalout;
% get variable out if run in job
if exist('evalout', 'var')
    fields = fieldnames(evalout);
    for i = 1:numel(fields)
        assignin('caller', fields{i}, evalout.(fields{i}));
    end
    clear('evalout');
end

%% %%%%%%%%%% Fine Process of FOV traces %%%%%%%%%%%%%%%
% initialize
n_fov = height(T_fovData);      % number of FOV
sevents = stim.segment(0.1);
seventRange = sevents.range();
Coord = cat(1, T_fovInfo.Coord{:});
T_fov = T_fovData;              % main FOV dataStore table
T_fov.ConstructName = T_fovInfo.Name;
T_fov.CoordWell = Coord(:, 2:3);
T_fov.CoordFOV = Coord(:, 4);
% load signal
traceRaw = sig.arbitrary.empty();
for i = 1:n_fov
    s = sig.arbitrary();
    s.setSignal(T_fov{i, 'trace_t'}{1}, T_fovData{i, 'trace'}{1});
    traceRaw(i, 1) = s;
end
T_fov.trace_t = [];
T_fov.trace = [];
T_fov.TraceRaw = traceRaw;

%% Aligmment and Detrend (You need to only run the section according to protocol)
% align using template 1p
traces = sig.arbitrary.empty();
delay = zeros(n_fov, 1);
% pb = ui.progressBar(n_fov);
parfor i = 1:n_fov
    s = T_fov.TraceRaw(i).copy;
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
    s.setSignal(s.x, Ffinal);
    traces(i, 1) = s;
%     pb.increment();         % increment progress
end
T_fov.Delay = delay;
T_fov.TraceDetrend = traces;

% align using excitation 1p
traces = sig.arbitrary.empty();
delay = zeros(n_fov, 1);
pb = ui.progressBar(n_fov);
for i = 1:n_fov
    s = T_fov.TraceRaw(i).copy;
    % time alignment
    delay(i) = img.temporalAlignment(s.y, s.x, 'method', 'excitation', 'figshow', false);
    s.lag(delay(i));
    % crop initial no excitation region
    s.crop(0.5);
    % detrend (remove response region)
    idleIdx = seventRange.t2v(s.x) == 0;
    Fd = img.detrend(s.y, 1, s.x, 'order', 3, 'figshow', false, ...
        'method', 'div', 'idleIdx', idleIdx);
    % normalize
    Ffinal = math.normbase(Fd, 1:100, 1);
    s.setSignal(s.x, Ffinal);
    traces(i, 1) = s;
    pb.increment();         % increment progress
end
T_fov.Delay = delay;
T_fov.TraceDetrend = traces;

% align using template with fast photobleaching 2p
traces = sig.arbitrary.empty();
sevent_t = seventRange.t_on(1:20) - 0.1;    % average 20 short stimulation to get alignment
stimavg = stim.select(1);
delay = zeros(n_fov, 1);
pb = ui.progressBar(n_fov);
for i = 1:n_fov
    s = T_fov.TraceRaw(i).copy;
    % detrend
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, 'method', 'div', 'order', 3);
    s.setSignal(s.x, Fd);
    % crop
    s.crop(0.01);
    % time averaging
    [Fda, T_sample] = img.eventAveraging(s.y, sevent_t, 1, s.x, 'dt', 0.002);
    % time alignment
    delay(i) = img.temporalAlignment(Fda, T_sample, 'method', 'template', ...
        'template', stimavg, 'figshow', false, 'delayRange', [-0.5, 0.5]);
    delay(i) = delay(i) + 0.3;  % offset by one period due to general delay
    s.lag(delay(i));
    % normalize
    Ffinal = math.normbase(s.y, 1:100, 1);
    s.setSignal(s.x, Ffinal);
    traces(i, 1) = s;
    pb.increment();         % increment progress
end
T_fov.Delay = delay;
T_fov.TraceDetrend = traces;

%% Optional Trace Conditioning
% !!!! If no additional step, you need to run the following line
T_fov.TraceFinal = T_fov.TraceDetrend;

% Event averaging
traces = sig.arbitrary.empty();
sevent_t1 = seventRange.t_on(1:3) - 0.2;
sevent_t2 = seventRange.t_on(4:6) - 0.2;
for i = 1:6
    traceIn = T_fov.TraceDetrend(i);
    [y1, x1] = img.eventAveraging(traceIn.y, sevent_t1, 1, traceIn.x, 'dt', 0.001);
    [y2, x2] = img.eventAveraging(traceIn.y, sevent_t2, 1, traceIn.x, 'dt', 0.001);
    s = sig.arbitrary();
    s.setSignal([x1;x2], [y1;y2]);
    traces(i, 1) = s;
end
T_fov.TraceFinal = traces;

% Event Cropping
T_fov.TraceFinal = T_fov.TraceDetrend.crop(4.9, 5.2);

% Event Time Shifting
T_fov.TraceFinal = T_fov.TraceFinal.lag(-3.5);

%% Resample all FOV Traces to use one common time point
% check if any roi is off too much
for i = 1:n_fov
    a = T_fov.TraceFinal(i);
    if a.x(1)>1
        disp(i)
    end
end
T_fov.TraceFinal = T_fov.TraceFinal.synchronize('dt', 0.001, 'method', 'intersect');

%% Pull From FOV to Construct and Wells
% %%%%%%%%% well summary table %%%%%%%%%%%%%%%%%
T_well = table();
% group FOVs by well
FOVgroups = findgroups(T_fov.CoordWell(:, 1), T_fov.CoordWell(:, 2));
% name of the well as construct
T_well.Name = splitapply(@unique, T_fov.ConstructName, FOVgroups);
% number of fov per well
T_well.n_fov = splitapply(@numel, T_fov.CoordFOV, FOVgroups);
% number of total pixels per well
T_well.Area = splitapply(@sum, T_fov.trace_n, FOVgroups);
% mean FOV trace per well
T_well.TraceMean = splitapply(@mean, T_fov.TraceFinal, T_fov.trace_n, FOVgroups);
% std of FOV traces per well
T_well.TraceStd = splitapply(@std, T_fov.TraceFinal, T_fov.trace_n, FOVgroups);

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

% %%% Construct summary table from FOV %%%%%%%
T_groupFromFOV = table();
% group FOVs by construct
FOVgroups = findgroups(T_fov.ConstructName);
% name of the construct
T_groupFromFOV.Name = splitapply(@unique, T_fov.ConstructName, FOVgroups);
% number of FOV per construct
T_groupFromFOV.n_fov = splitapply(@numel, T_fov.ConstructName, FOVgroups);
% number of total pixels per construct
T_groupFromFOV.Area = splitapply(@sum, T_fov.trace_n, FOVgroups);
% mean FOV trace per construct
T_groupFromFOV.TraceMean = splitapply(@mean, T_fov.TraceFinal, T_fov.trace_n, FOVgroups);
% std of FOV traces per construct
T_groupFromFOV.TraceStd = splitapply(@std, T_fov.TraceFinal, T_fov.trace_n, FOVgroups);


% FOVs stats per well -- plot mu and std traces on multi-panel figures
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'FOV Stats in Well');
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
    title(ax, sprintf('%s (n = %d)', T_well.Name(i), T_well.n_fov(i)));
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
    yf = [T_groupFromWell.TraceMean(i).y + T_groupFromFOV.TraceStd(i).y; 
          flipud(T_groupFromWell.TraceMean(i).y - T_groupFromFOV.TraceStd(i).y)] - 1;
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
title(t, 'FOV Stats in Construct');
t.Padding = 'none';
for i = 1:height(T_groupFromFOV)
    ax = nexttile(t);
    xf = [T_groupFromFOV.TraceMean(i).x; flipud(T_groupFromFOV.TraceMean(i).x)];
    yf = [T_groupFromFOV.TraceMean(i).y + T_groupFromFOV.TraceStd(i).y; 
          flipud(T_groupFromFOV.TraceMean(i).y - T_groupFromFOV.TraceStd(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_groupFromFOV.TraceMean(i).x, T_groupFromFOV.TraceMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4, 0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_groupFromFOV.Name(i), T_groupFromFOV.n_fov(i)));
end

%% Quantification
% response amplitude per construct
respEval = cell(height(T_groupFromFOV), 1);
for i = 1:height(T_groupFromFOV)
    respEval{i} = sig.evaluator.response.load(T_groupFromFOV.TraceMean(i).y, ...
        T_groupFromFOV.TraceMean(i).x, ...
        'event', seventRange, 'neighborhood', 0.2, 'extrema', true);
end
% plot
f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:height(T_groupFromFOV)
    ax = nexttile(t);
    respEval{i}.plot(ax);
    title(T_groupFromFOV.Name(i));
end
% please check the visualized evaluation quality before continue
% append response amplitude to table
n_stim = seventRange.n_pulse;   % number of stimulation events
for i = 1:n_stim
    dFF0 = zeros(height(T_groupFromFOV), 1);
    for j = 1:height(T_groupFromFOV)
        dFF0(j) = respEval{j}.event.dFF0(2 * i - 1);
    end
    T_groupFromFOV.("Stim_" + i) = dFF0;
end

% response amplitude per well
respEval = cell(height(T_well), 1);
pb = ui.progressBar(height(T_well));
for i = 1:height(T_well)
    respEval{i} = sig.evaluator.response.load(T_well.TraceMean(i).y, ...
        T_well.TraceMean(i).x, ...
        'event', seventRange, 'neighborhood', 0.2, 'extrema', true);
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

% brightness per FOV
GMean = zeros(n_fov, 1);
RMean = zeros(n_fov, 1);
for i = 1:n_fov
    RMap = T_fov{i, 'RMap'}{1};
    finalMask = T_fov{i, 'finalMask'}{1};
    % integrate a portion to minimize noise
    % you need to change the onset time specifically for each protocol
    % for excitation onset protocol, delay is needed
    delay = T_fov.Delay(i);
    GMean(i) = T_fov.TraceRaw(i).integral(0 - delay, 0.01 - delay)./0.01;
    % GMean(i) = T_fov.TraceRaw(i).integral(0, 0.01)./0.01;
    RMean(i) = mean(RMap(finalMask));
end
T_fov.GMean = GMean;
T_fov.RMean = RMean;
T_fov.GRRatio = T_fov.GMean./T_fov.RMean;

% photobleaching
T_fov.BleachingdFF0 = T_fov.TraceRaw.t2v(9)./T_fov.TraceRaw.t2v(0) - 1;
T_fov.AOC = T_fov.TraceRaw.integral(0, 10);

% pull data from FOV to well
% group FOVs by well
FOVgroups = findgroups(T_fov.CoordWell(:, 1), T_fov.CoordWell(:, 2));
% mean 
T_well.GRRatioMean = splitapply(@mean, T_fov.GRRatio, FOVgroups);
T_well.AOCMean = splitapply(@mean, T_fov.AOC, FOVgroups);
% std
T_well.GRRatioStd = splitapply(@std, T_fov.GRRatio, FOVgroups);
T_well.AOCStd = splitapply(@std, T_fov.AOC, FOVgroups);
% SEM
T_well.GRRatioSEM = splitapply(@(x) std(x)./sqrt(numel(x)), T_fov.GRRatio, FOVgroups);
T_well.AOCSEM = splitapply(@(x) std(x)./sqrt(numel(x)), T_fov.AOC, FOVgroups);

%% save result
% FOV
T_fov;
T_fovInfo;
% Well and well-level stats
T_well;
% Group and group-level stats
T_groupFromFOV;
T_groupFromWell;

save('result.mat', 'stim', 'T_fov', 'T_fovInfo', 'T_well', 'T_groupFromFOV', 'T_groupFromWell');
