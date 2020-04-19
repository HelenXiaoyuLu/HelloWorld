%% 2pFstim trace analysis and visual
% 
% Yueyang Gou, St-Pierre Lab, Apr. 2020
% 
% Load MATLAB that contains the following variables: 
% 
% T_fovData, T_fovInfo, stim，paths and conditions used for analysis in case
% need to be rerun



%% Fine Process of FOV traces
% Initialize (load and integrate all data into T_fov)

n_fov = height(T_fovData);      % number of FOV
sevents = stim.segment(0.1);
seventRange = sevents.range();
Coord = cat(1, T_fovInfo.Coord{:});

T_fov = T_fovData;              % main FOV dataStore table
T_fov.ConstructName = T_fovInfo.Name;
T_fov.wellName = T_fovInfo.wellName;
T_fov.CoordWell = Coord(:, 2:3);
T_fov.CoordFOV = Coord(:, 4);

% load signal
traceRaw = sig.arbitrary.empty();
for i = 1:n_fov
    s = sig.arbitrary();
    s.setSignal(T_fov{i, 'trace_t'}{1}, T_fovData{i, 'trace'}{1});
    traceRaw(i, 1) = s;
end

T_fov.TraceRaw = traceRaw;



 %% Align with long first and then use short stim
%long range [0,0.5] short range[-0.1,0.1]

traceraw = sig.arbitrary.empty();
tracenormlized= sig.arbitrary.empty();

longstimavg = stim.select(21);

shortsevent_t = seventRange.t_on(1:20) - 0.1;    % get the alignment timepoint
shortstimavg = stim.select(1);

Dshort=zeros(n_fov,1);
Dlong=zeros(n_fov,1);

for i = 1:n_fov
    i
    s_raw = T_fov.TraceRaw(i).copy();
    s=s_raw.copy();
    
    % detrend
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, 'method', 'div', 'order', 3);
    s.setSignal(s.x, Fd);
    
    
    s.crop(0.01);
    
    % align using long stimulation 

        % crop the signal
        s_long=s.crop(7.5);
        % time alignment
        Dlong(i) = img.temporalAlignment(s_long.y, s_long.x, 'method', 'template', ...
            'template', longstimavg, 'figshow', false, 'delayRange', [-0.5, 0.5]);
        %Align both raw and the detrened traces and save the lag
        s.lag(Dlong(i)); 
        s_raw.lag(Dlong(i));
       
    
    % align using short stimulation

        % time averaging
        [Fda, T_sample] = img.eventAveraging(s.y, shortsevent_t, 1, s.x, 'dt', 0.001); % Average different events
        % time alignment
        Dshort(i) = img.temporalAlignment(Fda, T_sample, 'method', 'template', ...
            'template', shortstimavg, 'figshow', false, 'delayRange', [-0.05,0.15]);
        %Align both raw and the detrened traces and save the lag
        
        s.lag(Dshort(i));
        s_raw.lag(Dshort(i));


        
    % normalize
    Ffinal = math.normbase(s.y, 1:20, 1); % Normalize the signal based on average of 1:100 frames
    s.setSignal(s.x, Ffinal);

    % save raw traces and normlized detrended signal
    traceraw(i, 1) = s_raw;
    tracenormlized(i, 1) = s;

      
end

 T_fov.TraceNormalized = tracenormlized;
 T_fov.TraceRaw = traceraw; 
 T_fov.delayshort=Dshort;
 T_fov.delaylong=Dlong;
 
 
%% Align using short stim with fast photobleaching 2p

traceraw = sig.arbitrary.empty();
tracenormlized= sig.arbitrary.empty();
sevent_t = seventRange.t_on(1:20) - 0.1;    % get the alignment timepoint
stimavg = stim.select(1);
for i = 1:n_fov
    i

    s_raw = T_fov.TraceRaw(i).copy;
    
    s=s_raw.copy();
    
    % detrend
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, 'method', 'div', 'order', 3);
    s.setSignal(s.x, Fd);
    % crop
    s.crop(0.01);
    % time averaging
    [Fda, T_sample] = img.eventAveraging(s.y, sevent_t, 1, s.x, 'dt', 0.002); % Average different events
    % time alignment
    D = img.temporalAlignment(Fda, T_sample, 'method', 'template', ...
        'template', stimavg, 'figshow', false, 'delayRange', [-0.5, 0.5]);
    s.lag(D);
    s_raw.lag(D);
    s.lag(0.3); % offset by one period due to general delay
    s_raw.lag(0.3);% offset by one period due to general delay
    
    % normalize
    Ffinal = math.normbase(s.y, 1:100, 1);
    s.setSignal(s.x, Ffinal);
    traceraw(i, 1) = s_raw;
    tracenormlized(i, 1) = s;
end

 T_fov.TraceNormalized = tracenormlized;
 T_fov.TraceRaw = traceraw; 


%% Align majorly using short but can also use long stimulation for sub
%long range [0,0.5] short range[-0.1,0.1]

traceraw = sig.arbitrary.empty();
tracenormlized= sig.arbitrary.empty();

longstimavg = stim.select(21);

shortsevent_t = seventRange.t_on(1:20) - 0.1;    % get the alignment timepoint
shortstimavg = stim.select(1);

Dshort=zeros(n_fov,1);
Dlong=zeros(n_fov,1);

for i = 1:n_fov
    i
    s_raw = T_fov.TraceRaw(i).copy();
    s=s_raw.copy();
    
    % detrend
    Fd = img.detrend(s.y, 1, s.x, 'figshow', false, 'method', 'div', 'order', 3);
    s.setSignal(s.x, Fd);
    
    
    s.crop(0.01);
    

    
    % align using short stimulation

        % time averaging
        [Fda, T_sample] = img.eventAveraging(s.y, shortsevent_t, 1, s.x, 'dt', 0.002); % Average different events
        % time alignment
        Dshort(i) = img.temporalAlignment(Fda, T_sample, 'method', 'template', ...
            'template', shortstimavg, 'figshow', false, 'delayRange', [-0.5,0.5]);
        %Align both raw and the detrened traces and save the lag
        
        s.lag(Dshort(i));
        s_raw.lag(Dshort(i));
        s.lag(0.3);
        s_raw.lag(0.3);

    % align using long stimulation 
        if Dshort(i)<=0.002
        % crop the signal
        s_long=s.crop(7.5);
        % time alignment
        Dlong(i) = img.temporalAlignment(s_long.y, s_long.x, 'method', 'template', ...
            'template', longstimavg, 'figshow', false, 'delayRange', [-0.15, 0.15]);
        %Align both raw and the detrened traces and save the lag
        s.lag(Dlong(i)); 
        s_raw.lag(Dlong(i));
        end
        
    % normalize
    Ffinal = math.normbase(s.y, 1:20, 1); % Normalize the signal based on average of 1:100 frames
    s.setSignal(s.x, Ffinal);

    % save raw traces and normlized detrended signal
    traceraw(i, 1) = s_raw;
    tracenormlized(i, 1) = s;

      
end

 T_fov.TraceNormalized = tracenormlized;
 T_fov.TraceRaw = traceraw; 
 T_fov.delayshort=Dshort;
 T_fov.delaylong=Dlong;
 

%% Crop averaged small stimulation AveSmallStimFinal and long stim LongStimFinal 

aveEvent=3:19;

traces = sig.arbitrary.empty();
sevent_t = seventRange.t_on(aveEvent) - 0.1;
for i = 1:n_fov
    traceIn = T_fov.TraceNormalized(i);
    [y, x] = img.eventAveraging(traceIn.y, sevent_t, 1, traceIn.x, 'dt', 0.001);
    s = sig.arbitrary();
    s.setSignal(x, y);
    traces(i, 1) = s;
end
T_fov.AveSmallStimFinal = traces;


T_fov.LongStimFinal = T_fov.TraceNormalized.crop(7.5, 9);


%% Resample all FOV Traces to use one common time point (Between differnet FOVs for easier alignment)
T_fov.TraceRaw = T_fov.TraceRaw.synchronize('dt', 0.002, 'method', 'intersect');
T_fov.TraceNormalized = T_fov.TraceNormalized.synchronize('dt', 0.002, 'method', 'intersect'); 
T_fov.AveSmallStimFinal = T_fov.AveSmallStimFinal.synchronize('dt', 0.002, 'method', 'intersect'); 
T_fov.LongStimFinal = T_fov.LongStimFinal.synchronize('dt', 0.002, 'method', 'intersect'); 


%% Pull From FOV to wells
% Well summary table

T_well = table();
% group FOVs by well
FOVgroups = findgroups(T_fov.CoordWell(:, 1), T_fov.CoordWell(:, 2));

% name of the well as construct
T_well.wellName = splitapply(@unique, T_fov.wellName, FOVgroups);
% name of the well as construct
T_well.ConstructName = splitapply(@unique, T_fov.ConstructName, FOVgroups);


% number of fov per well
T_well.n_fov = splitapply(@numel, T_fov.CoordFOV, FOVgroups);
% number of total pixels per well
T_well.Area = splitapply(@sum, T_fov.trace_n, FOVgroups);


% mean FOV TraceRaw per well
T_well.TraceRawMean = splitapply(@mean, T_fov.TraceRaw, T_fov.trace_n, FOVgroups);
% std of FOV TraceRaw per well
T_well.TraceRawStdfov = splitapply(@std, T_fov.TraceRaw, T_fov.trace_n, FOVgroups);

% mean FOV TraceNormalized per well
T_well.TraceNormalizedMean = splitapply(@mean, T_fov.TraceNormalized, T_fov.trace_n, FOVgroups);
% std of FOV TraceNormalized per well
T_well.TraceNormalizedStdfov = splitapply(@std, T_fov.TraceNormalized, T_fov.trace_n, FOVgroups);

% mean FOV small stimulation trace per well
T_well.AveSmallStimFinalMean = splitapply(@mean, T_fov.AveSmallStimFinal, T_fov.trace_n, FOVgroups);
% std of FOV small stimulation trace per well
T_well.AveSmallStimFinalStdfov = splitapply(@std, T_fov.AveSmallStimFinal, T_fov.trace_n, FOVgroups);

% mean FOV long stimulation trace per well
T_well.LongStimFinalMean = splitapply(@mean, T_fov.LongStimFinal, T_fov.trace_n, FOVgroups);
% std of FOV long stimulation trace per well
T_well.LongStimFinalStdfov = splitapply(@std, T_fov.LongStimFinal, T_fov.trace_n, FOVgroups);


%% brightness per FOV and per Well

% Calculate per FOV
GMean = zeros(n_fov, 1);
RMean = zeros(n_fov, 1);
GMean_final=zeros(n_fov, 1);
AveNormTrace=zeros(n_fov, 1);



for i = 1:n_fov
    GMap = T_fov{i, 'GMap'}{1};
    RMap = T_fov{i, 'RMap'}{1};
    finalMask = T_fov{i, 'finalMask'}{1};
    GMean(i) = mean(GMap(finalMask));
    RMean(i) = mean(RMap(finalMask));
    Gf=T_fov.trace{i,1};
    GMean_final(i) = mean(Gf(1,end-15,end));
    AveNormTrace(i)=mean(Gf,'all')/RMean(i);

end


T_fov.GMean = GMean;
T_fov.GMean_final = GMean_final;
T_fov.RMean = RMean;
T_fov.GRRatio = T_fov.GMean./T_fov.RMean;
T_fov.GbleachingRatio =1- T_fov.GMean_final./T_fov.GMean;
T_fov.AveNormTrace = AveNormTrace;

% Calculating for well data

weightmean=@(x1,x2) sum(x1.*x2)/sum(x2);
FOVgroups = findgroups(T_fov.CoordWell(:, 1), T_fov.CoordWell(:, 2));
fovtrace_n=T_fov.trace_n;



TgtBrightnessMean=zeros(height(T_well),1);
RefBrightnessMean=zeros(height(T_well),1);
GRratioMean=zeros(height(T_well),1);
TgtBrightnessFinalMean=zeros(height(T_well),1);
TgtBleachingRatioMean=zeros(height(T_well),1);
GRAveNormTrace=zeros(height(T_well),1);



for i=1:max(FOVgroups)
    TgtBrightnessMean(i)=weightmean(GMean(FOVgroups==i),fovtrace_n(FOVgroups==i));
    RefBrightnessMean(i)=weightmean(RMean(FOVgroups==i),fovtrace_n(FOVgroups==i));
    TgtBrightnessFinalMean(i)=weightmean(GMean_final(FOVgroups==i),fovtrace_n(FOVgroups==i));
    GRAveNormTrace(i)=weightmean(AveNormTrace(FOVgroups==i),fovtrace_n(FOVgroups==i));
end

% mean FOV TgtBrightness per well
T_well.TgtBrightnessMean = TgtBrightnessMean;
% mean FOV RefBrightness per well
T_well.RefBrightnessMean = RefBrightnessMean;
% mean FOV G/Rratio per well
T_well.GRratioMean = TgtBrightnessMean./RefBrightnessMean;
% mean FOV TgtFinalBrightness per well
T_well.TgtBrightnessFinalMean = TgtBrightnessFinalMean;
% mean FOV TgtBleachingRatio per well
T_well.TgtBleachingRatioMean = 1-TgtBrightnessFinalMean./TgtBrightnessMean;
% mean FOV AveTraceGRratio per well
T_well.AveTraceGRratioMean = GRAveNormTrace;


%% Well small stim trace plot -- FOVstd stat
% 


f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'Well short stim with FOV stat');
t.Padding = 'none';
for i = [13:36,61:84] 
    ax = nexttile(t);
    xf = [T_well.AveSmallStimFinalMean(i).x; flipud(T_well.AveSmallStimFinalMean(i).x)];
    yf = [T_well.AveSmallStimFinalMean(i).y + T_well.AveSmallStimFinalStdfov(i).y; 
          flipud(T_well.AveSmallStimFinalMean(i).y - T_well.AveSmallStimFinalStdfov(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well.AveSmallStimFinalMean(i).x, T_well.AveSmallStimFinalMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well.ConstructName(i), T_well.n_fov(i)));
end


%% Well Whole trace plot -- FOVstd stat


f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'Well whole trace with FOV stat');
t.Padding = 'none';
for i = [13:36,61:84] 
    ax = nexttile(t);
    xf = [T_well.TraceNormalizedMean(i).x; flipud(T_well.TraceNormalizedMean(i).x)];
    yf = [T_well.TraceNormalizedMean(i).y + T_well.TraceNormalizedStdfov(i).y; 
          flipud(T_well.TraceNormalizedMean(i).y - T_well.TraceNormalizedStdfov(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well.TraceNormalizedMean(i).x, T_well.TraceNormalizedMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well.ConstructName(i), T_well.n_fov(i)));
end


%% Well long trace plot -- FOVstd stat


f = figure();
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow');
title(t, 'Well long stim with FOV stat');
t.Padding = 'none';
for i = [13:36,61:84]
    ax = nexttile(t);
    xf = [T_well.LongStimFinalMean(i).x; flipud(T_well.LongStimFinalMean(i).x)];
    yf = [T_well.LongStimFinalMean(i).y + T_well.LongStimFinalStdfov(i).y; 
          flipud(T_well.LongStimFinalMean(i).y - T_well.LongStimFinalStdfov(i).y)] - 1;
    fill(ax, xf, yf, 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, T_well.LongStimFinalMean(i).x, T_well.LongStimFinalMean(i).y - 1, 'DisplayName', '\mu', 'LineWidth', 0.5);
    xlabel(ax, 'Time [s]');
    ylabel(ax, 'dF/F0');
    ylim(ax, [-0.4,0.4]);
    xlim(ax, [-Inf, Inf]);
    title(ax, sprintf('%s (n = %d)', T_well.ConstructName(i), T_well.n_fov(i)));
end


%% Plot response amplitude of short stim (averaged)

shortstim = stim.select(3);
shortstim = shortstim.segment(0.1);
shortstimeventRange = shortstim.range();
f=figure;
respEvalshort = cell(height(T_well), 1);
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = [13:36,61:84]
    respEvalshort{i} = sig.evaluator.response.load(T_well.AveSmallStimFinalMean(i).y, ...
        T_well.AveSmallStimFinalMean(i).x, ...
        'event', shortstimeventRange, 'extrema', true); % There used to be 'neighborhood',0.2
    ax = nexttile(t);
    respEvalshort{i}.plot(ax);
    title(T_well.ConstructName(i));
end
% please check the evaluation quality before continue


%% Quantify response amplitude (averaged)

shortstim = stim.select(3);
shortstim = shortstim.segment(0.1);
shortstimeventRange = shortstim.range();
respEvalshort = cell(height(T_well), 1);
for i = 1:height(T_well)
    respEvalshort{i} = sig.evaluator.response.load(T_well.AveSmallStimFinalMean(i).y, ...
        T_well.AveSmallStimFinalMean(i).x, ...
        'event', shortstimeventRange, 'extrema', true); % There used to be 'neighborhood',0.2
end
    % append response amplitude to table
    dFF0shortstim = zeros(height(T_well), 1);
    for j = 1:height(T_well)
        dFF0shortstim(j) = respEvalshort{j}.event.dFF0(1);
    end
    T_well.dFF0_shortStim = dFF0shortstim;


%% Plot response amplitude of long stim

longstim = stim.select(21);
longstim = longstim.segment(0.1);
longstimeventRange = longstim.range();
f=figure;
respEvallong = cell(height(T_well), 1);
f.Color = [1, 1, 1];
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = [13:36,61:84]
    respEvallong{i} = sig.evaluator.response.load(T_well.LongStimFinalMean(i).y, ...
        T_well.LongStimFinalMean(i).x, ...
        'event', longstimeventRange, 'extrema', true); % There used to be 'neighborhood',0.2
    ax = nexttile(t);
    respEvallong{i}.plot(ax);
    title(T_well.ConstructName(i));
end
% please check the evaluation quality before continue


%% Quantify response amplitude of long stim

longstim = stim.select(21);
longstim = longstim.segment(0.1);
longstimeventRange = longstim.range();
respEvallong = cell(height(T_well), 1);
for i = 1:height(T_well)
    respEvallong{i} = sig.evaluator.response.load(T_well.LongStimFinalMean(i).y, ...
        T_well.LongStimFinalMean(i).x, ...
        'event', longstimeventRange, 'extrema', true); % There used to be 'neighborhood',0.2
end

    % append response amplitude to table
    dFF0longstim = zeros(height(T_well), 1);
    for j = 1:height(T_well)
        dFF0longstim(j) = respEvallong{j}.event.dFF0(1);
    end
    T_well.dFF0_longStim = dFF0longstim;
    
    
%% Quantify photobleaching property

Common_t=T_well.TraceRawMean.range('method','intersect');
T_well.Bleachingratio_total = 1-T_well.TraceRawMean.t2v(Common_t.t_off)./T_well.TraceRawMean.t2v(Common_t.t_on+0.1);
T_well.Bleachingratio_2_1s = 1-T_well.TraceRawMean.t2v(Common_t.t_on+2.2)./T_well.TraceRawMean.t2v(Common_t.t_on+0.1);
T_well.TgtAUC_sig_NormaltoRandT = T_well.TraceRawMean.integral(Common_t.t_on, Common_t.t_off)./T_well.RefBrightnessMean./(Common_t.t_off-Common_t.t_on);


%% Output data to mat
savepath2ndmat='/Users/ericgou/Desktop/resultforgvei/2pFstim_20200202_GEVI_screening_cyOFPset_benchmarking_20X_objective_2ndmat.mat';
save(savepath2ndmat,'T_well','T_fov','stim','folderpath','folderV','folderI','refchselect','stimname','savepath','savepath2ndmat');


%% Output scalar result to excel file
Dataout= T_well;
Dataout = removevars(Dataout,{'TraceRawMean','TraceRawStdfov','TraceNormalizedMean','TraceNormalizedStdfov',...
    'AveSmallStimFinalMean','AveSmallStimFinalStdfov','LongStimFinalMean','LongStimFinalStdfov',});
Dataout= Dataout(:,[1 2 4 7 10 9 11 12 5 6 8 13 14 15 3]);
excelpath='/Users/ericgou/Desktop/resultforgvei/2pFstim_200126_GEVI_screening_cyOFPset_benchmarking_P1.xlsx';
writetable(Dataout,excelpath);

