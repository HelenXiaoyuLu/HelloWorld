% This is the example of sp. 
% Optimized for red single channel screening, 200Hz
% Xiaoyu, St-Pierre Lab, May. 2020
% xiaoyu.lu@rice.edu

clear
clc
% clear import
import sp.*

%% Phase 1: Import Files and Settings
lib = sp.libData('1P FieldStim');
lib.protocol.path = 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20190310-0318 Jun cp\20190308 Jun cp 200Hz';
lib.loadFrame('path', '200Hz', ...
            'f',1:8, ...
            'ncpu', Inf);
lib.loadFrame('path', 'brightness', ...
            'ch',1, ...
            'f',1:8, ...
            'asch', 2, ...
            'ncpu', Inf);
lib.loadFOV();
lib.regroup('method', 'perwell'); 

%%% define calibration
lib.protocol.getChannel(1).binSize = 2;
lib.protocol.getChannel(2).binSize = 2;

% define stimulation from file
lib.protocol.setStimulation('path', 'Dim2Bright Long-20180726-095811.ams4100', ...
                            'separation', 0.1);

% temporal alignment
lib.getChildren('f').syncChannel('ch', 1);

    
%% Phase 2: Image Processing
p = spcore.pipeline.pipeline('1P FieldStim');
p.Base = lib;
p.addBlock(...
    'Name', 'Load Pixel', ...
    'operation', '@(I) I', ...
    'source', {':', ':', ':'}, ...
    'input', {{'I', {'rawFrames', 'C', 1}}, {'I', {'rawFrames', 'C', 2}}}, ...
    'output', {'I', {'raw', 'C', '', '', 'T', '', '', 'Z', '', ''}});
p.addBlock(...
    'Name', 'Saturation Mask (per channel)', ...
    'operation', '@(I, channel) channel.getSaturationMask(I)', ...
    'source', {':', ':', ':'}, ...
    'input', {{'I', {'raw', 'C', 1}, 'tag', {'raw', 'C', 1}}, ...
              {'I', {'raw', 'C', 2}, 'tag', {'raw', 'C', 2}}}, ...
    'output', {'I', {'satMaskCh', 'C', '', ''}});
p.addBlock(...
    'Name', 'Saturation Mask', ...
    'operation', '@(I) any(I, 3)', ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'satMaskCh', 'C', [], 'sparse', false}}, ...
    'output', {'data', 'satMask'});
p.addBlock(...
    'Name', 'Camera Calibration and Binning', ...
    'operation', '@(I, channel) channel.applyCorrection(I)', ...
    'source', {':', ':', ':'}, ...
    'input', {{'I', {'raw', 'C', 1}, 'tag', {'raw', 'C', 1}}, ...
              {'I', {'raw', 'C', 2}, 'tag', {'raw', 'C', 2}}}, ...
    'output', {'I', {'rawCalib', 'C', '', '', 'T', '', '', 'Z', '', ''}});
p.addBlock(...
    'Name', 'Background Correction Ref Channel', ...
    'operation', '@(I, M) img.bgcrt(I, img.bglevel(img.nanmasking(I, ~M), 16))', ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'rawCalib', 'C', 2}, 'data', 'satMask'}, ...
    'output', {'I', {'bgcrted', 'C', '', '', 'T', '', '', 'Z', '', ''}});
p.addBlock(...
    'Name', 'Background Correction Target Channel', ...
    'operation', '@(I, Iref, M) img.bgcrt(I, img.bglevel(img.nanmasking(Iref, ~M), 16))', ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'rawCalib', 'C', 1}, 'I', {'rawCalib', 'C', 1, 'T', 1}, 'data', 'satMask'}, ...
    'output', {'I', {'bgcrted', 'C', '', '', 'T', '', '', 'Z', '', ''}});
p.addBlock(...
    'Name', 'Trace Time', ...
    'operation', '@(t) t', ...
    'source', {':', ':', ':'}, ...
    'input', {'tag', {'bgcrted', 'T', 'C', 1}}, ...
    'output', {'data', 'trace_t'});
p.addBlock(...
    'Name', 'Baseline', ...
    'operation', "@(I)  math.isbaseline(reshape(mean(squeeze(I), 1:2), [], 1), 'mode')", ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 1, 'T', ''}}, ...
    'output', {'data', 'baselinetf'});
p.addBlock(...
    'Name', 'Detrend', ...
    'operation', "@(I,bl,t) img.detrend(squeeze(I), 3, t, 'idleIdx', bl, 'figshow', false, 'order', 2)", ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 1, 'T', ''}, 'data', 'baselinetf','data', 'trace_t'}, ...
    'output', {'data', 'Idet'});
p.addBlock(...
    'Name', 'Detrend Reshape', ...
    'operation', "@(I) reshape(I, [], size(I, 3))", ...
    'source', {':', ':', ':'}, ...
    'input', {'data', 'Idet'}, ...
    'output', {'data', 'IdetLbg'});
p.addBlock(...
    'Name', 'Preliminary mask', ...
    'operation', "@(M1, M2, Iref) ~M1 & mean(M2, 3) > img.bglevel(img.nanmasking(Iref, ~M1), 16)", ...
    'source', {':', ':', ':'}, ...
    'input', {'data', 'satMask','data', 'Idet','I', {'rawCalib', 'C', 1, 'T', 1},}, ...
    'output', {'data', 'premask'});
p.addBlock(...
    'Name', 'Temporal Alignment', ...
    'operation', "@(It,t,mask,stim) img.temporalAlignment(mean(It(mask(:), 10:end),1), t(10:end), 'method', 'template', 'template', stim, 'figshow', false)",...
    'source', {':', ':', ':'}, ...
    'input', {'data', 'IdetLbg','data', 'trace_t', 'data', 'premask','lib.protocol.stimulation'}, ...
    'output', {'data', 'delay'});
p.addBlock(...
    'Name', 'Stimulation Resample', ...
    'operation', "@(stim,delay,t) fillmissing(resample(timeseries(stim.y,stim.x-delay),t).Data,'constant',0)",...
    'source', {':', ':', ':'}, ...
    'input', {'lib.protocol.stimulation','data', 'delay','data', 'trace_t'}, ...
    'output', {'data','RSstim'});
p.addBlock(...
    'Name', 'Correlation mask', ...
    'operation', "@(It, M, S) img.corrFilt(It(M(:), 10:end), 2, -S(10:700), 'groupSize', 20, 'selectRatio', [0.005, 0.5]))", ...
    'source', {':', ':', ':'}, ...
    'input', {'data', 'IdetLbg','data', 'premask','data','RSstim'}, ...
    'output', {'data', 'finalMask'});
p.addBlock(...
    'Name', 'Refch Map', ...
    'operation', '@(I) I', ...
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 2, 'T', 1}}, ...
    'output', {'data', 'RMap'});
p.addBlock(...
    'Name', 'TgtCh Map at Start', ...
    'operation', '@(I) mean(I, 4)', ... 
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 1, 'T', 1:5}}, ...
    'output', {'data', 'GMap'});
p.addBlock(...
    'Name', 'Trace', ...
    'operation', "@(I, M) img.groupfun(I, M, @nanmean, 'toCell', true, 'ignoreGroup', 0)", ... 
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 1}, 'data', 'finalMask'}, ...
    'output', {'data', 'trace', 'omit', [], 'data', 'trace_n'});
p.addBlock(...
    'Name', 'Trace Binned', ...
    'operation', "@(I, M, mask) img.groupfun(I, discretize(log2(M.*mask), [-Inf,4,5,6,7,8,9,10,11,12,13,14,15,16]), @nanmean,'ignoreGroup', 1)", ... 
    'source', {':', ':', ':'}, ...
    'input', {'I', {'bgcrted', 'C', 1}, 'data', 'GMap', 'data', 'finalMask'}, ...
    'output', {'data', 'trace_bin', 'data', 'trace_bin_g', 'data', 'trace_bin_n'});

% p.getChildren('t').plot();      % visualize pipeline
p.getChildren('t').independentTaskGroup{1}.plot();  % visualize first connected group
% execute tasks
p.getChildren('t').execute('ncpu', Inf);

%% Phase 3: Export and Save
% data is saved at p.dataStore
% You can genenrate report of ROI quantifications
f = lib.getChildren('f');
T_fovData = p.dataStore.report(f);

% generate FOV descriptor
T_fovInfo = table({f.Coord}', string({f.Name})', ...
    string({f.getParent().wellName})', ...
    arrayfun(@(w) w.plate.ID, f.getParent()), ...
    'VariableNames', ["Coord", "Name", "wellName", "plateID"]);

stim = lib.protocol.stimulation;

%% Visualize mask
n = 8;
figure()
hold on
for i = 1:n
m = T_fovData.('finalMask')(i);
tr = T_fovData.('trace')(i);
subplot(2,n,i)
imshow(m{:,:});
subplot(2,n,i+n)
plot(tr{:})
end