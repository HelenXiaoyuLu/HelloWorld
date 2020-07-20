function [T_group, T_well, T_fov] = Wellbase_1PPhotobleachingMonoF0offset(cfg)
    % Wellbase_1PPhotobleachingMono is a pipeline main function. 
    % This deals with photobleaching files and evaluate the traces during
    % defined frames ("tgtCh\frame": default = 0 = all frames), with a 
    % user-input threshold for brightness masking (default = 50). The
    % background correction and masking will be perform on the selected
    % frame ("tgtCh\sampleFrame": default = last frame)
    %
    % [T_group, T_well, T_fov] = Wellbase_1PPhotobleachingMonoF0offset(cfg)
    % INPUT
    %   cfg (optional): job configuration object. Default, test dataset is used. 
    %       See the source code for default values. 
    % OUTPUT
    %   T_group, T_well, T_fov: result tables. 
    % 
    % Xiaoyu Lu, St-Pierre Lab, June 2020
    % xiaoyu.lu@rice.edu
    
    arguments 
        cfg (1,1) db.jobConfig = db.jobConfig()
    end
    
    % Default Inputs
    cfg.default = db.jobConfig(...
        'path', 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\20200708 photoactivation train\Reversible photobleaching\30s x 50ms red',...
        'tgtCh\ch', 2,...
        'tgtCh\frame', 0 ,... % Exclude the unexpected signal at the beginning
        'tgtCh\sampleFrame', 0,... % appoint a frame for background correction
        'noiseLevel', 50); % Threshold for masking
    
    %% Phase 1: Import Files and Settings
    lib = sp.libData('1P Photobleaching');
    lib.protocol.path = cfg.getConfig('path');
    lib.loadFrame(...
        'path', '',...
        't',cfg.getConfig('tgtCh\frame'), ...
        'ch', cfg.getConfig('tgtCh\ch'));
    lib.loadFOV();
    lib.regroup('method', 'perwell');          % regroup
    T = cfg.getConfig('sampleFrame');
    noiseLvl = cfg.getConfig('noiseLevel');
    
%     %% Phase 1 alter: Import Files and Settings (default)
%     lib = sp.libData('1P Photobleaching');
%     lib.protocol.path = cfg.default.getConfig('path');
%     lib.loadFrame(...
%         'path', '',...
%         't',cfg.default.getConfig('tgtCh\frame'), ...
%         'ch', cfg.default.getConfig('tgtCh\ch'));
%     lib.loadFOV();
%     lib.regroup('method', 'perwell');          % regroup
%     T = cfg.default.getConfig('tgtCh\sampleFrame');
%     noiseLvl = cfg.default.getConfig('noiseLevel');
    
    %% Phase 2: Image Processing
    if T == 0 % Default: do background correction and masking on last frame
        T = lib.getChildren('f',1).nt;
    end
    p = spcore.pipeline.pipeline('1P Photobleaching');
    p.Base = lib;
    p.addBlock(...
        'Name', 'Load Pixel', ...
        'operation', '@(I) I', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'rawFrames', 'C', 1}}, ...
        'output', {'I', {'raw', 'C', '', '', 'T', '', '', 'Z', '', ''}});
    p.addBlock(...
        'Name', 'Saturation Mask', ...
        'operation', '@(I, channel) channel.getSaturationMask(I)', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'raw', 'C', 1}, 'tag', {'raw', 'C', 1}}, ...
        'output', {'data', 'satMask'});
    p.addBlock(...
        'Name', 'Background Correction', ... % using last frame
        'operation', '@(I, Iref, M) img.bgcrt(I, img.bglevel(img.nanmasking(Iref, ~M), 16))', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'raw', 'C', 1}, 'I', {'raw', 'C', 1, 'T', T}, 'data', 'satMask'}, ...
        'output', {'I', {'bgcrted', 'C', '', '', 'T', '', '', 'Z', '', ''}});
    p.addBlock(...
        'Name', 'Final Mask', ...
        'operation', '@(M1, I, noise) ~M1 & I > noise', ... 
        'source', {':', ':', ':'}, ...
        'input', {'data', 'satMask', 'I', {'bgcrted', 'C', 1, 'T', T}, 'literal', noiseLvl}, ...
        'output', {'data', 'finalMask'});
    p.addBlock(...
        'Name', 'Trace', ...
        'operation', "@(I, M) img.groupfun(I, M, @nanmean, 'toCell', true, 'ignoreGroup', 0)", ... 
        'source', {':', ':', ':'}, ...
        'input', {'I', {'bgcrted', 'C', 1}, 'data', 'finalMask'}, ...
        'output', {'data', 'trace', 'omit', [], 'data', 'Area'});
    p.addBlock(...
        'Name', 'Trace Time', ...
        'operation', '@(t) t', ...
        'source', {':', ':', ':'}, ...
        'input', {'tag', {'bgcrted', 'T', 'C', 1}}, ...
        'output', {'data', 'trace_t'});
    
    % execute tasks
    p.getChildren('t').execute('ncpu', Inf);
    
    %% Phase 3: Export Result
    % data is saved at p.dataStore
    % You can genenrate report of ROI quantifications
    f = lib.getChildren('f');
    T_fovData = p.dataStore.report(f);
    % generate FOV descriptor
    T_fov = table(cat(1, f.getParent().UID), ...
        cat(1, f.UID), ...
        arrayfun(@(fov) fov.frame(1).image(1).file.path, f), ...
        'VariableNames', ["Well", "FOV", "Path"]);
    T_fov = cat(2, T_fov, T_fovData);
    
    % create sig objects
    traceRaw = sig.arbitrary.empty();
    for i = 1:numel(f)
        s = sig.arbitrary();
        s.setSignal(T_fov.trace_t{i}, T_fov.trace{i});
        traceRaw(i, 1) = s;
    end
    T_fov.trace_t = [];
    T_fov.trace = [];
    T_fov.TraceRaw = traceRaw;
    
    % generate well descriptor
    w = lib.getChildren('w');
    T_well = table(cat(1, w.getParent().UID), ...
        cat(1, w.UID), ...
        string({w.wellName})', ...
        arrayfun(@(w) w.plate.ID, w), ...
        'VariableNames', ["Group", "Well", "wellName", "plateID"]);
    
    % generate group descriptor
    g = lib.getChildren('g');
    T_group = table(cat(1, g.UID), ...
        string({g.Name})', ...
        'VariableNames', ["Group", "Name"]);
    
    %% Phase 4: Trace Analysis
    % resample
    T_fov.TraceFinal = T_fov.TraceRaw.synchronize('dt', 0.5, 'method', 'intersect');
    % normalize
    for i = 1:height(T_fov)
        T_fov.TraceNorm(i) = T_fov.TraceFinal(i)./T_fov.TraceFinal(i).y(1);
    end
    innerjoinFcn = @(varargin) innerjoin(varargin{:});
    % Pull from FOV to Well
    T_fovWell = innerjoinFcn(T_fov, T_well,'Keys','Well');
    FOVwells = findgroups(T_fovWell.Group);
    T_well.n_FOV = splitapply(@numel, T_fov.Area, FOVwells);
    % number of total pixels per construct
    T_well.Area = splitapply(@sum, T_fov.Area, FOVwells);
    % mean ROI trace per construct
    T_well.TraceNormMean = splitapply(@mean, T_fov.TraceNorm, FOVwells);
    % std of ROI traces per construct
    T_well.TraceNormStd = splitapply(@std, T_fov.TraceNorm, FOVwells);
    
    % Pull from FOV to Construct
    T_fovGroup = innerjoinFcn(T_fovWell, T_group,'Keys','Group');
    FOVgroups = findgroups(T_fovGroup.Group);
    
    % number of ROI per construct
    T_group.n_FOV = splitapply(@numel, T_fov.Area, FOVgroups);
    % number of total pixels per construct
    T_group.Area = splitapply(@sum, T_fov.Area, FOVgroups);
    % mean ROI trace per construct
    T_group.TraceNormMean = splitapply(@mean, T_fov.TraceNorm, FOVgroups);
    % std of ROI traces per construct
    T_group.TraceNormStd = splitapply(@std, T_fov.TraceNorm, FOVgroups);
    
end