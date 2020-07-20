function [T_group, T_well, T_fov, T_roi] = Singlecell_1PPAManual(cfg)
    % Singlecell_1PFieldStimManual is a pipeline main function. 
    % [T_group, T_well, T_fov, T_roi] = Singlecell_1PFieldStimManual(cfg)
    % INPUT
    %   cfg (optional): job configuration object. Default, test dataset is used. 
    %       See the source code for default values. 
    % OUTPUT
    %   T_group, T_well, T_fov, T_roi: result tables. 
    % 
    % Zhuohe Liu, St-Pierre Lab, June 2020
    % harry.liu@rice.edu
    
    arguments
        cfg (1,1) db.jobConfig = db.jobConfig()
    end
    
    % Default Inputs
    cfg.default = db.jobConfig(...
        'path', 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20200628 photoactivation', ...
        'tgtCh\path', 'exp50ms', ...
        'segment\path', 'exp50ms\maskTiff', ...
        'segment\matchingRule', '<filename>*', ...
        'stim', 'pAtriggerMIX-20200628-163515.ams4100');
    
    %% Phase 1: Import Files and Settings
    lib = sp.libData('1P Photoactivation Manual');
    lib.protocol.path = cfg.getConfig('path');
    lib.loadFrame(...
        'path', cfg.getConfig('tgtCh\path'), ...
        'ch', [1,2], ...
        'ncpu', Inf);
    lib.loadFOV();
    lib.regroup('method', 'perwell');      % regroup

    % define stimulation from file
    lib.protocol.setStimulation('path', cfg.getConfig('stim'), ...
                                'separation', 0.1);
    %% Phase 1 alter: Import Files and Settings default
    lib = sp.libData('1P Photoactivation Manual');
    lib.protocol.path = cfg.default.getConfig('path');
    lib.loadFrame(...
        'path', cfg.default.getConfig('tgtCh\path'), ...
        'ch', [1,2], ...
        'ncpu', Inf);
    lib.loadFOV();
    lib.regroup('method', 'perwell');      % regroup

    % define stimulation from file
    lib.protocol.setStimulation('path', cfg.default.getConfig('stim'), ...
                                'separation', 0.1);    
    % temporal alignment
%     lib.getChildren('f').syncChannel('ch', 1);
    
    lib.protocol.setSegmenter(...
        'ID', 1, ...
        'method','import', ...
        'importPath', cfg.default.getConfig('segment\path'), ...
        'matchingRule', cfg.default.getConfig('segment\matchingRule'));
    
    %% Phase 2: Image Processing
    p = spcore.pipeline.pipeline('1P Photoactivation Manual');
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
        'Name', 'Background Correction Green Channel', ...
        'operation', '@(I, Iref, M) img.bgcrt(I, img.bglevel(img.nanmasking(Iref, ~M), 16))', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'rawCalib', 'C', 1}, 'I', {'rawCalib', 'C', 1, 'T', 10}, 'data', 'satMask'}, ...
        'output', {'I', {'bgcrted', 'C', '', '', 'T', '', '', 'Z', '', ''}});
    p.addBlock(...
        'Name', 'Background Correction Red Channel', ...
        'operation', '@(I, Iref, M) img.bgcrt(I, img.bglevel(img.nanmasking(Iref, ~M), 16))', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'rawCalib', 'C', 2}, 'I', {'rawCalib', 'C', 2, 'T', 10}, 'data', 'satMask'}, ...
        'output', {'I', {'bgcrted', 'C', '', '', 'T', '', '', 'Z', '', ''}});
    p.addBlock(...
        'Name', 'Import Mask', ...
        'operation', "@(f, segmenter) segmenter(1).segment(f)", ...
        'source', {':', ':', ':'}, ...
        'input', {'source', '', 'literal', lib.protocol.segmenter}, ...
        'output', {'I', {'mask', 'C', 1, 1}});
    p.addBlock(...
        'Name', 'Populate ROI', ...
        'operation', "@(M, M2) sp.roiData(M .* uint16(~M2))", ... % remove saturation pixels
        'source', {':', ':', ':'}, ...
        'input', {'I', {'mask', 'C', 1}, 'data', 'satMask'}, ...
        'output', {'roi', {'index', []}});
    p.addBlock(...
        'Name', 'Exampler Frame Green', ...
        'operation', '@(I) I', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'bgcrted', 'C', 1, 'T', 10}}, ...
        'output', {'data', 'GMap'});
    p.addBlock(...
        'Name', 'Exampler Frame Red', ...
        'operation', '@(I) I', ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'bgcrted', 'C', 2, 'T', 10}}, ...
        'output', {'data', 'RMap'});
    p.addBlock(...
        'Name', 'ROI Trace Green', ...
        'operation', "@(I, r) img.groupfun(I, r, @nanmean, 1:2, 'toCell', true)", ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'bgcrted', 'C', 1}, 'roi', 'index'}, ...
        'output', {'roi', 'trace_green', 'omit', [], 'roi', 'Area'});
    p.addBlock(...
        'Name', 'ROI Trace Red', ...
        'operation', "@(I, r) img.groupfun(I, r, @nanmean, 1:2, 'toCell', true)", ...
        'source', {':', ':', ':'}, ...
        'input', {'I', {'bgcrted', 'C', 2}, 'roi', 'index'}, ...
        'output', {'roi', 'trace_red', 'omit', [], 'roi', 'Area'});
    p.addBlock(...
        'Name', 'ROI Trace Time', ...
        'operation', "@(t, ~) t", ...
        'source', {':', ':', ':'}, ...
        'input', {'tag', {'bgcrted', 'T', 'C', 1}, 'roi', 'index'}, ... % include ROI index for dependency
        'output', {'roi', 'trace_t'});
    
    % execute tasks
    p.getChildren('t').execute('ncpu', 1);
    
    %% Phase 3: Export and Save
    % data is saved at p.dataStore
    % You can genenrate report of ROI quantifications
    r = lib.getChildren('r');
    T_roiData = p.dataStore.report(r);
    
    % generate ROI descriptor
    T_roi = table(cat(1, r.getParent().UID), 'VariableNames', "FOV");
    T_roi.Area = T_roiData.Area;
    roinames = cell(height(T_roiData),1);
    roiwellnames = cell(height(T_roiData),1);
    tracesred = sig.arbitrary.empty();
    tracesgreen = sig.arbitrary.empty();
    for i = 1:height(T_roiData)
        roinames{i} = r(i).getParent('w').Name;
        roiwellnames{i} = r(i).getParent('w').wellName;
        sr = sig.arbitrary();
        sr.setSignal(T_roiData.trace_t{i}, T_roiData.trace_red{i});
        sg = sig.arbitrary();
        sg.setSignal(T_roiData.trace_t{i} ,T_roiData.trace_green{i});
        tracesred(i, 1) = sr;
        tracesgreen(i,1) = sg;
    end
    T_roi.Name = roinames;
    T_roi.wellName = roiwellnames;
    T_roi.trace_red_bgcrt= tracesred;
    T_roi.trace_green_bgcrt= tracesgreen;

    % generate FOV descriptor
    f = lib.getChildren('f');
    fovgroups = findgroups(T_roi.FOV);
    T_fov = table(cat(1, f.getParent().UID), ...
        cat(1, f.UID), ...
        string({f.Name})', ...
        arrayfun(@(fov) fov.frame(1).image(1).file.path, f), ...
        'VariableNames', ["Well", "FOV", "Name", "Path"]);
    T_fov.trace_red_bgcrt = splitapply(@mean, T_roi.trace_red_bgcrt, T_roi.Area, fovgroups);
    T_fov.trace_green_bgcrt = splitapply(@mean, T_roi.trace_green_bgcrt, T_roi.Area, fovgroups);
    
    % generate well descriptor
    w = lib.getChildren('w');
    wellgroups = findgroups(strcat(T_roi.Name,T_roi.wellName));
    T_well = table(cat(1, w.getParent().UID), ...
        cat(1, w.UID), ...
        string({w.Name})', ...
        string({w.wellName})', ...
        arrayfun(@(w) w.plate.ID, w), ...
        'VariableNames', ["Group", "Well", "Name", "wellName", "plateID"]);
    T_well.trace_red_bgcrt = splitapply(@mean, T_roi.trace_red_bgcrt, T_roi.Area, wellgroups);
    T_well.trace_green_bgcrt = splitapply(@mean, T_roi.trace_green_bgcrt, T_roi.Area, wellgroups);
    
    % generate group descriptor
    g = lib.getChildren('g');
    groupgroups = findgroups(T_roi.Name);
    T_group = table(cat(1, g.UID), ...
        string({g.Name})', ...
        'VariableNames', ["Group", "Name"]);
    T_group.trace_red_bgcrt = splitapply(@mean, T_roi.trace_red_bgcrt, T_roi.Area, groupgroups);
    T_group.trace_green_bgcrt = splitapply(@mean, T_roi.trace_green_bgcrt, T_roi.Area, groupgroups);
    
end