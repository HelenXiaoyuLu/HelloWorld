function [T_group, T_well, T_fov, T_power] = FOVbase_2PSpectra_local(cfg)
% Wellbase_2PSpectra is a pipeline function that rename captured 2p images 
% imcompatible with file mananger and plot spectra
    
    arguments
        cfg (1,1) db.jobConfig = db.jobConfig()
    end
    
% Default Inputs
    cfg.default = db.jobConfig(...
        'path', 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200904 ramping power', ...
        'tgtCh\path', 'Spectra 10nm spacing 700-1080', ...
        'lookup', '20200905Lookup_spectra',...
        'powerspectra','20200904_peakvaluest_jobs_100mwrange.mat',...
        'brightnessthreshold',90,...
        'rename','false');
    
%% Phase 1: Load Raw Data and Rename
    rename = cfg.getConfig('rename');
    if eval(rename)       
        lib = spcore.library('Rename');
        lib.protocol.path = cfg.getConfig('path');
        T_lookup = readtable(fullfile(cfg.getConfig('path'),...
            cfg.getConfig('lookup')),'ReadRowNames',true,...
            'TextType','string','ReadVariablenames',true);
        lib.loadFrame(...
            'path',  cfg.getConfig('tgtCh\path'), ...
            'ncpu', Inf);  
        lib.loadFOV();

        % Rename original files 
        n_well =  numel(lib.getChildren('w'));
        for i = 1:n_well
            wi = lib.getChildren('w',i);
            wellname = wi.wellName();
            OrignalName = split(wi.Name(),'_Seq'); % for well name ignore frame
            T_lookup.OriginalName(wellname) = string(OrignalName{1});
        end

        % Get Original flist
        T_lookup = rmmissing(T_lookup);
        groupping = findgroups(T_lookup.Name);
        for i = 1:max(groupping) % for each group
            wellidx = find(groupping == i);
            for j = 1:numel(wellidx) % for each well 
                flist = dir(strcat(cfg.getConfig('path'),...
                    '\',cfg.getConfig('tgtCh\path'),'\',...
                    T_lookup.OriginalName(wellidx(j)),'*.nd2'));
                newname = T_lookup.Name(wellidx(j));
                for k = 1:numel(flist)
                    oldPath = flist(k).folder;
                    oldfName = flist(k).name;
                    newfName = strcat(newname,'_',num2str(j),'-',num2str(k),'.nd2');
                    oldf = fullfile(oldPath,oldfName);
                    newf = fullfile(oldPath,newfName);
                    movefile(oldf, newf, 'f');
                end
            end
        end
    else
        disp('No renaming is needed');
    end        
    
%% Phase 2: Load Raw Data again
    lib = spcore.library('Spectra');
    lib.protocol.path = cfg.getConfig('path');
    lib.loadFrame(...
        'path',  cfg.getConfig('tgtCh\path'), ...
        'ncpu', Inf);  
    lib.loadFOV();
    lib.regroup('method','name');
      % define calibration
    lib.protocol.getChannel(1).binSize = 1;
    brgate = cfg.getConfig('brightnessthreshold');
    T = load(fullfile(cfg.getConfig('path'),cfg.getConfig('powerspectra')));
    T_powerspectra = T.(string(fieldnames(T(1))));
    T_power = table();
    if isa(T_powerspectra,"table") % Use the power lookup table directly if it's table
        if size(T_powerspectra,2) == 2
            T_power = T_powerspectra;
        else
            warning('Dimensions of input power lookup table is not valid')
        end
    elseif isa(T_powerspectra,"double") 
        if size(T_powerspectra,1) == 2
            T_power.Wavelength = T_powerspectra(1,:)';
            T_power.("Power at Sample Plane") =  T_powerspectra(2,:)';
        elseif size(T_powerspectra,2) == 2
            T_power.Wavelength = T_powerspectra(:,1);
            T_power.("Power at Sample Plane") =  T_powerspectra(:,2);
        else
            warning('Dimensions of input power lookup table is not valid')
        end   
    else
        warning('Format of input power lookup table is not valid')
    end
    % empty table will be saved in output for invalid cases
        
%% Phase 3: image processing 
    p = spcore.pipeline('Spectra');
    p.Base = lib;
    p.addBlock(...
        'Name', 'Load Pixel', ...
        'operation', '@(I) I', ...
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "rawFrames", 'C', 1}, 'T', ':'}}, ...                 
        'output', {'frame', {'T', ':', 'C', ':', 'L', "raw"}}); 
    p.addBlock(...
        'Name', 'Saturation Mask', ...
        'operation', '@(I, channel) channel.getSaturationMask(I)', ...
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "raw", 'C', 1}, 'T', ':'}, 'literal', lib.protocol.getChannel(1)}, ...
        'output', {'frame', {'L', "satMaskCh", 'C', ':', 'T', NaN}}); 
    p.addBlock(...
        'Name', 'Saturation Mask', ...
        'operation', '@(I) any(I, 3)', ...
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "satMaskCh"}, 'C', ':'}}, ...
        'output', {'data', 'satMask'});   
    p.addBlock(...
        'Name', 'Camera Calibration and Binning', ...
        'operation', '@(I, channel) channel.applyCorrection(I)', ...
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "raw", 'C', 1}, 'T', ':'}, 'literal', lib.protocol.getChannel(1)}, ...
        'output', {'frame', {'T', ':', 'C', ':', 'L', "rawCalib"}});
    p.addBlock(...
        'Name', 'Background Correction Target Channel per frame', ...
        'operation', '@(I, M) img.bgcrt(I, img.bglevel(img.nanmasking(I, ~M), 12))', ...
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "rawCalib", 'C', 1}, 'T', ':'}, ...
                  'data', 'satMask'}, ...
        'output', {'frame', {'T', ':', 'C', ':', 'L', "bgcrted"}}); 
    p.addBlock(...
        'Name', 'TgtCh Map', ...
        'operation', '@(I) max(I, [], 3)', ... % Get the brightest frame
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "bgcrted", 'C', 1}, 'T', '@(x) x(:)'}}, ...
        'output', {'data', 'GMap'});      
    p.addBlock(...
        'Name', 'Final Mask', ...
        'operation', '@(M1, M2, threshold) ~M1 & M2 > threshold', ... % set a hard threshold 
        'source', {':', ':', ':'}, ...
        'input', {'data', 'satMask', 'data', 'GMap', 'literal', brgate}, ...
        'output', {'data', 'finalMask'});
    p.addBlock(...
        'Name', 'Trace', ...
        'operation', "@(I, M) img.groupfun(I, M, @nanmean, 'toCell', true, 'ignoreGroup', 0)", ... 
        'source', {':', ':', ':'}, ...
        'input', {'frame', {{'L', "bgcrted", 'C', 1}, 'T', ':'}, 'data', 'finalMask'}, ...
        'output', {'data', 'trace', 'omit', [], 'data', 'Area'});    
    p.getChildren('t').execute('ncpu', Inf);
    dataStore = p.dataStore;
    
    %% Phase 4: Export Result
    % Generate FOV descriptor
    f = lib.getChildren('f');
    T_fovData = dataStore.report(f);
    T_fov = table(cat(1, f.getParent().getParent().UID), ...
        cat(1, f.getParent().UID), ...
        cat(1, f.UID), ...
        arrayfun(@(fov) fov.frame.values(1).file.path, f), ...
        string({f.getParent().Name}'),...
        'VariableNames', ["Group","Well", "FOV", "Path", "Group Name"]);
    T_fov = cat(2, T_fov, T_fovData);
    
    % Generate well descriptor
    w = lib.getChildren('w');
    T_well = table(cat(1, w.getParent().UID), ...
        cat(1, w.UID), ...
        string({w.wellName})', ...
        arrayfun(@(w) w.plate.ID, w), ...
        'VariableNames', ["Group", "Well", "Well Name", "plateID"]);  
    
    % generate group descriptor
    g = lib.getChildren('g');
    T_group = table(cat(1, g.UID), ...
        string({g.Name})', ...
        'VariableNames', ["Group", "Name"]);   
 
      %% Phase 5: Stats
    % Reformat tracedata
    T_fov.TraceFinal = [T_fov.trace{:}]';        
    weightedMean = @(x, w) sum(x.*w, 1)./sum(w, 1);
    weightedStd = @(x, w) std(x, w, 1);
    
    % Pull from FOV to Well
    grouping = findgroups(T_fov.Well);
    T_fov2Well = table();
    T_fov2Well.Well = splitapply(@unique, T_fov.Well, grouping);
    T_fov2Well.("Group Name") = splitapply(@unique, T_fov.("Group Name"), grouping);
    % number of FOV per well
    T_fov2Well.n_FOV = splitapply(@numel, T_fov.Area, grouping);
    % number of total pixels per well
    T_fov2Well.Area = splitapply(@sum, T_fov.Area, grouping);
    % mean FOV trace per well
    T_fov2Well.TraceMean = splitapply(weightedMean, T_fov.TraceFinal(:,:), T_fov.Area, grouping);
    % std of FOV traces per well
    T_fov2Well.TraceStd = splitapply(weightedStd, T_fov.TraceFinal(:,:), T_fov.Area(:), grouping);   
    % Merge
    T_well = innerjoin(T_fov2Well, T_well, 'Keys', 'Well');
        
    % Pull from FOV to Group (spectra is done at FOV level)
    grouping = findgroups(T_fov.Group);
    T_fov2Group = table();
    T_fov2Group.Group = splitapply(@unique, T_fov.Group, grouping);
    % Name of group 
    T_fov2Group.("Group Name") = splitapply(@unique, T_fov.("Group Name"), grouping);
    % number of FOV per group
    T_fov2Group.n_FOV = splitapply(@numel, T_fov.Area, grouping);
    % number of total pixels per group
    T_fov2Group.Area = splitapply(@sum, T_fov.Area, grouping);
    % mean FOV trace per group
    T_fov2Group.TraceMean = splitapply(weightedMean, T_fov.TraceFinal(:,:), T_fov.Area, grouping);
    % std of FOV traces per group
    T_fov2Group.TraceStd = splitapply(weightedStd, T_fov.TraceFinal(:,:), T_fov.Area(:), grouping);   
    % Merge
    T_group = innerjoin(T_fov2Group, T_group, 'Keys', 'Group');
    
end