function [mask, ROIs, savePath] = cellSel(I, p)
    % CELLSEL displays an image and let user to draw ROI on it and 
    % returns a mask contains all ROIs. 
    % [mask, ROIs, savePath] = ui.cellSel(I, ...)
    % INPUT
    %   I (optional): image to be drawn on. Default: use examplar image. 
    %   figTitle (optional name-value pair): title to be prompt, 
    %       default: 'ROI Selection'.
    %   merge (optional name-value pair): method to merge ROIs, 
    %       'none': ROIs will be saved in separate images, use this if
    %           there are overlapping ROIs. 
    %       'label' (default): ROI will be stored as a uint16 label map. 
    %       'binary': ROI will be stored as a logical (uint8) binary map,
    %           touching ROIs will likely be merged. 
    %   initROI (optional name-value pair): initial ROIs points in cell 
    %       array, default: {}, i.e. none. 
    %   timeOut (optional name-value pair): time out seconds to close the 
    %       UI, default: 0, no time out. 
    % OUTPUT
    %   mask: the logical mask drawn by the user. 
    %   ROIs: cell array of ROI boundary points. 
    %   savePath: path of the exported matfile. 
    % 
    % Zhuohe Liu, harry.liu@rice.edu
    % Kevin Colbert, kevin.colbert@bcm.edu
    % St-Pierre Lab, May 2020
    
    arguments
        I {mustBeNonempty, mustBeNumeric} = imread('cameraman.tif');
        p.figTitle (1,1) string = "ROI Selection"
        p.merge (1,1) string {mustBeMember(p.merge, ...
            ["none", "label", "binary"])} = "label"
        p.initROI cell = {}
        p.timeOut (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    end
    
    warning('ui.cellSel will be deprecated soon. Use ui.segmenter instead. ');
    
    figTitle = p.figTitle;
    merge = p.merge;
    initROI = p.initROI;
    timeOut = p.timeOut;
    
    f = imtool(I, [min(I(:)), max(I(:))], 'InitialMagnification', 300);
    ax = imgca(f);
    f.Position =  [92 558 1680 420];
    imcontrast(f);
    f.UserData = struct('fin', false, ...
                        'savePath', '', ...
                        't', 0, ...
                        'timeOut', timeOut, ...
                        'merge', merge);
    f.CloseRequestFcn = @closeFigure;
    
    % add menu
    mMenu = uimenu('Parent', f, ...
        'Text', '&ROI Selection', ...
        'Tag', 'Edit');
    uimenu('Parent', mMenu, ...
        'Text', '&Add ROI', ...
        'Tag', 'AddROI', ...
        'Accelerator', 'i', ...
        'MenuSelectedFcn', @addROIMenuCallback);
    mMerge = uimenu('Parent', mMenu, ...
        'Text', '&Merge Method', ...
        'Tag', 'Merge');
    uimenu('Parent', mMerge, ...
        'Text', '&None', ...
        'Tag', 'MergeNone', ...
        'UserData', 'none', ...
        'MenuSelectedFcn', @mergeMenuCallback);
    uimenu('Parent', mMerge, ...
        'Text', '&Label', ...
        'Tag', 'MergeLabel', ...
        'UserData', 'label', ...
        'MenuSelectedFcn', @mergeMenuCallback);
    uimenu('Parent', mMerge, ...
        'Text', '&Binary', ...
        'Tag', 'MergeBinary', ...
        'UserData', 'binary', ...
        'MenuSelectedFcn', @mergeMenuCallback);
    currentMergeMenu = findobj(mMerge, 'UserData', merge);
    currentMergeMenu.Checked = true;
    uimenu('Parent', mMenu, ...
        'Text', '&Export ROIs...', ...
        'Accelerator', 'e', ...
        'Separator', true, ...
        'MenuSelectedFcn', @exportROIMenuCallback);
    uimenu('Parent', mMenu, ...
        'Text', '&Apply ROIs and Quit', ...
        'Accelerator', 'q', ...
        'MenuSelectedFcn', @(~, event) closeFigure(f, event));
    
    % load initial ROIs
    for i = 1:numel(initROI)
        h_roi = images.roi.Freehand('Parent', ax, ...
            'Position', initROI{i}, ...
            'Label', num2str(i), ...
            'Color', 'y', ...
            'UserData', i);
        h_roi.addlistener('MovingROI', @ROImoved);  % stop time out timer
        h_roi.addlistener('AddingWaypoint', @ROImoved);
        h_roi.addlistener('RemovingWaypoint', @ROImoved);
        h_roi.addlistener('DeletingROI', @ROImoved);
    end
    while f.UserData.fin == false
        if f.UserData.timeOut > 0
            f.Name = sprintf('%s (time out %d s)', figTitle, ...
                f.UserData.timeOut - f.UserData.t);
            pause(1);        
            f.UserData.t = f.UserData.t + 1;
            if f.UserData.timeOut - f.UserData.t <= 0 && ...
                    f.UserData.timeOut > 0
                f.UserData.fin = true;
            end
        else
            f.Name = figTitle;
            waitfor(f, 'UserData');
        end
    end
    [mask, ROIs] = getMask(f);
    savePath = f.UserData.savePath;
    delete(f);
end

function mergeMenuCallback(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    siblingMenus = src.Parent.Children;
    [siblingMenus.Checked] = deal(false);
    src.Checked = true;
    f.UserData.merge = src.UserData;
end

function [mask, ROIs] = getMask(f)
    ax = imgca(f);
    hROI = findobj(ax.Children, '-class', 'images.roi.Freehand');
    I = findobj(ax.Children, '-class', 'matlab.graphics.primitive.Image');
    w = size(I.CData, 2);
    h = size(I.CData, 1);
    ROIs = {};
    mask = {};
    if isempty(hROI)
        mask = false(h, w);
    else
        [~, idx] = sort([hROI.UserData], 'ascend'); % sort ROI by index
        hROI = hROI(idx);
    end
    for i = 1:numel(hROI)
        poly = hROI(i).Position;
        ROIs{i} = poly;                                     %#ok<AGROW>
        mask{i} = poly2mask(poly(:,1), poly(:,2), h, w);    %#ok<AGROW>
    end
    switch f.UserData.merge
        case 'none'
            mask = cat(3, mask{:});
        case 'binary'
            masktol = false(h, w);
            for i = 1:length(hROI)
                if any(masktol(:) & mask{i}(:))  % detect overlap
                    warning('Overlapping ROIs detected when merging ROIs. ');
                end
                masktol = masktol | mask{i};
            end
            mask = masktol;
        otherwise % 'label'
            masktol = zeros(h, w, 'uint16');
            for i = 1:length(hROI)
                if any(masktol(:) & mask{i}(:))  % detect overlap
                    warning('Overlapping ROIs detected when merging ROIs. ');
                end
                masktol(mask{i}) = uint16(i);
            end
            mask = masktol;
    end
end

function closeFigure(src, ~)
    qans = questdlg('Do you want to apply ROI and quit?', ...
        'uicellsel', 'Yes', 'Cancel', 'Cancel');
    if strcmp(qans, 'Yes')
        src.UserData.fin = true;
    end
end

function ROImoved(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    f.UserData.timeOut = 0;   % stop time out
end

function exportROIMenuCallback(src, ~)
    % Menu: Export ROIs
    f = ancestor(src, 'matlab.ui.Figure');
    ax = findobj(f.Children, '-class', 'matlab.graphics.axis.Axes');
    hROI = findobj(ax.Children, '-class', 'images.roi.Freehand');
    f.UserData.timeOut = 0;   % stop time out
    if isempty(hROI)
        msgbox('No ROI is selected! ', 'ROI Selection', 'warn');
        return;
    end
    [filename, pathname] = uiputfile( ...
        {'*.mat', 'MAT Files (*.mat)';
         '*.tif;*.tiff', 'TIF Files (*.tif;*.tiff)';
         '*.*',  'All Files (*.*)'}, ...
        'Save ROI', pwd);
    if isequal(filename, 0) || isequal(pathname, 0)
        return
    end
    f.UserData.savePath = fullfile(pathname, filename);
    [~, ~, fExt] = fileparts(f.UserData.savePath);
    [mask, ROIs] = getMask(f);
    switch lower(fExt)
        case {'.mat'}
            roistruct = struct('mask', mask, ...
                       'ROIs', {ROIs}, ...
                       'merge', f.UserData.merge);
            save(f.UserData.savePath, '-struct', 'roistruct', '-v7.3');
        case {'.tif', 'tiff'}
            switch f.UserData.merge
                case 'binary'
                    imwrite(mask, f.UserData.savePath);
                case 'none'
                    for i = 1:size(mask, 3)
                        imwrite(mask(:,:,i), f.UserData.savePath, ...
                            'WriteMode', 'append');
                    end
                otherwise % 'label'
                    imwrite(mask, f.UserData.savePath);
            end
        otherwise
            error('Unsupported export file extension: %s. ', fExt);
    end
    fprintf('ROI saved at: %s\n', f.UserData.savePath);
end

function addROIMenuCallback(src, ~)
    % Menu: Add ROI
    f = ancestor(src, 'matlab.ui.Figure');
    ax = findobj(f.Children, '-class', 'matlab.graphics.axis.Axes');
    f.UserData.timeOut = 0;   % stop time out
    hROI = findobj(ax.Children, '-class', 'images.roi.Freehand');
    % get new ROI ID
    if isempty(hROI)
        ID = 1;
    else
        ID = setdiff(1:(numel(hROI) + 1), [hROI.UserData]);
        ID = ID(1);
    end
    h = drawfreehand(ax, 'Label', num2str(ID), ...
        'LabelVisible', 'hover', ...
        'UserData', ID);
%    h = drawassisted(ax, 'Label', num2str(idx));
    if isempty(h.Position)  % cancel or empty
        h.delete();
        return
    end
    uimenu('Parent', h.ContextMenu, ...
        'Text', '&Rename ROI', ...
        'Tag', 'RenameROI', ...
        'MenuSelectedFcn', {@renameROI, h});
end

function renameROI(src, ~, h)
    f = ancestor(src, 'matlab.ui.Figure');
    me = h;
    ax = findobj(f.Children, '-class', 'matlab.graphics.axis.Axes');
    f.UserData.timeOut = 0;
    hROI = findobj(ax.Children, '-class', 'images.roi.Freehand');
    
    % dialog input box
    prompt = {'Enter new index for current ROI:'};
    dlgtitle = 'Input';
    dims = [1 35];
    answer = inputdlg(prompt, dlgtitle, dims);
    if isempty(answer)  % abort or empty
        return
    else
        newLabel = answer{1};
        newidx = str2double(newLabel);
        try
            validateattributes(newidx, {'numeric'}, {'nonnan', 'integer', 'positive', '<=', 65535});
        catch
            warning('Your input must be an integer that is from 1 to 65535. ');
            return
        end
    end
    if me.UserData ~= newidx && ismember(newidx, [hROI.UserData])
        warning('You are attempting to rename your ROI to an index already taken, please try again. ');
        return
    end
    me.UserData = newidx;
    me.Label = newLabel;
end