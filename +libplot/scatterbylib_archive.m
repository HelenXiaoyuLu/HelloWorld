% scatterbylib do a scatter plot for x = property1 (pp1) and y = property2 (pp2)
% USAGE
%
% INPUT 
% libTable: table of library data 
% 
function scatterbylib(libTable, pp1, pp2, varargin)
    %% parse input
    c = colormap(lines);
    p = inputParser;
    p.addRequired('libTable', @istable);
    p.addRequired('pp1', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addRequired('pp2', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('pp3', '', @(x) any(strcmp(x, libTable.Properties.VariableNames))); % add a third property -> size of dots
    p.addParameter('norm2ctrl',{}, @(n) validateattributes(n,{'cell'},{'size',[1,1]})) % assign the control for normalization 
    p.addParameter('controls', {}, @iscell) % assign the control for scatter
    p.addParameter('selectLib', {}, @iscell) % plot specified libraries
    p.addParameter('colormap', c, @isnumeric);
    p.addParameter('mskth', 10, @isnumeric);
    p.addParameter('figshow', true, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));
    p.addParameter('datatip', {}, @iscell) % plot specified libraries    
    p.parse(libTable,pp1,pp2,varargin{:})
    pp3 = p.Results.pp3;    
    controls = p.Results.controls;  
    selectLib = p.Results.selectLib;  
    figshow = p.Results.figshow; 
    cmap = p.Results.colormap;
    norm2ctrl = p.Results.norm2ctrl;
    mskth =  p.Results.mskth;
    datatip = p.Results.datatip;
    % Standardize table content
    if ~any(strcmp('pPlateWell', libTable.Properties.VariableNames))
        libTable = libplot.tableparse(libTable,mskth);
    end
    
    %% Decide properties to plot
    LibPos = unique(libTable.plibName);
    % Default: plot no control
    ctrlsel = [];
    if ~isempty(controls)  
        for i = 1:length(controls)
            ctrlsel = [ctrlsel,strmatch(controls{i}, LibPos,'exact')] ;
        end
    end
    % Default: plot all libraries
    libsel = [];
    if ~isempty(selectLib)
        for i = 1:length(selectLib)
            libsel = [libsel,reshape(strmatch(selectLib{i}, LibPos),1,[])] ;
        end
    else
        libsel = 1:length(LibPos);
    end
    libsel = setdiff(libsel,ctrlsel);
    pp1col = find(strcmp(pp1,libTable.Properties.VariableNames));
    pp2col = find(strcmp(pp2,libTable.Properties.VariableNames));
    libcol = find(strcmp('plibName',libTable.Properties.VariableNames));
    pwcol = find(strcmp('pPlateWell',libTable.Properties.VariableNames));
    if ~isempty(pp3)
        pp3col = find(strcmp(pp3,libTable.Properties.VariableNames));
    end
    
    %% Scatter libraries
    figure()
    hold on
    for l = 1:length(libsel)
        librange = find(cellfun(@(x)~isempty(strfind(x,LibPos{libsel(l)})), libTable.plibName));
        if iscell(libTable{librange,pp1col})
            p1 = cell2mat(libTable{librange,pp1col});
        else
            p1 = libTable{librange,pp1col};
        end
        if iscell(libTable{librange,pp2col})
            p2 = cell2mat(libTable{librange,pp2col});
        else
            p2 = libTable{librange,pp2col};
        end
        s = scatter(p1,p2,[],cmap,'filled','MarkerFaceAlpha',0.5);
        s.DataTipTemplate.DataTipRows(1).Label = pp1;
        s.DataTipTemplate.DataTipRows(2).Label = pp2;
        row = dataTipTextRow('Library',libTable{librange,libcol});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('PlateWell',libTable{librange,pwcol});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        if ~isempty(datatip)
            for i = 1:length(datatip)
                tipcol = find(strcmp(datatip{i},libTable.Properties.VariableNames));
                row = dataTipTextRow(datatip{i},libTable{librange,tipcol});
                s.DataTipTemplate.DataTipRows(end+1) = row;
            end
        end
    end
    
    % Scatter controls    
    for c = 1:length(ctrlsel)
        ctrlrange = find(cellfun(@(x)~isempty(strfind(x,LibPos{ctrlsel(c)})), libTable.plibName));
        if iscell(libTable{ctrlrange,pp1col})
            p1 = cell2mat(libTable{ctrlrange,pp1col});
        else
            p1 = libTable{ctrlrange,pp1col};
        end
        if iscell(libTable{ctrlrange,pp2col})
            p2 = cell2mat(libTable{ctrlrange,pp2col});
        else
            p2 = libTable{ctrlrange,pp2col};
        end
        s = scatter(p1,p2,70,'rx','linewidth',2);
        s.DataTipTemplate.DataTipRows(1).Label = pp1;
        s.DataTipTemplate.DataTipRows(2).Label = pp2;
        row = dataTipTextRow('Library',libTable{ctrlrange,libcol});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('PlateWell',libTable{ctrlrange,pwcol});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        if ~isempty(datatip)
            for i = 1:length(datatip)
                tipcol = find(strcmp(datatip{i},libTable.Properties.VariableNames));
                row = dataTipTextRow(datatip{i},libTable{ctrlrange,tipcol});
                s.DataTipTemplate.DataTipRows(end+1) = row;
            end
        end
    end
%     xmin = min(libTable{:,pp1col});
%     xmax = max(libTable(:,pp1col));
    % x = xmin:0.1:v(xmax+0.5);
    % DI1 = 1./sqrt(x);
    % DI2 = 2./sqrt(x);
    % plot(x,DI1,'--k')
    % plot(x,DI2,'-.k')
    g = legend(LibPos([reshape(libsel,1,[]),reshape(ctrlsel,1,[])]),'location','eastoutside')
    g.ItemHitFcn = @ui.legendClick;
    xlabel(pp1)
    ylabel(pp2)
end


