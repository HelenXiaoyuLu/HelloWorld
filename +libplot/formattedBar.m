% Plot formatted errorbars from selected rows and column of a table
% USAGE:
%   [b,T_well_stats_sorted] = libplot.formattedBar(Table, 'dFF0_longStim_Mean', ...
%   'dFF0_shortStim_Mean', 'pp1std','dFF0_longStim_STD','pp2std','dFF0_shortStim_STD', ...
%   'sort_on', 'dFF0_shortStim_Mean', 'Name','Legal Name','ylabel', ...
%   'Response Amplitude \DeltaF/F0','norm2ctrl','ASAP1', 'savePath',savefigpath);
%
% INPUT:
%   REQUIRED: libTable, property 1 (pp1) which belongs to one of the variable names
%   OPTIONAL: property 2 (pp2) which belongs to one of the variable names
%   PARAMETERS (Name-value pairs):
%       pp1std/pp2std: standard deviation to be plotted as errorbar
%       sort_on: sort on which column, default: pp1
%       sort_direction: 'ascend','descend' or false. default: ascend
%       norm2ctrl: normalize bar values to control row. default: no control
%       Name: title of the column with standardized names. default: 'Name'
%       selectNames/selectRows: specify the rows to be plotted, either as
%       cell or arrary. Only one criteria is needed. default: plot all rows
%       xlabel: plot xlabel or no. default: true
%       ylabel: specify ylable. default: pp1
%       colormap: 3xn matrix. n >= number of properties. default: lines
%       masked_on: filter the table by which column. default: 'Area_Mean'
%       mskth: a number specify the masking threshold. default is empty.
%       datatip: properties & values to be added into the datatip
%       savePath: save to the path with specified extension at dpi 600
%
% OUTPUT: 
%   hBar: a handle to the barplot
%   libTableSorted: sorted table based on specified criteria
% Xiaoyu updated on May 20th, 2020

function [hBar,libTableSorted] = formattedBar(libTable,pp1,varargin)
    c = colormap(lines);
    p = inputParser;
    p.addRequired('libTable', @istable);
    p.addRequired('pp1', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addOptional('pp2', '', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('pp1std', '', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('pp2std', '', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('sort_on', pp1, @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('sort_direction','ascend', @ischar);    
    p.addParameter('norm2ctrl', '',  @ischar)% assign the control for normalization 
    p.addParameter('Name', 'Name', @(x) any(strcmp(x, libTable.Properties.VariableNames))); % specify the name column
    p.addParameter('selectNames', {}, @iscell) % plot specified libraries
    p.addParameter('selectRows', [], @isnumeric) % plot specified libraries
    p.addParameter('xlabel', true, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}))
    p.addParameter('ylabel', pp1, @ischar) 
    p.addParameter('colormap', c, @isnumeric);
    p.addParameter('masked_on', 'Area_Mean', @ischar);
    p.addParameter('mskth', [], @isnumeric);
    p.addParameter('savePath', '', @ischar);
    p.addParameter('datatip', {}, @iscell) % properties in labels    
    p.parse(libTable,pp1,varargin{:})
    pp2 = p.Results.pp2;  
    pp1std = p.Results.pp1std;
    pp2std = p.Results.pp2std;
    sort_on = p.Results.sort_on;
    sort_direction = p.Results.sort_direction;
    norm2ctrl = p.Results.norm2ctrl;
    Name = p.Results.Name;
    selectNames = p.Results.selectNames;  
    selectRows =  p.Results.selectRows;
    xlbl = p.Results.xlabel;
    ylbl = p.Results.ylabel; 
    cmap = p.Results.colormap;
    masked_on =  p.Results.masked_on;
    mskth =  p.Results.mskth;
    datatip = p.Results.datatip;
    savePath = p.Results.savePath; 
    
    libTable.Properties.RowNames = libTable.(Name);
    
    % Remove unwanted rows
    if ~isempty(selectRows)
        libTable = libTable(selectRows,:);
    elseif ~isempty(selectNames)
        libTable = libTable(selectNames,:);
    end
    
    % Sort table (default: sort by pp1)
    if sort_direction
        libTableSorted = sortrows(libTable,sort_on,sort_direction,'MissingPlacement','last');
    end
    
    % Mask table by specified column (default: no masking)
    if ~isempty(mskth)
        libTableSorted = libplot.tableparse(libTableSorted,'masked_on',masked_on,'mskth',mskth,'parseName', false);
    end
    
    % Normalized to ctrl (default: no normalization)
    ctrlval = [1,1]; 
    if ~isempty(norm2ctrl)
        normctrlsel = strmatch(norm2ctrl, libTableSorted.(Name),'exact');
        Tnormctrl = libTableSorted(normctrlsel,:);
        ctrlval(1) = nanmean(Tnormctrl.(pp1));
        if ~isempty(pp2)
            ctrlval(2) = nanmean(Tnormctrl.(pp2));
        end
    end
    
    % Plot w/ errorbar
    switch pp2
        case '' % 1-D            
            figure() % Relative Brightness
            hold on
            NameCell = cellstr(libTableSorted.(Name));
            X = 1:height(libTableSorted);
            hBar = bar(X,libTableSorted.(pp1)./ctrlval(1),'BarWidth',0.8,'LineWidth', 1);
            if ~isempty(pp1std)
                errorbar(hBar(1).XEndPoints,hBar(1).YData,(libTableSorted.(pp1std))./ctrlval(1),'.k');
            end
            if xlbl
                set(gca, 'XTickLabel', NameCell, 'XTick',1:numel(NameCell));
            end
            ylabel(ylbl);
            ax = gca;
            ax.XTickLabelRotation=90;
            ax.XAxis.TickLength = [0 0];
            hBar(1).FaceColor = cmap(1,:);	
            addTip(libTableSorted, datatip, pp1, pp2, hBar);
            set(gcf, 'Position',  [110, 100, 600, 400]);
            libplot.saveEXT(savePath);
            
        otherwise
            figure()
            hold on
            NameCell = cellstr(libTableSorted.(Name));
            X = 1:height(libTableSorted);
            hBar = bar(X,[(libTableSorted.(pp1))./ctrlval(1),(libTableSorted.(pp2))./ctrlval(2)],'grouped','BarWidth',1,'LineWidth', 1);
            if ~isempty(pp1std)
                errorbar(hBar(1).XEndPoints,hBar(1).YData,(libTableSorted.(pp1std))./ctrlval(1),'.k');
            end
            if ~isempty(pp2std)
                errorbar(hBar(2).XEndPoints,hBar(2).YData,(libTableSorted.(pp2std))./ctrlval(2),'.k');
            end
            if xlbl
                set(gca, 'XTickLabel',NameCell, 'XTick',1:numel(NameCell));
            end
            ylabel(ylbl);
            legend(pp1,pp2,'location','best','interpreter','none');
            ax = gca;
            ax.XTickLabelRotation = 90;
            ax.XAxis.TickLength = [0 0];
            hBar(1).FaceColor = cmap(1,:);	
            hBar(2).FaceColor = cmap(2,:);	           
            addTip(libTableSorted, datatip, pp1, pp2, hBar);            
            set(gcf, 'Position',  [710, 100, 600, 400]);
            libplot.saveEXT(savePath);
    end
end

% Add tips
function addTip(libTable, datatip, pp1, pp2, h)
    if ~isempty(datatip)
    h(1).DataTipTemplate.DataTipRows(1).Label = 'Rank';
    h(1).DataTipTemplate.DataTipRows(2).Label = pp1;
    for i = 1:length(datatip)
        tipcol = find(strcmp(datatip{i},libTable.Properties.VariableNames));
        row = dataTipTextRow(datatip{i},libTable{:,tipcol});
        h(1).DataTipTemplate.DataTipRows(end+1) = row;
    end
        if ~isempty(pp2)      
            h(2).DataTipTemplate.DataTipRows(1).Label = 'Rank';
            h(2).DataTipTemplate.DataTipRows(2).Label = pp2;
            for i = 1:length(datatip)
                tipcol = find(strcmp(datatip{i},libTable.Properties.VariableNames));
                row = dataTipTextRow(datatip{i},libTable{:,tipcol});
                h(2).DataTipTemplate.DataTipRows(end+1) = row;
            end
        end
    end
end

