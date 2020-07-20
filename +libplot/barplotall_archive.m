function [T,meanctrl] = barplotall(libTable,pp1,varargin)
    %% parse input
    c = colormap(lines);
    p = inputParser;
    p.addRequired('libTable', @istable);
    p.addRequired('pp1', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addOptional('pp2', '', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addParameter('norm2ctrl',{}, @(n) validateattributes(n,{'cell'},{'size',[1,1]})) % assign the control for normalization 
    p.addParameter('controls', {}, @iscell) % assign the control for scatter
    p.addParameter('selectLib', {}, @iscell) % plot specified libraries
    p.addParameter('colormap', c, @isnumeric);
    p.addParameter('masked_on', 'masked', @ischar);
    p.addParameter('mskth', 10, @isnumeric);
    p.addParameter('FaceAlpha', 1, @isnumeric);
    p.addParameter('sorted','ascend',@isstr); % 'descend' ; 'false';
    p.addParameter('figshow', true, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));
    p.addParameter('plotctrl', true, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));  
    p.addParameter('datatip', {}, @iscell) % plot specified libraries        
    p.parse(libTable,pp1,varargin{:})
    pp2 = p.Results.pp2;
    controls = p.Results.controls;  
    selectLib = p.Results.selectLib;  
    figshow = p.Results.figshow; 
    cmap = p.Results.colormap;
    norm2ctrl = p.Results.norm2ctrl;
    masked_on = p.Results.masked_on;
    mskth =  p.Results.mskth;
    sorted = p.Results.sorted;
    plotctrl = p.Results.plotctrl;
    pp2 = p.Results.pp2;
    FaceAlpha = p.Results.FaceAlpha;
    datatip = p.Results.datatip;
    % Standardize table content
    if ~any(strcmp('pPlateWell', libTable.Properties.VariableNames))
        libTable = libplot.tableparse(libTable,'masked_on',masked_on,'mskth',mskth);
    end    
    %% Decide properties to plot    
    LibPos = unique(libTable.plibName);
    % Default: plot all libraries
    libsel = [];
    if ~isempty(selectLib)
        for i = 1:length(selectLib)
            libsel = [libsel;strmatch(selectLib{i}, libTable.plibName)] ;
        end
    else
        libsel = 1:height(libTable);
    end
    f1 = figure();
    hold on
    pp1col = find(strcmp(pp1,libTable.Properties.VariableNames));
    T = libTable(libsel,:);
    if ~strcmp(sorted,'false')
        T = sortrows(T,pp1,sorted);
    end
    b1 = bar(T{:,pp1col},'FaceColor',cmap(1,:),'FaceAlpha',FaceAlpha);
    b1.DataTipTemplate.DataTipRows(1).Label = 'rankOrder';
    b1.DataTipTemplate.DataTipRows(2).Label = pp1;
    xlabel('Number of construct screened')
    ylabel(pp1)
    % Default: plot no control
    meanctrl = table;
    meanctrl.Name = reshape(controls,[],1);
    pltctrl(libTable,controls,pp1,pp1col,libsel,plotctrl,meanctrl);  
    
    if ~isempty(pp2)
        pp2col = find(strcmp(pp2,libTable.Properties.VariableNames));
        b2 = bar(T{:,pp2col},'FaceColor',cmap(2,:),'FaceAlpha',FaceAlpha);
        b2.DataTipTemplate.DataTipRows(1).Label = 'rankOrder';
        b2.DataTipTemplate.DataTipRows(2).Label = pp2;
        pltctrl(libTable,controls,pp2,pp2col,libsel,plotctrl,meanctrl);
    end
    
    if ~isempty(datatip)
        for i = 1:length(datatip)
            tipcol = find(strcmp(datatip{i},T.Properties.VariableNames));
            row = dataTipTextRow(datatip{i},T{:,tipcol});
            b1.DataTipTemplate.DataTipRows(end+1) = row;
            if ~isempty(pp2)
                b2.DataTipTemplate.DataTipRows(end+1) = row;
            end
        end
    end
    
    if ~isempty(norm2ctrl)
        normctrlsel = strmatch(norm2ctrl{1}, libTable.plibName,'exact') ;
        Tnormctrl = libTable(normctrlsel,:);
        meanctrlppl = nanmean(Tnormctrl{:,pp1col});
        T.normpp1 = T{:,pp1col}/meanctrlppl;
        libTable.normpp1 = libTable{:,pp1col}/meanctrlppl;
        f2 = figure();
        hold on
        b1 = bar(T.normpp1,'FaceColor',cmap(1,:),'FaceAlpha',FaceAlpha);
        b1.DataTipTemplate.DataTipRows(1).Label = 'rankOrder';
        b1.DataTipTemplate.DataTipRows(2).Label = pp1;
        xlabel('Number of construct screened')
        ylabel(['Normalized ',pp1]) 
        normpp1col = find(strcmp('normpp1',T.Properties.VariableNames));
%         meannormctrl = table;
%         meannormctrl.Name = reshape(controls,[],1);
        pltctrl(libTable,controls,pp1,normpp1col,libsel,plotctrl,meanctrl);
        if ~isempty(pp2)
            T.normpp2 = T{:,pp2col}/meanctrlppl;
            libTable.normpp2 = libTable{:,pp2col}/meanctrlppl;
            b2 = bar(T.normpp2,'FaceColor',cmap(2,:),'FaceAlpha',FaceAlpha);
            b2.DataTipTemplate.DataTipRows(1).Label = 'rankOrder';
            b2.DataTipTemplate.DataTipRows(2).Label = pp2;
            normpp2col = find(strcmp('normpp2',T.Properties.VariableNames));
            pltctrl(libTable,controls,pp2,normpp2col,libsel,plotctrl,meanctrl);
        end
        if ~isempty(datatip)
            for i = 1:length(datatip)
                tipcol = find(strcmp(datatip{i},libTable.Properties.VariableNames));
                row = dataTipTextRow(datatip{i},T{:,tipcol});
                b1.DataTipTemplate.DataTipRows(end+1) = row;
                if ~isempty(pp2)
                    b2.DataTipTemplate.DataTipRows(end+1) = row;
                end
            end
        end
    end   
end

function pltctrl(libTable,controls,pp,ppcol,libsel,plotctrl,meanctrl)
    if ~isempty(controls)  
        ctrlsel = [];
        Varname = strcat(pp,'_',num2str(ppcol));
        meanctrl{:,end+1} = zeros(size(controls))';
        meanctrl.Properties.VariableNames{end} = Varname;
        for i = 1:length(controls)
            ctrlsel = [ctrlsel;strmatch(controls{i}, libTable.plibName,'exact')] ;
            Tctrl = libTable(ctrlsel,:);
            meanctrl{i,end} = mean(Tctrl{:,ppcol});
            if plotctrl
                stimName = split(pp,{'(',')'});
                if length(stimName)>1;
                    stimName = stimName{2};
                else
                    stimName = stimName{1};
                end
                libplot.pltHorizon(mean(Tctrl{:,ppcol}),length(libsel),{['\leftarrow',controls{i}],stimName})
            end
        end
    end 
end
