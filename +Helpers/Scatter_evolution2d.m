function [f, T_plot] = Scatter_evolution2d(j)
%   Plot the evolutionary path based on group level stats, 
%   on a 2D scatter plot, with user-defined axis
%   DEMO
%   [f, T_plot] = Scatter_2Devolution2d(j)
%   INPUT
%       j: a compatible job.
%   OUTPUT
%       f: figure handle. 
%       T_plot: data used for the plot. If this output is used, figure 
%           plotting will be suppressed.
% 
%   Xiaoyu Lu, xiaoyu.lu@rice.edu
%   St-Pierre Lab, Sep. 2020

    arguments
        j (1,1) struct
    end
    
    %% Preprocesssing
    T_plot =table();
    
    % brightness
    [~,T_ploti] = Wellbase_2PFieldStim.Bar_Brightness(j);
    T_plot.Name = T_ploti.Name;
    T_plot.Properties.RowNames = T_plot.Name;
    T_plot.N = T_ploti.N;
    T_plot.("Brightness Mean") = T_ploti.Mean;
    T_plot.("Brightness STD") = T_ploti.Std;
    
    % photostability
    [~,T_ploti] = Wellbase_2PFieldStim.Bar_Photostability(j);
    T_ploti.Properties.RowNames = T_ploti.Name;
    T_plot.("Photostability Mean") = T_ploti.Mean(T_plot.Properties.RowNames);
    T_plot.("Photostability STD") = T_ploti.Std(T_plot.Properties.RowNames);
    
    % spike
    [~,T_ploti] = Wellbase_2PFieldStim.Bar_TraceShort(j);
    T_ploti.Properties.RowNames = T_ploti.Name;
    T_plot.("-dF/F0 Short Mean") = - T_ploti.Mean(T_plot.Properties.RowNames);
    T_plot.("dF/F0 Short STD") = T_ploti.Std(T_plot.Properties.RowNames);    
    
    % Derived dimensions
    T_plot.("Detectability Index") = T_plot.("-dF/F0 Short Mean").* sqrt(T_plot.("Brightness Mean"));
    T_plot.("AUC unnormalized") = T_plot.("Photostability Mean").* T_plot.("Brightness Mean");
    isonevolpath = zeros(height(T_plot),1);
    for i = 1:height(T_plot)
        parsename = split(T_plot.Name(i),'-cyOFP');
        [~, ~, p] = isevol(parsename{1});
        T_plot.Protein(i) = p;
        if p.n == 427
            isonevolpath(i) = 1;
        end
    end
    % T_plot.("Is on evol path") = isonevolpath;
    T_plot = T_plot(find(isonevolpath),:);
    
    %% Plotting
    f = figure('Color', [1, 1, 1]);
    f.UserData.PlotData = T_plot;   % bind data
%     db.watermark(f, j);
    t = tiledlayout(f, 'flow', 'Padding', 'none');
    ax = nexttile(t);
    hold(ax, 'on');
    % Default: plot photostability vs dFF0, brightness as size
    cmap = lines;
    pp1 = "-dF/F0 Short Mean"; % x value
    pp2 = "Brightness Mean"; % y value
    pp3 = "Photostability Mean"; % size
    T_plot.(strcat(pp3," Normalize")) = normalize(T_plot.(pp3),'range')*50 + 20;
    s = gobjects(height(T_plot), 1);
    for i = 1:height(T_plot)
        s(i) = scatter(ax,T_plot.(pp1)(i), T_plot.(pp2)(i), ...
            T_plot.(strcat(pp3," Normalize"))(i), cmap(i,:),'filled',...
            'DisplayName',T_plot.Name(i),'UserData',i);
    end
    lgd = legend(s, 'Location', 'eastoutside','Interpreter','none');
    ax.XAxis.Label.String = pp1;
    ax.YAxis.Label.String = pp2;
    
    % Plot countour
    operation = "X.*sqrt(Y)";
    xvals = T_plot.(pp1);
    yvals = T_plot.(pp2);
    x = 0.9*min(xvals):0.1:1.1*max(xvals);
    y = 0.9*min(yvals):0.1:1.1*max(yvals);
    [ax,C,h] = plotcontour(x,y,'ax',ax,'operation',operation,...
        'note',operation);
    h.UserData = [pp1,pp2,operation];
    
    for i = 1:height(T_plot)
        proteins(i).Header = char(T_plot.Name(i));
        proteins(i).Sequence = char(T_plot.Protein(i).Sequence);
    end
JC_distances = seqpdist(proteins,'method','jukes-cantor','alphabet','AA', ...
                           'indels','pairwise-delete','squareform',true);
JC_distances(JC_distances == 0) = nan;
[mindist, minidx] = min(JC_distances,[],1,'omitnan');
findParent = T_plot.Name(minidx)
   %% Add selection to change axis
    mData = uimenu('Parent', f, ...
        'Tag', 'Data', ...
        'Text', 'Set Axis');
    uimenu('Parent', mData, ...
        'Text', 'Set X-axis', ...
        'Tag', 'SetX', ...
        'MenuSelectedFcn', @setX);
   uimenu('Parent', mData, ...
        'Text', 'Set Y-axis', ...
        'Tag', 'SetY', ...
        'MenuSelectedFcn', @setY); 
    
end

function [isinDB, isderiveDB, p] = isevol(proteinname)
    pdb = readtable(fullfile(pwd,'Helper\+mb\ProteinDB.csv'),...
        'ReadVariableNames',false,...
        'ReadRowNames',true,...
        'TextType','string');
    searchpdb = @(x1, x2) {contains(x1,x2)};
    cmppdb = @(x1, x2) {strcmp(x1,x2)};
    G = findgroups(pdb.Properties.RowNames);
    pdbmatch = cell2mat(splitapply(cmppdb, repmat(proteinname,[height(pdb),1]), string(pdb.Properties.RowNames), G));
    pdbcontain = cell2mat(splitapply(searchpdb, repmat(proteinname,[height(pdb),1]), string(pdb.Properties.RowNames), G));
    isinDB = any(pdbmatch);
    isderiveDB = any(pdbcontain);
    if isinDB  %in Database: assign sequence directly
        p = mb.protein.parseNotation(strcat('[',proteinname,']'));
        p.Sequence = pdb.Var1(pdbcontain);
    elseif ~isinDB && isderiveDB % derive sequence from Database
        precursorname = string(pdb.Properties.RowNames(pdbcontain));
        mutations = split(proteinname,precursorname);
        standardname = strcat('[',precursorname,']',mutations(2));
        islegal = mb.protein.standardizeNotation(standardname,'method',"syntax").Success;
        if islegal
            p = mb.protein.parseNotation(standardname);
            precursorseq = pdb.Var1(pdbcontain);
            ptemp = mb.protein('Name', p.Precursor.Name, 'Sequence', precursorseq);
            for k = 1:height(p.notationTable)
                ptemp = ptemp.addMutagenesis(p.notationTable.Position(k),p.notationTable.Destination(k));
            end
            p = ptemp;
            p.Name = strcat('[',precursorname,']',mutations(2));
        else
           p = mb.protein();
        end
    else
        p = mb.protein();
    end
end

function [ax,C,h] = plotcontour(x,y,p)
    arguments 
        x (1,:) double
        y (1,:) double
        p.ax (1,1) handle = nexttile(tiledlayout(figure(), 'flow', 'Padding', 'none'))
        p.operation (1,1) string = "X.*Y"
        p.LevelList (1,:) double = []
        p.note (1,:) string = ""
    end
    ax = p.ax;
    operation = p.operation;
    note = p.note;
    LevelList = p.LevelList;
    hold(ax, 'on');
    [X,Y] = meshgrid(x,y);
    Z = eval(operation);
    switch note
        case "" % default: text numbers)
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z,'ShowText','on');
            else
                [C,h] = contour(X,Y,Z,'ShowText','on','LevelList',LevelList);
            end
        case "off"
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z);
            else
                [C,h] = contour(X,Y,Z,'LevelList',LevelList);
            end
        otherwise
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z,'ShowText','on');     
                h.DisplayName = note;
            else
                [C,h] = contour(X,Y,Z,'ShowText','on','LevelList',LevelList);
                h.DisplayName = note;
            end
    end        
end

function setX(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    ax = f.CurrentAxes;
    T_plot = f.UserData.PlotData;
    ngroup = height(T_plot);
    xlist = string(T_plot.Properties.VariableNames);
    isvariablenum = zeros(1,length(xlist));
    for i = 1:length(xlist)
        if isnumeric(T_plot.(xlist(i))(1))
            if ~contains(xlist(i),'STD')
                isvariablenum(i) = 1;
            end
        end
    end
    namelist = xlist(find(isvariablenum));    
    [idx, tf] = listdlg('ListString', namelist, ...
        'Name', 'Variables', ...
        'SelectionMode', 'single', ...
        'ffs', 0, 'ListSize', [250, 160], ...
        'PromptString', 'Please select X-axis');
    if tf == 0  % cancel
        return
    end
    s = findobj(f, '-class', 'matlab.graphics.chart.primitive.Scatter');
    [~,sidx] = sort(cell2mat({s.UserData}),'ascend');
    s = s(sidx);
    pp1 = namelist(idx);
    for i = 1 : ngroup 
        s(i).XData = T_plot.(pp1)(i);
    end
    h = findobj(f, '-class', 'matlab.graphics.chart.primitive.Contour');
    operation = h(1).UserData(3);
    pp2 = h.UserData(2);
    xvals = T_plot.(pp1);
    yvals = T_plot.(pp2);
    x = 0.9*min(xvals):0.1:1.1*max(xvals);
    y = 0.9*min(yvals):0.1:1.1*max(yvals);
    [~,~,h1] = plotcontour(x,y,'ax',ax,'operation',operation,...
        'note',operation);
    h1.UserData = [pp1,pp2,operation];
    ax.XAxis.Label.String = pp1;
    delete(h)
end

function setY(src, ~)
    f = ancestor(src, 'matlab.ui.Figure');
    ax = f.CurrentAxes;
    T_plot = f.UserData.PlotData;
    ngroup = height(T_plot);
    ylist = string(T_plot.Properties.VariableNames);
    isvariablenum = zeros(1,length(ylist));
    for i = 1:length(ylist)
        if isnumeric(T_plot.(ylist(i))(1))
            if ~contains(ylist(i),'STD')
                isvariablenum(i) = 1;
            end
        end
    end
    namelist = ylist(find(isvariablenum));    
    [idx, tf] = listdlg('ListString', namelist, ...
        'Name', 'Variables', ...
        'SelectionMode', 'single', ...
        'ffs', 0, 'ListSize', [250, 160], ...
        'PromptString', 'Please select Y-axis');
    if tf == 0  % cancel
        return
    end
    s = findobj(f, '-class', 'matlab.graphics.chart.primitive.Scatter');
    [~,sidx] = sort(cell2mat({s.UserData}),'ascend');
    s = s(sidx);
    pp2 = namelist(idx);
    for i = 1 : ngroup 
        s(i).YData = T_plot.(pp2)(i);
    end
    h = findobj(f, '-class', 'matlab.graphics.chart.primitive.Contour');
    operation = h(1).UserData(3);
    pp1 = h.UserData(1);
    xvals = T_plot.(pp1);
    yvals = T_plot.(pp2);
    x = 0.9*min(xvals):0.1:1.1*max(xvals);
    y = 0.9*min(yvals):0.1:1.1*max(yvals);
    [~,~,h1] = plotcontour(x,y,'ax',ax,'operation',operation,...
        'note',operation);
    h1.UserData = [pp1,pp2,operation];
    ax.YAxis.Label.String = pp2;
    delete(h)
end