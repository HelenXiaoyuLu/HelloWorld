%% Load plate reader data
clear
clc
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20210408 JEDI pH Spectra and Isosbestic point\Isosbestic point'; 
fName1 = '20210408 1P spectra Xiaoyu 2nd trial default export - merged.xlsx';
fLookUp = '20210408Lookup_GEVIIsosbesticpoint';

TEx1 = readtable(fullfile(path,fName1),'Range','B47:AL233','ReadVariableNames',true);
% TEx2 = readtable(fullfile(path,fName2),'Range','B47:AL98','ReadVariableNames',true);
TEm1 = readtable(fullfile(path,fName1),'Range','B237:AM398','ReadVariableNames',true);
% TEm2 = readtable(fullfile(path,fName2),'Range','B102:AL166','ReadVariableNames',true);
lookUp = readtable(fullfile(path,fLookUp),'Range','A:B','ReadVariableNames',true);
lookUp.Properties.RowNames = lookUp.Well;
refBuffer = 'External_6pt55';
refCell = 'Kir';
movmeanwindow = 5;

%% Re-arrange data
TEx = TEx1;
TEm = TEm1;

TEx = fillmissing(TEx, 'constant', 0);
TEm = fillmissing(TEm, 'constant', 0);

ExScan1 = TEx.Wavelength';
EmScan1 = TEm.Wavelength';
Ex = table(); 
Ex.Well = TEx.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Ex)
    Ex.Name{i} = lookUp.Name{Ex.Well{i}};
    s = sig.arbitrary();
    s.setSignal('x',ExScan1, 'y', movmean(TEx{:,Ex.Well{i}}', movmeanwindow));    
    traces(i, 1) = s;
end
Ex.rawTrace = traces;

Em = table();
Em.Well = TEm.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Em)
    Em.Name{i} = lookUp.Name{Em.Well{i}};
    s = sig.arbitrary();
    s.setSignal('x', EmScan1, 'y', TEm{:,Em.Well{i}}');    
    traces(i, 1) = s;
end
Em.rawTrace = traces;

%% Group data
parsegevi = @(s) extractBefore(s, '_');
bgcrt = @(tr1,tr2) tr1.y - tr2.y;
gEx = table();
Exgroups = findgroups(Ex.Name);
gEx.Name = splitapply(@unique, Ex.Name, Exgroups);
gEx.Properties.RowNames = gEx.Name;
gEx.n_well = splitapply(@numel, Ex.Name, Exgroups);
gEx.rawTrace = splitapply(@mean, Ex.rawTrace,  Exgroups);
gEx.rawTraceStd = splitapply(@std, Ex.rawTrace,  Exgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
ctcrttraces = sig.arbitrary.empty();
ExPeakWL = zeros(height(gEx), 1);
for i = 1:height(gEx)
    s = sig.arbitrary();
    s.setSignal('x', ExScan1, 'y', gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary();
    s1.setSignal('x', ExScan1, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    s2 = sig.arbitrary();
    name = split(gEx.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2.setSignal('x', ExScan1, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(refCell).y,'range'));
    else
        s2.setSignal('x', ExScan1, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(name{1}).y,'range'));
    end
    [~, idx] = max(s2.y);
    ExPeakWL(i) = s2.x(idx);
    ctcrttraces(i, 1) = s2;
end
gEx.('background corrected traces') = traces;
gEx.('background corrected traces normalized') = normtraces;
gEx.('celltype corrected Traces normalized') = ctcrttraces;
gEx.("Excitation Peak") = ExPeakWL;
gEx.("GEVI") = string(cat(2, arrayfun(parsegevi, gEx.Name)));

gEm = table();
Emgroups = findgroups(Em.Name);
gEm.Name = splitapply(@unique, Em.Name, Emgroups);
gEm.Properties.RowNames = gEm.Name;
gEm.n_well = splitapply(@numel, Em.Name, Emgroups);
gEm.rawTrace = splitapply(@mean, Em.rawTrace,  Emgroups);
gEm.rawTraceStd = splitapply(@std, Em.rawTrace,  Emgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
ctcrttraces = sig.arbitrary.empty();
EmPeakWL = zeros(height(gEm), 1);
for i = 1:height(gEm)
    s = sig.arbitrary();
    s.setSignal('x', EmScan1, 'y', gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary();
    s1.setSignal('x', EmScan1, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    s2 = sig.arbitrary();
    name = split(gEm.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2.setSignal('x', EmScan1, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(refCell).y,'range'));
    else
        s2.setSignal('x', EmScan1, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(name{1}).y,'range'));
    end
    [~, idx] = max(s2.y);
    EmPeakWL(i) = s2.x(idx);
    ctcrttraces(i, 1) = s2;  
end
gEm.('background corrected traces') = traces;
gEm.('background corrected traces normalized') = normtraces;
gEm.('celltype corrected Traces normalized') = ctcrttraces;
gEm.("Emission Peak") = EmPeakWL;
gEm.("GEVI") = string(cat(2, arrayfun(parsegevi, gEm.Name)));

%% Figure in general
cmap = colormap(lines);
Variants = lookUp.Name(~contains(lookUp.Name, 'External') & ~contains(lookUp.Name, 'lyn-dpOPT'));
colNames  = {'rawTrace', 'background corrected traces','background corrected traces normalized',...
    'celltype corrected Traces normalized'};

for tracetypes = 1:length(colNames)
    f = figure('Color', [1, 1, 1]);
    t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
    
    tracetype = string(colNames(tracetypes));
    for i = 1:numel(Variants) 
        ax = nexttile(t);
        hold(ax, 'on');
        Variant = Variants{i};
        trEx = gEx.(tracetype)(Variant);
        trEm = gEm.(tracetype)(Variant);
        plot(trEx.x, trEx.y, 'color',cmap(6,:), 'DisplayName', Variant);
        plot(gEx.("Excitation Peak")(Variant), trEx.x2y(gEx.("Excitation Peak")(Variant)), '+r');
        plot(trEm.x, trEm.y, 'color',cmap(5,:), 'DisplayName', Variant);
        plot(gEm.("Emission Peak")(Variant), trEm.x2y(gEm.("Emission Peak")(Variant)), '+r');
        title(Variant);
    end
    
    sgtitle(tracetype);
end

%% Group by pH
gExrmvExternal = gEx(~contains(gEx.GEVI, 'External'), :);
gExrmvExternal = gExrmvExternal(~(gExrmvExternal.GEVI == ""),:);
T_ex = table();
GEVIs = gExrmvExternal.GEVI;
pH = ["6pt55", "7pt97"];
celltype = ["HEK", "Kir"];
grouping = findgroups(GEVIs);
T_ex.Name = splitapply(@unique, GEVIs, grouping);
for i = 1:numel(pH)
    for j = 1:numel(celltype)
        condition = "pH" + pH(i) + "_" + celltype(j);
        conditionidx = contains(gExrmvExternal.Name, pH(i)) & contains(gExrmvExternal.Name, celltype(j));
        T_ex.(condition) = gExrmvExternal(conditionidx,:);
    end
end
T_ex.Properties.RowNames = T_ex.Name;

gEmrmvExternal = gEm(~contains(gEm.GEVI, 'External'), :);
T_em = table();
T_em.Name = splitapply(@unique, GEVIs, grouping);
for i = 1:numel(pH)
    for j = 1:numel(celltype)
        condition = "pH" + pH(i) + "_" + celltype(j);
        conditionidx = contains(gEmrmvExternal.Name, pH(i)) & contains(gEmrmvExternal.Name, celltype(j));
        T_em.(condition) = gEmrmvExternal(conditionidx,:);
    end
end
T_em.Properties.RowNames = T_em.Name;


%% Plot compare pH
close all
GEVIsel = ["EGFP-CAAX", "lyn-EGFP", "JEDI-1P", "JEDI-2P"];
ls = {'-', '--'};
TraceType = "celltype corrected Traces normalized";
for g = 1:numel(GEVIsel)
    f = figure('Color', [1, 1, 1]);
    t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
    for j = 1:numel(celltype)
        ax = nexttile(t);
        hold(ax, 'on');
        for i = 1:numel(pH)
            condition = "pH" + pH(i) + "_" + celltype(j);
            T_crop = T_ex(GEVIsel(g),:);
            plot(ax, T_crop.(condition).(TraceType)(1).x, ...
                T_crop.(condition).(TraceType)(1).y, ...
                'DisplayName', 'Excitation ' + pH(i), 'color',cmap(6,:), 'LineStyle', ls{i});
            T_crop = T_em(GEVIsel(g),:);
            plot(ax, T_crop.(condition).(TraceType)(1).x, ...
                T_crop.(condition).(TraceType)(1).y, ...
                'DisplayName', 'Emission ' + pH(i), 'color',cmap(5,:), 'LineStyle', ls{i})
            title(condition, 'Interpreter', 'latex');
            legend();
        end
    end
    sgtitle(GEVIsel(g));
end
%%

save(fullfile(path,'1Pspectra_from_platereader.mat'),'gEm','gEx','Ex1','Em1','lookUp');
