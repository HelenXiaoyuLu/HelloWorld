%% Load plate reader data
clear
clc
path = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\1P spectra\20210406 JEDI pH Spectra'; 
fName1 = '2021_04_06_GEVI Spectra II_1nm steps - merged.xlsx';
fName2 = '2021_04_06_GEVI Spectra.xlsx';
fLookUp = '20210406Lookup_GEVIpH';
TEx1 = readtable(fullfile(path,fName1),'Range','B47:AM203','ReadVariableNames',true);
TEx2 = readtable(fullfile(path,fName2),'Range','B47:AM98','ReadVariableNames',true);
% TEx3 = [readtable(fullfile(path,fName3),'Range','B47:B188','ReadVariableNames',true),...
%     readtable(fullfile(path,fName3),'Range','G47:L188','ReadVariableNames',true)];
TEm1 = readtable(fullfile(path,fName1),'Range','B207:AM398','ReadVariableNames',true);
TEm2 = readtable(fullfile(path,fName2),'Range','B102:AM166','ReadVariableNames',true);
% TEm3 = [readtable(fullfile(path,fName3),'Range','B192:B256','ReadVariableNames',true),...
%     readtable(fullfile(path,fName3),'Range','G192:L256','ReadVariableNames',true)];
lookUp = readtable(fullfile(path,fLookUp),'Range','A:B','ReadVariableNames',true);
lookUp.Properties.RowNames = lookUp.Well;
refBuffer = 'External_6pt55';
refCell = 'Kir';
movmeanwindow = 10;

%% Re-arrange data
TEx = TEx1;
TEm = TEm1;

TEx = fillmissing(TEx, 'constant', 0);
TEm = fillmissing(TEm, 'constant', 0);

ExScan1 = TEx.Wavelength';
EmScan1 = TEm.Wavelength';
Ex1 = table(); 
Ex1.Well = TEx.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Ex1)
    Ex1.Name{i} = lookUp.Name{Ex1.Well{i}};
    s = sig.arbitrary();
    s.setSignal('x',ExScan1, 'y', TEx{:,Ex1.Well{i}}');    
    traces(i, 1) = s;
end
Ex1.rawTrace = traces;
Em1 = table();
Em1.Well = TEm.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Em1)
    Em1.Name{i} = lookUp.Name{Em1.Well{i}};
    s = sig.arbitrary();
    s.setSignal('x', EmScan1, 'y', TEm{:,Em1.Well{i}}');    
    traces(i, 1) = s;
end
Em1.rawTrace = traces;

%% Group data
bgcrt = @(tr1,tr2) tr1.y - tr2.y;
gEx = table();
Exgroups = findgroups(Ex1.Name);
gEx.Name = splitapply(@unique, Ex1.Name, Exgroups);
gEx.Properties.RowNames = gEx.Name;
gEx.n_well = splitapply(@numel, Ex1.Name, Exgroups);
gEx.rawTrace = splitapply(@mean, Ex1.rawTrace,  Exgroups);
gEx.rawTraceStd = splitapply(@std, Ex1.rawTrace,  Exgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
normRawtraces = 'rawTrace normalized';
ctcrttraces = sig.arbitrary.empty();
ctsmoothtraces = sig.arbitrary.empty();
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
    smoothex = movmean(s2.y, movmeanwindow);
    s3 = sig.arbitrary();
    s3.setSignal('x', ExScan1, 'y', smoothex);
    [~, idx] = max(smoothex);
    ExPeakWL(i) = s2.x(idx);
    ctcrttraces(i, 1) = s2;
    ctsmoothtraces(i, 1) = s3;
end
gEx.('background corrected traces') = traces;
gEx.('background corrected traces normalized') = normtraces;
gEx.('celltype corrected Traces normalized') = ctcrttraces;
gEx.("celltype corrected Traces normalized movemean = " + movmeanwindow) = ctsmoothtraces;
gEx.("Excitation Peak") = ExPeakWL;

gEm = table();
Emgroups = findgroups(Em1.Name);
gEm.Name = splitapply(@unique, Em1.Name, Emgroups);
gEm.Properties.RowNames = gEm.Name;
gEm.n_well = splitapply(@numel, Em1.Name, Emgroups);
gEm.rawTrace = splitapply(@mean, Em1.rawTrace,  Emgroups);
gEm.rawTraceStd = splitapply(@std, Em1.rawTrace,  Emgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
ctcrttraces = sig.arbitrary.empty();
ctsmoothtraces = sig.arbitrary.empty();
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
    smoothem = movmean(s2.y, movmeanwindow);
    [~, idx] = max(smoothem);
    s3 = sig.arbitrary();
    s3.setSignal('x', EmScan1, 'y', smoothem);
    EmPeakWL(i) = s2.x(idx);
    ctcrttraces(i, 1) = s2;  
    ctsmoothtraces(i, 1) = s3;
end
gEm.('background corrected traces') = traces;
gEm.('background corrected traces normalized') = normtraces;
gEm.('celltype corrected Traces normalized') = ctcrttraces;
gEm.("celltype corrected Traces normalized movemean = " + movmeanwindow) = ctsmoothtraces;
gEm.("Emission Peak") = EmPeakWL;

%% Save
save(fullfile(path,'1Pspectra_from_platereader.mat'),'gEm','gEx','Ex1','Em1','lookUp');

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
        title(Variant, 'Interpreter', 'latex');
    end
    
    sgtitle(tracetype);
end

%% Group by pH
T_ex = table();
GEVIs = ["EGFP-CAAX", "dpOPT-CAAX", "lyn-EGFP", "lyn-dpOPT",  "JEDI-1P", "JEDI-2P"];
pH = ["6pt55", "7pt97"];
celltype = ["HEK", "Kir"];
T_ex.Name = GEVIs';
for i = 1:numel(pH)
    for j = 1:numel(celltype)
        condition = "pH" + pH(i) + "_" + celltype(j);
        conditionidx = contains(gEx.Name, pH(i)) & contains(gEx.Name, celltype(j));
        T_ex.(condition) = gEx(conditionidx,:);
    end
end
T_ex.Properties.RowNames = T_ex.Name;

T_em = table();
T_em.Name = GEVIs';
for i = 1:numel(pH)
    for j = 1:numel(celltype)
        condition = "pH" + pH(i) + "_" + celltype(j);
        conditionidx = contains(gEm.Name, pH(i)) & contains(gEm.Name, celltype(j));
        T_em.(condition) = gEm(conditionidx,:);
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
            title(celltype(j), 'Interpreter', 'latex');
            xlabel('wavelength (nm)')
            ylabel('Normalize Ex/Em')
            legend();
        end
    end
    sgtitle(GEVIsel(g));
end
