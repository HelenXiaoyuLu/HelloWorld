%% Load plate reader data
clear
clc
path = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\1P spectra\20201104_1Pspectra'; 
fName = '20201104 1P Spectra 1nm spacing GEVI.xlsx';
fLookUp = '20201104Lookup_JEDIs';
TEx = [readtable(fullfile(path,fName),'Range','B47:J193','ReadVariableNames',true)];
TEm = [readtable(fullfile(path,fName),'Range','B197:J393','ReadVariableNames',true)];
lookUp = readtable(fullfile(path,fLookUp),'Range','A1:B9','ReadVariableNames',true);
lookUp.Properties.RowNames = lookUp.Well;
refBuffer = 'External';
refCell = 'Kir';

%% Re-arrange data
ExScan = TEx.Wavelength';
EmScan = TEm.Wavelength'; 
Ex = table(); 
Ex.Well = TEx.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Ex)
    Ex.Name{i} = lookUp.Name{Ex.Well{i}};
    s = sig.arbitrary('x', ExScan, 'y', TEx{:,Ex.Well{i}}');    
    traces(i, 1) = s;
end
Ex.rawTrace = traces;

Em = table();
Em.Well = TEm.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Em)
    Em.Name{i} = lookUp.Name{Em.Well{i}};
    s = sig.arbitrary('x', EmScan,'y', TEm{:,Em.Well{i}}');    
    traces(i, 1) = s;
end
Em.rawTrace = traces;

%% Group data
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
for i = 1:height(gEx)
    s = sig.arbitrary('x', ExScan, 'y', gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary('x', ExScan, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    name = split(gEx.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2 = sig.arbitrary('x', ExScan, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(refCell).y,'range'));
    else
        s2 = sig.arbitrary('x', ExScan, 'y', normalize(gEx.rawTrace(i).y-gEx.rawTrace(name{1}).y,'range'));
    end
    ctcrttraces(i, 1) = s2;
end
gEx.('background corrected traces') = traces;
gEx.('background corrected traces normalized') = normtraces;
gEx.('celltype corrected Traces normalized') = ctcrttraces;

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
for i = 1:height(gEm)
    s = sig.arbitrary('x', EmScan, 'y', gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary('x', EmScan, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    name = split(gEm.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2 = sig.arbitrary('x', EmScan, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(refCell).y,'range'));
    else
        s2 = sig.arbitrary('x', EmScan, 'y', normalize(gEm.rawTrace(i).y-gEm.rawTrace(name{1}).y,'range'));
    end
    ctcrttraces(i, 1) = s2;    
end
gEm.('background corrected traces') = traces;
gEm.('background corrected traces normalized') = normtraces;
gEm.('celltype corrected Traces normalized') = ctcrttraces;

%% Save
save(fullfile(path,'20201104 1P Spectra 1nm spacing GEVI.mat'),'gEm','gEx','Ex','Em','lookUp');

%%
cmap = colormap(lines);
colNames = ["background corrected traces";...
    "background corrected traces normalized";...
    "celltype corrected Traces normalized"];
sel = ["EGFP"; "EGFP-CAAX"; "ASAP2s"; "JEDI-1P"; "JEDI-2P"; "Kir"];
for j = 1:length(colNames)
    f = figure();
    t = tiledlayout(f, 'flow', 'Padding', 'none');
    colName = colNames{j};
    t.Title.String = colName;
    for i = 1:length(sel)
        ax = nexttile(t);
        hold(ax, 'on');
        chrom = sel(i);
        if chrom ~= refCell;
    %         plot(gEx.(colName)('Kir_EGFP').x,gEx.(colName)('Kir_EGFP').y,'b-')
            l(1) = plot(gEx.(colName)("HEK_" + chrom).x,gEx.(colName)("HEK_" + chrom).y, ...
                'color', cmap(6,:), 'DisplayName', "Ex-HEK");
    %         plot(gEm.(colName)('Kir_EGFP').x,gEm.(colName)('Kir_EGFP').y,'g-')
            l(2) = plot(gEm.(colName)("HEK_" + chrom).x,gEm.(colName)("HEK_" + chrom).y,...
                'color',cmap(5,:), 'DisplayName', "Em-HEK");
        else
                %         plot(gEx.(colName)('Kir_EGFP').x,gEx.(colName)('Kir_EGFP').y,'b-')
            l(1) = plot(gEx.(colName)(chrom).x,gEx.(colName)(chrom).y, ...
                'color', cmap(6,:), 'DisplayName', "Ex-HEK");
    %         plot(gEm.(colName)('Kir_EGFP').x,gEm.(colName)('Kir_EGFP').y,'g-')
            l(2) = plot(gEm.(colName)(chrom).x,gEm.(colName)(chrom).y,...
                'color',cmap(5,:), 'DisplayName', "Em-HEK");
        end
        ax.Title.String = chrom;
        legend(l)
        xlabel(ax, 'Wavelength (nm)');
        ylabel(ax, 'Fluorescent Intensity (a.u.)');
    end
end
%% merge traces
f = figure();
cmap = lines;
t = tiledlayout(f, 'flow', 'Padding', 'none');
colName = colNames{3};
t.Title.String = colName;
% Excitation
ax = nexttile(t);
hold(ax, 'on');
for i = 1:length(sel)-1
    chrom = sel(i);
    l(i) = plot(gEx.(colName)("HEK_" + chrom).x,gEx.(colName)("HEK_" + chrom).y, ...
                'color', cmap(i,:), 'DisplayName', chrom);
end
legend(l);
ax.Title.String = "Excitation Spectra";
xlabel(ax, 'Wavelength (nm)');
ylabel(ax, 'Fluorescent Intensity (a.u.)');

% Emission
ax = nexttile(t);
hold(ax, 'on');
for i = 1:length(sel)-1
    chrom = sel(i);
    l(i) = plot(gEm.(colName)("HEK_" + chrom).x,gEm.(colName)("HEK_" + chrom).y, ...
                'color', cmap(i,:), 'DisplayName', chrom);
end
legend(l);
ax.Title.String = "Excitation Spectra";
xlabel(ax, 'Wavelength (nm)');
ylabel(ax, 'Normalized Fluorescence');
