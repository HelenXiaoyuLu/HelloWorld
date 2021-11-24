%% Load plate reader data
clear
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20210915 JEDI-1P One-photon spectra'; 
fName1 = '20210908-15 JEDI Spectra Summary';
T{1} = readtable(fullfile(path,fName1), 'Sheet', 'Excitation (pH7.4)',  'ReadVariableNames', true, 'ReadRowNames', true);
T{2} = readtable(fullfile(path,fName1), 'Sheet', 'Emission (pH7.4)',  'ReadVariableNames', true, 'ReadRowNames', true);
T{3} = readtable(fullfile(path,fName1), 'Sheet', 'Excitation (pH8.1)',  'ReadVariableNames', true, 'ReadRowNames', true);
T{4} = readtable(fullfile(path,fName1), 'Sheet', 'Emission (pH8.1)',  'ReadVariableNames', true, 'ReadRowNames', true);
clc

%% Re-arrange data'
Tg = {};
expression = "ASAP2s|ASAP3|JEDI_1P|EGFP_CAAX";
Condition = ["Excitation (pH7.4)", "Emission (pH7.4)", "Excitation (pH8.1)", "Emission (pH8.1)"];
for i = 1:numel(T)
    wavlen = str2num(cell2mat(T{i}.Properties.RowNames));
    for j = 1:width(T{i})
        s(j) = sig.arbitrary('x', wavlen, 'y', T{i}{:, j});   
    end
    variantnames = regexp(T{i}.Properties.VariableNames, expression, 'match');
    variantnames = string(cat(2, variantnames{:}));
    grouping = findgroups(variantnames);
    Tg{i} = table();
    Tg{i}.GroupName = splitapply(@unique, variantnames', grouping');
    Tg{i}.(Condition(i)) = table(splitapply(@mean, s', grouping'),...
        splitapply(@std, s', grouping'), ...
        splitapply(@numel, s', grouping'), ...
        'VariableNames', {'Mean', 'Std', 'N'});
    clearvars s
    Tg{i}.Properties.RowNames = Tg{i}.GroupName;
end
T_group = outerjoin(Tg{1}, Tg{2}, 'Keys', 'GroupName', 'LeftVariables', ...
    Tg{1}.Properties.VariableNames, 'RightVariables', 'Emission (pH7.4)');
T_group = outerjoin(T_group, Tg{3}, 'Keys', 'GroupName', 'LeftVariables',...
    T_group.Properties.VariableNames, 'RightVariables', 'Excitation (pH8.1)');
T_group = outerjoin(T_group, Tg{4}, 'Keys', 'GroupName', 'LeftVariables',...
    T_group.Properties.VariableNames, 'RightVariables', 'Emission (pH8.1)');
T_group.Properties.RowNames = T_group.GroupName;

%% Plot spectra per variant
cmap = lines;
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
% Plot pH = 7.4    
for i = 1:height(T_group)
    ax = nexttile(t);
    hold(ax, 'on');
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH7.4)"}.N;
    tr_ex = T_group{varname, "Excitation (pH7.4)"}.Mean;
    plot(ax, tr_ex.x, tr_ex.y, 'Color', cmap(6, :), ...
        'DisplayName', "Excitation", 'LineWidth', 1.5);
    tr_em = T_group{varname, "Emission (pH7.4)"}.Mean;
    plot(ax, tr_em.x, tr_em.y, 'Color', cmap(5, :), ...
        'DisplayName', "Emission", 'LineWidth', 1.5)
    ax.XLim = [350, 650];
    ax.YLim = [0, 1.1];
    ax.XLabel.String = "Wavelength (nm)";
    ax.YLabel.String = "Ex/Em normalized";
    ax.Title.String = sprintf('%s (N = %d)', varname, groupn);
    ax.Title.Interpreter = "Latex";
    legend('Location', 'northeast')
end

f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
% Plot pH = 8.1    
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH8.1)"}.N;
    if ~isnan(groupn)
        ax = nexttile(t);
        hold(ax, 'on');
        tr_ex = T_group{varname, "Excitation (pH8.1)"}.Mean;
        plot(ax, tr_ex.x, tr_ex.y, 'Color', cmap(6, :), ...
            'DisplayName', "Excitation", 'LineWidth', 1.5);
        tr_em = T_group{varname, "Emission (pH8.1)"}.Mean;
        plot(ax, tr_em.x, tr_em.y, 'Color', cmap(5, :), ...
            'DisplayName', "Emission", 'LineWidth', 1.5)
        ax.XLim = [350, 650];
        ax.YLim = [0, 1.1];
        ax.XLabel.String = "Wavelength (nm)";
        ax.YLabel.String = "Ex/Em normalized";
        ax.Title.String = sprintf('%s (N = %d)', varname, groupn);
        ax.Title.Interpreter = "Latex";
        legend('Location', 'northeast')
    end
end

%% Plot spectra per variant (overlay pH)
cmap = lines;
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
% Plot pH = 7.4    
for i = 1:height(T_group)
    ax = nexttile(t);
    hold(ax, 'on');
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH7.4)"}.N;
    tr_ex = T_group{varname, "Excitation (pH7.4)"}.Mean;
    plot(ax, tr_ex.x, tr_ex.y, 'Color', cmap(6, :), ...
        'DisplayName', "Excitation (pH7.4)", 'LineWidth', 1.5);
    tr_em = T_group{varname, "Emission (pH7.4)"}.Mean;
    plot(ax, tr_em.x, tr_em.y, 'Color', cmap(5, :), ...
        'DisplayName', "Emission (pH7.4)", 'LineWidth', 1.5)
    tr_ex2 = T_group{varname, "Excitation (pH8.1)"}.Mean;
        plot(ax, tr_ex2.x, tr_ex2.y, 'Color', cmap(6, :), ...
            'DisplayName', "Excitation (pH8.1)", 'LineWidth', 1.5, 'LineStyle', '--');
    tr_em2 = T_group{varname, "Emission (pH8.1)"}.Mean;
        plot(ax, tr_em2.x, tr_em2.y, 'Color', cmap(5, :), ...
            'DisplayName', "Emission (pH8.1)", 'LineWidth', 1.5, 'LineStyle', '--')
    ax.XLim = [350, 650];
    ax.YLim = [0, 1.1];
    ax.XLabel.String = "Wavelength (nm)";
    ax.YLabel.String = "Ex/Em normalized";
    ax.Title.String = sprintf('%s', varname);
    ax.Title.Interpreter = "Latex";
    legend('Location', 'northeast')
end

%% Plot spectra per condition
cmap = lines;
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');
% Plot pH = 7.4 excitation  
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH7.4)"}.N;
    tr_ex = T_group{varname, "Excitation (pH7.4)"}.Mean;
    plot(ax, tr_ex.x, tr_ex.y, 'Color', cmap(i, :), ...
        'DisplayName', sprintf('%s (N = %d)', varname, groupn), 'LineWidth', 1.5);
end
ax.XLim = [350, 550];
ax.YLim = [0, 1.1];
ax.XLabel.String = "Wavelength (nm)";
ax.YLabel.String = "Excitation normalized";
ax.Title.String = "Excitation";
legend('Location', 'eastoutside', 'Interpreter', "Latex");

ax = nexttile(t);
hold(ax, 'on');
% Plot pH = 7.4 emision  
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH7.4)"}.N;
    tr_em = T_group{varname, "Emission (pH7.4)"}.Mean;
    plot(ax, tr_em.x, tr_em.y, 'Color', cmap(i, :), ...
        'DisplayName', sprintf('%s (N = %d)', varname, groupn), 'LineWidth', 1.5)
end
ax.XLim = [450, 650];
ax.YLim = [0, 1.1];
ax.XLabel.String = "Wavelength (nm)";
ax.YLabel.String = "Emission normalized";
ax.Title.String = "Emission";
legend('Location', 'eastoutside', 'Interpreter', "Latex");

f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');
% Plot pH = 8.1 excitation  
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH8.1)"}.N;
    tr_ex = T_group{varname, "Excitation (pH8.1)"}.Mean;
    plot(ax, tr_ex.x, tr_ex.y, 'Color', cmap(i, :), ...
        'DisplayName', sprintf('%s (N = %d)', varname, groupn), 'LineWidth', 1.5);
end
ax.XLim = [350, 550];
ax.YLim = [0, 1.1];
ax.XLabel.String = "Wavelength (nm)";
ax.YLabel.String = "Excitation normalized";
ax.Title.String = "Excitation";
legend('Location', 'eastoutside', 'Interpreter', "Latex");

ax = nexttile(t);
hold(ax, 'on');
% Plot pH = 8.1 emision  
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH8.1)"}.N;
    tr_em = T_group{varname, "Emission (pH8.1)"}.Mean;
    plot(ax, tr_em.x, tr_em.y, 'Color', cmap(i, :), ...
        'DisplayName', sprintf('%s (N = %d)', varname, groupn), 'LineWidth', 1.5)
end
ax.XLim = [450, 650];
ax.YLim = [0, 1.1];
ax.XLabel.String = "Wavelength (nm)";
ax.YLabel.String = "Emission normalized";
ax.Title.String = "Emission";
legend('Location', 'eastoutside', 'Interpreter', "Latex");
%% Load filter
T_filt = readtable('\\Stpierrelab7910\e\Images\Xiaoyu\20210915 JEDI-1P One-photon spectra\Chroma 525-50', 'ReadVariableNames', true);
T_JEDI2P = readtable('D:\OneDrive - Baylor College of Medicine\JEDI-2P\Francois St-Pierre_In Vitro Data and Figures\Source\1P Spectra Summary', 'Sheet', 'Emission', 'ReadVariableNames', true, 'ReadRowNames', true, );
Tg2 = {};
expression = "JEDI_2P|EGFP";
wavlen = str2num(cell2mat(T_JEDI2P.Properties.RowNames));
for j = 1:width(Tg2)
    s(j) = sig.arbitrary('x', wavlen, 'y', T{i}{:, j});   
end
variantnames = regexp(T{i}.Properties.VariableNames, expression, 'match');
variantnames = string(cat(2, variantnames{:}));
grouping = findgroups(variantnames);
Tg2{i} = table();
Tg2{i}.GroupName = splitapply(@unique, variantnames', grouping');
Tg2{i}.(Condition(i)) = table(splitapply(@mean, s', grouping'),...
    splitapply(@std, s', grouping'), ...
    splitapply(@numel, s', grouping'), ...
    'VariableNames', {'Mean', 'Std', 'N'});
clearvars s
Tg2{i}.Properties.RowNames = Tg{i}.GroupName;


f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');
selvar = ["ASAP2s", "ASAP3", "JEDI-2P"]
for i = 1:height(T_group)
    varname = T_group.GroupName(i);
    groupn = T_group{varname, "Emission (pH7.4)"}.N;
    tr_em = T_group{varname, "Emission (pH7.4)"}.Mean;
    plot(ax, tr_em.x, tr_em.y, 'Color', cmap(i, :), ...
        'DisplayName', sprintf('%s (N = %d)', varname, groupn), 'LineWidth', 1.5)
end

%%
save(fullfile(path,'1Pspectra_from_platereader.mat'),'gEm','gEx','Ex','Em', 'T_ex', 'T_em', 'lookUp');
