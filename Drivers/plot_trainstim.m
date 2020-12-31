%% load data
clear
clc
path = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201215 TrainStim test';
fname = 'Wellbase_1pFieldStimTrigger_output_T';
load(fullfile(path, fname))
T_well = T.Well;
T_group = T.Group;
sevents = T.Library.Stimulation;
seventRange = sevents.segment(0.1).range();
seventdx = seventRange.dxEvent;

%% Regroup data
% Pull from Well to Group
T_well2Group = table();
[T_well2Group.Group_UID, ~, grouping] = unique(T_well.Well_Parent, 'rows', 'stable');
nstim = sum(contains(T_well.Properties.VariableNames, 'dFF0 Stim'));
for i = 1:nstim
    T_well2Group.("dFF0 Stim " + i + " Mean") = splitapply(@mean, T_well.("dFF0 Stim "+ i + "(G)"), grouping);
    T_well2Group.("dFF0 Stim " + i + " Std") = splitapply(@std, T_well.("dFF0 Stim "+ i + "(G)"), grouping);
end
% Merge
T_group = innerjoin( T_group, T_well2Group,'Keys', 'Group_UID');
T_group.GEVI = splitapply(@parseGevi, T_group.Group_Name(:), (1:height(T_group))');
T_group.Voltage = splitapply(@parseVolt, T_group.Group_Name(:), (1:height(T_group))');
dff0idx = contains(T_group.Properties.VariableNames, 'dFF0 Stim') & contains(T_group.Properties.VariableNames, 'Mean');
dff0stdidx = contains(T_group.Properties.VariableNames, 'dFF0 Stim') & contains(T_group.Properties.VariableNames, 'Std');
T_group = fillmissing(T_group, 'constant', 0, 'DataVariables',@isnumeric);   
% Grouping the groups by sensor
grouping = findgroups(T_group.GEVI);

%% For the same sensor, check dFF0 incresing with stimulation time
V = 60;
gevis = ["ASAP1-cyOFP", "ASAP2s-cyOFP", "ASAP3-cyOFP", "ASAP2s-H152E-cyOFP",...
    "JEDI-1P-cyOFP", "JEDI-2P-cyOFP", "ASAP1-N124V-R406K-cyOFP", "ArcLightA242-cyOFP"];
f = figure('Name', V + "V", 'Color', [1, 1, 1]);
f.UserData.Voltage = V;
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:max(grouping)
    ax = nexttile(t);
    hold(ax, 'on');
    for idx = find(T_group.GEVI == gevis(i) & T_group.Voltage == V)
        errorbar(seventdx*1000, -T_group{idx, dff0idx}, T_group{idx, dff0stdidx});
        ax.XAxis.Limits = [0, 105];
        ax.YAxis.Limits = [-0.05, 0.3];
        ax.XAxis.Label.String = "Duration (ms)";
        ax.YAxis.Label.String = "-dF/F0";
        ax.Title.String = gevis(i);
    end
end

%% For the same sensor & duration, check dFF0 incresing with voltage increasing
S = 2;
Sspec = round(seventdx(S)*1000, 1) + "ms";
gevis = ["ASAP1-cyOFP", "ASAP2s-cyOFP", "ASAP3-cyOFP", "ASAP2s-H152E-cyOFP",...
    "JEDI-1P-cyOFP", "JEDI-2P-cyOFP", "ASAP1-N124V-R406K-cyOFP", "ArcLightA242-cyOFP"];
f = figure('Name', Sspec, 'Color', [1, 1, 1]);
f.UserData.Voltage = S;
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:max(grouping)
    ax = nexttile(t);
    hold(ax, 'on');
    for idx = find(T_group.GEVI == gevis(i))
        [~, ord] = sort(T_group.Voltage(idx),'ascend');
        idx = idx(ord);
        errorbar(T_group.Voltage(idx), -T_group.("dFF0 Stim " + i + " Mean")(idx), ...
            T_group.("dFF0 Stim " + i + " Std")(idx), 'LineStyle', '-');
        ax.XAxis.Limits = [0, 105];
        ax.YAxis.Limits = [-0.05, 0.3];
        ax.XAxis.Label.String = "Voltage (V)";
        ax.YAxis.Label.String = "-dF/F0";
        ax.Title.String = gevis(i);
    end
end

%% For the same sensor & voltage, plot FWHM
V = 60;
gevis = ["ASAP1-cyOFP", "ASAP2s-cyOFP", "ASAP3-cyOFP", "ASAP2s-H152E-cyOFP",...
    "JEDI-1P-cyOFP", "JEDI-2P-cyOFP", "ASAP1-N124V-R406K-cyOFP", "ArcLightA242-cyOFP"];
f = figure('Name', V + "V", 'Color', [1, 1, 1]);
f.UserData.Voltage = V;
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:max(grouping)
    ax = nexttile(t);
    hold(ax, 'on');
    for idx = find(T_group.GEVI == gevis(i) & T_group.Voltage == V)
        l = gobjects(1, numel(stimrange));
        for stimidx = 1:numel(stimrange)
            sidx = stimrange(stimidx);
            s = T_group.TraceWellMean(idx).copy();
            s.crop(seventRange.xOn(sidx)-0.2, seventRange.xOff(sidx)+ 0.4);   
            s = sig.arbitrary('x', s.x, 'y', math.normbase(s.y, 300:395, 1));
            Sroot = @(s) (s.y - 1)./min(s.y - 1)-0.5;
            Sroot.roots()
            ax.XAxis.Limits = [-0.1, 0.4];
            ax.YAxis.Limits = [-0.1, 1.05];
            ax.XAxis.Label.String = "Time (s)";
            ax.YAxis.Label.String = "-dF/F0";
            ax.Title.String = gevis(i);            
        end
        if ~isempty(lnote)
            legend(l,lnote)
        else
            legend(l)
        end
    end
end

%% For the same sensor & voltage, plot overlapping spike
V = 60;
stimrange = flip([1:9]);
% lnote = ["100ms", "20X1.5ms @ 200Hz", "10X1.5ms @ 100Hz"];
lnote = [];
cmap = jet(numel(stimrange));
% cmap = lines;
gevis = ["ASAP1-cyOFP", "ASAP2s-cyOFP", "ASAP3-cyOFP", "ASAP2s-H152E-cyOFP",...
    "JEDI-1P-cyOFP", "JEDI-2P-cyOFP", "ASAP1-N124V-R406K-cyOFP", "ArcLightA242-cyOFP"];
ReNormTrace = @(s, dff0) (s - 1)./(max(abs(s.y  - 1)))*sign(dff0);
f = figure('Name', V + "V", 'Color', [1, 1, 1]);
f.UserData.Voltage = V;
t = tiledlayout(f, 'flow', 'Padding', 'none');
for i = 1:max(grouping)
    ax = nexttile(t);
    hold(ax, 'on');
    for idx = find(T_group.GEVI == gevis(i) & T_group.Voltage == V)
        l = gobjects(1, numel(stimrange));
        for stimidx = 1:numel(stimrange)
            sidx = stimrange(stimidx);
            s = T_group.TraceWellMean(idx).copy();
            s.crop(seventRange.xOn(sidx)-0.4, seventRange.xOff(sidx)+ 0.4);   
            s = sig.arbitrary('x', s.x, 'y', math.normbase(s.y, 300:395, 1));
            l(stimidx) = plot(s.x - seventRange.xOn(sidx), (s.y - 1)./min(s.y - 1), 'Color', cmap(stimidx, :), ...
                'DisplayName' , round(seventdx(sidx)*1000,1) + " ms");
%             l(stimidx) = plot(s.x - seventRange.xOn(sidx), -(s.y - 1), 'Color', cmap(stimidx, :), ...
%                 'DisplayName', round(seventdx(sidx)*1000,1) + " ms");
            ax.XAxis.Limits = [-0.1, 0.4];
            ax.YAxis.Limits = [-0.1, 1.05];
            ax.XAxis.Label.String = "Time (s)";
            ax.YAxis.Label.String = "-dF/F0";
            ax.Title.String = gevis(i);            
        end
        if ~isempty(lnote)
            legend(l,lnote)
        else
            legend(l)
        end
    end
end

%%
function volt = parseVolt(name_volt)
    names = split(name_volt,{'_Train','V'});
    volt = str2double(names{end-1});
end

function gevi = parseGevi(name_volt)
    names = split(name_volt,{'_Train'});
    gevi = string(names{1});
end
