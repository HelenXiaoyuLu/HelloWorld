load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\2P photobleaching power test\20190925 2P photobleaching power scan v20200728.mat')

%% Regroup the table
for i = 1:height(T_fov)
    Configuration = split(T_fov.Path{i},{'\','_','.','PW','HV'});
    T_fov.Name{i} = Configuration{10};
    T_fov.Power(i) = str2num(Configuration{12});
    T_fov.HV(i) = str2num(Configuration{14});
    T_fov.PlateWell{i} = Configuration{15};
    T_fov.Settings{i} = strcat('Power: ',Configuration{12},' HV: ',Configuration{14});
    T_fov.AOCNorm(i) = T_fov.TraceNorm(i).integral(0, max(T_fov.TraceNorm(i).x));
    T_fov.AOCFinal(i) = T_fov.TraceFinal(i).integral(0, max(T_fov.TraceFinal(i).x));
    T_fov.BleachRatio(i) = T_fov.TraceNorm(i).y(1)-T_fov.TraceNorm(i).y(end);
    T_fov.RemainedRatio(i) = T_fov.TraceNorm(i).y(end);
end

conditionGroups = findgroups(strcat(T_fov.Name, T_fov.Settings));
T_fovCondition = table();
T_fovCondition.Name = splitapply(@unique, T_fov.Name, conditionGroups);
T_fovCondition.Area = splitapply(@mean, T_fov.Area, conditionGroups);
T_fovCondition.Settings = splitapply(@unique, T_fov.Settings, conditionGroups);
T_fovCondition.TraceNorm = splitapply(@mean, T_fov.TraceNorm, T_fov.Area, conditionGroups);
T_fovCondition.TraceFinal = splitapply(@mean, T_fov.TraceFinal, T_fov.Area, conditionGroups);
T_fovCondition.AOCNorm =  splitapply(@mean, T_fov.AOCNorm, conditionGroups);
T_fovCondition.AOCFinal = T_fovCondition.TraceFinal.integral(0,max(T_fovCondition.TraceFinal(1).x));
T_fovCondition.BleachRatio = splitapply(@mean,  T_fov.BleachRatio, conditionGroups);
T_fovCondition.RemainedRatio = splitapply(@mean,  T_fov.RemainedRatio, conditionGroups);

%% Plot per variant
c = hsv(10);
wellGroups = findgroups(T_fovCondition.Name);
nVariants = max(wellGroups);
figure(1)
clf
for i = 1:max(wellGroups)
    T_sub = T_fovCondition(find(wellGroups ==i),:);
    subplot(4,3,i)
    hold on
    for j = 1:height(T_sub)
        plot(T_sub.TraceNorm(j).x,T_sub.TraceNorm(j).y, 'color',c(j,:));
    end
    if i == max(wellGroups)
        legend(T_sub.Settings)
    end
    xlim([0 120])
    title(T_sub.Name{1})
    xlabel('t [s]')
    ylabel('Normalized fluorescence')
    box on
end

%% Plot per condition
c = hsv(12);
wellGroups = findgroups(T_fovCondition.Settings);
nVariants = max(wellGroups);
figure(2)
clf
for i = 1:max(wellGroups)
    T_sub = T_fovCondition(find(wellGroups ==i),:);
    subplot(3,3,i)
    hold on
    for j = 1:height(T_sub)
        plot(T_sub.TraceNorm(j).x,T_sub.TraceNorm(j).y, 'color',c(j,:));
    end
    if i == max(wellGroups)
        legend(T_sub.Name)
    end
    xlim([0 120])
    title(T_sub.Settings{1})
    xlabel('t [s]')
    ylabel('Normalized fluorescence')
    box on
end
figure(3)
clf
for i = 1:max(wellGroups)
    T_sub = T_fovCondition(find(wellGroups ==i),:);
    subplot(3,3,i)
    hold on
    for j = 1:height(T_sub)
        plot(T_sub.TraceFinal(j).x,T_sub.TraceFinal(j).y,'color',c(j,:));
    end
    if i == max(wellGroups)
        legend(T_sub.Name)
    end
    xlim([0 120])
    title(T_sub.Settings{1})
    xlabel('t [s]')
    ylabel('fluorescence [a.u.]')
    box on
end
%% Bar plot ranking
T_sub = T_fovCondition(find(strcmp('Power:20 HV:20',T_fovCondition.Settings)),:);
libplot.formattedBar(T_sub, 'AOCNorm','AOCFinal','sort_on', 'AOCNorm', 'Name','Name','ylabel', ...
   'AUC normalized to ASAP1','norm2ctrl','ASAP1')
libplot.formattedBar(T_sub, 'RemainedRatio','sort_on', 'RemainedRatio','Name','Name','ylabel', ...
   '1-\DeltaF/F0','norm2ctrl','ASAP1')

%% Correlation plot re-arrange
T_sub = T_fovCondition(find(strcmp('Power:20 HV:20',T_fovCondition.Settings)),:);
for i = 1:height(T_sub)
    T_sub.LegalName(i) = string(libplot.fzmatch(T_sub.Name{i}));
    T_sub.Properties.RowNames(i) = T_sub.LegalName(i);
end
for i = 1:height(T_well_stats)
    T_well_stats.LegalName(i) = string(libplot.fzmatch(T_well_stats.Name{i}));
end
T_sub = sortrows(T_sub, 'LegalName', 'ascend');
T_well_stats.Properties.RowNames = T_well_stats.LegalName;
T_well_stats_sub =  T_well_stats(T_sub.LegalName,:)
T_well_stats_sub = sortrows(T_well_stats_sub, 'LegalName', 'ascend');

%% Correlation plot
figure(4)
clf
hold on
pp1 = 'BleachRatio';
pp2 = 'Bleachingratio_totalMean';
for i = 1:height(T_sub)
scatter(T_sub.(pp1)(i),T_well_stats_sub.(pp2)(i),[],c(i,:),'filled')
end
legend(T_sub.LegalName)
[k, R2, pVal] = libplot.lnrFitting(T_sub.(pp1),T_well_stats_sub.(pp2),false,true,4)

%%
for i = 1:height(T_well)
end