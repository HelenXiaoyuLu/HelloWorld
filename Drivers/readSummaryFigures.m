%% 20190315 1P Fstim
clear
clc
T_well_1PFstim = readtable('D:\OneDrive - rice.edu\Paper\Data\Screenings\SummarizeFigures.xlsx',...
    'Sheet','190315 1P Fstim','Range','A1:CH37','ReadVariableNames',true);
T_well_1PFstim_508 = readtable('D:\OneDrive - rice.edu\Paper\Data\Screenings\SummarizeFigures.xlsx',...
    'Sheet','190315 1P Fstim','Range','A55:CS67','ReadVariableNames',true);

%% re-arrange
for i = 1:height(T_well_1PFstim)
    variantName = split(T_well_1PFstim.Name{i},'(P');
    T_well_1PFstim.StandardName(i) = libplot.fzmatch(variantName{1}); %Standardize the name
    T_well_1PFstim.('dFF0 1.5ms')(i) = mean([T_well_1PFstim.('dFF0OfMean_C1S1_')(i),...
        T_well_1PFstim.('dFF0OfMean_C1S2_')(i), T_well_1PFstim.('dFF0OfMean_C1S3_')(i)]);
    T_well_1PFstim.('dFF0 2.5ms')(i) = mean([T_well_1PFstim.('dFF0OfMean_C1S4_')(i),...
        T_well_1PFstim.('dFF0OfMean_C1S5_')(i), T_well_1PFstim.('dFF0OfMean_C1S6_')(i)]);    
end
for i = 1:height(T_well_1PFstim_508)
    variantName = split(T_well_1PFstim_508.Name{i},'(P');
    T_well_1PFstim_508.StandardName(i) = libplot.fzmatch(variantName{1}); %Standardize the name
    T_well_1PFstim_508.('dFF0 1.5ms')(i) = mean([T_well_1PFstim_508.('dFF0OfMean_C1S1_')(i),...
        T_well_1PFstim_508.('dFF0OfMean_C1S2_')(i), T_well_1PFstim_508.('dFF0OfMean_C1S3_')(i)]);
    T_well_1PFstim_508.('dFF0 2.5ms')(i) = mean([T_well_1PFstim_508.('dFF0OfMean_C1S4_')(i),...
        T_well_1PFstim_508.('dFF0OfMean_C1S5_')(i), T_well_1PFstim_508.('dFF0OfMean_C1S6_')(i)]);    
end
T_well_1PFstim = [T_well_1PFstim; T_well_1PFstim_508(:,T_well_1PFstim.Properties.VariableNames)];

wellgroups = findgroups(T_well_1PFstim.StandardName);
T_group_1PFstim = table();
T_group_1PFstim.Name = splitapply(@unique, T_well_1PFstim.StandardName, wel lgroups);
T_group_1PFstim.('Mean Relative Brightness (C2/C3)') = splitapply(@mean, T_well_1PFstim.('MeanRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('Median Relative Brightness (C2/C3)') = splitapply(@mean, T_well_1PFstim.('MedianRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('Mode Relative Brightness (C2/C3)') = splitapply(@mean, T_well_1PFstim.('ModeRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('dFF0 1.5ms') = splitapply(@mean, T_well_1PFstim.('dFF0 1.5ms'), wellgroups);
T_group_1PFstim.('dFF0 2.5ms') = splitapply(@mean, T_well_1PFstim.('dFF0 2.5ms'), wellgroups);
T_group_1PFstim.('dFF0 200Hz') = splitapply(@mean, T_well_1PFstim.('dFF0OfMean_C1S7_'), wellgroups);
T_group_1PFstim.('dFF0 100Hz') = splitapply(@mean, T_well_1PFstim.('dFF0OfMean_C1S8_'), wellgroups);
T_group_1PFstim.('dFF0 20ms') = splitapply(@mean, T_well_1PFstim.('dFF0OfMean_C1S9_'), wellgroups);

T_group_1PFstim.('Mean Relative Brightness (C2/C3) STD') = splitapply(@std, T_well_1PFstim.('MeanRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('Median Relative Brightness (C2/C3) STD') = splitapply(@std, T_well_1PFstim.('MedianRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('Mode Relative Brightness (C2/C3) STD') = splitapply(@std, T_well_1PFstim.('ModeRelativeBrightness_C2_C3_'), wellgroups);
T_group_1PFstim.('dFF0 1.5ms STD') = splitapply(@std, T_well_1PFstim.('dFF0 1.5ms'), wellgroups);
T_group_1PFstim.('dFF0 2.5ms STD') = splitapply(@std, T_well_1PFstim.('dFF0 2.5ms'), wellgroups);
T_group_1PFstim.('dFF0 200Hz STD') = splitapply(@std, T_well_1PFstim.('dFF0OfMean_C1S7_'), wellgroups);
T_group_1PFstim.('dFF0 100Hz STD') = splitapply(@std, T_well_1PFstim.('dFF0OfMean_C1S8_'), wellgroups);
T_group_1PFstim.('dFF0 20ms STD') = splitapply(@std, T_well_1PFstim.('dFF0OfMean_C1S9_'), wellgroups);

save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190315 benchmarking\20190315 1P Fstim 1000Hz.mat',...
    'T_group_1PFstim','T_well_1PFstim','T_well_1PFstim_508');
c = lines;
figure(1)
clf
hold on
errorbar(T_group_1PFstim.('Median Relative Brightness (C2/C3)'),...
    T_group_1PFstim.('Mean Relative Brightness (C2/C3)'),...
    T_group_1PFstim.('Mean Relative Brightness (C2/C3) STD'),...
    T_group_1PFstim.('Mean Relative Brightness (C2/C3) STD'),...
    T_group_1PFstim.('Median Relative Brightness (C2/C3) STD'),...
    T_group_1PFstim.('Median Relative Brightness (C2/C3) STD'),'.')
scatter(T_group_1PFstim.('Median Relative Brightness (C2/C3)')',...
    T_group_1PFstim.('Mean Relative Brightness (C2/C3)')',...
    ones(height(T_group_1PFstim),1)*20,c(1:height(T_group_1PFstim),:),'o','filled')
for i = 1:height()
legend(S,T_group_1PFstim.Name)

%% 20190315 1P photobleaching 
clear
clc
bleaching1PMean = readtable('D:\OneDrive - rice.edu\Paper\Data\Screenings\SummarizeFigures.xlsx',...
    'Sheet','190315 1P Photobleaching','Range','AH1:CC612','ReadVariableNames',true);
bleaching1PMedian = readtable('D:\OneDrive - rice.edu\Paper\Data\Screenings\SummarizeFigures.xlsx',...
    'Sheet','190315 1P Photobleaching','Range','CD1:DY612','ReadVariableNames',true);
time = readtable('D:\OneDrive - rice.edu\Paper\Data\Screenings\SummarizeFigures.xlsx',...
    'Sheet','190315 1P Photobleaching','Range','DZ1:FU612','ReadVariableNames',true);
deltat = 0.5; % For synchronization

% re-arrange
T_well_1Pbleaching = table();
for i = 1:width(bleaching1PMean)
    variantName = bleaching1PMean.Properties.VariableNames{i};
    T_well_1Pbleaching.Name(i) = libplot.fzmatch(variantName); %Standardize the name
    s = sig.arbitrary();
    s.setSignal(time.(variantName)(:)', bleaching1PMean{:,i}');
    T_well_1Pbleaching.traceMean(i, 1) = s;
    s2 = sig.arbitrary();
    s2.setSignal(time.(variantName)(:)', bleaching1PMedian.(variantName)(:)');
    T_well_1Pbleaching.traceMedian(i, 1) = s2;
end
T_well_1Pbleaching.traceMeanFinal = T_well_1Pbleaching.traceMean.synchronize('dt', deltat, 'method', 'intersect');
T_well_1Pbleaching.traceMedianFinal = T_well_1Pbleaching.traceMedian.synchronize('dt', deltat, 'method', 'intersect');

wellgroups = findgroups(T_well_1Pbleaching.Name);
T_group_1Pbleaching = table();
T_group_1Pbleaching.Name = splitapply(@unique, T_well_1Pbleaching.Name, wellgroups);
T_group_1Pbleaching.meanofMean = splitapply(@mean, T_well_1Pbleaching.traceMeanFinal, wellgroups);
T_group_1Pbleaching.stdofMean = splitapply(@std, T_well_1Pbleaching.traceMeanFinal, wellgroups);
T_group_1Pbleaching.meanofMedian = splitapply(@mean, T_well_1Pbleaching.traceMedianFinal, wellgroups);
T_group_1Pbleaching.stdofMedian = splitapply(@std, T_well_1Pbleaching.traceMedianFinal, wellgroups);

save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190315 benchmarking\1P photobleaching traces summary.mat',...
    'T_group_1Pbleaching','T_well_1Pbleaching');

% Plot figures
figure(1)
clf
hold on
c = lines;
for i = 1: height(T_group_1Pbleaching)
    t = T_group_1Pbleaching.meanofMean(i).x';
    meany = T_group_1Pbleaching.meanofMean(i).y';
    stdy = T_group_1Pbleaching.stdofMean(i).y';
    patch = fill([t,fliplr(t)], [meany - stdy, fliplr(meany+stdy)], c(i,:),...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');    
    h(i) = plot(T_group_1Pbleaching.meanofMean(i).x, T_group_1Pbleaching.meanofMean(i).y,'color', c(i,:));
end
legend(h,T_group_1Pbleaching.Name)
xlabel('t [s]')
ylabel('Normalized mean brightness')

figure(2)
clf
hold on
c = lines;
for i = 1: height(T_group_1Pbleaching)
    t = T_group_1Pbleaching.meanofMedian(i).x';
    meany = T_group_1Pbleaching.meanofMedian(i).y';
    stdy = T_group_1Pbleaching.stdofMedian(i).y';
    patch = fill([t,fliplr(t)], [meany - stdy, fliplr(meany+stdy)], c(i,:),...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');    
    h(i) = plot(T_group_1Pbleaching.meanofMedian(i).x, T_group_1Pbleaching.meanofMedian(i).y,'color', c(i,:));
end
legend(h,T_group_1Pbleaching.Name)
xlabel('t [s]')
ylabel('Normalized median brightness')