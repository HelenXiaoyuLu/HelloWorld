%% plot ephys result (archive) & correlation plot of ephys (1P) vs Fstim (2P)
% Xiaoyu Lu, xiaoyu.lu@rice.edu 
% St-Pierre Lab, Aug. 2020

%% Load data 
clear
dirp = '%TEST%\plot_ephys_vs_2pFstim';
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Ephys\plot_ephys_vs_2pFstim';
load(fullfile(dirp,'ThreeTest_Single_AP_wASAP1'));
load(fullfile(dirp,'ThreeTest_Train_AP_wASAP1'));
load(fullfile(dirp,'ThreeTest_Transition_curve_wASAP1'));
load(fullfile(dirp,'AP_4msFWHM_100mV_10kHz_TRIMMED'));
load(fullfile(dirp,'20200126 2pFstim plate1'));
T_group = T_out.T_group;
T_well = T_out.T_well;
CtrlFP = 'cyOFP';

%% Single Action Potential Trace Plot 
figure(1) % Show average trace from all
clf
for i = 1:height(Single_AP)-1
    subplot(2,3,i)
    hold on
    trace_t = cell2mat(Single_AP(i,7).Variables);
    trace_mean = cell2mat(Single_AP(i,3).Variables);
    trace_std = cell2mat(Single_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,'k','LineStyle','none','FaceAlpha','0.2');
    plot(trace_t, trace_mean,'Color',[0.39,0.83,0.07]);
    title(Single_AP(i,1).Variables)
    axis([-0.05 0.05 0.65 1.05])
    axis square
    xlabel('Time/s')
    ylabel('Normalized Fluorescence')
end

figure(2) % Show stimulation and selected traces
clf
subplot(2,1,1) % Command waveform
T_exd = [0:0.01:0.05];
Baseline_exd = c002_IN_0(1)*ones(size(T_exd));
T_full_singleAP = [T_exd, T_exd(end)+c001_Time];
Stim_full_singleAP = [Baseline_exd, c002_IN_0];
plot(T_full_singleAP-0.05,1000*Stim_full_singleAP,'Color',[0.20,0.46,0.02])
axis([-0.05 0.05 -80 60])
axis square
xlabel('Time (s)')
ylabel('Command Voltage (mV)')
set(gcf, 'Position',  [10, 10, 400, 800])

subplot(2,1,2) % Show selected traces
hold on
ColorStore = [0.00,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
for i = 1:3
    trace_t = cell2mat(Single_AP(i,7).Variables);
    trace_mean = 1-cell2mat(Single_AP(i,3).Variables);
    trace_std = cell2mat(Single_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,ColorStore(i,:),'LineStyle','none','FaceAlpha','0.2');
    set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h(i)= plot(trace_t, trace_mean,'Color',ColorStore(i,:),'LineWidth',1);
end
axis([-0.05 0.05 -0.05 0.3])
axis square
xlabel('Time (s)')
ylabel('Normalized Fluorescence')
legend(Single_AP.Name(1:3),'location','best')

%% Train Artificial Action Potential Trace Plot 
figure(3) % Show average trace from all
clf
for i = 1:height(Train_AP)-1
    subplot(2,3,i)
    hold on
    trace_t = cell2mat(Train_AP(i,7).Variables);
    trace_mean = cell2mat(Train_AP(i,3).Variables);
    trace_std = cell2mat(Train_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,'k','LineStyle','none','FaceAlpha','0.2');
    plot(trace_t, trace_mean,'Color',[0.39,0.83,0.07])
    title(Train_AP(i,1).Variables)
    axis([-0.05 0.15 0.6 1.05])
    axis square
    xlabel('Time/s')
    ylabel('Normalized Fluorescence')
end

figure(4) % Show stimulation waveform and selected traces
clf
subplot(2,1,1) % Command voltage
T_exd = [0:0.01:0.05];
Baseline_exd = c002_IN_0(1)*ones(size(T_exd));
T_full_trainAP = T_exd;
Stim_full_trainAP = Baseline_exd;
for i = 1:10
    T_full_trainAP = [T_full_trainAP, T_full_trainAP(end)+c001_Time(1:100)];
    Stim_full_trainAP = [Stim_full_trainAP, c002_IN_0(1:100)];
end
T_full_trainAP = [T_full_trainAP,T_full_trainAP(end)+T_exd];
Stim_full_trainAP = [Stim_full_trainAP,Baseline_exd];
plot(T_full_trainAP-0.05,1000*Stim_full_trainAP,'Color',[0.20,0.46,0.02])
axis([-0.05 0.15 -80 60])
axis square
xlabel('Time (s)')
ylabel('Command Voltage (mV)')
set(gcf, 'Position',  [510, 10, 400, 800])

subplot(2,1,2) % selected traces
hold on
ColorStore = [0.00,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
for i = 1:3
    trace_t = cell2mat(Train_AP(i,7).Variables);
    trace_mean = 1-cell2mat(Train_AP(i,3).Variables);
    trace_std = cell2mat(Train_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h1 = fill([trace_t, fliplr(trace_t)],Shade,ColorStore(i,:),'LineStyle','none','FaceAlpha','0.2');
    set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h(i) = plot(trace_t, trace_mean,'Color',ColorStore(i,:),'LineWidth',1);   
end
axis([-0.05 0.15 -0.05 0.3])
axis square
xlabel('Time (s)')
ylabel('Normalized Fluorescence')
legend(Train_AP.Name(1:3),'location','best')
uistack(h(2),'top')
uistack(h(1),'top')

%% Regroup Fstim Result 
% Selected variants to be included in the fitting curve
sel = ["JEDI-1P","JEDI-2P","ASAP1","ASAP2s","ASAP3","ASAP2s-H152E"];

% Create a temporary table to pull well2group stats for dF/F0
weightedMean = @(x, w) sum(x.*w, 'all')./sum(w, 'all');
grouping = findgroups(T_well.Group);
T_group_temp = table();
T_group_temp.Group = splitapply(@unique, T_well.Group, grouping);
T_group_temp.("dFF0 Short Stim") = splitapply(weightedMean, T_well.("dFF0 Short Stim"), T_well.Area, grouping);
T_group_temp.("dFF0 Short Stim STD") = splitapply(@std, T_well.("dFF0 Short Stim"), T_well.Area, grouping);
T_group_temp.("dFF0 Long Stim") = splitapply(weightedMean, T_well.("dFF0 Long Stim"), T_well.Area, grouping);
T_group_temp.("dFF0 Long Stim STD") = splitapply(@std, T_well.("dFF0 Long Stim"), T_well.Area, grouping);
T_group = innerjoin(T_group_temp, T_group, 'Keys', 'Group');
% Label each row
for i = 1:height(T_group)
    groupname = split(T_group.Name(i),['-',CtrlFP]);
    T_group.Properties.RowNames(i) = groupname(1);
end

%% Compare Ephys with Fstim (Voltage Clamp V.S. Long Stim)
f = figure(6); 
clf
ax = axes(f);
hold (ax,'on')
EphysdFF0Mean = zeros(1,numel(sel));
EphysdFF0Std = zeros(1,numel(sel));
FstimdFF0Mean = zeros(1,numel(sel));
FstimdFF0Std = zeros(1,numel(sel));
for i = 1:numel(sel)  
    EphysdFF0Mean(i) = Transition_curve.V30Mean(sel(i));
    EphysdFF0Std(i) = Transition_curve.V30Std(sel(i));
    FstimdFF0Mean(i) = T_group.("dFF0 Long Stim")(sel(i));
    FstimdFF0Std(i) = T_group.("dFF0 Long Stim STD")(sel(i));
    errorbar(-EphysdFF0Mean(i),-FstimdFF0Mean(i),FstimdFF0Std(i),...
        FstimdFF0Std(i),EphysdFF0Std(i),EphysdFF0Std(i),'.','MarkerSize',30)
end
x = -EphysdFF0Mean;
y = -FstimdFF0Mean;
box on
[k, R2] = libplot.lnrFitting(x,y,false,true,ax,'polyfit');
xlabel({'Ephys VClamp = +30mV','-\DeltaF/F0'})
ylabel({'Fstim = 100ms@100Hz','-\DeltaF/F0'})
axis([0 inf 0 inf])
legend(sel,'location','best')
axis square
title(['Ephys vs Fstim (',CtrlFP,')'])
set(gcf, 'Position',  [50, 10, 500, 500])

%% Compare Ephys with Fstim (Train AP V.S. Long Stim)
f = figure(7);
clf
ax = axes(f);
hold (ax,'on')
EphysdFF0Mean = zeros(1,numel(sel));
EphysdFF0Std = zeros(1,numel(sel));
FstimdFF0Mean = zeros(1,numel(sel));
FstimdFF0Std = zeros(1,numel(sel));
for i = 1:numel(sel)  
    EphysdFF0Mean(i) = Train_AP.SpikeMean(sel(i));
    EphysdFF0Std(i) = Train_AP.SpikeStd(sel(i));
    FstimdFF0Mean(i) = T_group.("dFF0 Long Stim")(sel(i));
    FstimdFF0Std(i) = T_group.("dFF0 Long Stim STD")(sel(i));
    errorbar(-EphysdFF0Mean(i),-FstimdFF0Mean(i),FstimdFF0Std(i),...
        FstimdFF0Std(i),EphysdFF0Std(i),EphysdFF0Std(i),'.','MarkerSize',30)
end
x = -EphysdFF0Mean;
y = -FstimdFF0Mean;
box on
[k, R2] = libplot.lnrFitting(x,y,false,true,ax,'polyfit');
xlabel({'Ephys Train AP','-\DeltaF/F0'})
ylabel({'Fstim = 100ms@100Hz','-\DeltaF/F0'})
legend(sel,'location','best')
axis square
title(['Ephys vs Fstim (',CtrlFP,')'])
set(gcf, 'Position',  [50, 10, 500, 500])

%% Compare Ephys with Fstim (Single AP V.S. Short Stim)
f = figure(8);
clf
ax = axes(f);
hold (ax,'on')
EphysdFF0Mean = zeros(1,numel(sel));
EphysdFF0Std = zeros(1,numel(sel));
FstimdFF0Mean = zeros(1,numel(sel));
FstimdFF0Std = zeros(1,numel(sel));
for i = 1:numel(sel)  
    EphysdFF0Mean(i) = Single_AP.SpikeMean(sel(i));
    EphysdFF0Std(i) = Single_AP.SpikeStd(sel(i));
    FstimdFF0Mean(i) = T_group.("dFF0 Short Stim")(sel(i));
    FstimdFF0Std(i) = T_group.("dFF0 Short Stim STD")(sel(i));
    errorbar(-EphysdFF0Mean(i),-FstimdFF0Mean(i),FstimdFF0Std(i),...
        FstimdFF0Std(i),EphysdFF0Std(i),EphysdFF0Std(i),'.','MarkerSize',30)
end
x = -EphysdFF0Mean;
y = -FstimdFF0Mean;
box on
[k, R2] = libplot.lnrFitting(x,y,false,true,ax,'polyfit');
xlabel({'Ephys Single AP','-\DeltaF/F0'})
ylabel({'Fstim = 1ms','-\DeltaF/F0'})
axis([0 inf 0 inf])
legend(sel,'location','best')
axis square
title(['Ephys vs Fstim (',CtrlFP,')'])
set(gcf, 'Position',  [50, 10, 500, 500])