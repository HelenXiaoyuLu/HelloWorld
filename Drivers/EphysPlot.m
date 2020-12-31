% Ephys plot
DataPath = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\Paper_GEVI_summary_figures_bcm.xlsx';

%% Transition Curve data processing
Name = {'ASAP2s','ASAP3','JEDI-2P','JEDI-1P',...
    'ASAP2s-H152E','ASAP2s-H152E-Q397H'};
ASAP2s_TC_raw = readmatrix(DataPath,'Sheet','4','Range','C3:F10');
JEDI1P_TC_raw = readmatrix(DataPath,'Sheet','8','Range','C3:F10');
JEDI2P_TC_raw = readmatrix(DataPath,'Sheet','12','Range','C3:F8');
ASAP3_TC_raw = readmatrix(DataPath,'Sheet','34','Range','C3:F9');
GV3_TC_raw = readmatrix(DataPath,'Sheet','38','Range','C3:F7');
GV75_TC_raw = readmatrix(DataPath,'Sheet','42','Range','C3:F10');
RawData = {ASAP2s_TC_raw,ASAP3_TC_raw,JEDI2P_TC_raw,JEDI1P_TC_raw,...
    GV3_TC_raw,GV75_TC_raw};
for j = 1:length(RawData) 
    RawMat = RawData{i};
    V_100Mean{j} = mean(RawMat(:,1));
    V_100Std{j} = std(RawMat(:,1)); 
    V_70Mean{j} = mean(RawMat(:,2));
    V_70Std{j} = std(RawMat(:,2));     
    V_40Mean{j} = mean(RawMat(:,3));
    V_40Std{j} = std(RawMat(:,3)); 
    V30Mean{j} = mean(RawMat(:,4));
    V30Std{j} = std(RawMat(:,4));
end
Transition_curve_raw = struct('Name',Name,'RawData',RawData,'V_100Mean',V_100Mean,...
    'V_100Std',V_100Std,'V_70Mean',V_70Mean,'V_70Std',V_70Std,'V_40Mean',V_40Mean,'V_40Std',V_40Std,...
    'V30Mean',V30Mean,'V30Std',V30Std);
Transition_curve = struct2table(Transition_curve_raw);

%% Transition Curve data plot (Figure 1)
figure(1) % Compare all variants
hold on
% Transition_curve_T = cell2table(table2cell(Transition_curve)','RowNames',...
%     Transition_curve.Properties.VariableNames,...
%     'VariableNames',Transition_curve.Name)
Voltages = [-100,-70,-40,30];
for i = 1:height(Transition_curve)
    for j = 3:2:width(Transition_curve)
        dFF0((j-1)/2) = Transition_curve(i,j).Variables;
        err((j-1)/2) = Transition_curve(i,j+1).Variables;
    end
    errorbar(Voltages,dFF0,err,'.-')
end
legend(Transition_curve.Name,'location','best')
xlabel('Command Voltage (mV)')
ylabel({'Steady-state response','\DeltaF/F0'})
title('Ephys 3-test Voltage Clamp')
% save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\ThreeTest_Transition_curve.mat','Transition_curve');
% h = get(gca,'Children')
% h(1).delete()

%% Single Artificial Action Potential
ASAP2s_AP_raw = readmatrix(DataPath,'Sheet','5','Range','D2:DL17');
JEDI1P_AP_raw = readmatrix(DataPath,'Sheet','9','Range','D2:DL17');
JEDI2P_AP_raw = readmatrix(DataPath,'Sheet','13','Range','D2:DL13');
ASAP3_AP_raw = readmatrix(DataPath,'Sheet','35','Range','D2:DL15');
GV3_AP_raw = readmatrix(DataPath,'Sheet','39','Range','D2:DL11');
GV75_AP_raw = readmatrix(DataPath,'Sheet','43','Range','D2:DL17');
Raw_AP_Data = {ASAP2s_AP_raw,ASAP3_AP_raw,JEDI2P_AP_raw,JEDI1P_AP_raw,...
    GV3_AP_raw,GV75_AP_raw};
for j = 1:length(Raw_AP_Data) 
    RawMat = Raw_AP_Data{j};
    TraceMean{j} = rmmissing(mean(RawMat(2:2:end,:),1));
    TraceStd{j} = rmmissing(std(RawMat(2:2:end,:),1)); 
    SpikeMean{j} = mean(min(RawMat(2:2:end,:),[],2));
    SpikeStd{j} = std(min(RawMat(2:2:end,:),[],2));     
    Time{j} = RawMat(1,1:length(TraceMean{j}));
end
Single_AP_raw = struct('Name',Name,'RawData',Raw_AP_Data,'TraceMean',TraceMean,...
    'TraceStd', TraceStd, 'SpikeMean', SpikeMean,'SpikeStd',SpikeStd,'Time',Time);
Single_AP = struct2table(Single_AP_raw);
% save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\ThreeTest_Single_AP.mat','Single_AP');

%% Single Artificial Action Potential Plot(Figure 2-3)
sel = ["ASAP2s", "ASAP3", "JEDI-1P"];
figure(2) % Show average trace from all
clf
for i = 1:height(Single_AP)
    subplot(2,6,i)
    hold on
    trace_t = cell2mat(Single_AP(i,7).Variables);
    trace_mean = cell2mat(Single_AP(i,3).Variables);
    trace_std = cell2mat(Single_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,'k','LineStyle','none','FaceAlpha','0.2');
    plot(trace_t, trace_mean,'Color',[0.39,0.83,0.07])
    title(Single_AP(i,1).Variables)
    axis([-0.05 0.05 0.65 1.05])
    axis square
    xlabel('Time/s')
    ylabel('Normalized Fluorescence')
end

figure(3)
clf
subplot(2,1,2)
hold on
ColorStore = [0.00,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
for i = 1:3
    trace_t = cell2mat(Single_AP.Time(sel(i)));
    trace_mean = 1-cell2mat(Single_AP.TraceMean(sel(i)));
    trace_std = cell2mat(Single_AP.TraceStd(sel(i)));
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,ColorStore(i,:),'LineStyle','none','FaceAlpha','0.2');
    set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h(i) = plot(trace_t, trace_mean,'Color',ColorStore(i,:),'LineWidth',1)
end
axis([-0.05 0.05 -0.05 0.4])
axis square
xlabel('Time (s)')
ylabel('Normalized Fluorescence')
legend(sel,'location','best')
%uistack(h(2),'top')
%uistack(h(1),'top')

subplot(2,1,1)
load('D:\OneDrive - Rice University\Francois\Paper\JEDI-2P\Figures\Ephys\AP_4msFWHM_100mV_10kHz_TRIMMED.mat');
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

%% Single Artificial Action Potential Trace plot optional (Figure ())
figure()
for i = 1:size(Single_AP,1)
    subplot(2,6,i+6)
    hold on
    trace_t = cell2mat(Single_AP(i,7).Variables);
    trace_mean = cell2mat(Single_AP(i,3).Variables);
    trace_std = cell2mat(Single_AP(i,4).Variables);
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h = fill([trace_t, fliplr(trace_t)],Shade,'k','LineStyle','none','FaceAlpha','0.2');
    plot(trace_t, trace_mean,'Color',[0.39,0.83,0.07])
    title(Single_AP(i,1).Variables)
    axis([-0.05 0.05 0.65 1.05])
    axis square
    xlabel('Time/s')
    ylabel('Normalized Fluorescence')
end

%% Artificial Action Potential 100Hz
ASAP2s_trainAP_raw = readmatrix(DataPath,'Sheet','6','Range','D2:SE17');
JEDI1P_trainAP_raw = readmatrix(DataPath,'Sheet','10','Range','D2:SE17');
JEDI2P_trainAP_raw = readmatrix(DataPath,'Sheet','14','Range','D2:SE13');
ASAP3_trainAP_raw = readmatrix(DataPath,'Sheet','36','Range','D2:SE15');
GV3_trainAP_raw = readmatrix(DataPath,'Sheet','40','Range','D2:SE11');
GV75_trainAP_raw = readmatrix(DataPath,'Sheet','44','Range','D2:SE17');
Raw_trainAP_Data = {ASAP2s_trainAP_raw,ASAP3_trainAP_raw,JEDI2P_trainAP_raw,...
    JEDI1P_trainAP_raw,GV3_trainAP_raw,GV75_trainAP_raw};
for j = 1:length(Raw_trainAP_Data) 
    RawMat = Raw_trainAP_Data{j};
    TraceMean{j} = rmmissing(mean(RawMat(2:2:end,:),1));
    TraceStd{j} = rmmissing(std(RawMat(2:2:end,:),1)); 
    SpikeMean{j} = mean(min(RawMat(2:2:end,:),[],2));
    SpikeStd{j} = std(min(RawMat(2:2:end,:),[],2));     
    Time{j} = RawMat(1,1:length(TraceMean{j}));
end
Train_AP_raw = struct('Name',Name,'RawData',Raw_trainAP_Data,'TraceMean',TraceMean,...
    'TraceStd', TraceStd, 'SpikeMean', SpikeMean,'SpikeStd',SpikeStd,'Time',Time);
Train_AP = struct2table(Train_AP_raw);
save('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\ThreeTest_Train_AP.mat','Train_AP');

%% Train Artificial Action Potential Plot (Figure 4-5)
figure(4) % Show average trace from all
clf
for i = 1:size(Train_AP,1)
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

figure(5)
clf
subplot(2,1,2)
hold on
ColorStore = [0.00,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
for i = 1:3
    trace_t = cell2mat(Train_AP.Time(sel(i)));
    trace_mean = 1-cell2mat(Train_AP.TraceMean(sel(i)));
    trace_std = cell2mat(Train_AP.TraceStd(sel(i)));
    Shade = [trace_mean+trace_std, fliplr(trace_mean-trace_std)];
    h1 = fill([trace_t, fliplr(trace_t)],Shade,ColorStore(i,:),'LineStyle','none','FaceAlpha','0.2');
    set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h(i) = plot(trace_t, trace_mean,'Color',ColorStore(i,:),'LineWidth',1)   
end
axis([-0.05 0.15 -0.05 0.4])
axis square
xlabel('Time (s)')
ylabel('Normalized Fluorescence')
legend(sel,'location','best')
uistack(h(2),'top')
uistack(h(1),'top')

subplot(2,1,1)
 load('D:\OneDrive - Rice University\Francois\Paper\JEDI-2P\Figures\Ephys\AP_4msFWHM_100mV_10kHz_TRIMMED.mat');
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

%% Compare Ephys with Fstim result
% Ephys2Fstim = [4,6,12,11,5,1];
% Ephys2Fstim = [6,7,12,16,4,5];
Ephys2Fstim = [4,7,13,12,5,1];
sel = [1:5,7];
LongStim = 'Stim5';
CtrlFP = 'cyOFP';
figure(6) % Compare Ephys with Fstim
clf
hold on

EphysdFF0Mean = [];
EphysdFF0Std = [];
FstimLongdFF0Mean = [];
FstimLongdFF0Std = [];
for i = 1:length(sel) %height(Transition_curve)   
    EphysdFF0Mean(i) = Transition_curve(sel(i),9).Variables;
    EphysdFF0Std(i) = Transition_curve(sel(i),10).Variables;
    FstimLongdFF0Mean(i) = T_well_stats(Ephys2Fstim(i),strcat(LongStim,'Mean')).Variables;
    FstimLongdFF0Std(i) = T_well_stats(Ephys2Fstim(i),strcat(LongStim,'Std')).Variables;
    errorbar(-EphysdFF0Mean(i),-FstimLongdFF0Mean(i),FstimLongdFF0Std(i),FstimLongdFF0Std(i),EphysdFF0Std(i),EphysdFF0Std(i),'.','MarkerSize',30)
end
x = -EphysdFF0Mean;
y = -FstimLongdFF0Mean;
[k, R2] = libplot.lnrFitting(x,y,false,true,6);
xlabel({'Ephys VClamp = +30mV','-\DeltaF/F0'})
ylabel({'Fstim = 100ms@100Hz','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis([0 inf 0 inf])
axis square
legend(Transition_curve.Name(sel),'location','best')
set(gcf, 'Position',  [50, 10, 500, 500])

figure(7)
clf
hold on
EphysAPdFF0Mean = [];
EphysAPdFF0Std = [];
FstimShortdFF0Mean = [];
FstimShortdFF0Std = [];
for i = 1:length(sel) %1:size(Single_AP,1)
    EphysAPdFF0Mean(i) = 1-Single_AP(sel(i),5).Variables;
    EphysAPdFF0Std(i) = Single_AP(sel(i),6).Variables;
    FstimShortdFF0Mean(i) = T_well_stats(Ephys2Fstim(i),'StimShortMean').Variables;
    FstimShortdFF0Std(i) = T_well_stats(Ephys2Fstim(i),'StimShortStd').Variables;
    errorbar(EphysAPdFF0Mean(i),-FstimShortdFF0Mean(i),FstimShortdFF0Std(i),FstimShortdFF0Std(i),EphysAPdFF0Std(i),EphysAPdFF0Std(i),'.','MarkerSize',30)
end
x = EphysAPdFF0Mean;
y = -FstimShortdFF0Mean;
[k, R2] = libplot.lnrFitting(x,y,false,true,7);
xlabel({'Ephys Single AP','-\DeltaF/F0'})
ylabel({'Fstim = 1ms ','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis([0 inf 0 inf])
axis square
legend(Single_AP.Name(sel),'location','best')
set(gcf, 'Position',  [550, 10, 500, 500])

figure(8)
clf
hold on
EphysTrainAPdFF0Mean = [];
EphysTrainAPdFF0Std = [];
FstimLongdFF0Mean = [];
FstimLongdFF0Std = [];
for i = 1:length(sel) %size(Train_AP,1)
    EphysTrainAPdFF0Mean(i) = 1-Train_AP(sel(i),5).Variables;
    EphysTrainAPdFF0Std(i) = Train_AP(sel(i),6).Variables;
    FstimLongdFF0Mean(i) = T_well_stats(Ephys2Fstim(i),strcat(LongStim,'Mean')).Variables;
    FstimLongdFF0Std(i) = T_well_stats(Ephys2Fstim(i),strcat(LongStim,'Std')).Variables;
    errorbar(EphysTrainAPdFF0Mean(i),-FstimLongdFF0Mean(i),FstimLongdFF0Std(i),FstimLongdFF0Std(i),EphysTrainAPdFF0Std(i),EphysTrainAPdFF0Std(i),'.','MarkerSize',30)
end
x = EphysTrainAPdFF0Mean;
y = -FstimLongdFF0Mean;
p = fitlm(x,y,'Intercept',false);
[k, R2] = libplot.lnrFitting(x,y,false,true,8);
xlabel({'Ephys train AP','-\DeltaF/F0'})
ylabel({'Fstim = 100ms 100Hz','-\DeltaF/F0'})
title(['Ephys vs Fstim (',CtrlFP,')'])
axis square
axis([0 inf 0 inf])
legend(Train_AP.Name(sel),'location','best')
set(gcf, 'Position',  [1550, 10, 500, 500])


