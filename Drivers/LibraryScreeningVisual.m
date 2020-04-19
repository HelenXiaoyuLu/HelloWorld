% Library visualization
%% Load multiple libraries from Eric
clear
dirPath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\2pLibraries\*.mat';
files = dir(dirPath);
Libraries = struct;
for i=1:length(files)
    L = load(fullfile(dirPath,files(i).name));
    Libraries(i).Data = L;
end
 %% load library from Wellbase_RefineVisual

% load('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20200223 repeated screening\200223 repeat screening for GV75and GV79_plate1_1ptrigger_result_ver200409.mat');
% T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3','Stim_4'}}, 2); 
% T_well.stim_shortstd = std(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3','Stim_4'}}, 1, 2); 
T_well.stim_shortmean = mean(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3'}}, 2); 
T_well.stim_shortstd = std(T_well{:, {'Stim_1', 'Stim_2' ,'Stim_3'}}, 1, 2); 
T_well.stim_midmean = mean(T_well{:, {'Stim_4', 'Stim_5' ,'Stim_6'}}, 2); 
T_well.stim_midstd = std(T_well{:, {'Stim_4', 'Stim_5' ,'Stim_6'}}, 1, 2); 
T_well.Kineticsmean = T_well.stim_shortmean./; 
T_well.Kineticsstd = std(T_well{:, {'Stim_4', 'Stim_5' ,'Stim_6'}}, 1, 2); 

%% Barplot: Short & Long stim
SortOn = 'Stim_5'; % Input Sorting Criteria 
Direction = 'descend';
T_well_sorted = sortrows(T_well,SortOn,Direction);
ctrlRank = find(T_well_sorted.Name == 'ASAP2'); % Get control

figure()
clf
subplot(3,1,1)
hold on
X = 1:height(T_well);
hBarLong = bar(X,-T_well_sorted.Stim_5,'BarWidth',1,'LineWidth', 1,'FaceAlpha',1);
hBarShort = bar(X,-T_well_sorted.stim_shortmean,'BarWidth',0.8,'LineWidth', 1,'FaceAlpha',1);
hBarShort(1).FaceColor = [0.8500, 0.3250, 0.0980];	
hBarLong(1).FaceColor = [0.9290, 0.6940, 0.1250];
axis([0 100 0 inf])
title({'20190223 Plate 1 1P Fstim','well-level stat New Pipeline'})
ctrlY = -T_well_sorted(ctrlRank,[10,17]).Variables;
t = text(ctrlRank,ctrlY(1),{'\bf{\downarrow}'},'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r','FontSize',16);
plot(X,ones(size(X,2))*ctrlY(1),'--k')
plot(X,ones(size(X,2))*ctrlY(2),'--k')
text(X(end)+1,ctrlY(1),{'\leftarrow ASAP2s', '100 ms'}','HorizontalAlignment','left','VerticalAlignment', 'baseline')
text(X(end)+1,ctrlY(2),{'\leftarrow ASAP2s', '1 ms'}','HorizontalAlignment','left','VerticalAlignment', 'baseline')
legend('100 ms','1 ms','location','best')
xlabel('NO. of Variant Screened')
ylabel({'Response Amplitude','-\DeltaF/F0'})

subplot(3,1,2)
hold on
X = 1:height(T_well);
hBar = bar(X,T_well_sorted.GRRatioMean,'BarWidth',1,'LineWidth', 1,'FaceAlpha',1);
hBar.FaceColor = [0.4660, 0.6740, 0.1880];	
ctrlY = T_well_sorted(ctrlRank,11).Variables;
plot(X,ones(size(X))*ctrlY,'--k')
text(ctrlRank,ctrlY(1),{'\bf{\downarrow}'},'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r','FontSize',16);
text(X(end)+1,ctrlY,'\bf{\leftarrow ASAP2s}','HorizontalAlignment','left')% 
axis([0 100 0 inf])
xlabel('NO. of Variant Screened')
ylabel('Relative Brightness (G/R)')
title({'20190223 Plate 1 GR Ratio','well-level stat New Pipeline'})

subplot(3,1,3)
hold on
hBar = bar(X,T_well_sorted.AOCMean,'BarWidth',1,'LineWidth', 1,'FaceAlpha',0.5);
hBar.FaceColor = [0.4660, 0.6740, 0.1880];	
ctrlY = T_well_sorted(ctrlRank,12).Variables;
plot(X,ones(size(X))*ctrlY,'--k')
text(X(end)+1,ctrlY,'\bf{\leftarrow ASAP2s}','HorizontalAlignment','left')
text(ctrlRank,ctrlY,{'\bf{\downarrow}'},'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r','FontSize',16)
axis([0 100 0 inf])
xlabel('NO. of Variant Screened')
ylabel('Photostability (Area-Under-Curve)')
title({'20190223 Plate 1 Photostability','well-level stat New Pipeline'})
set(gcf, 'Position',  [10, 10, 650, 950])
%% N-D plot
ShortStim = nan(96,length(Libraries));
LongStim = nan(96,length(Libraries));
GRratio = nan(96,length(Libraries));
LibName = cell(1,length(Libraries));
skip = [];
for i = 1:length(Libraries)
    if isfield(Libraries(i).Data,'T_well')
        S = Libraries(i).Data.T_well.dFF0_shortStim;
        L = Libraries(i).Data.T_well.dFF0_longStim;
        GR = Libraries(i).Data.T_well.GRratioMean;
        ShortStim(1:length(L),i) = S;
        LongStim(1:length(L),i) = L;
        GRratio(1:length(L),i) = GR;
        LibName{i} = files(i).name;
    else
        skip = [skip,i];
    end
end
figure()
hold on
plot3(LongStim,ShortStim./LongStim,GRratio,'.','MarkerSize',20)
zlim([0,2])
axis on
grid on
xlabel('dF/F0 Long Stim')
ylabel('kinetics (short/long)')
zlabel('GR ratio Mean')

%% Violin plot
figure()
hold on
Y = nan(96,length(Libraries));
LibName = cell(1,length(Libraries));
skip = [];
for i = 1:length(Libraries)
    if isfield(Libraries(i).Data,'T_well')
        L = Libraries(i).Data.T_well.dFF0_shortStim;
        Y(1:length(L),i) = L;
        LibName{i} = files(i).name;
    else
        skip = [skip,i];
    end
end
notEmpty = 1:length(Libraries);
notEmpty(skip) = [];
cmap = colormap(lines);
[h,L,MX,MED]=violin(Y(:,notEmpty),'facecolor',cmap(1:length(notEmpty),:),'mc','k-','medc','k--');
% scatter(1:7,Y,'filled','MarkerFaceAlpha',0.5)
set(gca, 'XTickLabel',LibName(notEmpty), 'XTick',1:numel(notEmpty))
ylabel('dFF0_shortStim')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];

figure(4)
clf
hold on
Y = nan(96,length(Libraries));
LibName = cell(1,length(Libraries));
skip = [];
for i = 1:length(Libraries)
    if isfield(Libraries(i).Data,'T_well')
        L = Libraries(i).Data.T_well.dFF0_shortStim;
        Y(1:length(L),i) = L;
        LibName{i} = files(i).name;
    else
        skip = [skip,i];
    end
end
notEmpty = 1:length(Libraries);
notEmpty(skip) = [];
cmap = colormap(lines);
figure(2)
clf
hold on
violin(LongStim(:,notEmpty))
figure(3)
clf
hold on
[h,L,MX,MED]=dualViolin(LongStim(:,notEmpty),ShortStim(:,notEmpty),'facecolor',cmap(1:length(notEmpty),:),'mc','k-','medc','k--');
% scatter(1:7,Y,'filled','MarkerFaceAlpha',0.5)
set(gca, 'XTickLabel',LibName(notEmpty), 'XTick',1:numel(notEmpty))
ylabel('dFF0_shortStim')
ax = gca;
ax.XTickLabelRotation=90;
ax.XAxis.TickLength = [0 0];
