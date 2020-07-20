%% Read excel from google drive
% Get spread sheet
clear
clc
JEDImCherry = table();
JEDImCherry.sheetName = {'1pFstim_1khz';'1pbleaching';'2pbleaching'};
JEDImCherry.docid = {'17m5iO7dJM1xuijDYHSNsL6hGuOHNOqEpKPnMD-04Hs4';...
    '17m5iO7dJM1xuijDYHSNsL6hGuOHNOqEpKPnMD-04Hs4';...
    '17m5iO7dJM1xuijDYHSNsL6hGuOHNOqEpKPnMD-04Hs4'};
JEDImCherry.gid = {'716650162';'1486321860';'1604302427'};
JEDIcyOFP = table();
% JEDIcyOFP.sheetName = {'JEDI1P-cyOFP library with V2 analysis pipeline';...
%     '1pPhotobleaching JEDI1P-cyOFP'};
% JEDIcyOFP.docid = {'16IL72dATUdLvuhA-hzg8qzaxR8266RBePQxEDbqEoos';...
%     '16IL72dATUdLvuhA-hzg8qzaxR8266RBePQxEDbqEoos'};
% JEDIcyOFP.gid = {'445171161';'1909282833'};
JEDIcyOFP.sheetName = {'JEDI1P-cyOFP library with V2 analysis pipeline';...
    '1pPhotobleaching JEDI1P-cyOFP'};
JEDIcyOFP.docid = {'1ZayP7wc3aNayzZa4sFXPuHHed_LpWJqnb15syBKjVLo';...
    '1ZayP7wc3aNayzZa4sFXPuHHed_LpWJqnb15syBKjVLo'};
JEDIcyOFP.gid = {'0';'1991604789'};

%% Load data from google sheet (mCherry)
for i = 1:height(JEDImCherry)
    result = GetGoogleSpreadsheet(JEDImCherry.docid{i},JEDImCherry.gid{i});
    JEDImCherry.evalout{i} = numcell(result);
end

%% Load data from google sheet (cyOFP)
for i = 1:height(JEDIcyOFP)
    result = GetGoogleSpreadsheet(JEDIcyOFP.docid{i},JEDIcyOFP.gid{i});
    JEDIcyOFP.evalout{i} = numcell(result);
end

%% Parse mCherry libraries & 
% scatter plot by library
FstimC = libplot.tableparse(JEDImCherry.evalout{1});
Bleaching1PC = libplot.tableparse(JEDImCherry.evalout{2},'parseName',true);
Bleaching1PC(ismember(Bleaching1PC.Coord(:,5),1),:)=[];
Bleaching1PCstat = avgtable(Bleaching1PC);
Bleaching1PCstat.('1-dFF0 median') = 1 + Bleaching1PCstat.('dF/F0 of median');
Bleaching1PCstat.('1-dFF0 mean') = 1 + Bleaching1PCstat.('dF/F0 of mean');
Bleaching1PCstat = libplot.normalize2ctrl(Bleaching1PCstat,'1-dFF0 median','JEDI2P');
Bleaching1PCstat = libplot.normalize2ctrl(Bleaching1PCstat,'1-dFF0 mean','JEDI2P');
Bleaching1PCstat(:,7:10) = [];

Bleaching2PC = libplot.tableparse(JEDImCherry.evalout{3},'parseName',true);
Bleaching2PC(ismember(Bleaching2PC.Coord(:,5),1),:)=[];
Bleaching2PCstat = avgtable(Bleaching2PC);
Bleaching2PCstat.('1-dFF0 median') = 1 + Bleaching2PCstat.('dF/F0 of median');
Bleaching2PCstat.('1-dFF0 mean') = 1 + Bleaching2PCstat.('dF/F0 of mean');
Bleaching2PCstat = libplot.normalize2ctrl(Bleaching2PCstat,'1-dFF0 median','JEDI2P');
Bleaching2PCstat = libplot.normalize2ctrl(Bleaching2PCstat,'1-dFF0 mean','JEDI2P');
Bleaching2PCstat(:,7:10) = [];
%
revdFF0S1 = -FstimC{:,find(strcmp('dFF0 of Mean (C1 S1)',FstimC.Properties.VariableNames))};
revdFF0S2 = -FstimC{:,find(strcmp('dFF0 of Mean (C1 S2)',FstimC.Properties.VariableNames))};
revdFF0S3 = -FstimC{:,find(strcmp('dFF0 of Mean (C1 S3)',FstimC.Properties.VariableNames))};
revdFF0S4 = -FstimC{:,find(strcmp('dFF0 of Mean (C1 S4)',FstimC.Properties.VariableNames))};
revdFF0S5 = -FstimC{:,find(strcmp('dFF0 of Mean (C1 S5)',FstimC.Properties.VariableNames))};
GRratio = FstimC{:,find(strcmp('Ratio of Mean Brightness (C2/C3)',FstimC.Properties.VariableNames))};
FstimC{:,end+1} = mean([revdFF0S1,revdFF0S2,revdFF0S3,revdFF0S4],2);
FstimC.Properties.VariableNames{end} = '-dFF0 of Mean (60V 1ms)';
FstimC{:,end+1} = revdFF0S5;
FstimC.Properties.VariableNames{end} = '-dFF0 of Mean (30V 100ms)'; % 2.5ms@100Hz
FstimC{:,end+1} = sqrt(GRratio).* abs(mean([revdFF0S1,revdFF0S2,revdFF0S3,revdFF0S4],2));
FstimC.Properties.VariableNames{end} = 'Detectability index (60V 1ms)';
FstimC{:,end+1} = sqrt(GRratio).* abs(revdFF0S5);
FstimC.Properties.VariableNames{end} = 'Detectability index (30V 100ms)';
FstimC.Properties.VariableNames{end} = 'Detectability index (30V 100ms)';
FstimC = libplot.normalize2ctrl(FstimC,'-dFF0 of Mean (60V 1ms)','JEDI2P');
FstimC = libplot.normalize2ctrl(FstimC,'-dFF0 of Mean (30V 100ms)','JEDI2P');
FstimC = libplot.normalize2ctrl(FstimC,'Ratio of Mean Brightness (C2/C3)','JEDI2P');
FstimC = libplot.normalize2ctrl(FstimC,'Detectability index (60V 1ms)','JEDI2P');
FstimC = libplot.normalize2ctrl(FstimC,'Detectability index (30V 100ms)','JEDI2P');

Bleaching1Pstat = Bleaching1PCstat;
Bleaching2Pstat = Bleaching2PCstat;
Fstim = FstimC;
Screening = outerjoin(Bleaching1Pstat,Bleaching2Pstat,'LeftKeys','Name','RightKeys','Name','MergeKeys',true);
ScreeningC = outerjoin(FstimC,Screening,'LeftKeys','Name','RightKeys','Name','MergeKeys',true);
refcol = find(strcmp('Ratio of Mean Brightness (C2/C3)',ScreeningC.Properties.VariableNames));
ScreeningC(isnan(ScreeningC{:,refcol}),:)=[];


%% Plot mCherry library screening
libplot.scatterbylib(ScreeningC,'Ratio of Mean Brightness (C2/C3)',...
    '-dFF0 of Mean (60V 1ms)','controls',{'JEDI2P'},'selectLib',{'JEDI2P'},...
    'datatip',{'-dFF0 of Mean (30V 100ms)','dF/F0 median_Bleaching1Pstat','dF/F0 median_Bleaching2Pstat'});
libplot.scatterbylib(ScreeningC,'Detectability index (30V 100ms)','Detectability index (60V 1ms)',...
    'controls',{'JEDI2P'},'selectLib',{'JEDI2P'},'datatip',...
    {'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)',...
    '-dFF0 of Mean (30V 100ms)','dF/F0 median_Bleaching1Pstat','dF/F0 median_Bleaching2Pstat'});

% barplot mCherry library screening
close all
cmap = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
[TJ2,meanJ2]=libplot.barplotall(ScreeningC,'-dFF0 of Mean (30V 100ms)','-dFF0 of Mean (60V 1ms)','norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap,'FaceAlpha',0.9);
% [TJ1,meanJ1]=libplot.barplotall(Fstim,'-dFF0 of Mean (30V 100ms)','-dFF0 of Mean (60V 1ms)','controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'colormap',cmap,'FaceAlpha',0.9);
cmap = [0.4660, 0.6740, 0.1880];
libplot.barplotall(TJ2,'dF/F0 median_Bleaching2Pstat','norm2ctrl',{'JEDI2P'},...
    'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap,...
    'sorted','ascend','FaceAlpha',0.9,'datatip',{'Name','Ratio of Mean Brightness (C2/C3)',...
    '-dFF0 of Mean (60V 1ms)','-dFF0 of Mean (30V 100ms)','dF/F0 median_Bleaching1Pstat'},'plotctrl',false);
libplot.barplotall(TJ2,'dF/F0 median_Bleaching1Pstat','norm2ctrl',{'JEDI2P'},...
    'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap,...
    'sorted','ascend','FaceAlpha',0.9,'datatip',{'Name','Ratio of Mean Brightness (C2/C3)',...
    '-dFF0 of Mean (60V 1ms)','-dFF0 of Mean (30V 100ms)','dF/F0 median_Bleaching1Pstat'},'plotctrl',false);
libplot.barplotall(TJ2,'Ratio of Mean Brightness (C2/C3)','norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap,'FaceAlpha',0.9);
libplot.barplotall(TJ2,'Ratio of Mean Brightness (C2/C3)','norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap,'sorted','false','FaceAlpha',0.9);

%%
    normctrlsel = strmatch('JEDI2P', FstimO.plibName,'exact') ;
    Tnormctrl = FstimO(normctrlsel,:);
    meanctrlppl = nanmean(Tnormctrl{:,15});
    meanctrlpp2 = nanmean(Tnormctrl{:,16});
    meanctrlpp3 = nanmean(Tnormctrl{:,17});
    meanctrlpp4 = nanmean(Tnormctrl{:,18});
    
%% Parse cyOFP libraries & 
% scatter plot by library
close all

FstimO = libplot.tableparse(JEDIcyOFP.evalout{1});
Bleaching1PO = libplot.tableparse(JEDIcyOFP.evalout{2},'parseName',true);
Bleaching1PO(ismember(Bleaching1PO.Coord(:,5),1),:)=[];
Bleaching1POstat = avgtable(Bleaching1PO);
Bleaching1POstat.('1-dFF0 median') = 1 + Bleaching1POstat.('dF/F0 of median');
Bleaching1POstat.('1-dFF0 mean') = 1 + Bleaching1POstat.('dF/F0 of mean');
Bleaching1POstat = libplot.normalize2ctrl(Bleaching1POstat,'1-dFF0 median','JEDI2P');
Bleaching1POstat = libplot.normalize2ctrl(Bleaching1POstat,'1-dFF0 mean','JEDI2P');
Bleaching1POstat(:,10:13) = [];
%
revdFF0S1 = -FstimO{:,find(strcmp('dFF0 of Mean (C1 S1)',FstimO.Properties.VariableNames))};
revdFF0S2 = -FstimO{:,find(strcmp('dFF0 of Mean (C1 S2)',FstimO.Properties.VariableNames))};
revdFF0S3 = -FstimO{:,find(strcmp('dFF0 of Mean (C1 S3)',FstimO.Properties.VariableNames))};
revdFF0S4 = -FstimO{:,find(strcmp('dFF0 of Mean (C1 S4)',FstimO.Properties.VariableNames))};
revdFF0S5 = -FstimO{:,find(strcmp('dFF0 of Mean (C1 S5)',FstimO.Properties.VariableNames))};
GRratio = FstimO{:,find(strcmp('Ratio of Mean Brightness (C2/C3)',FstimO.Properties.VariableNames))};
FstimO{:,end+1} = mean([revdFF0S1,revdFF0S2,revdFF0S3,revdFF0S4],2);
FstimO.Properties.VariableNames{end} = '-dFF0 of Mean (60V 1ms)';
FstimO{:,end+1} = revdFF0S5;
FstimO.Properties.VariableNames{end} = '-dFF0 of Mean (30V 100ms)'; % 2.5ms@100Hz
FstimO{:,end+1} = sqrt(GRratio).* abs(mean([revdFF0S1,revdFF0S2,revdFF0S3,revdFF0S4],2));
FstimO.Properties.VariableNames{end} = 'Detectability index (60V 1ms)';
FstimO{:,end+1} = sqrt(GRratio).* abs(revdFF0S5);
FstimO.Properties.VariableNames{end} = 'Detectability index (30V 100ms)';
FstimO = libplot.normalize2ctrl(FstimO,'-dFF0 of Mean (60V 1ms)','JEDI2P');
FstimO = libplot.normalize2ctrl(FstimO,'-dFF0 of Mean (30V 100ms)','JEDI2P');
FstimO = libplot.normalize2ctrl(FstimO,'Ratio of Mean Brightness (C2/C3)','JEDI2P');
FstimO = libplot.normalize2ctrl(FstimO,'Detectability index (60V 1ms)','JEDI2P');
FstimO = libplot.normalize2ctrl(FstimO,'Detectability index (30V 100ms)','JEDI2P');

extrappName = {'Plate name','Ratio of Mean Brightness (C1/C3)','Variant Name','FOV Name','ROI Name'};
for i = 1:length(extrappName)
    extracol = find(strcmp(extrappName{i},FstimO.Properties.VariableNames));
    FstimO(:,extracol) = [];
    extracol = find(strcmp(extrappName{i},Bleaching1POstat.Properties.VariableNames));
    Bleaching1POstat(:,extracol) = [];
end
%
Bleaching1Pstat = Bleaching1POstat;
Bleaching2Pstat = table();
for i = 1:width(Bleaching2PCstat)
    if strmatch(Bleaching2PCstat.Properties.VariableNames{i},'Name','exact')
        Bleaching2Pstat.Name = Bleaching1POstat.Name;
    elseif iscell(Bleaching2PCstat{2,i})
        Bleaching2Pstat.(i) = cell(height(Bleaching1POstat),1);
        Bleaching2Pstat.Properties.VariableNames{i} = Bleaching2PCstat.Properties.VariableNames{i};
    elseif isnumeric(Bleaching2PCstat{2,i})
        Bleaching2Pstat.(i) = nan(height(Bleaching1POstat),1);
        Bleaching2Pstat.Properties.VariableNames{i} = Bleaching2PCstat.Properties.VariableNames{i};
    end
end

Fstim = FstimO;
Screening = outerjoin(Bleaching1Pstat,Bleaching2Pstat,'LeftKeys','Name','RightKeys','Name','MergeKeys',true);
ScreeningO = outerjoin(Fstim,Screening,'LeftKeys','Name','RightKeys','Name','MergeKeys',true);
refcolO = find(strcmp('Ratio of Mean Brightness (C2/C3)',ScreeningO.Properties.VariableNames));
ScreeningO(isnan(ScreeningO{:,refcolO}),:)=[];


%% plot cyOFP library screening
libplot.scatterbylib(FstimO,'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)','controls',{'JEDI1P'},'selectLib',{'JEDI1P'});
libplot.scatterbylib(FstimO,'Detectability index (30V 100ms)','Detectability index (60V 1ms)','controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'datatip',{'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)','-dFF0 of Mean (30V 100ms)'});
libplot.scatterbylib(FstimO,'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)','controls',{'JEDI1P'},'selectLib',{'JEDI1P'});

% barplot cyOFP library screening
cmap = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
[TJ1O,meanJ2]=libplot.barplotall(FstimO,'-dFF0 of Mean (30V 100ms)','-dFF0 of Mean (60V 1ms)','norm2ctrl',{'JEDI1P'},'controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'colormap',cmap,'FaceAlpha',0.9);
% [TJ1,meanJ1]=libplot.barplotall(Fstim,'-dFF0 of Mean (30V 100ms)','-dFF0 of Mean (60V 1ms)','controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'colormap',cmap,'FaceAlpha',0.9);
cmap = [0.4660, 0.6740, 0.1880];
libplot.barplotall(TJ1O,'Ratio of Mean Brightness (C2/C3)','norm2ctrl',{'JEDI1P'},'controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'colormap',cmap,'FaceAlpha',0.9);
libplot.barplotall(TJ1O,'Ratio of Mean Brightness (C2/C3)','norm2ctrl',{'JEDI1P'},'controls',{'JEDI1P'},'selectLib',{'JEDI1P'},'colormap',cmap,'sorted','false','FaceAlpha',0.9);
LibPos = unique(FstimAll.plibName);

%% Merge cyOFP and mCherry
close all
rm = [];
for i = 1:length(ScreeningC.Properties.VariableNames)
    diff = strmatch(ScreeningC.Properties.VariableNames{i},ScreeningO.Properties.VariableNames,'exact');
    if isempty(diff)
        rm = [rm,i];
    end
end
ScreeningC(:,rm) = [];
rm = [];
for i = 1:length(ScreeningO.Properties.VariableNames)
    diff = strmatch(ScreeningO.Properties.VariableNames{i},ScreeningC.Properties.VariableNames,'exact');
    if isempty(diff)
        rm = [rm,i];
    end
end 
ScreeningO(:,rm) = [];

dtS = {'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)',...
    '-dFF0 of Mean (30V 100ms)','dF/F0 of median_Bleaching1Pstat',...
    'dF/F0 of median_Bleaching2Pstat'};
dtB = {'Name','Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)',...
    '-dFF0 of Mean (30V 100ms)','dF/F0 of median_Bleaching1Pstat',...
    'dF/F0 of median_Bleaching2Pstat'};
FstimAll = [FstimC;FstimO];
ScreeningAll = [ScreeningC;ScreeningO];

%%
close all
libplot.scatterbylib(ScreeningAll,'Detectability index (30V 100ms)Norm2JEDI2P',...
    'Detectability index (60V 1ms)Norm2JEDI2P','controls',{'JEDI2P'},'norm2ctrl',...
    {'JEDI2P'},'datatip',dtS);
cmapSC = [0, 0.4470, 0.7410];
linkdata on
brush on
libplot.scatterbylib(ScreeningAll,'Detectability index (30V 100ms)Norm2JEDI2P',...
    'Detectability index (60V 1ms)Norm2JEDI2P','controls',{'JEDI2P'},'norm2ctrl',...
    {'JEDI2P'},'datatip',dtS,'colormap',cmapSC);
C = get(gca).Children;
C(1).DisplayName = 'JEDI-2P';
C(2).DisplayName = 'Libraries';
legend([C(1),C(2)],'location','best')
xlabel({'Detectability index (30V 100ms)';'normalized to JEDI-2P'})
ylabel({'Detectability index (60V 1ms)';'normalized to JEDI-2P'})
axis square
set(gcf, 'Position',  [110, 100, 300, 300])
linkdata on
brush on

cmap1 = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
% Fstim bar overlap plot ranked
[TJ1all,meanJ2,~,F]=libplot.barplotall(ScreeningAll,'-dFF0 of Mean (30V 100ms)Norm2JEDI2P',...
    '-dFF0 of Mean (60V 1ms)Norm2JEDI2P','norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},...
    'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap1,'FaceAlpha',0.9,'datatip',dtB);

% 100ms bar plot ranked
[TJ1all,meanJ2,~,F]=libplot.barplotall(ScreeningAll,'-dFF0 of Mean (30V 100ms)Norm2JEDI2P',...
    'norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},...
    'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap1(1,:),'FaceAlpha',0.9,'datatip',dtB);
ylim([-0.5,inf])
xlabel('Number of construct screened')
ylabel({'\DeltaFF0 (30V 100ms)';'normalized to JEDI-2P'})
set(gcf, 'Position',  [110, 100, 800, 300])

% 1ms bar plot ranked
[TJ1all,meanJ2,~,F]=libplot.barplotall(ScreeningAll,...
    '-dFF0 of Mean (60V 1ms)Norm2JEDI2P','norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},...
    'selectLib',{'JEDI1P','JEDI2P'},'colormap',cmap1(2,:),'FaceAlpha',0.9,'datatip',dtB);
ylim([-0.5,inf])
set(gcf, 'Position',  [110, 100, 1200, 400])
ylabel({'\DeltaFF0 (60V 1ms)';'normalized to JEDI-2P'})

% GRratio ranked
cmap2 = [0.4660, 0.6740, 0.1880];
[~,~,~,gGRSorted] = libplot.barplotall(TJ1all,'Ratio of Mean Brightness (C2/C3)Norm2JEDI2P',...
    'norm2ctrl',{'JEDI2P'},'controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},...
    'colormap',cmap2,'FaceAlpha',0.9,'datatip',dtB);
ylabel({'Relative brightness (G/R)';'normalized to JEDI-2P'})

libplot.barplotall(TJ1all,'Ratio of Mean Brightness (C2/C3)','norm2ctrl',...
    {'JEDI2P'},'controls',{'JEDI2P'},'colormap',cmap2,'sorted','false','FaceAlpha',0.9,'datatip',dtB);
libplot.barplotall(ScreeningAll,'1-dFF0 medianNorm2JEDI2P_Bleaching1Pstat','norm2ctrl',...
    {'JEDI2P'},'controls',{'JEDI2P'},'colormap',cmap2,'FaceAlpha',0.9,'plotctrl',true,'datatip',dtB);
xlim([0 510])
set(gcf, 'Position',  [110, 100, 700, 300])
 plot(0:1:510, ones(1,511),'--k')
ylabel({'1P photostability';'normalized to JEDI-2P'})


libplot.barplotall(ScreeningAll,'dF/F0 of median_Bleaching2Pstat','norm2ctrl',...
    {'JEDI2P'},'controls',{'JEDI2P'},'colormap',cmap2,'FaceAlpha',0.9,'plotctrl',false,'datatip',dtB);


ylim([-0.5,inf])
set(gcf, 'Position',  [110, 100, 1200, 400])

LibPos = unique(FstimO.plibName);

libplot.scatterbylib(FstimAll,'-dFF0 of Mean (60V 1ms)Norm2JEDI1P','Detectability index (60V 1ms)Norm2JEDI1P','controls',{'JEDI2P'},'selectLib',{'JEDI1P','JEDI2P'},'datatip',{'Ratio of Mean Brightness (C2/C3)','-dFF0 of Mean (60V 1ms)','-dFF0 of Mean (30V 100ms)'});

%% function collection
function tcell = numcell(c)
    tcell = table();
    dataNum = str2double(c);
    for i = 1:size(dataNum,2)
        if ~isnan(dataNum(2,i))
            T = table(dataNum(2:end,i),'VariableNames',c(1,i));
            tcell  =[tcell,T];      
        elseif strcmp(c{1,i},'Coord')
            coordC = [];
            for j = 2:size(c,1)
                coordC(j-1,:) = str2num(c{j,i});                
            end
            T = table(coordC,'VariableNames',c(1,i));
            tcell  =[tcell,T];    
        else
            T = table(c(2:end,i),'VariableNames',c(1,i));
            tcell  =[tcell,T];      
        end
        tcell.Properties.VariableNames{end} = c{1,i};
    end
end

function avgT = avgtable(inputT)
    avgT = table();
    G = findgroups(inputT.Name);
    for i = 1:width(inputT)
        varName = inputT.Properties.VariableNames{i};
        if isnumeric(inputT{1,i}) && ~strcmp(varName,'Coord')
            tempT = varfun(@mean,inputT,'InputVariables',varName,...
           'GroupingVariables','Name');
            avgT = [avgT,tempT(:,end)];
            avgT.Properties.VariableNames{end} = varName;
        elseif iscell(inputT{1,i})
            tempT = varfun(@unique,inputT,'InputVariables',varName,...
                   'GroupingVariables','Name');
            avgT = [avgT,tempT(:,end)];
            avgT.Properties.VariableNames{end} = varName;
        end
    end
end






