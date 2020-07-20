%% Read excel from google drive
% Get spread sheet
JEDIup = table();
JEDIup.sheetName = {'1pFstim_1khz'};
JEDIup.docid = {'1he6FBO_qr0Luy96mhRdDinxZkMBjf8q47ct3GHb8RxM'}
JEDIup.gid = {'0'};

%% Load data from google sheet (d2B)
for i = 1:height(JEDIup)
    result = GetGoogleSpreadsheet(JEDIup.docid{i},JEDIup.gid{i});
    JEDIup.evalout{i} = numcell(result);
end

%%
FstimUP = libplot.tableparse(JEDIup.evalout{1});
GRratio = FstimUP{:,find(strcmp('Ratio of Mean Brightness (C2/C3)',FstimUP.Properties.VariableNames))};
FstimUP.('Detectability index (20V 1.5ms)') = sqrt(GRratio).* abs(FstimUP{:,...
    find(strcmp('dFF0 of Mean (C1 S1)',FstimUP.Properties.VariableNames))});
FstimUP.('Detectability index (20V 20ms)') = sqrt(GRratio).* abs(FstimUP{:,...
    find(strcmp('dFF0 of Mean (C1 S2)',FstimUP.Properties.VariableNames))});
FstimUP.('Detectability index (20V 100ms)') = sqrt(GRratio).* abs(FstimUP{:,...
    find(strcmp('dFF0 of Mean (C1 S3)',FstimUP.Properties.VariableNames))});
FstimUP = libplot.normalize2ctrl(FstimUP,'dFF0 of Mean (C1 S1)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'dFF0 of Mean (C1 S2)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'dFF0 of Mean (C1 S3)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'Ratio of Mean Brightness (C2/C3)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'Detectability index (20V 1.5ms)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'Detectability index (20V 20ms)','JEDIalpha');
FstimUP = libplot.normalize2ctrl(FstimUP,'Detectability index (20V 100ms)','JEDIalpha');
rmRows = find(FstimUP.('dFF0 of Mean (C1 S3)')<0);
FsimUPfilt = FstimUP;
FsimUPfilt(rmRows,:)=[];
rmRows = find(FsimUPfilt.('masked')<50);
FsimUPfilt(rmRows,:)=[];

%% Plot
% Multicolor 2D DI scatter
libplot.scatterbylib(FsimUPfilt,'Detectability index (20V 100ms)','Detectability index (20V 1.5ms)',...
    'controls',{'JEDIalpha'},'mskth', 40,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
xlabel({'Detectability index (20V 100ms)','normalized to JEDI-alpha'})
ylabel({'Detectability index (20V 1.5ms)','normalized to JEDI-alpha'})
axis square
set(gcf, 'Position',  [110, 100, 500, 500])

% Multicolor 2D response scatter
libT = libplot.scatterbylib(FsimUPfilt,'dFF0 of Mean (C1 S3)','dFF0 of Mean (C1 S1)',...
    'controls',{'JEDIalpha'},'mskth', 40,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
xlabel({'dFF0 of Mean (20V 100ms)','normalized to JEDI-alpha'})
ylabel({'dFF0 of Mean (20V 1.5ms)','normalized to JEDI-alpha'})
axis square
set(gcf, 'Position',  [110, 100, 500, 500])

cmapSC = [0, 0.4470, 0.7410];
libplot.scatterbylib(FsimUPfilt,'Detectability index (20V 100ms)','Detectability index (20V 1.5ms)',...
    'controls',{'JEDIalpha'},'mskth', 0,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',},'colormap',cmapSC);
C = get(gca).Children;
C(1).DisplayName = 'JEDI-alpha';
C(2).DisplayName = 'Libraries';
legend([C(1),C(2)],'location','best')
xlabel({'Detectability index (20V 100ms)','normalized to JEDI-alpha'})
ylabel({'Detectability index (20V 1.5ms)','normalized to JEDI-alpha'})
axis square
set(gcf, 'Position',  [110, 100, 500, 500])

% Plot
libplot.scatterbylib(FsimUPfilt,'dFF0 of Mean (C1 S3)','dFF0 of Mean (C1 S1)',...
    'controls',{'JEDIalpha'},'datatip',...
    {'Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});

libplot.scatterbylib(FsimUPfilt,'Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S3)',...
    'controls',{'JEDIalpha'},'datatip',...
    {'Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});

cmap1 = [0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
libplot.barplotall(FsimUPfilt,'dFF0 of Mean (C1 S1)',...
    'controls',{'JEDIalpha'},'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'Name','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)'},'colormap',cmap1(2,:));
xlabel('Number of construct screened')
ylabel({'\DeltaF/F0 (20V 1.5ms)';'normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 800, 300])
ylim([-0.5,inf])

libplot.barplotall(FsimUPfilt,'dFF0 of Mean (C1 S3)',...
    'controls',{'JEDIalpha'},'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'Name','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)'},'colormap',cmap1(1,:));
xlabel('Number of construct screened')
ylabel({'\DeltaFF0 (20V 100ms)';'normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 800, 300])
ylim([-0.5,inf])

cmap2 = [0.4660, 0.6740, 0.1880];
libplot.barplotall(FsimUPfilt,'Ratio of Mean Brightness (C2/C3)',...
    'controls',{'JEDIalpha'},'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'Name','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)'},'colormap',cmap2(1,:));
xlabel('Number of construct screened')
ylabel({'\DeltaFF0 (20V 100ms)';'normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 800, 300])
ylim([-0.5,inf])
%%
sel = {'JEDIalpha','JEDIalpha-A147','JEDIalpha-A148','JEDIalpha-D395','JEDIalpha-F152',...
    'JEDIalpha-F391','JEDIalpha-N150','JEDIalpha-N392','JEDIalpha-Q369','JEDIalpha-S393',...
    'JEDIalpha-T394','JEDIalpha-V151','JEDIalpha-Y149','JEDIbeta'};
libplot.scatterbylib1d(FsimUPfilt,'dFF0 of Mean (C1 S2)',...
    'mskth', 40,'selectLib',sel,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
ylabel({'\DeltaFF0 (20V 20ms)','normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 500, 500])
set(gca, 'XTickLabel',sel, 'XTick',1:numel(sel));
ax = gca;
ax.XTickLabelRotation = 45;


libplot.scatterbylib1d(FsimUPfilt,'dFF0 of Mean (C1 S1)',...
    'mskth', 40,'selectLib',sel,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
ylabel({'\DeltaFF0 (20V 1.5ms)','normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 500, 500])
set(gca, 'XTickLabel',sel, 'XTick',1:numel(sel));
ax = gca;
ax.XTickLabelRotation = 45;

libplot.scatterbylib1d(FsimUPfilt,'dFF0 of Mean (C1 S3)',...
    'mskth', 40,'selectLib',sel,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
ylabel({'\DeltaFF0 (20V 100ms)','normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 500, 500])
set(gca, 'XTickLabel',sel, 'XTick',1:numel(sel));
ax = gca;
ax.XTickLabelRotation = 45;

libplot.scatterbylib1d(FsimUPfilt,'Ratio of Mean Brightness (C2/C3)',...
    'mskth', 40,'selectLib',sel,'norm2ctrl',{'JEDIalpha'},'datatip',...
    {'masked','Ratio of Mean Brightness (C2/C3)','dFF0 of Mean (C1 S1)',...
    'dFF0 of Mean (C1 S2)','dFF0 of Mean (C1 S3)',...
    'dF/F0 median_Bleaching1Pstat','Detectability index (20V 1.5ms)',...
    'Detectability index (20V 20ms)','Detectability index (20V 100ms)',});
ylabel({'\DeltaFF0 (20V 100ms)','normalized to JEDI-alpha'})
set(gcf, 'Position',  [110, 100, 500, 500])
set(gca, 'XTickLabel',sel, 'XTick',1:numel(sel));
ax = gca;
ax.XTickLabelRotation = 45;

%%

    normctrlsel = strmatch('JEDIalpha', FstimUP.plibName,'exact') ;
    Tnormctrl = FstimUP(normctrlsel,:);
    meanctrlppl = nanmean(Tnormctrl{:,7});
    meanctrlpp2 = nanmean(Tnormctrl{:,8});
    meanctrlpp3 = nanmean(Tnormctrl{:,9});
    meanctrlpp4 = nanmean(Tnormctrl{:,10});
 

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






