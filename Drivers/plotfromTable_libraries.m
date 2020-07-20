%% Read and combine data from tables
clc;
clear;
savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\201912_photostablelibrary';
saveName = '201912_photostablelibrary_cyOFP';
fNames1P = '1pFstim*.mat';
fNames2P = '2pFstim*.mat';

% Load results table 1P and 2P
fList1P = dir(fullfile(savepath,fNames1P));
lib1P = struct;
T_well_1P = table;
for i = 1:length(fList1P)
    lib1P(i).T_well_all = load(fullfile(fList1P(i).folder,fList1P(i).name),'T_well_all').T_well_all;
    lib1P(i).T_well_all.Plate(:) = {strcat('P',num2str(i))};
    T_well_1P = [T_well_1P;lib1P(i).T_well_all]; 
end

fList2P = dir(fullfile(savepath,fNames2P));
lib2P = struct;
T_well_2P = table;
for i = 1:length(fList2P)
    lib2P(i).T_well_all = load(fullfile(fList2P(i).folder,fList2P(i).name),'T_well').T_well;
    lib2P(i).T_well_all.Plate(:) = {strcat('P',num2str(i))};
    T_well_2P = [T_well_2P;lib2P(i).T_well_all]; 
end

LookUp = readtable(fullfile(savepath,'Lookup.xlsx'));

%% Add kinetics and group rows by name 
selpp1P = {'Legal Name','wellName_T_well','ConstructName','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','dFF0_shortStim','dFF0_longStim','Bleachingratio_total','TgtAUC_sig_NormaltoRandT',...
'Areableach','TgtBrightnessMeanBleach','RefBrightnessMeanBleach','GRratio_pixelMeanBleach',...
'GRratio_meanMeanBleach','Bleachingratio_totalBleach','TgtAUCNormaltoRandTBleach',...
'Kinetics','Plate','DI_short','DI_long','TraceRawMean','TraceRawStdfov','TraceNormalizedMean',...
'TraceNormalizedStdfov','TraceMeanbleach','TraceStdbleach','AveSmallStimFinalMean','AveSmallStimFinalStdfov','LongStimFinalMean','LongStimFinalStdfov'};
selpp2P = {'Legal Name','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','TgtBrightnessFinalMean','TgtBleachingRatioMean','AveTraceGRratioMean',...
'dFF0_shortStim','dFF0_longStim','Bleachingratio_total','Bleachingratio_2_1s',...
'TgtAUC_sig_NormaltoRandT','Kinetics','DI_short','DI_long','TraceRawMean','TraceRawStdfov','TraceNormalizedMean',...
'TraceNormalizedStdfov','AveSmallStimFinalMean','AveSmallStimFinalStdfov','LongStimFinalMean','LongStimFinalStdfov'};
for i = 1:height(T_well_1P)
    idx = find(strcmp(T_well_1P.ConstructName{i},LookUp.ConstructName));
    if ~isempty(idx)
        T_well_1P.ConstructName{i} = LookUp.StandardName{idx};
    end
end
for i = 1:height(T_well_2P)
    idx = find(strcmp(T_well_2P.ConstructName{i},LookUp.ConstructName));
    if ~isempty(idx)
        T_well_2P.ConstructName{i} = LookUp.StandardName{idx};
    end
end
T_well_1P.Kinetics = T_well_1P.dFF0_shortStim./T_well_1P.dFF0_longStim;
T_well_1P.DI_short = abs(T_well_1P.dFF0_shortStim) .* sqrt(T_well_1P.GRratioMean);
T_well_1P.DI_long = abs(T_well_1P.dFF0_longStim) .* sqrt(T_well_1P.GRratioMean);
T_well_1P.('Legal Name') = strcat(T_well_1P.ConstructName,' (',T_well_1P.Plate,T_well_1P.wellName_T_well,')'); 
T_well_2P.Kinetics = T_well_2P.dFF0_shortStim./T_well_2P.dFF0_longStim;
T_well_2P.DI_short = abs(T_well_2P.dFF0_shortStim) .* sqrt(T_well_2P.GRratioMean);
T_well_2P.DI_long = abs(T_well_2P.dFF0_longStim) .* sqrt(T_well_2P.GRratioMean);
T_well_2P.('Legal Name') = strcat(T_well_2P.ConstructName,' (',T_well_2P.Plate,T_well_2P.wellName,')'); 
OnePhoton = T_well_1P(:,selpp1P);
TwoPhoton = T_well_2P(:,selpp2P);
Screening = outerjoin(OnePhoton,TwoPhoton,'LeftKeys','Legal Name','RightKeys','Legal Name','MergeKeys',true);
Screening = sortrows(Screening,'Legal Name','ascend');
Screening = libplot.tableparse(Screening,'parse_on','Legal Name','mask_on','Area_OnePhoton');
Screening.plibName = [];
Screening.Properties.VariableNames{find(strcmp(Screening.Properties.VariableNames,'ConstructName'))} = 'plibName';

save(fullfile(savepath,strcat(saveName,'.mat')),'Screening');

%% Scatter by library
close all
Screening = rmmissing(Screening,'DataVariables',{'plibName'});
namegroups = splitapply(@unique, Screening.plibName, findgroups(Screening.plibName));
nangroups = find(contains(Screening.('Legal Name'),'cyOFP'));
sel = reshape(cellstr(namegroups(setdiff(1:length(namegroups),nangroups))),1,[]);
dtS = {'Area_OnePhoton','GRratioMean_OnePhoton','dFF0_shortStim_OnePhoton',...
    'dFF0_longStim_OnePhoton','Bleachingratio_total_OnePhoton','GRratioMean_TwoPhoton',...
    'dFF0_shortStim_TwoPhoton','dFF0_longStim_TwoPhoton','Bleachingratio_total_TwoPhoton'};
cmap = colormap(colorcube(30));
% 1P DI scatter
libplot.scatterbylib(Screening,'DI_long_OnePhoton',...
    'DI_short_OnePhoton','controls',{'JEDI-1P'},'selectLib',sel,'norm2ctrl',...
    {'JEDI-1P'},'datatip',dtS,'colormap',cmap);
xlabel({'1P Detectability Index 100 ms','normalized to ASAP2s'})
ylabel({'1P Detectability Index 1 ms','normalized to ASAP2s'})
axis square

% 1P dFF0 scatter
libplot.scatterbylib(Screening,'dFF0_longStim_OnePhoton',...
    'dFF0_shortStim_OnePhoton','controls',{'JEDI-1P'},'selectLib',sel,'datatip',dtS,'colormap',cmap);
xlabel({'1P dFF0 100 ms','normalized to ASAP2s'})
ylabel({'1P dFF0 1 ms','normalized to ASAP2s'})
axis square

% 2P DI scatter
libplot.scatterbylib(Screening,'DI_long_TwoPhoton',...
    'DI_short_TwoPhoton','controls',{'JEDI-1P'},'selectLib',sel,...
    'norm2ctrl',{'JEDI-1P'},'datatip',dtS,'colormap',cmap);
xlabel({'2P Detectability Index 100 ms','normalized to ASAP2s'})
ylabel({'2P Detectability Index 1 ms','normalized to ASAP2s'})
axis square

% DI-AUC scatter
libplot.scatterbylib(Screening,'DI_short_OnePhoton',...
    'TgtAUC_sig_NormaltoRandT_OnePhoton','controls',{'JEDI-1P'},'selectLib',sel,'norm2ctrl',...
    {'JEDI-1P'},'datatip',dtS,'colormap',cmap);
xlabel({'1P Detectability Index 1 ms','normalized to ASAP2s'})
ylabel({'TgtAUC_sig_NormaltoRandT_OnePhoton','normalized to ASAP2s'})
axis square

libplot.scatterbylib(Screening,'GRratioMean_OnePhoton',...
    'Bleachingratio_total_OnePhoton','controls',{'ASAP2s','JEDI-1P'},'selectLib',...
    sel,'norm2ctrl',{'JEDI-1P'},'datatip',dtS,'colormap',cmap);
libplot.scatterbylib(Screening,'GRratioMean_OnePhoton',...
    'Bleachingratio_total_OnePhoton','controls',{'ASAP2s','JEDI-1P'},'selectLib',...
    sel,'norm2ctrl',{'JEDI-1P'},'datatip',dtS,'colormap',cmap);
%% 1D scatter by lib
[~,ng] = libplot.scatterbylib1d(Screening,'dFF0_longStim_OnePhoton',...
    'controls',{'JEDI-1P'},'selectLib',sel,'datatip',dtS);
xlbls = sel(ng);
set(gcf, 'Position',  [110, 100, 500, 500])
set(gca, 'XTickLabel',xlbls, 'XTick',1:numel(sel));
ax = gca;
ax.XTickLabelRotation = 45;