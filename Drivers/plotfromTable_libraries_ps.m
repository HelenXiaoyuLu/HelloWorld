%% Read and combine data from tables
clc;
clear;
savepath = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\tempMat20200510\201912_photostablelibrary';
saveName = '201912_photostablelibrary';
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

%% Add kinetics and group rows by name 

selpp1P = {'Legal Name','wellName_T_well','ConstructName','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','dFF0_shortStim','dFF0_longStim','Bleachingratio_total','TgtAUC_sig_NormaltoRandT',...
'Areableach','TgtBrightnessMeanBleach','RefBrightnessMeanBleach','GRratio_pixelMeanBleach',...
'GRratio_meanMeanBleach','Bleachingratio_totalBleach','TgtAUCNormaltoRandTBleach',...
'Kinetics','Plate','DI_short','DI_long','TraceRawMean','TraceRawStdfov','TraceNormalizedMean',...
'TraceNormalizedStdfov','AveSmallStimFinalMean','AveSmallStimFinalStdfov','LongStimFinalMean','LongStimFinalStdfov'};
selpp2P = {'Legal Name','Area','TgtBrightnessMean','RefBrightnessMean',...
'GRratioMean','TgtBrightnessFinalMean','TgtBleachingRatioMean','AveTraceGRratioMean',...
'dFF0_shortStim','dFF0_longStim','Bleachingratio_total','Bleachingratio_2_1s',...
'TgtAUC_sig_NormaltoRandT','Kinetics','DI_short','DI_long','TraceRawMean','TraceRawStdfov','TraceNormalizedMean',...
'TraceNormalizedStdfov','AveSmallStimFinalMean','AveSmallStimFinalStdfov','LongStimFinalMean','LongStimFinalStdfov'};
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
Screening = rmmissing(Screening,'DataVariables',{'Legal Name'});
Screening = libplot.tableparse(Screening,'parse_on','Legal Name','mask_on','Area_OnePhoton');
Screening.plibName = [];
Screening.Properties.VariableNames{find(strcmp(Screening.Properties.VariableNames,'ConstructName'))} = 'plibName';
Screening.Properties.RowNames = Screening.('Legal Name');
save(fullfile(savepath,strcat(saveName,'_screening_withTrace.mat')),'Screening');

%% Scatter by library
Screening = rmmissing(Screening,'DataVariables',{'plibName'});

dtS = {'Area_OnePhoton','GRratioMean_OnePhoton','dFF0_shortStim_OnePhoton',...
    'dFF0_longStim_OnePhoton','Bleachingratio_total_OnePhoton','GRratioMean_TwoPhoton',...
    'dFF0_shortStim_TwoPhoton','dFF0_longStim_TwoPhoton','Bleachingratio_total_TwoPhoton'};
cmap = colormap(colorcube(30));
libplot.scatterbylib(Screening,'DI_long_OnePhoton',...
    'DI_short_OnePhoton','controls',{'JEDI-1P'},'datatip',dtS,'colormap',cmap);
libplot.scatterbylib(Screening,'DI_long_TwoPhoton',...
    'DI_short_TwoPhoton','controls',{'JEDI-1P'},'datatip',dtS,'colormap',cmap);
libplot.scatterbylib(Screening,'GRratioMean_OnePhoton',...
    'Bleachingratio_total_OnePhoton','controls',{'JEDI1P'},'datatip',dtS,'colormap',cmap);
libplot.scatterbylib(Screening,'GRratioMean_TwoPhoton',...
    'Bleachingratio_total_TwoPhoton','controls',{'JEDI1P'},'datatip',dtS,'colormap',cmap);


%% Plot single traces
