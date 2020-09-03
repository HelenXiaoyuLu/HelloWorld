%% Load plate reader data
clear
clc
path = 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Spectrum\20200720_1Pspectra'; 
fName1 = '2020_07_20_Gevi_spectra_trial3_2_EGFPwcontrol.xlsx';
fName2 = '2020_07_20_Gevi_spectra_trial3_1_HEKcells_gain100.xlsx';
fName3 = '2020_07_20_Gevi_spectra_trial3_4_2_Kir21cells_gain120_step1.xlsx';
fLookUp = '20200720Lookup_JEDIs';
TEx1 = readtable(fullfile(path,fName1),'Range','B47:H94','ReadVariableNames',true);
TEx2 = [readtable(fullfile(path,fName2),'Range','B47:F94','ReadVariableNames',true),...
    readtable(fullfile(path,fName2),'Range','K47:L94','ReadVariableNames',true)];
TEx3 = [readtable(fullfile(path,fName3),'Range','B47:B188','ReadVariableNames',true),...
    readtable(fullfile(path,fName3),'Range','G47:L188','ReadVariableNames',true)];
TEm1 = readtable(fullfile(path,fName1),'Range','B98:H162','ReadVariableNames',true);
TEm2 = [readtable(fullfile(path,fName2),'Range','B98:F162','ReadVariableNames',true),...
    readtable(fullfile(path,fName2),'Range','K98:L162','ReadVariableNames',true)];
TEm3 = [readtable(fullfile(path,fName3),'Range','B192:B256','ReadVariableNames',true),...
    readtable(fullfile(path,fName3),'Range','G192:L256','ReadVariableNames',true)];
lookUp = readtable(fullfile(path,fLookUp),'Range','A1:B16','ReadVariableNames',true);
lookUp.Properties.RowNames = lookUp.Well;
refBuffer = 'External';
refCell = 'HEK';

%% Re-arrange data
ExScan1 = TEx2.Wavelength';
EmScan1 = TEm2.Wavelength'; 
Ex1 = table(); 
Ex1.Well = TEx2.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Ex1)
    Ex1.Name{i} = lookUp.Name{Ex1.Well{i}};
    s = sig.arbitrary();
    s.setSignal(ExScan1, TEx2{:,Ex1.Well{i}}');    
    traces(i, 1) = s;
end
Ex1.rawTrace = traces;

Em1 = table();
Em1.Well = TEm2.Properties.VariableNames(2:end)';
traces = sig.arbitrary.empty();
for i = 1:height(Em1)
    Em1.Name{i} = lookUp.Name{Em1.Well{i}};
    s = sig.arbitrary();
    s.setSignal(EmScan1, TEm2{:,Em1.Well{i}}');    
    traces(i, 1) = s;
end
Em1.rawTrace = traces;

%% Group data
bgcrt = @(tr1,tr2) tr1.y - tr2.y;
gEx = table();
Exgroups = findgroups(Ex1.Name);
gEx.Name = splitapply(@unique, Ex1.Name, Exgroups);
gEx.Properties.RowNames = gEx.Name;
gEx.n_well = splitapply(@numel, Ex1.Name, Exgroups);
gEx.rawTrace = splitapply(@mean, Ex1.rawTrace,  Exgroups);
gEx.rawTraceStd = splitapply(@std, Ex1.rawTrace,  Exgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
ctcrttraces = sig.arbitrary.empty();
for i = 1:height(gEx)
    s = sig.arbitrary();
    s.setSignal(ExScan1, gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary();
    s1.setSignal(ExScan1, normalize(gEx.rawTrace(i).y-gEx.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    s2 = sig.arbitrary();
    name = split(gEx.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2.setSignal(ExScan1, normalize(gEx.rawTrace(i).y-gEx.rawTrace(refCell).y,'range'));
    else
        s2.setSignal(ExScan1, normalize(gEx.rawTrace(i).y-gEx.rawTrace(name{1}).y,'range'));
    end
    ctcrttraces(i, 1) = s2;
end
gEx.('background corrected traces') = traces;
gEx.('background corrected traces normalized') = normtraces;
gEx.('celltype corrected Traces normalized') = ctcrttraces;

gEm = table();
Emgroups = findgroups(Em1.Name);
gEm.Name = splitapply(@unique, Em1.Name, Emgroups);
gEm.Properties.RowNames = gEm.Name;
gEm.n_well = splitapply(@numel, Em1.Name, Emgroups);
gEm.rawTrace = splitapply(@mean, Em1.rawTrace,  Emgroups);
gEm.rawTraceStd = splitapply(@std, Em1.rawTrace,  Emgroups);
traces = sig.arbitrary.empty();
normtraces = sig.arbitrary.empty();
ctcrttraces = sig.arbitrary.empty();
for i = 1:height(gEm)
    s = sig.arbitrary();
    s.setSignal(EmScan, gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y);    
    traces(i, 1) = s;
    s1 = sig.arbitrary();
    s1.setSignal(EmScan, normalize(gEm.rawTrace(i).y-gEm.rawTrace(refBuffer).y,'range')); 
    normtraces(i, 1) = s1;
    s2 = sig.arbitrary();
    name = split(gEm.Name{i},'_');
    if sum(strcmp(name{1},lookUp.Name)) == 0
        s2.setSignal(EmScan, normalize(gEm.rawTrace(i).y-gEm.rawTrace(refCell).y,'range'));
    else
        s2.setSignal(EmScan, normalize(gEm.rawTrace(i).y-gEm.rawTrace(name{1}).y,'range'));
    end
    ctcrttraces(i, 1) = s2;    
end
gEm.('background corrected traces') = traces;
gEm.('background corrected traces normalized') = normtraces;
gEm.('celltype corrected Traces normalized') = ctcrttraces;

%% Save
save(fullfile(path,'1Pspectra_from_platereader.mat'),'gEm','gEx','Ex','Em','lookUp');

%%
d1 = 1; d2 = 4;
cmap = colormap(lines);
colNames = {'background corrected traces','background corrected traces normalized','celltype corrected Traces normalized'};
for i = 1:3
colName = colNames{i};
figure(i)
sgtitle(colName)
% subplot(d1,d2,1)
% hold on
% plot(gEx.(colName)('Kir_EGFP').x,gEx.(colName)('Kir_EGFP').y,'b-')
% plot(gEx.(colName)('HEK_EGFP').x,gEx.(colName)('HEK_EGFP').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_EGFP').x,gEm.(colName)('Kir_EGFP').y,'g-')
% plot(gEm.(colName)('HEK_EGFP').x,gEm.(colName)('HEK_EGFP').y,'color',cmap(5,:))
% title('EGFP')
% legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')
% 
% subplot(d1,d2,2)
% hold on
% plot(gEx.(colName)('Kir_EGFP-CAAX').x,gEx.(colName)('Kir_EGFP-CAAX').y,'b-')
% plot(gEx.(colName)('HEK_EGFP-CAAX').x,gEx.(colName)('HEK_EGFP-CAAX').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_EGFP-CAAX').x,gEm.(colName)('Kir_EGFP-CAAX').y,'g-')
% plot(gEm.(colName)('HEK_EGFP-CAAX').x,gEm.(colName)('HEK_EGFP-CAAX').y,'color',cmap(5,:))
% title('EGFP-CAAX')
% legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')

subplot(d1,d2,1)
hold on
% plot(gEx.(colName)('Kir_ASAP1').x,gEx.(colName)('Kir_ASAP1').y,'b-')
plot(gEx.(colName)('HEK_ASAP1').x,gEx.(colName)('HEK_ASAP1').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_ASAP1').x,gEm.(colName)('Kir_ASAP1').y,'g-')
plot(gEm.(colName)('HEK_ASAP1').x,gEm.(colName)('HEK_ASAP1').y,'color',cmap(5,:))
title('ASAP1')
legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')

subplot(d1,d2,2)
hold on
% plot(gEx.(colName)('Kir_ASAP2s').x,gEx.(colName)('Kir_ASAP2s').y,'b-')
plot(gEx.(colName)('HEK_ASAP2s').x,gEx.(colName)('HEK_ASAP2s').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_ASAP2s').x,gEm.(colName)('Kir_ASAP2s').y,'g-')
plot(gEm.(colName)('HEK_ASAP2s').x,gEm.(colName)('HEK_ASAP2s').y,'color',cmap(5,:))
title('ASAP2s')
legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')

subplot(d1,d2,3)
hold on
% plot(gEx.(colName)('Kir_JEDI-1P').x,gEx.(colName)('Kir_JEDI-1P').y,'b-')
plot(gEx.(colName)('HEK_JEDI-1P').x,gEx.(colName)('HEK_JEDI-1P').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_JEDI-1P').x,gEm.(colName)('Kir_JEDI-1P').y,'g-')
plot(gEm.(colName)('HEK_JEDI-1P').x,gEm.(colName)('HEK_JEDI-1P').y,'color',cmap(5,:))
title('JEDI-1P')
legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')

subplot(d1,d2,4)
hold on
% plot(gEx.(colName)('Kir_JEDI-2P').x,gEx.(colName)('Kir_JEDI-2P').y,'b-')
plot(gEx.(colName)('HEK_JEDI-2P').x,gEx.(colName)('HEK_JEDI-2P').y,'color',cmap(6,:))
% plot(gEm.(colName)('Kir_JEDI-2P').x,gEm.(colName)('Kir_JEDI-2P').y,'g-')
plot(gEm.(colName)('HEK_JEDI-2P').x,gEm.(colName)('HEK_JEDI-2P').y,'color',cmap(5,:))
title('JEDI-2P')
legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')
set(gcf, 'Position',  [50, 10, 1000, 600])

end

%%

figure()
hold on
plot(gEx.(colName)('Kir').x,gEx.(colName)('Kir').y,'b-')
plot(gEx.(colName)('HEK').x,gEx.(colName)('HEK').y,'color',cmap(6,:))
plot(gEm.(colName)('Kir').x,gEm.(colName)('Kir').y,'g-')
plot(gEm.(colName)('HEK').x,gEm.(colName)('HEK').y,'color',cmap(5,:))
legend('Ex-Kir','Ex-HEK','Em-Kir','Em-HEK')
