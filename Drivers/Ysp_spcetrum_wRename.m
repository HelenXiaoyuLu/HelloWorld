clear
%clear import
import ysp.*
import spcore.*

    cfg.default = db.jobConfig(...
        'path', 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200904 ramping power', ...
        'tgtCh\path', 'Power ramping 920nm', ...
        'lookup', '*.xlsx');
%% Rename Files 2P
lib = libData('0129Rename');
Lookup = readcell('','FileType','spreadsheet');
Jobs = 'HV20 Merge - Copy';
fNames = '\20200128_195856_004_NDExp_*.nd2';
Parent = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200128 Spectra\';
lib.loadFrame(...
    'path',strcat(Jobs,...
    );
wells = cell(3,length(lib.getChildren('w')));

for i = 1:length(lib.getChildren('w'))
    wells{1,i} = lib.getChildren('w',i).wellName();
    wells{2,i} = Lookup{find(strcmp(Lookup, wells{1,i})),2};
    S = lib.getChildren('w',i).Name();
    wells{3,i} = S(1:end-8);
end
% Get Original flist
% loop for wells
updateCount = 1;
updateList = cell(0, 2); 
for i = 1:length(lib.getChildren('w'))
    groupNames = strcat('\20200128_195856_004_NDExp_Well',sprintf('%04d',i-1),'_Point0000*.nd2');
    subDir = strcat(Parent,Jobs,groupNames);
    dirtable = dir(subDir);
    subcount = 1;
    for j = 1:length(dirtable)
        [oldPath, oldfName, fext] = fileparts(strcat(dirtable(j).folder,'\',dirtable(j).name));
        varName = wells{2,find(strcmp(wells(3,:), oldfName(1:end-8)))};
        newfName = strcat(varName,'_',num2str(i),'-',num2str(subcount));
        subcount = subcount + 1;
        updateList{updateCount, 1} = fullfile(oldPath, [newfName, fext]);
        updateList{updateCount, 2} = fullfile(oldPath, [oldfName, fext]);
        movefile(updateList{updateCount, 2}, updateList{updateCount, 1}, 'f');
        updateCount = updateCount + 1;
    end
    % Now generate the tiff file for 'maximum value'
    subDirStruct = dir(strcat(oldPath,'\',varName,'_',num2str(i),'*.nd2'));
    stackdSat(subDirStruct,false);
    maxStack(subDirStruct,false);
end
% Alternative: for renaming only
% for i = 1:length(lib.getChildren('f'))
%     dirtable = dir(jobDir);
%     subcount = 1;
%     f = lib.getChildren('f',i);
%     for j = 1:length(dirtable)
%         [oldPath, oldfName, fext] = fileparts(f.frame(j).path);
%         varName = wells{2,1};
%         newfName = strcat(varName,'_',num2str(i),'-',num2str(subcount));
%         subcount = subcount + 1;
%         updateList{updateCount, 1} = fullfile(oldPath, [newfName, fext]);
%         updateList{updateCount, 2} = fullfile(oldPath, [oldfName, fext]);
%         movefile(updateList{updateCount, 2}, updateList{updateCount, 1}, 'f');
%         updateCount = updateCount + 1;
%     end
% end
%% After renaming
% Reload dataset
clear
%clear import
import ysp.*
import spcore.*
lib = ysp.libData('GEVI set');
Jobs = 'HV20 Merge';
fNames = '\*.nd2';
Parent = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200128 Spectra';
lib.protocol.path = fullfile(Parent,Jobs);
lib.loadSet('path', '');
lib.regroup('method','name');      % regroup wells into constructs
numFrame = lib.getChildren('f',1).nt;
% Segmentation method1: Discard saturated pixels and generate mask for each fov
% lib.getChildren('fov').setSegmentFun(ysp.fovData.foregroundSegment('ch', 1, 't', 1:numFrame)); 
% Segmentation method2: Deploy manual segmentation
lib.getChildren('fov').setSegmentFun(fovData.importSegment(...
    'path', 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200128 Spectra\HV20 Merge\ExSat\maskTiff', ...
    'matchingRule', '<filename>*', ...
    'reader', @imread));
% Set background correction 
% channel = lib.protocol.getChannel(1);
% channel.setBackground('method', 'custom','bg',50);
lib.getChildren('fov').segment();               % perform segmentation
lib.getChildren('fov').populateROI();           % create ROIs

%% Define Metric
workspace = lib.protocol.workspace; % workspace
workspace.deleteMetric()
workspace.addMetric('roi r');
workspace.addMetric("F0(r) = r.getMeanIntensity('ch', 1,'t',1)", ...
    struct('Info', 'Absolute intensity', ...
           'Unit', 'a.u.'));
% workspace.addMetric("F0(r) = r.getMeanIntensity('ch', 1,'t',1:g.getChildren('r',1).nt)-r.getParent().getChildren('r',1).getMeanIntensity('ch',1,'t',1:g.getChildren('r',1).nt)", ...
%     struct('Info', 'Absolute intensity', ...
%            'Unit', 'a.u.'));
workspace.addMetric("group g");
% workspace.addMetric("meanB(g) = reshape(nanmean(g.getChildren('r',1: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt), 1), [], 1)");
% workspace.addMetric("stdB(g) = reshape(nanstd(g.getChildren('r',2: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt), 0, 1), [], 1)");
% workspace.addMetric("stdE(g) = reshape(nanstd(g.getChildren('r',2: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt), 0, 1), [], 1)./ sqrt(-1+length(g.getChildren('r')))");
workspace.addMetric("meanB(g) = reshape(nanmean(g.getChildren('r',2: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt)-g.getChildren('r',1).getMeanIntensity('ch',1,'t',1:g.getChildren('r',1).nt), 1), [], 1)");
workspace.addMetric("stdB(g) = reshape(nanstd(g.getChildren('r',2: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt)-g.getChildren('r',1).getMeanIntensity('ch',1,'t',1:g.getChildren('r',1).nt), 0, 1), [], 1)");
workspace.addMetric("stdE(g) = reshape(nanstd(g.getChildren('r',2: numel(g.getChildren('r'))).getMeanIntensity('ch', 1, 't', 1:g.getChildren('r',1).nt)-g.getChildren('r',1).getMeanIntensity('ch',1,'t',1:g.getChildren('r',1).nt), 0, 1), [], 1)./ sqrt(length(g.getChildren('r')))");
metricsToCompute = workspace{'meanB','stdB','stdE'};
% metricsToCompute = workspace{'meanB'};

results = struct();
g = lib.getChildren('g');
for i = 1:numel(g)
    for j = 1:numel(metricsToCompute)
        results(i).(metricsToCompute(j).Name) = metricsToCompute(j).eval(g(i));
    end
end

lib.visualize()
Variants = lib.getChildren('g',6).getChildren('w',1).getChildren('r',4);
Variants = lib.getChildren('r');
Variants.plotIntensity('method','frame')
%% Plot Mono Direction
% load the Citrine/EGFP spectra from fpbase
mCitrinesp = readmatrix('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\fpbase_spectra_mCitrine');
EGFPsp = readmatrix('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\fpbase_spectra_EGFP');
% load the actual wavelength-power relationships as a matrix
load('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\20200215_peakvaluest_jobs_700-20-1080.mat'); 
% Derive the wavelength-power relationships for this specific library
% waveLens = [900,700:20:1080,900];
powers = zeros(1, length(waveLens));
for i = 1:numFrame
    powers(i) = peakvaluest(2,find(peakvaluest(1,:)== waveLens(i)));
end
figure(1)
plot(peakvaluest(1,:),peakvaluest(2,:))
title('Power vs wavelength')
ylabel('Power / mW')
xlabel('WaveLength / nm')
varsel = [11:16];
legends = cell(1,length(varsel));
% EGFPspsparse = peakvaluest;
% for i = 1:length(peakvaluest)
%     EGFPspsparse(2,i) = EGFPsp(find(EGFPsp== EGFPspsparse(1,i)),4);
% end

% figure 1: Plot directly (no background correction, no ref, no normalize, no power correction)
figure()
hold on
for i = varsel %1:numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    errorbar(waveLens(2:end-1),(results(i).meanB(2:end-1))',(results(i).stdE(2:end-1)'));
end
legend(legends(varsel))
%title({'no background correction';'no power correction';'no normalization'})
xlabel('wavelength')
ylabel('fluorescence (a.u.)')
title({'background correction = TRUE ';'no power correction';'no normalization'})

% figure 2: no background correction, no ref, no normalize, power correction
figure()
hold on
for i = varsel %numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    plot(waveLens(2:end-1),(results(i).meanB(2:end-1))'./powers(2:end-1));
end
legend(legends(varsel))
%title({'no background correction';'power correction = TRUE ';'no normalization'})
xlabel('wavelength')
ylabel('fluorescence (a.u.)')
title({'background correction = TRUE ';'power correction = TRUE ';'no normalization'})

% figure 3: Normalize all spectra value to its maximum 
figure()
hold on
for i = varsel
    legends{i} = lib.getChildren('g',i).Name();
    plot(waveLens(2:end-1),normalize((results(i).meanB(2:end-1))'./powers(2:end-1),'range'));
end
plot(EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),1),EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),4),'Color',[0.4660, 0.6740, 0.1880])
plot(EGFPsp(find(EGFPsp==700:20:1080),1),EGFPsp(find(EGFPsp==700:20:1080),4),'Color',[0.4660, 0.6740, 0.1880])
legend([legends(varsel),'EGFP(fpbase)']);

%title({'no background correction';'power correction = TRUE ';'normalization = TRUE'})
xlabel('wavelength')
ylabel('Normalized fluorescence')
title({'background correction = TRUE ';'power correction = TRUE ';'normalization = TRUE'})
% Save library
save('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\Results\HV20_spectrum_lib.mat', 'lib', '-v7.3');       

%% Define Dual Direction
% load the EGFP spectra from fpbase
EGFPsp = readmatrix('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\fpbase_spectra_EGFP');
mCitrinesp = readmatrix('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\fpbase_spectra_mCitrine');
% load the actual wavelength-power relationships as a matrix
load('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200227 Spectrum\20200227_JOBs_peakvaluest.mat'); 
% Derive the wavelength-power relationships for this specific library
precision = 20; % spacing between each two sampled wavelengths
ControlWaveLength = 900; % if control images is taken 
% waveLens = [900,700:precision:1080,900, flip(700:precision:1080), 900]; %700-1080-700_20nm
waveLens = [900,flip(700:precision:1080),900, 700:precision:1080, 900]; %700-1080-700_20nm
ControlIdx = find(waveLens(1,:) == ControlWaveLength);
ForSt = ControlIdx(1) + 1;
ForEnd = ControlIdx(3) - 1;
RevSt = ControlIdx(3) + 1;
RevEnd = ControlIdx(5) - 1;
powers = zeros(1, numFrame);
for i = 1:numFrame
    powers(i) = peakvaluest(2,find(peakvaluest(1,:)== waveLens(i)));
end
varsel = [1:6];
legends = cell(1,length(varsel));
%% Generate figures
cmap = lines;
figure(2) 
clf
% subplot 1: w/ background correction, w/o power correction, w/o normalize, 
subplot(1,3,1)
hold on
for i = varsel %1:numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    errorbar(waveLens(ForSt:ForEnd),(results(i).meanB(ForSt:ForEnd))',(results(i).stdE(ForSt:ForEnd)'),'Color',cmap(i+1,:),'LineStyle','-');
end
for i = varsel %numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    errorbar(waveLens(RevSt:RevEnd),(results(i).meanB(RevSt:RevEnd))',(results(i).stdE(RevSt:RevEnd)'),'Color',cmap(i+1,:),'LineStyle','--');
end
legend(legends(varsel),'location','best')
%title({'no background correction';'no power correction';'no normalization'})
xlabel('2P laser wavelength /nm')
ylabel('fluorescence (a.u.)')
title({'background correction = TRUE ';'power correction = FALSE';'normalization = FALSE'})
axis([700 1080 0 inf])

% subplot 2: w/ background correction, w/ power correction, w/o normalize, 
figure(2)
subplot(1,3,2)
hold on
for i = varsel %numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    plot(waveLens(ForSt:ForEnd),(results(i).meanB(ForSt:ForEnd))'./(powers(ForSt:ForEnd).^2),'Color',cmap(i+1,:));  
end
for i = varsel %numel(g)
    legends{i} = lib.getChildren('g',i).Name();
    plot(waveLens(RevSt:RevEnd),(results(i).meanB(RevSt:RevEnd))'./(powers(RevSt:RevEnd).^2),'LineStyle','--','Color',cmap(i+1,:));
end
legend(legends(varsel),'location','best')
xlabel('2P laser wavelength /nm')
ylabel('fluorescence (a.u.)')
title({'background correction = TRUE ';'power correction = TRUE ';'normalization = FALSE'})
axis([700 1080 0 inf])

% subplot 3: w/ background correction, w/ power correction, w normalize
figure(2)
subplot(1,3,3)
hold on
for i = varsel
    legends{i} = lib.getChildren('g',i).Name();
    pCintenseFor = (results(i).meanB(ForSt:ForEnd))'./(powers(ForSt:ForEnd).^2);
    plot(waveLens(ForSt:ForEnd),pCintenseFor./max(pCintenseFor),'LineStyle','-','Color',cmap(i+1,:));
end
% plot(mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),1),mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),4),'Color',[0.4660, 0.6740, 0.1880])
% plot(EGFPspsparse(1,:),EGFPspsparse(2,:),'Color',[0.4660, 0.6740, 0.1880])
plot(EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),1),EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),4),'Color',[0.4660, 0.6740, 0.1880])
for i = varsel
    legends{i} = lib.getChildren('g',i).Name();
    pCintenseRev = (results(i).meanB(RevSt:RevEnd))'./(powers(RevSt:RevEnd).^2);
    plot(waveLens(RevSt:RevEnd),normalize((results(i).meanB(RevSt:RevEnd))'./(powers(RevSt:RevEnd).^2),'range'),'LineStyle','--','Color',cmap(i+1,:));
end
legend([legends(varsel),'EGFP (fpbase)'],'location','best')
axis([700 1080 0 1])
xlabel('2P laser wavelength /nm')
ylabel('Normalized fluorescence')
title({'background correction = TRUE ';'power correction = TRUE ';'normalization = TRUE'})
set(gcf,'Position',[100 100 1800 400])

% figure 4 Average back and forth
figure(3)
clf
hold on
% for i = varsel
%     legends{i} = lib.getChildren('g',i).Name();
%     plot(waveLens(2:end-1),normalize((results(i).meanB(2:end-1))'./(powers(2:end-1).^2),'range'),'+','Color',cmap(i+11,:));
% end
%plot(mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),1),normalize(mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),4),'range'),'Color',cmap(1,:))
for i = varsel
    legends{i} = lib.getChildren('g',i).Name();
    Intense = zeros(1,length(peakvaluest));
    for j = 1:length(peakvaluest)
        Intense(j) = mean((results(i).meanB(find(waveLens == peakvaluest(1,j)))));
    end
    powerCorrect  = Intense./(peakvaluest(2,:).^2);
    plot(peakvaluest(1,:),powerCorrect./max(powerCorrect),'LineStyle','-','Color',cmap(i+1,:));
end
% plot(EGFPspsparse(1,:),EGFPspsparse(2,:),'Color',[0.4660, 0.6740, 0.1880])
plot(EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),1),EGFPsp(find(EGFPsp==700):find(EGFPsp==1080),4),'Color',[0.4660, 0.6740, 0.1880])
% plot(mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),1),mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),4),'Color',[0.4660, 0.6740, 0.1880])
legend([legends(varsel),'mCitrine (fpbase)'],'location','best')
title({'background correction = TRUE';'power correction = TRUE ';'normalization = TRUE'})
xlabel('2P laser wavelength /nm')
ylabel('Normalized fluorescence')
axis([700 1080 0 1])
% % no normalization
% figure()
% hold on
% for i = varsel
%     legends{i} = lib.getChildren('g',i).Name();
%     plot(waveLens(2:end-1),(results(i).meanB(2:end-1))'./(powers(2:end-1).^2),'+','Color',cmap(i+12,:));
% end
% %plot(mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),1),21*mCitrinesp(find(mCitrinesp==700):find(mCitrinesp==1080),4),'Color',cmap(2,:))
% for i = varsel
%     legends{i} = lib.getChildren('g',i).Name();
%     Intense = zeros(1,length(peakvaluest));
%     for j = 1:length(peakvaluest)
%         Intense(j) = mean((results(i).meanB(find(waveLens == peakvaluest(1,j)))))
%     end
%     plot(peakvaluest(1,:),Intense./peakvaluest(2,:).^2,'LineStyle','-','Color',cmap(i+14,:));
% end
% legend([legends(varsel),'mCitrine (fpbase)'])
% legend('mGold','mGold average','mVenus','mVenus average','mCitrine (fpbase)')
% title({'no background correction';'power correction = TRUE ';'normalization = TRUE'})
% xlabel('wavelength')
% ylabel('Normalized fluorescence')
%%
save('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200227 Spectrum\Results\20200227_Jobs_mGold_700-1080-700_spectrum_lib.mat', 'lib', '-v7.3');       
load('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\Results\mGold700-1080-700_spectrum_ws.mat')

%% Merge Existing Figures
f2=hgload('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Results\20200306 Manual mGold 10nM norm2range10mW.fig');
f1=hgload('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Results\20200306 Manual mGold mVenus 10nM norm2range10mW average only.fig');
figure()
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);
copyobj(allchild(get(f1,'CurrentAxes')),h(1));
copyobj(allchild(get(f2,'CurrentAxes')),h(2));
subplot(1,2,1)
legend('mGold (datapoints)','mCitrine (fpbase)', 'mGold (average)','location','best')
xlabel('2P laser wavelength /nm')
ylabel('Normalized fluorescence')
title('mGold vs mCitrine')
subplot(1,2,2)
Children = get(gca,'children');
legend('mGold','mCitrine (fpbase)', 'mGold average','location','best')
xlabel('2P laser wavelength /nm')
ylabel('Normalized fluorescence')
title('mGold vs mVenus')
% Children = get(gca,'children')
% delete(Children(3));
set(gcf,'Position',[100 100 1000 400])
sgtitle('background correction = TRUE, power correction = TRUE, normalization = TRUE')