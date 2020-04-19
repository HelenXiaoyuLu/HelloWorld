clear
% clear import
import ysp.*
import spcore.*

 

%% Phase 1: Load Raw Data
lib = libData('MSP Experiment');
Jobs = '700-1080-700';
fNames = '\20200215_225808_778_NDExp*.nd2';
jobDir = strcat('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200215 mGold spectrum\',Jobs,fNames);
lib.loadSet(...
    'path',jobDir, ...
    'plateFormat', 96);
lib.regroup();      % regroup wells into constructs
lib.getChildren('fov').setSegmentFun(ysp.fovData.foregroundSegment('ch', 1, 't', 1:43));

% lib.getChildren('fov').setSegmentFun(fovData.importSegment(...
%     'path', 'D:\OneDrive - rice.edu\Images\SingleCellMSP\20190227_Benchmarking_plate1_1P_Brightness\maskTiff', ...
%     'matchingRule', '<filename>*', ...
%     'reader', @imread));

lib.calibration(1).setBackground('method', 'custom', 'bg', 50);
lib.getChildren('fov').segment('ncpu', 1);
lib.getChildren('fov').populateROI();

%% Phase 3: Define Scores
%lib.metrics.delete()
lib.metrics.addMetric("roi r");
lib.metrics.addMetric("group g");
lib.metrics.addMetric("refF0(r) = r.getMeanIntensity('ch', 1, 't', 1)", ...
    struct('Info', 'Absolute intensity', ...
           'Unit', 'a.u.', ...
           'bound', [0, 4000], ...
           'bw', 1));
lib.metrics.addMetric("refF0(r) = r.getMeanIntensity('ch', 1, 't', 1)", ...
metricsToCompute = lib.metrics{'refF0'};
resutls = metricsToCompute(1).eval(lib.getChildren('g',1));

%% Phase 4 plot
Variants = lib.getChildren('g');
%plotROI(Variants, 'scores', {'refF0'});
%lib.visualize();
% plot raw image
%lib.getChildren('f', 1).plotImg('ch', 1, 't', 1);
% plot mask
%lib.getChildren('f', 1).plotMask();
r = lib.getChildren('roi');
%f = Variants.getChildren('fov',1);
r.plotIntensity('method','frame')
waveLens = [900,flipud( 700:20:1080 ), 900,700:20:1080,900];
Powers = peakvaluest;
%[10.088 10.473 10.405 10.875 10.336 10.04 9.9558 10.639 10.294 10.719 10.454 10.088 9.8246 9.5222 9.721 9.7423 9.9879 10.446 10.802 10.874 10.652 10.088];
Intensities = zeros(12,22);
legends = cell(1,length(Variants.getChildren('roi')));
for i = 1:length(Variants.getChildren('roi'))
    legends{i} = Variants.getChildren('roi',i).Parent.Name();
    Intensities(i,:) = Variants.getChildren('roi',i).meanIntensity();
end
figure()
hold on 
plot(waveLens(23:end-1), Intensities(23:end-1)./(peakvaluest(2:end-1)-0.5))
title("Normalized to maximum")
legend(legends)
%plotROI(Variants(1),'scores', {'refF0'});

%%
VariantsName = {{'ASAP1'},{'ASAP2s'},{'JEDI-1P'},{'JEDI-2P'},{'EGFP'},{'EGFP-CAAX'}};
VarNameCat = cat(1, VariantsName{:});
IntesWell = [mean(----ntensities([1,6],:),1);mean(Intensities([2,5],:),1);mean(Intensities([9,10],:),1);mean(Intensities([8,11],:),1);mean(Intensities([3,4],:),1);mean(Intensities([7,12],:),1)];
% findobj(H).isvalid
figure()
hold on 
plot(waveLens(2:21), IntesWell(:,2:21)./Powers(2:21))
legend(VarNameCat)
title(strcat(Jobs,' Averaged'));

%% load egfp spectrum
[EGFPspectra,~,~] = xlsread('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\fpbase_spectra_EGFP.csv')
figure(5)
hold on
plot(EGFPspectra(401:781,1),140*EGFPspectra(401:781,4))