%% Init
close all
clear
clc
WaveLengths = 700:10:1080;
PowerPcts = 0:10:100;
PowerReads = zeros(length(WaveLengths),11);
WaveLengths20 = 700:20:1080;
%%
clc
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200928 Power adjust 960\Power spectra\20200928_afteradjust_PFSoff_100mWrange_18mm_920-700-10-1080-920_2*.csv';
flist = dir(dirp);
T_data = table();
for i = 1:length(flist)
    T_partial = readtable(fullfile(flist(i).folder,flist(i).name));
    T_data = [T_data; T_partial];
end
value_num = T_data.Power_W_.*10^3;

%%
[~,data2] = xlsread("D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Power Config\20200305_manual_850-990nm.csv");
% wavelength2 = strsplit(data2{4,:},';');
% wavelength2 = str2num(wavelength2{2});
datavalue2 = string(data2(8:end-20));
value2 = split(datavalue2,';');
value_num2 = double(value2(:,3))*10^3;

%%
[~,data3] = xlsread("D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Power Config\20200305_manual_1000-1080nm.csv");
% wavelength3 = strsplit(data3{4,:},';');
% wavelength3 = str2num(wavelength3{2});
datavalue3 = string(data3(8:end-20));
value3 = split(datavalue3,';');
value_num3 = double(value3(:,3))*10^3;

%%
[~,data4] = xlsread("D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\20200227_test1.csv");
wavelength4 = strsplit(data4{4,:},';');
wavelength4 = str2num(wavelength4{2});
datavalue4 = string(data4(8:end-20));
%datavalue4 = string(data4(8:102200));
value4 = split(datavalue4,';');
value_num4 = double(value4(:,3))*10^3;

%%
%PowerSpec = cat(1,value_num, value_num2, value_num3);
PowerSpec = fillmissing(value_num,'constant',0);
figure(2)
clf
hold on
plot(1:length(PowerSpec),PowerSpec)
xlabel("sample#")
ylabel("Power/mW")

%% 
PowerSpecDiff = diff(PowerSpec);
yDiff = PowerSpec;
more1 = abs(PowerSpec) > 1;
yDiff(more1) = 1;
less1 = abs(PowerSpec) < 1;
yDiff(less1) = 0;
yDiffDiff = diff(yDiff);
yTurnOFF = find(yDiffDiff == -1);
yTurnON = find(yDiffDiff == 1);

%% Get read
peakvaluest = zeros(1,length(yTurnOFF));
peakLoc = zeros(1,length(yTurnOFF));
for i = 1:length(yTurnOFF)
     peakvaluest(i) = max(PowerSpec(yTurnON(i):yTurnOFF(i)));
     yDiff(yTurnON(i):yTurnOFF(i)) = yDiff(yTurnON(i):yTurnOFF(i))*peakvaluest(i);
end
figure(2)
hold on
ylabel("Power/mW")
plot(yDiff,'LineWidth',2);
title('power spectra')

%% compare method
figure(3)
clf
hold on
peakvalm = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\peakvaluest_measure.mat').peakvaluest;
plot(WaveLengths20,peakvalm)
plot(WaveLengths,peakvaluest(2:end-1))
ylabel("Power/mW")
xlabel("Wavelength/nm")
title('manual vs jobs')
legend('Manual','Jobs')
%%
figure(4)
clf
hold on
peakval1 = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Power Config\20200305_peakvaluest_jobs_700-10-1080.mat').peakvaluest;
peakval2 = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\PowerSpectrum\20200305 SpectraFinal\Power Config\20200305_peakvaluest_manual_700-10-1080.mat').peakvaluest;
peakval3 = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\PowerSpectrum\20200227 Spectrum\peakvaluest_test1_0227.mat').peakvaluest;
peakval4 = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\PowerSpectrum\20200904 Spectra power ramp\20200904_peakvaluest_jobs_10mwrange.mat').peakvaluest;
peakval5 = load('D:\OneDrive - Rice University\Francois\Paper\Data\Power measurements\PowerSpectrum\20200904 Spectra power ramp\20200904_peakvaluest_jobs_100mwrange.mat').peakvaluest;
peakval6 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200917 Power measurement\20200917_peakvaluest_jobs_100mwrange_700-1080_10nm').peakvaluest;
peakval7 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200917 Power measurement\20200917_peakvaluest_jobs_100mwrange_1080-700_10nm').peakvaluest;
peakval8 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200918 Power ramping during warm up\20200918_peakvaluest_jobs_autorange_10nm').peakvaluest;
peakval9 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200928 Power adjust 960\Power spectra\20200928_peakvaluest_jobs_range100mW_700-1080_10nm').peakvaluest;
peakval10 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200928 Power adjust 960\Power spectra\20200928_peakvaluest_jobs_range100mW_1080-700_10nm').peakvaluest;
% plot(peakval1(1,:),peakval1(2,:)) % 0305 jobs
% plot(peakval2(1,:),peakval2(2,:)) % 0305 manual
% plot(WaveLengths20,peakval3(2:end-1)) % 0227 jobs
% plot(peakval4(1,:),peakval4(2,:)) % 0904 jobs 10mw
% plot(peakval5(1,:),peakval5(2,:)) % 0904 jobs 100mw
% plot(peakval6(1,:),peakval6(2,:)) % 0917 jobs 100mw
% plot(peakval7(1,:),peakval7(2,:)) % 0917 jobs 100mw
% plot(peakval8(1,:),peakval8(2,:)) % 0917 jobs 100mw
plot(peakval9(1,:),peakval9(2,:)) % 0928 jobs 100mw
plot(peakval10(1,:),peakval10(2,:)) % 0928 jobs 100mw
plot(WaveLengths,peakvaluest(2:end-1)) % Current

ylabel("Power/mW")
xlabel("Wavelength/nm")
title('manual vs jobs')
legend('700-1080','1080-700',...
    'current')
axis([700 1080 9 15 ])
%'0305 jobs','0305 manual','0227 jobs','0904 jobs 10mw',...
% '0917 jobs 100mW 1080-700','0918 jobs autorange',0904 jobs 100mw','0917 jobs 100mW 700-1080',
%% legend('Manual','JOBs','location','best')
peakvaluest = [WaveLengths; peakvaluest(2:end-1)];
save(fullfile(flist(1).folder,'20200928_peakvaluest_jobs_range100mW_700-1080_10nm_2.mat'),'peakvaluest')
%%
[~,data900] = xlsread("20200123 PowerScan 900nm125 126 135 128 125 126.csv");
datavalue900 = string(data900(8:end-20));
value900 = split(datavalue900,';');
value_num900 = double(value900(:,3))*10^3;
figure()
plot(1:length(value_num900),value_num900)
xlabel("sample#")
ylabel("Power/mW")
Square900 = zeros(1,length(value_num900));
more2 = abs(value_num900) > 0.1;
Square900(more2) = 1;
less2 = abs(value_num900) < 0.1;
Square900(less2) = 0;
Diff900 = diff(Square900);
yTurnOFF900 = find(Diff900 == -1);
yTurnON900 = find(Diff900 == 1);
hold on
peakvalues900 = zeros(1,length(yTurnOFF900));
peakLoc900 = zeros(1,length(yTurnOFF900));
for i = 1:length(yTurnOFF900)
     peakvalues900(i) = max(value_num900(yTurnON900(i)+200:yTurnOFF900(i)-200));
     Square900(yTurnON900(i):yTurnOFF900(i)) = Square900(yTurnON900(i):yTurnOFF900(i))*peakvalues900(i);
end
plot(Square900,'LineWidth',1.5);
title({"re-adjust 900nm:", "Power = 1.25, 1.26, 1.35, 1.28, 1.25, 1.26"})
xlim([0,1E+5])

%%  
T_collection = load(fullfile(pwd,'Visual\+FOVbase_2PSpectra\T_power_collection.mat'))
