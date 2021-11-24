%% Power ramp at wavelength = X nm
% Oct 5th, 2021

%% Load files
clear
dirp = '\\Stpierrelab7910\e\Images\Xiaoyu\20211021 KCR1 spectra';
fname = "20211021_KCR_PowerAtSamplePlane_1040nmPowerRamp_dwell2pt7_PFSoff_S170_64x64_postadjust";
flist = dir(fullfile(dirp,fname + ".csv"));
[~, idx] = sort({flist.date});
flist = flist(idx);
T_data = table();
for i = 1:numel(flist)
    T_partial = readtable(fullfile(flist(i).folder,flist(i).name));
    T_data = [T_data; T_partial]; 
end
PowerSpec = T_data.Power_W_.*10^3;
timeline = T_data.TimeOfDay_hh_mm_ss_fff_;

%% Plot raw data
figure(1)
clf
hold on
plot(timeline,PowerSpec)
xlabel("Time of the day")
ylabel("Power [mW]")

%% Find on-off time and get read
yDiff = PowerSpec;
movmeanwin = 5;
more1 = abs(movmean(PowerSpec, movmeanwin)) > 0.3;
yDiff(more1) = 1;
less1 = abs(movmean(PowerSpec, movmeanwin)) < 0.3;
yDiff(less1) = 0;
yDiffDiff = diff(yDiff);
yTurnOFF = find(yDiffDiff == -1);
yTurnON = find(yDiffDiff == 1);
peakvaluest = zeros(1,length(yTurnOFF));
peakLoc = zeros(1,length(yTurnOFF));
for i = 1:length(yTurnOFF)
     peakvaluest(i) = mean(PowerSpec(yTurnON(i):yTurnOFF(i)));
     yDiff(yTurnON(i):yTurnOFF(i)) = yDiff(yTurnON(i):yTurnOFF(i))*peakvaluest(i);
end
figure(), plot(yDiff)
%% Overlay
figure(1)
hold on
plot(timeline, yDiff,'LineWidth',2);
title('power spectra')

%% Plot againts power
figure(2)
clf 
hold on
x = [1, 2, 4, 6, 8, 12, 16, 20];
plot(tp,peakvaluest(2:end-1),'+-')
% ylim([0,max(peakvaluest)])
xlabel("Percent of power")
ylabel({'1000nm excitation';'power at sample plane [mW]'})

%% Plot againts wavelength
figure(2)
clf 
hold on
x = 800:40:1080;
plot(x,peakvaluest(2:end-1),'+-')
% ylim([0,max(peakvaluest)])
xlabel("Percent of power")
ylabel({'power at sample plane [mW]'})
ylim([0, inf])

%% Save table Spectra
T_power = table();
T_power.Wavelength = x';
T_power.Power = peakvaluest(2:end-1)';
powerarray = [920, 800:40:1080, 920; peakvaluest];
save(fullfile(dirp, fname + "summary"), 'T_power', 'powerarray')

%% Save table power ramp
T_power = table();
T_power.PowerPercent = x;
T_power.Power = peakvaluest(2:end-1);
powerarray = [8, tp, 8; peakvaluest];
save(fullfile(dirp, fname + "summary"))