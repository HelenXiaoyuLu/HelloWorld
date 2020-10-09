%% Power ramp up 
% Sep 18th, 2020

%% Load files
clear
dirp = 'C:\Users\13801\OneDrive\Documents\Thorlabs\Optical Power Monitor';
flist = dir(fullfile(dirp,'20200918_Power_ramping_during_warm_up*.csv'));
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
more1 = abs(PowerSpec) > 1;
yDiff(more1) = 1;
less1 = abs(PowerSpec) < 1;
yDiff(less1) = 0;
yDiffDiff = diff(yDiff);
yTurnOFF = find(yDiffDiff == -1);
yTurnON = find(yDiffDiff == 1);
peakvaluest = zeros(1,length(yTurnOFF));
peakLoc = zeros(1,length(yTurnOFF));
for i = 1:length(yTurnOFF)
     peakvaluest(i) = max(PowerSpec(yTurnON(i):yTurnOFF(i)));
     yDiff(yTurnON(i):yTurnOFF(i)) = yDiff(yTurnON(i):yTurnOFF(i))*peakvaluest(i);
end

%% Overlay
figure(1)
hold on
plot(timeline, yDiff,'LineWidth',2);
title('power spectra')

%% Plot againts time
figure(2)
clf 
hold on
tp = timeline(yTurnON);
plot(tp,peakvaluest,'+-')
% ylim([0,max(peakvaluest)])
xlabel("Time of the day")
ylabel({'920nm excitation';'power at sample plane [mW]'})
