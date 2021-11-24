close all
clear
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\KCR2-mCherry (was KCR1 instead of KCR1)';
atfpath = '20211005_KCR1-mCherry_HEK293A_UltraII_1000ex_40x_23C_powerramp_001.atf';
T = readtable(fullfile(path, atfpath), 'FileType', 'text', ...
        'VariableNamingRule', 'preserve');
path_command = '\\Stpierrelab7910\e\Images\Xiaoyu\20211004 KCR preliminary test\Vcmd 20211004190634.csv';
T_cmd = readtable(path_command);
figure(), plot(T.("Signals="), T.Im_prime);
xlabel("Time (s)");
ylabel("Current (pA)")
xlim([10, 500])
title("920nm -> 800nm : 40nm: 1080nm -> 920nm")
title("Power ramp at 1000 nm")
load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\20211005_KCR_PowerAtSamplePlane_1000nmPowerRamp_dwell2pt7_PFSoff_S170_64x64summary')

%% find events
ydiff = diff(T_cmd.CommandVoltage_a_u__);
y1 = find(ydiff > 4);
y2 = find(ydiff < -4);
x1 = T_cmd.Time_s_(y1);
x2 = T_cmd.Time_s_(y2);
deltaT = mode(diff(x1(1:end-1)));
eventdur = mode(x2 - x1);
figure(), plot(T.("Signals=")-25.2558, T.Im_prime);
hold on, plot(T_cmd.Time_s_-33.5125, T_cmd.CommandVoltage_a_u__*10)

%% Find events alternative
t_raw = T.("Signals=");
Ismooth = movmean(T.Im_prime, 500);
[~, locs] = findpeaks(Ismooth, "MinPeakDistance", 100000, 'MinPeakHeight', 30)
event_on = t_raw(locs);

%% Chop
percentpower = [8, 1, 2, 4, 6, 8, 12, 16, 20, 8];
actualpower = powerarray(2, :);
x0 = event_on(1);
n_event = 10;
movwin = 100;
s_full = sig.arbitrary('x', T.("Signals="), 'y', T.Im_prime);
stim = sig.square('xOn', 0, 'dxEvent', 1.8, 'yEvent', 1);
xshifts = []; % manual
for i = 1:n_event
    eventstart = event_on(i) -1.5;
    eventend = event_on(i) + 3;
    s(i) = s_full.crop(eventstart, eventend);
    s(i).xShift(-eventstart);
    [~, peakloc] = max(diff(movmean(s(i).y, movwin)));
    s(i).xShift(-s(i).x(peakloc));
%         s(i) = s(i).crop(-0.5, 2);
end
figure(), hold on,
cmap = jet(10);
for i = 2:9
    plot(s(i).x, s(i).y, 'Color', cmap(i, :), 'DisplayName', percentpower(i) + "%", 'LineWidth', 1.5)
end
xlabel('Time (s)')
ylabel('Current (pA)')
legend('Location', 'eastoutside')
xlim([-1, 2])


%% normalize
baseidx = 1:500;
movwin = 100;
slp = zeros(1, numel(s));
for i = 1:numel(s)
    ssmth(i) = sig.arbitrary('x', s(i).x, 'y', movmean(s(i).y, movwin));
    snorm(i) = sig.arbitrary('x', ssmth(i).x, 'y', ssmth(i).y - mean(ssmth(i).y(baseidx)));
    slp(i) = max(diff(ssmth(i).y));
end
figure(), hold on,
for i = 2:9
    plot(ssmth(i).x, ssmth(i).y, 'Color', cmap(i, :), 'DisplayName', percentpower(i) + "%", 'LineWidth', 1.5)
end
xlabel('Time (s)')
ylabel('Current (pA)')
legend('Location', 'eastoutside')
xlim([-1, 2])

%% Plot power ramp
figure()
plot(actualpower, slp);
figure('Color', [1,1,1]), hold on
plot(actualpower(2:end-1), slp(2:end-1)./max(slp), 'LineWidth', 1.5);
scatter(actualpower, slp./max(slp), 'LineWidth', 1.5);
xlabel("Power (mW)")
ylabel("Slope normalized")

%% Plot spectra
figure()
plot([920, 800:40:1080, 920], slp);

figure('Color', [1,1,1]), hold on
plot([800:40:1080], slp(2:end-1)./max(slp), 'LineWidth', 1.5);
scatter([920, 800:40:1080, 920], slp./max(slp), 'LineWidth', 1.5);
xlabel("Wavelength (nm)")
ylabel("Slope normalized")

%%
% wl = 800:40:1080;
wl = actualpower(2:end-1);
T_sum = table();
T_sum.PowerPercent = wl';
T_sum.Slope = slp(2:end-1)';
T_sum.Raw = s(2:end-1)';
T_sum.Smooth = ssmth(2:end-1)';
T_sum.Norm = snorm(2:end-1)';
T_sum.Power = T_power.Power';

%%
save(fullfile('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results',...
    '20211005_KCR1-mCherry_HEK293A_UltraII_1000ex_40x_23C_powerramp_001'), 'T_sum')

%% 
Tsum(1) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results\20211005_KCR1-mCherry_HEK293A_UltraII_40x_23C_2pt5mW_Spectra40nmspacing_001');
Tsum(2) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results\20211005_KCR1-mCherry_HEK293A_UltraII_40x_23C_2pt5mW_Spectra40nmspacing_002');
Tsum(3) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results\20211005_KCR1-mCherry_HEK293A_UltraII_40x_23C_2pt5mW_Spectra40nmspacing_003');

%% 
figure(), hold on
plot(T_power.Wavelength, T_power.Power, 'LineWidth', 1.5)
xlabel("Wavelength (nm)")
ylabel("Power (mW)")
ylim([0, 3])
%%
figure('Color', [1,1,1]), hold on
plot(Tsum(1).T_sum.Wavelength, Tsum(1).T_sum.Slope./max(Tsum(1).T_sum.Slope))
plot(Tsum(2).T_sum.Wavelength, Tsum(2).T_sum.Slope./max(Tsum(2).T_sum.Slope))
plot(Tsum(3).T_sum.Wavelength, Tsum(3).T_sum.Slope./max(Tsum(3).T_sum.Slope))
xlabel("Wavelength (nm)")
ylabel("Slope normalized")

figure('Color', [1,1,1]), hold on
meanslope = mean([Tsum(1).T_sum.Slope./max(Tsum(1).T_sum.Slope),Tsum(2).T_sum.Slope./max(Tsum(2).T_sum.Slope), Tsum(3).T_sum.Slope./max(Tsum(3).T_sum.Slope)], 2);
stdslope = std([Tsum(1).T_sum.Slope./max(Tsum(1).T_sum.Slope),Tsum(2).T_sum.Slope./max(Tsum(2).T_sum.Slope), Tsum(3).T_sum.Slope./max(Tsum(3).T_sum.Slope)]');
plot(Tsum(1).T_sum.Wavelength, meanslope)
errorbar(Tsum(1).T_sum.Wavelength',meanslope, stdslope)
xlabel("Wavelength (nm)")
ylabel("Slope normalized")

%% 
Tsumramp(1) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results\20211005_KCR1-mCherry_HEK293A_UltraII_1000ex_40x_23C_powerramp_001');
Tsumramp(2) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results\20211005_KCR1-mCherry_HEK293A_UltraII_1000ex_40x_23C_powerramp_002');
figure('Color', [1,1,1]), hold on
plot(Tsumramp(1).T_sum.PowerPercent, Tsumramp(1).T_sum.Slope./max(Tsumramp(1).T_sum.Slope))
plot(Tsumramp(2).T_sum.PowerPercent, Tsumramp(2).T_sum.Slope./max(Tsumramp(2).T_sum.Slope))
% plot(Tsum(3).T_sum.Wavelength, Tsum(3).T_sum.Slope./max(Tsum(3).T_sum.Slope))
xlabel("Power (mW)")
ylabel("Slope normalized")

figure('Color', [1,1,1]), hold on
meanslope = mean([Tsumramp(1).T_sum.Slope./max(Tsumramp(1).T_sum.Slope),Tsumramp(2).T_sum.Slope./max(Tsumramp(2).T_sum.Slope)], 2);
stdslope = std([Tsumramp(1).T_sum.Slope./max(Tsumramp(1).T_sum.Slope),Tsumramp(2).T_sum.Slope./max(Tsumramp(2).T_sum.Slope)]');
plot(Tsumramp(1).T_sum.PowerPercent, meanslope)
errorbar(Tsumramp(1).T_sum.PowerPercent',meanslope, stdslope)
xlabel("Power (mW)")
ylabel("Slope normalized")