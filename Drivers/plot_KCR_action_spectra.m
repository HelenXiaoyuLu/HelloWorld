%% Load data
close all
clear
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20211018 KCR2 spectra\KCR2';
fname = "20211018_KCR2-mCherry_HEK293A_UltraII_40x_23C_005";
ttlrecorded = true;
atfpath = fname + ".atf";
T = readtable(fullfile(path, atfpath), 'FileType', 'text', ...
        'VariableNamingRule', 'preserve');
if ttlrecorded
    cmdpath = fname + ".csv";
    path_command = fullfile(path, cmdpath);
    T_cmd = readtable(path_command);
end 
% Load power lookup table
load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\20211005_KCR_PowerAtSamplePlane_2pt5mW_dwell2pt7_PFSoff_S170_64x64_postadjustsummary')
fname = "20211005_KCR1-mCherry_HEK293A_UltraII_40x_23C_2pt5mW_Spectra40nmspacing_003";

%% Plot raw data and command
cmap = lines;
I_raw = T.Im_prime;
t_raw = T.("Signals=");
dt = 0.0001;
for i = 1:(width(T)-1)/2-1
    I_raw = [I_raw; T.("Im_prime_" + i)]; 
    t_raw = [t_raw; t_raw(end) + T.("Signals=") + dt];
end
if ttlrecorded
    scmd = sig.arbitrary('x', T_cmd.Time_s_, 'y', T_cmd.CommandVoltage_a_u__);
    D = img.temporalAlignment(I_raw, t_raw, 'method', 'arbitrary', ...
        'template', scmd, 'figshow', true);
    delay = rem(D, 51);
    figure(), hold on;
    plot(scmd.x, scmd.y, 'DisplayName', "Light ON", 'Color', cmap(1,:));
    plot(t_raw + delay, I_raw./max(I_raw)*5, 'DisplayName', "Light ON", 'Color', cmap(2,:));
    xlabel("Time (s)");
    ylabel("Current (pA)")
    title("960nm -> 8       00nm : 40nm: 1080nm -> 960nm")
    legend('Location', 'southoutside')
    % Find events
    event_on = scmd.x(diff(scmd.y)>2)-delay;
else
    figure(), hold on;
    plot(t_raw, I_raw./max(I_raw)*5, 'DisplayName', "Light ON", 'Color', cmap(2,:));
    xlabel("Time (s)");
    ylabel("Current (pA)")
    title("960nm -> 800nm : 40nm: 1080nm -> 960nm")
    legend('Location', 'southoutside');     
    Ismooth = movmean(I_raw, 500);
    findpeaks(Ismooth, "MinPeakDistance", 100000, 'MinPeakHeight', 20, 'MinPeakProminence', 15, 'MinPeakWidth', 0.5)
    [~, locs] = findpeaks(Ismooth, "MinPeakDistance", 100000, 'MinPeakHeight', 20, 'MinPeakProminence', 15);
    event_on = t_raw(locs);
end

%% Chop
spectra = [920, 800:40:1080, 920];
actualpower = powerarray(2, :);
x0 = event_on(1);
n_event = 10;
movwin = 100;
s_full = sig.arbitrary('x', t_raw, 'y', I_raw);
% stim = sig.square('xOn', 0, 'dxEvent', 1.8, 'yEvent', 1);
xshifts = []; % manual
for i = 1:n_event
    eventstart = event_on(i)-1.5;
    eventend = event_on(i) + 4;
    s(i) = s_full.crop(eventstart, eventend);
    s(i).xShift(-event_on(i));
    [~, peakloc] = max(diff(movmean(s(i).y, movwin)));
    s(i).xShift(-s(i).x(peakloc));
    s(i) = s(i).crop(-2, 2);
end
figure(), hold on,
cmap = jet(10);
for i = 2:9
    plot(s(i).x, s(i).y, 'Color', cmap(i, :), 'DisplayName', spectra(i) + "nm", 'LineWidth', 1.5)
end
xlabel('Time (s)')
ylabel('Current (pA)')
legend('Location', 'eastoutside')
xlim([-inf, inf])


%% normalize
baseidx = 1:5000;
movwin = 10;
derivmax = zeros(1, numel(s));
eventslp = zeros(1, numel(s));
eventmax = zeros(1, numel(s));
eventplateau = zeros(1, numel(s));
figure(), hold on, 
r = cell(1, numel(s));
rflip = cell(1, numel(s));
for i = 1:numel(s)
    ssmth(i) = sig.arbitrary('x', s(i).x, 'y', movmean(s(i).y, movwin));
    snorm(i) = sig.arbitrary('x', ssmth(i).x, 'y', ssmth(i).y - mean(ssmth(i).y(baseidx)));
    snormflip(i) = sig.arbitrary('x', ssmth(i).x, 'y', flipud(ssmth(i).y - mean(ssmth(i).y(baseidx))));
    derivmax(i) = max(diff(movmean(s(i).y, 500)));
    r{i} = sig.evaluator.response.load(snorm(i).y, snorm(i).x, 'event', sig.square('xOn', 0, 'dxEvent', 1));
    rflip{i} = sig.evaluator.response.load(snormflip(i).y, snormflip(i).x, 'event', sig.square('xOn', 0, 'dxEvent', 1));
    eventmax(i) = r{i}.event.Ff(1);
    eventslp(i) = snorm(i).x2y(r{i}.event.t_f(1))/(r{i}.event.t_f(1) - r{i}.event.t_s(1));
    eventplateau(i) = snormflip(i).x2y(rflip{i}.event.t_f(1));
end
figure(), hold on,
for i = 2:9
    plot(snorm(i).x, snorm(i).y, 'Color', cmap(i, :), 'DisplayName', spectra(i) + "nm", 'LineWidth', 1.5);
end
xlabel('Time (s)')
ylabel('Current (pA)')
legend('Location', 'eastoutside')
title("Current smoothed with a movmean window of " + movwin)
xlim([-inf, inf])

%% Plot spectra
f = figure('Color', [1,1,1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on')
% plot(spectra, slp);
plot(spectra(2:end-1), derivmax(2:end-1)./max(derivmax(2:end-1)), 'LineWidth', 1.5);
scatter(spectra, derivmax./max(derivmax(2:end-1)), 'LineWidth', 1.5);
xlabel("Wavelength (nm)")
ylabel("Max derivative normalized")
title("Max derivative (movmean = 500)")

ax = nexttile(t);
hold(ax, 'on')
% plot(spectra, eventslp);
plot(spectra(2:end-1), eventslp(2:end-1)./max(eventslp(2:end-1)), 'LineWidth', 1.5);
scatter(spectra, eventslp./max(eventslp(2:end-1)), 'LineWidth', 1.5);
xlabel("Wavelength (nm)")
ylabel("Slope normalized")
title("Slope")

ax = nexttile(t);
hold(ax, 'on')
% plot(spectra, eventmax);
plot(spectra(2:end-1), eventmax(2:end-1)./max(eventmax(2:end-1)), 'LineWidth', 1.5);
scatter(spectra, eventmax./max(eventmax(2:end-1)), 'LineWidth', 1.5);
xlabel("Wavelength (nm)")
ylabel("Maximum current normalized")
title("Maximum current")

ax = nexttile(t);
hold(ax, 'on')
% plot(spectra, eventplateau);
plot(spectra(2:end-1), eventplateau(2:end-1)./max(eventplateau(2:end-1)), 'LineWidth', 1.5);
scatter(spectra, eventplateau./max(eventplateau(2:end-1)), 'LineWidth', 1.5);
xlabel("Wavelength (nm)")
ylabel("Plateau current normalized")
title("Plateau current")

%%
T_sum = table();
T_sum.Wavelength = spectra(2:end-1)';
T_sum.("Max derivative") = derivmax(2:end-1)';
T_sum.("Slope") = eventslp(2:end-1)';
T_sum.("Maximum current") = eventmax(2:end-1)';
T_sum.("Plateau current") = eventplateau(2:end-1)';
T_sum.("Trace Raw") = s(2:end-1)';
T_sum.("Trace Smooth") = ssmth(2:end-1)';
T_sum.("Trace Norm") = snorm(2:end-1)';
T_sum.("Trace Power") = T_power.Power';

T_raw = table();
T_raw.Wavelength = powerarray(1, :)';
T_raw.("Max derivative") = derivmax';
T_raw.("Slope") = eventslp';
T_raw.("Maximum current") = eventmax';
T_raw.("Plateau current") = eventplateau';
T_raw.("Trace Raw") = s';
T_raw.("Trace Smooth") = ssmth';
T_raw.("Trace Norm") = snorm';
T_raw.("Trace Power") = powerarray(2,:)';

%%
save(fullfile('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\Results',...
    fname + "_optimized"), 'T_sum', 'T_raw')

%% 
Tsummary(1) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211018 KCR spectra\Results\20211018_KCR2-mCherry_HEK293A_UltraII_40x_23C_002_optimized');
Tsummary(2) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211018 KCR spectra\Results\20211018_KCR2-mCherry_HEK293A_UltraII_40x_23C_003_optimized');
Tsummary(3) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211018 KCR spectra\Results\20211018_KCR2-mCherry_HEK293A_UltraII_40x_23C_004_optimized');
Tsummary(4) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211018 KCR spectra\Results\20211018_KCR2-mCherry_HEK293A_UltraII_40x_23C_005_optimized');

%% 
figure(), hold on
plot(T_power.Wavelength, T_power.Power, 'LineWidth', 1.5)
xlabel("Wavelength (nm)")
ylabel("Power (mW)")

%%
Method = "Plateau current";
figure('Color', [1,1,1]), hold on
wavelength = Tsummary(1).T_sum.Wavelength;
valsnorm = zeros(numel(Tsummary), numel(wavelength));
for i = 1:numel(Tsummary)
    valsnorm(i, :) = Tsummary(i).T_sum.(Method)./max(Tsummary(i).T_sum.(Method));
end
meanslope = mean(valsnorm, 1);
stdslope = std(valsnorm);
errorbar(wavelength, meanslope, stdslope, 'Color', 'k', 'LineWidth', 1.0, 'DisplayName', 'Mean \pm STD')

for i = 1:numel(Tsummary)
    scatter(wavelength, Tsummary(i).T_sum.(Method)./max(Tsummary(i).T_sum.(Method)), ...
        80, '+', 'LineWidth', 1.5, 'DisplayName', sprintf('trial %d', i))
end

xlabel("Wavelength (nm)")
ylabel(Method + " normalized")
legend('location', 'eastoutside');
title('KCR2')
