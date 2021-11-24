close all
clear
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\KCR2-mCherry (was KCR1 instead of KCR1)';
atfpath = '20211005_KCR1-mCherry_HEK293A_UltraII_40x_23C_2pt5mW_Spectra40nmspacing_001.atf';
T = readtable(fullfile(path, atfpath), 'FileType', 'text', ...
        'VariableNamingRule', 'preserve');
path_command = '\\Stpierrelab7910\e\Images\Xiaoyu\20211004 KCR preliminary test\Vcmd 20211004190634.csv';
T_cmd = readtable(path_command);
figure(), subplot(2,1,1), plot(T.("Signals="), T.Im_prime);
subplot(2,1,2), plot(T_cmd.Time_s_, T_cmd.CommandVoltage_a_u__)

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
[~, locs] = findpeaks(Ismooth, "MinPeakDistance", 100000, 'MinPeakHeight', 40);
event_on = t_raw(locs);

%% Chop
x0 = 21;
n_event = 10;
movwin = 100;
t_full = T.("Signals=");
s_full = T.Im_prime;
stim = sig.square('xOn', 0, 'dxEvent', 1.8, 'yEvent', 1);
xshifts = []; % manual
for i = 1:n_event
    [~, eventstart] = min(abs(t_full - x0));
    [~, eventend] = min(abs(t_full - x0 - deltaT));
    s(i) = sig.arbitrary('x', t_full(eventstart:eventend), 'y', s_full(eventstart:eventend));
    s(i).xShift(-x0);
    [~, peakloc] = max(movmean(s(i).y, movwin));
    s(i).xShift(-s(i).x(peakloc));
        s(i) = s(i).crop(-0.5, 2);
    x0 = x0 + deltaT;
end
s.plot()

%% normalize
baseidx = 1:500;
movwin = 100;
slp = zeros(1, numel(s));
for i = 1:numel(s)
    ssmth(i) = sig.arbitrary('x', s(i).x, 'y', movmean(s(i).y, movwin));
    snorm(i) = sig.arbitrary('x', ssmth(i).x, 'y', ssmth(i).y - mean(ssmth(i).y(baseidx)));
    slp(i) = max(diff(ssmth(i).y));
end
snorm.plot()
ssmth.plot()
%% 
figure()
plot(slp./max(slp));

%%
figure()
wl = 800:40:1080;
plot(wl, slp(2:end-1)./max(slp(2:end-1)));

%%
load('\\Stpierrelab7910\e\Images\Xiaoyu\20211005 KCR power ramp\20211005_KCR_PowerAtSamplePlane_2pt5mW_dwell2pt7_PFSoff_S170_64x64_postadjustsummary')
wl = 800:40:1080;
T_sum = table();
T_sum.Wavelength = wl';
T_sum.Slope = slp(2:end-1)';
T_sum.Raw = s(2:end-1)';
T_sum.Smooth = ssmth(2:end-1)';
T_sum.Norm = snorm(2:end-1)';
T_sum.Power = T_power.Power;
T_sum.Corrected = T_sum.Raw./(T_power.Power.^2);

%%
save('\\Stpierrelab7910\e\Images\Xiaoyu\20211004 KCR preliminary test\KCR1_2', 'T_sum')

%% 
Tsum(1) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211004 KCR preliminary test\KCR2_2');
Tsum(2) = load('\\Stpierrelab7910\e\Images\Xiaoyu\20211004 KCR preliminary test\KCR2_1');

%%
figure('Color', [1,1,1]), hold on
plot(Tsum(1).T_sum.Wavelength, Tsum(1).T_sum.Slope./max(Tsum(1).T_sum.Slope))
plot(Tsum(2).T_sum.Wavelength, Tsum(2).T_sum.Slope./max(Tsum(2).T_sum.Slope))
xlabel("Wavelength (nm)")
ylabel("Slope normalized")