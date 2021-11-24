%% Load recorded current
clear
close all
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20210129 Model cell test';
fname = 'sp10kHz_bessel4kHz_001';
fname2 = 'sp20kHz_bessel10kHz_001';
cmap = lines;
T_temp1 = readtable(fullfile(path, strcat(fname, '.atf')),'FileType','text');
T_temp2 = readtable(fullfile(path, strcat(fname2, '.atf')),'FileType','text');
temp1 = sig.arbitrary('x', T_temp1.Signals_, 'y', T_temp1.Vm_sec);
temp2 = sig.arbitrary('x', T_temp2.Signals_, 'y', T_temp2.Vm_sec);
temp3 = sig.arbitrary('x', T_temp1.Signals_, 'y', T_temp1.Im_prime);
temp4 = sig.arbitrary('x', T_temp2.Signals_, 'y', T_temp2.Im_prime);

f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');
plot(ax, temp1.x, temp1.y, 'Color', cmap(1,:),'DisplayName', 'Sample @ 10kHz');
plot(ax, temp2.x, temp2.y, 'Color', cmap(2,:),'DisplayName', 'Sample @ 20kHz');
ax.XAxis.Limits = [0, 0.05];
ax.Title.String = "Command voltage";
ax.XAxis.Label.String = "Time (s)";
ax.YAxis.Label.String = "Command voltage (scaled)/mV";
legend(ax);
ax = nexttile(t);
hold(ax, 'on');
plot(ax, temp3.x, temp3.y, 'Color', cmap(1,:),'DisplayName', 'Sample @ 10kHz');
plot(ax, temp4.x, temp4.y, 'Color', cmap(2,:),'DisplayName', 'Sample @ 20kHz');
ax.XAxis.Limits = [0, 0.05];
ax.Title.String = "Current";
ax.XAxis.Label.String = "Time (s)";
ax.YAxis.Label.String = "Current/pA";
legend(ax);
f.Position(3:4) = [700, 300];

%%
path = '\\Stpierrelab7910\e\Images\Xiaoyu\20210129 Model cell test';
fname = 'sp10kHz_bessel4kHz_modelcell';
fname2 = 'sp20kHz_bessel4kHz_modelcell';

T_temp1 = readtable(fullfile(path, strcat(fname, '.csv')));
T_temp2 = readtable(fullfile(path, strcat(fname2, '.csv')));
temp1 = sig.arbitrary('x', T_temp1.Time_s_, 'y', T_temp1.CommandVoltage_a_u__);
temp2 = sig.arbitrary('x', T_temp2.Time_s_, 'y', T_temp2.CommandVoltage_a_u__);
temp1 = temp1.crop(15.9,16.1);
temp2 = temp2.crop(15.95,16.15);
temp1.xShift(0.0524);
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');
plot(ax, temp1.x, temp1.y, 'Color', cmap(1,:),'DisplayName', 'Sample @ 10kHz');
plot(ax, temp2.x, temp2.y, 'Color', cmap(2,:),'DisplayName', 'Sample @ 20kHz');
ax.Title.String = "Current";
ax.XAxis.Label.String = "Time (s)";
ax.YAxis.Label.String = "Current (non-scaled)/pA";
legend(ax);
f.Position(3:4) = [350, 300];