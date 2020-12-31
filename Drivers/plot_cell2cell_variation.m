%% Single cell variation
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking\plate2';
j0.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask foreground no thresholding'));
j1.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask foreground'));
j2.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask background'));

%% Load T_roi
j = j0;
T_roi = j.output.T_roi;
T_fov = j.output.T_fov;
T_roi = innerjoin(T_roi, T_fov, ...
    'LeftKeys', 'ROI_Parent',...
    'RightKeys','FOV_UID',...
    'RightVariables', 'FOV_Name');
grouping = findgroups(T_roi.FOV_Name(:,1));
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
t.Title.String = 'Channel intensity distribution';
for groupn = 1:max(grouping)
    ax = nexttile(t);
    hold(ax,'on')
    groupidx = find(grouping == groupn);
    s = scatter(ax, T_roi.R(groupidx), T_roi.G(groupidx), 5, cmap(5,:),'filled',...
        'DisplayName','corrFilt gating v2');  
    ax.Title.String = T_roi.FOV_Name(groupidx(1));
    ax.XAxis.Label.String = 'Red Intensity';
    ax.YAxis.Label.String = 'Green Intensity';
    ax.XAxis.Limits = [0,4096];
    ax.YAxis.Limits = [0,4096];
end

%%
dirp ='D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking\plate2\maskTiffilastik';
img = imread(fullfile(dirp, 'JEDI-1P-cyOFP_P2A3_1-2.tiff'));
figure(),imshow(double(img).*(img == 1))