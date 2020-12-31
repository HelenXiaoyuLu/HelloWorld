clear;
clc;
maindir = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking\plate2';
uncent = 'maskTiffilastik_strict_uncentainty';
mask = 'maskTiffilastik_strict';
dir_uncent = dir(fullfile(maindir, uncent, '*.tiff'));
dir_mask = dir(fullfile(maindir, mask, '*.tiff'));

%% Merge files
for i = 1 : length(dir_mask)
    fname = dir_mask(i).name;
    msk = imread(fullfile(dir_mask(i).folder, fname));
    uncertainty = imread(fullfile(dir_uncent(i).folder, fname));
    mskfg = (msk == 1) & (uncentainty < 0.05); % "membrane"
    mskbg = (msk == 2) & (uncentainty < 0.01); % "background"
    mskclp = (msk == 3) & (uncentainty < 0.05); % "clump"
    imwrite(mskfg,fullfile(maindir, 'fgmask', fname));
    imwrite(mskbg,fullfile(maindir, 'bgmask', fname));
    imwrite(mskclp,fullfile(maindir, 'clpmask', fname));
end

%%    