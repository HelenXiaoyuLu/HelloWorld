%% Load image
dirp = '\\Stpierrelab7910\e\Images\Xiaoyu\20210507 VADER truncation\VADER Localization';
flist = dir(fullfile(dirp, '*.nd2'));
i = 1;
reader = bfGetReader(fullfile(dirp, flist(1).name));
Iraw = io.nd2.read(reader, 't', 1, 'ch', 0, 'dimOrder', 'TCSZ', 'type', 'uint16');   
I = squeeze(Iraw); % X*Y*Ch
fileMeta = io.nd2.getOptics(reader);
bitDepth = fileMeta.BitDepth;

%% Col