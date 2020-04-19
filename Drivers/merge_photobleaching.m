%% Batch merge
% Get a ref channel for all. I took a mCherry channel from a control
% All channels will be normalized to this LUT 
close all
clear
dirp = 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate2_1P_Brightness_photobleaching\';
tracelist = dir(strcat(dirp,'Trace','\*.nd2'));
brightnesslist = dir(strcat(dirp,'Brightness','\*.nd2'))
imgTar = nd2.read(strcat(tracelist(1).folder,'\',tracelist(1).name));
imgTar = imgTar(:,:,1);
imgRef = nd2.read(strcat(brightnesslist(1).folder,'\',brightnesslist(1).name));
imgTar16 = uint16(imgTar-1);
imgRef16 = uint16(imgRef-1);
refC = imadjust(imgRef16);
Anorm = imhistmatch(imgTar16,refC,65536); 
Bnorm = imhistmatch(imgRef16,refC,65536); % match green LUT to red LUT
M = 0.5*(im2double(Anorm)+im2double(Bnorm));
%%
% Batch merge
dirpath = 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate2_1P_Brightness_photobleaching';
mkdir(strcat(dirpath,'\Merge_adj_redo'));
numfiles = length(tracelist);
for k = 1:numfiles 
        imgTar = nd2.read(strcat(tracelist(k).folder,'\',tracelist(k).name));
        imgTar = imgTar(:,:,1);
        imgRef = nd2.read(strcat(brightnesslist(k).folder,'\',brightnesslist(k).name));
        imgTar16 = uint16(imgTar-1);
        imgRef16 = uint16(imgRef-1);
        Anorm = imhistmatch(imgTar16,refC,65536); 
        Bnorm = imhistmatch(imgRef16,refC,65536); % match green LUT to red LUT
        M = 0.5*(im2double(Anorm)+im2double(Bnorm));
        fname = strsplit(tracelist(k).name,'.');
        imwrite(M,strcat(dirpath,'\Merge_adj_redo\',fname{1},'.tiff'));
end
%%