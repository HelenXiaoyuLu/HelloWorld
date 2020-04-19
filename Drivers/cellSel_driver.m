%% cellSel ui driver
% Get merged files
clear
dirp = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200128 Spectra\HV20 Merge\ExSat';
flist = dir(strcat(dirp,'\*.tiff'));
mkdir(strcat(dirp,'\maskTiff'));
mkdir(strcat(dirp,'\maskMat'));
merge = 'label';
% ui ROI select
for k = 1:length(flist)
    I = imread(strcat(flist(k).folder,'\',flist(k).name));
    [~,ImgName,~] = fileparts(flist(k).name)
    disp(strcat('k = ',num2str(k),'; ImgName = ',ImgName));
    [mask, ROIs, ~] = ui.cellSel(I);  
    % Save tiff
    Savep_tif = strcat(dirp,'\maskTiff\',ImgName,'.tiff');
    imwrite(mask,Savep_tif)
    % Save mat
    Savep_mat = strcat(dirp,'\maskMat\',ImgName,'.mat')
    save(Savep_mat,'mask','merge','ROIs', '-v7.3');
    closereq
end

%% Add ROI on saved mat
% ui ROI select
clear
% load established mat file
dirp = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200215 mGold spectrum\700-1080-700\ExSat';
flist = dir(strcat(dirp,'\*.tiff'));
tifflist = dir(strcat(dirp,'\maskMat\*.mat'));
mkdir(strcat(dirp,'\newmaskTiff'));
mkdir(strcat(dirp,'\newmaskMat'));
merge = 'label';
% ui ROI select
for k = 1:length(flist)
    I = imread(strcat(flist(k).folder,'\',flist(k).name));
    load(fullfile(tifflist(k).folder, tifflist(k).name)); 
    [~,ImgName,~] = fileparts(flist(k).name);
    disp(strcat('k = ',num2str(k),'; ImgName = ',ImgName));
    [mask, ROIs, ~] =  ui.cellSel(I,'initROI', ROIs)
    % Save tiff
    Savep_tif = strcat(dirp,'\newmaskTiff\',ImgName,'.tiff');
    imwrite(mask,Savep_tif)
    % Save mat
    Savep_mat = strcat(dirp,'\newmaskMat\',ImgName,'.mat');
    save(Savep_mat,'mask','merge','ROIs', '-v7.3');
    closereq
end