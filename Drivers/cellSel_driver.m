%% cellSel ui driver
% Get merged files
clear
dirp = 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20200628 photoactivation\exp50ms';
flist = dir(strcat(dirp,'\Tiff\*.tiff'));
mkdir(strcat(dirp,'\maskTiff'));
mkdir(strcat(dirp,'\maskMat'));
merge = 'label';
%% ui ROI select
for k = 2:length(flist)
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
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200918 Power ramping during warm up';
flist = dir(strcat(dirp,'\maskTiff\*.tiff'));
tifflist = dir(strcat(dirp,'\maskMat\*.mat'));
mkdir(strcat(dirp,'\newmaskTiff'));
mkdir(strcat(dirp,'\newmaskMat'));
merge = 'label';
%% ui ROI select
for k = 1:1
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