%% Generate manual mask for mono channel
% saturated pixel will be discarded at this step

%% Phase 1 mCherry: discard saturated pixels, match histogram and save as tiff for later masking
% Set path: parental and subfolders
% ImgsDir: list of all imgs to be processed and converted
clear 
Parent = 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking';
subJobs = '\20190227_benchmarking_Merge';
subfNames = '\*.nd2';
ImgsDir = dir(strcat(Parent,subJobs,subfNames));
refName = 'ASAP2s_P1D1_1-6.nd2';
maxVal = 4095;
refImg = uint16(max(squeeze(io.nd2.read(fullfile(ImgsDir(1).folder,refName),'t',0,'ch',1)),[],3)-1);
trdSat(ImgsDir,1,maxVal,refImg,false)

%% Phase 1 cyOFP Mono channel: discard saturated pixels, match histogram and save as tiff for later masking
% Set path: parental and subfolders
% ImgsDir: list of all imgs to be processed and converted
clear 
Parent = 'E:\Images\Eric\repeat_result_GEVI\200126 Fstim platform characterization\plate1';
subJobs = '\2p Fstim video';
subfNames = '\*.nd2';
ImgsDir = dir(strcat(Parent,subJobs,subfNames));
refName = 'ASAP1-cyOFP_P1E3_1-4.nd2';
maxVal = 4095;
refImg = max(squeeze(io.nd2.read(fullfile(ImgsDir(1).folder,refName),'t',0,'ch',1)),[],3);
refSat = find(refImg == maxVal);
refImg(refSat) = min(refImg(:));
refImg = uint16(refImg-1);
trdSat(ImgsDir,1,maxVal,refImg,false)

%% Phase 1 cyOFP: discard saturated pixels, match histogram and save as tiff for later masking
% % Set path: parental and subfolders
% % ImgsDir: list of all imgs to be processed and converted
% clear 
% Parent = 'E:\Images\Eric\repeat_result_GEVI\20200202 Fstim platform characterization 10x vs 20x\20x plate2';
% tgtFolder = '\2p Fstim video';
% refFolder = '\2p image';
% subfNames = '\*.nd2';
% tgtDir = dir(strcat(Parent,tgtFolder,subfNames));
% refDir = dir(strcat(Parent,refFolder,subfNames));
% refName = 'ASAP1-cyOFP.1_P2F3_1-1.nd2';
% maxVal = 4095;
% refImg = max(squeeze(io.nd2.read(fullfile(tgtDir(1).folder,refName),'t',0,'ch',1)),[],3);
% refSat = find(refImg == maxVal);
% refImg(refSat) = min(refImg(:));
% refImg = uint16(refImg-1);
% trdSatDual(tgtDir,refDir,maxVal,refImg,true)
% 
%% Phase 1 no histmatch
dirp = "D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201008 double blind\Benchmarking repeat\2pFstim data";
tgtlist = dir(fullfile(dirp, "\video\*.nd2"));
reflist = dir(fullfile(dirp, "\image\*.nd2"));
mkdir(strcat(dirp,'\2p 440hz DualCh ExSat'));
maxVal = 4095;
for i = 1:numel(tgtlist)
    tgtTr = squeeze(io.nd2.read(fullfile(tgtlist(i).folder,tgtlist(i).name),'t',1:20,'ch',1)); 
    refTr = squeeze(io.nd2.read(fullfile(reflist(i).folder,reflist(i).name),'t',1:20,'ch',1)); 
    tgtMax = max(tgtTr, [], 3);
    refMax = max(refTr, [], 3);
    idxSat = (tgtMax == maxVal) & (refMax == maxVal);
    tgtMap = mean(tgtTr, 3);
    tgtMap(idxSat) = 0;
    tgtdSat16 = uint16(tgtMap-1);
    fname = strsplit(tgtlist(i).name,'.');
    imwrite(tgtdSat16,strcat(dirp,'\2p 440hz DualCh ExSat\',fname{1},'.tiff'));
end

%% Draw mask using ui.cellSel
% dirp = strcat(Parent, subJobs);
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201008 double blind\Benchmarking repeat\2pFstim data';
flist = dir(strcat(dirp,'\2p 440hz DualCh ExSat\*.tiff'));
mkdir(strcat(dirp,'\fgmaskTiff'));
% mkdir(strcat(dirp,'\maskMat'));
merge = 'label';
%%
% ui ROI select
for k = 1:length(flist)
    I = imread(strcat(flist(k).folder,'\',flist(k).name));
    [~,ImgName,~] = fileparts(flist(k).name)
    disp(strcat('k = ',num2str(k),'; ImgName = ',ImgName));
    [mask, ROIs, ~] = ui.cellSel(I);  
    % Save tiff
    Savep_tif = strcat(dirp,'\fgmaskTiff\',ImgName,'.tiff');
    imwrite(mask,Savep_tif)
%     % Save mat
%     Savep_mat = strcat(dirp,'\maskMat\',ImgName,'.mat')
%     save(Savep_mat,'mask','merge','ROIs', '-v7.3');
    closereq
end
close('force')
%% Updated ui.segmenter
% I = imread(strcat(flist(1).folder,'\',flist(1).name));
% s = ui.segmenter('I', I);
for k = 56:length(flist)
    I = imread(strcat(flist(k).folder,'\',flist(k).name));
    [~,ImgName,~] = fileparts(flist(k).name)
    disp(strcat('k = ',num2str(k),'; ImgName = ',ImgName));
    s = ui.segmenter('I', I);
    waitfor(s, 'finalized', true);
    Savep_tif = strcat(dirp,'\maskTiff\',ImgName,'.tiff');
    s.export(Savep_tif);
    Savep_mat = strcat(dirp,'\maskMat\',ImgName,'.mat');
    s.export(Savep_mat);
    s.clear();
    delete(s);
end


%% Mono channel trace discard saturated pixel & match histogram
function trdSat(flist,tgtch,maxVal,refC,figureOn)
    mkdir(strcat(flist(1).folder,'\ExSat'));
    numfiles = length(flist);
    for i = 1:numfiles
        tr = squeeze(io.nd2.read(fullfile(flist(i).folder,flist(i).name),'t',0,'ch',tgtch)); 
        imgMax = max(tr,[],3);
        idxSat = find(imgMax == maxVal);
        imgMax(idxSat) = 0; % Replace the saturated pixel with a 0
        imgdSat16 = imhistmatch(uint16(imgMax-1),refC,maxVal);
        fname = strsplit(flist(i).name,'.nd2');
        imwrite(imgdSat16,strcat(flist(i).folder,'\ExSat\',fname{1},'.tiff'));   
    end   
    if figureOn
        figure()
        imshow(imgdSat16,[])
    end
end

%% Dual channel trace discard saturated pixel & match histogram
function trdSatDual(tgtlist,reflist,maxVal,refC,figureOn)
    mkdir(strcat(tgtlist(1).folder,'\ExSat'));
    numfiles = length(tgtlist);
    for i = 1:numfiles
        tgtTr = squeeze(io.nd2.read(fullfile(tgtlist(i).folder,tgtlist(i).name),'t',0,'ch',1)); 
        refTr = squeeze(io.nd2.read(fullfile(reflist(i).folder,reflist(i).name),'t',0,'ch',2)); 
        tgtMax = max(tgtTr,[],3);
        refMax = max(refTr,[],3);
        idxSat = find(tgtMax == maxVal) && find(refMax == maxVal);
        tgtMax(idxSat) = 0; % Replace the saturated pixel with a 0
        refMax(idxSat) = 0; % Replace the saturated pixel with a 0
        tgtdSat16 = imhistmatch(uint16(tgtMax-1),refC,maxVal);
        refdSat16 = imhistmatch(uint16(refMax-1),refC,maxVal);
        M = 0.5*(im2double(tgtdSat16)+im2double(refdSat16));
        fname = strsplit(flist(i).name,'.');
        imwrite(M,strcat(flist(i).folder,'\ExSat\',fname{1},'.tiff'));
    end

end