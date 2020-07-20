%% discard saturated pixels, match histogram and save as tiff --> later masking
% Set path: parental and subfolders
% ImgsDir: list of all imgs to be processed and converted
Parent = 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20200628 photoactivation';
subJobs = '\exp50msTest';
subfNames = '\*.nd2';
ImgsDir = dir(strcat(Parent,subJobs,subfNames));
trdSat(ImgsDir,2,65535,true)

%%
function trdSat(flist,tgtch,maxVal,figureOn)
    mkdir(strcat(flist(1).folder,'\ExSat'));
    numfiles = length(flist);
    for i = 1:numfiles
        tr = squeeze(io.nd2.read(fullfile(flist(i).folder,flist(i).name),'t',0,'ch',tgtch)); 
        imgMax = max(tr,[],3);
        idxSat = find(imgMax == maxVal);
        imgMax(idxSat) = min(imgMax(:)); % Replace the saturated pixel with a bcground value
        imgdSat16 = uint16(imgMax-1);
        fname = strsplit(flist(i).name,'.');
        imwrite(imgdSat16,strcat(flist(i).folder,'\ExSat\',fname{1},'_dSat.tiff'));   
    end   
    if figureOn
        figure()
        imshow(imgdSat16,[])
    end
end