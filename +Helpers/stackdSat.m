% This function take a stack of nd2 images on the same fov and generate a tiff
% file out of it, excluding all saturated pixels.
function stackdSat(flist,figureOn)
    mkdir(strcat(flist(1).folder,'\ExSat'));
    numfiles = length(flist);
    imgMax = nd2.read(strcat(flist(1).folder,'\',flist(1).name)); 
    for k = 2:numfiles 
        img = nd2.read(strcat(flist(k).folder,'\',flist(k).name)); 
        imgMax = max(img,imgMax);
    end
    idxSat = find(imgMax == 4095);
    imgdSat = imgMax;
    imgdSat(idxSat) = min(imgdSat(:));
%     refC = uint16(ref-1);
%     refC = imadjust(refC);
%     imgdSatNorm = imhistmatch(imgdSat,refC,65536); 
    fname = strsplit(flist(1).name,'.');
    imgdSat16 = uint16(imgdSat-1);
    imwrite(imgdSat16,strcat(flist(1).folder,'\ExSat\',fname{1},'_dSat.tiff'));   
%     imwrite(imgdSatNorm,strcat(flist(1).folder,'\ExSat\',fname{1},'_dSatNorm.tiff'));   
    if figureOn
        figure()
%         subplot(1,2,1)
        imshow(imgdSat16,[])
%         subplot(1,2,2)
%         imshow(imgdSatNorm)
    end
end
