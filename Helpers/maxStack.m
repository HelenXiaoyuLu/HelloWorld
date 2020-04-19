% This function take a stack of nd2 images on the same fov and generate a tiff
% file out of it, keeping all saturated pixels.
function maxStack(flist,figureOn)
    mkdir(strcat(flist(1).folder,'\maxStack'));
    numfiles = length(flist);
    imgMax = nd2.read(strcat(flist(1).folder,'\',flist(1).name)); 
    for k = 2:numfiles 
        img = nd2.read(strcat(flist(k).folder,'\',flist(k).name)); 
        imgMax = max(img,imgMax);
    end
%     refC = uint16(ref-1);
%     refC = imadjust(refC);
%     imgdSatNorm = imhistmatch(imgdSat,refC,65536); 
    fname = strsplit(flist(1).name,'.');
    imgStackMax16 = uint16(imgMax-1);
    imwrite(imgStackMax16,strcat(flist(1).folder,'\maxStack\',fname{1},'_maxStack.tiff'));      
    if figureOn
        figure()
        imshow(imgStackMax16,[])
    end
end
