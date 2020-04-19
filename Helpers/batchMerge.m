%% function batch merge 
% dirpath = 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness';
   
function batchMerge(dirpath,refC)
    mkdir(strcat(dirpath,'\Merge_adj'))
    flist = dir(strcat(dirpath,'\*.nd2'));
    numfiles = length(flist);
    for k = 1:numfiles 
        img = nd2.read(strcat(flist(k).folder,'\',flist(k).name)); 
        img16 = uint16(img-1);
        A = img16(:,:,1); % green channel
        B = img16(:,:,2); % red channel
        Anorm = imhistmatch(A,refC,65536); 
        Bnorm = imhistmatch(B,refC,65536); % match green LUT to red LUT
        M = 0.5*(im2double(Anorm)+im2double(Bnorm));
%         M = im2int16(M1)
%         C = rgb2gray(imfuse(img(:,:,1),img(:,:,2),'Scaling','independent'));
        fname = strsplit(flist(k).name,'.');
        imwrite(M,strcat(flist(k).folder,'\Merge_adj\',fname{1},'.tiff'));
    end
end
