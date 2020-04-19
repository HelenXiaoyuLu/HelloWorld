clear
% clear import
import ysp.*
import spcore.*

%% Phase 1: Load Raw Data
lib = libData('MSP Experiment');
lib.protocol.path = 'D:\OneDrive - rice.edu\Francois\Platfrom\SingleCell\20200304 Fstim single cell mock co transfection\Fstim 1000Hz\Sel';
imgName = '20200304_120636_842_Well0000_Point0000_Seq0001_afterPA.nd2'
lib.loadSet(...
    'path', imgName, ...
    'ch', 1);   % B vs G

lib.regroup();      % regroup wells into constructs
lib.getChildren('fov').setSegmentFun(ysp.fovData.foregroundSegment('ch', 1, 't', 1));

lib.getChildren('fov').setSegmentFun(ysp.fovData.mammalianSegment('ch', 1, 't', 1));
lib.getChildren('fov').segment();
lib.getChildren('fov').populateROI();
lib.visualize()

%%
close all
I = nd2.read(fullfile(lib.protocol.path,imgName).path);
figure()
imshow(I,[])

mask = mammalianSegmentHelper(I,4500)
figure()
imshow(mask,[])

function mask = mammalianSegmentHelper(I,thres)
    % Threshold image - adaptive threshold
    idx = find(I>thres);
    mask = zeros(size(I));
    mask(idx)=I(idx);
    mask = reshape(mask,size(I));
    figure()
    imshow(mask)
    % Dilate mask with disk
    radius = 3;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    mask = imdilate(mask, se);
    
    % Erode mask with disk
    radius = 5;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    mask = imerode(mask, se);
    
    % Fill holes
    mask = imfill(mask, 'holes');
    
    % Clear borders
    mask = imclearborder(mask);
    
    mask = uint16(bwlabel(cleanSegmentation(mask, 100, 500, 0.95)));
end