I = imread('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200227 Spectrum\1080-700-1080 BFP LP out\maxStack\HEK_EGFP_5-1_maxStack.tiff');
figure(1)
hold on
subplot(2,2,1)
imshow(I,[])

subplot(2,2,2)
BW1 = edge(I,'sobel');
imshow(BW1,[])


subplot(2,2,3)
BW2 = edge(I,'canny',0.1);
imshow(BW2,[])

se = offsetstrel('ball',3,3);
I2 = imerode(I,se);
imshow(I2,[])
level = graythresh(I2);
BW = im2bw(I2,level);
D = -bwdist(~BW);
D(~BW) = -Inf;
L = watershed(D);
imshow(L,[])

L = watershed(I);
imgcp = I;
imgcp(L == 0) = 0;
imshow(imgcp,[])