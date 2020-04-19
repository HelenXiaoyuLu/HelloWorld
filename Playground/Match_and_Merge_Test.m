%% Match LUT and merge dual channel as gray scale
img = nd2.read('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness\JEDI-1P_P1G1_1-1.nd2');
img16 = uint16(img-1);
% A:green B:red
A = img16(:,:,1);
B = img16(:,:,2);
figure(1)
subplot(6,1,1)
imshow(A,[min(img16(:)),max(img16(:))])
title('ch1=green')
subplot(6,1,2)
histogram(A(:),100)
title('ch1=green LUT')
subplot(6,1,3)
imshow(B,[min(img16(:)),max(img16(:))])
title('ch2=red')
subplot(6,1,4)
histogram(B(:),100)
title('ch1=red LUT')
subplot(6,1,5)
J = imhistmatch(A,B,65536);
imshow(J,[min(img16(:)),max(img16(:))])
title('ch1 after correction')
subplot(6,1,6)
histogram(J(:),100)
title('ch1 LUT after correction')

figure()
subplot(3,1,1)
C1 = J+B;
imshow(C1,[])
title('Manual correction')

%imfuse can only output unit8
subplot(3,1,2)
C2 = rgb2gray(imfuse(img(:,:,1),img(:,:,2),'Scaling','independent'));
imshow(C2,[])
title('imfuse correction')

subplot(3,1,3)
imshow(C2-C1,[])
title('show difference')

cdiff = C1-C2;
sum(cdiff(:))
sum(C1(:))
sum(C2(:))

%% imfuse can also normalize the LUT
figure()
x = imfuse(img(:,:,1),img(:,:,2),'Scaling','independent','ColorChannels','green-magenta');
imshow(x)

figure()
subplot(2,1,1)
imshow(C1)
subplot(2,1,2)
C1adj = imadjust(C1);
imshow(C1adj)


%% compare LUT of imadjust
img1 = nd2.read('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness\JEDI-1P_P1G1_1-1.nd2');
img2 = nd2.read('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness\ASAP1_P1C1_1-2.nd2');
img3 = nd2.read('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness\Bongwoori-P6_P1C4_1-3.nd2');
C1 = rgb2gray(imfuse(img1(:,:,1),img1(:,:,2),'Scaling','independent'));
C2 = rgb2gray(imfuse(img2(:,:,1),img2(:,:,2),'Scaling','independent'));
C3 = rgb2gray(imfuse(img3(:,:,1),img3(:,:,2),'Scaling','independent'));
%%
figure(1)
subplot(6,1,1)
imshow(C1)
title('JEDI-1P')
subplot(6,1,2)
histogram(C1)
title('JEDI-1P LUT')
subplot(6,1,3)
imshow(C2) 
title('ASAP1')
subplot(6,1,4)
histogram(C2)
title('ASAP1 LUT')
subplot(6,1,5)
imshow(C3)
title('Bongwoori')
subplot(6,1,6)
histogram(C3)
title('Bongwoori LUT')

%% adjust image
figure(2)
subplot(4,1,1)
imshow(imadjust(C1),[])
subplot(4,1,2)
BW1 = im2bw(imadjust(C1),graythresh(imadjust(C1))); 
imshow(BW1)
C11 = C1*2
subplot(4,1,3)
imshow(C11,[])
subplot(4,1,4)
BW1 = im2bw(C11,graythresh(C11));
imshow(BW1)