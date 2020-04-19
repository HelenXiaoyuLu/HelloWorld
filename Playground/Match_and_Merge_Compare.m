%% Batch merge
% Get a ref channel for all. I took a mCherry channel from a control
% All channels will be normalized to this LUT 
close all
clear
img = nd2.read('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness\JEDI-1P_P1G1_1-1.nd2');
img16 = uint16(img-1);
A = img16(:,:,1); % ch1:green
B = img16(:,:,2); % ch2:red
Cadj = imadjust(B); % adjust range to boost contrast

% Call batch merge
dirp = 'D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\Masking\20190227_Benchmarking_plate1_1P_Brightness';
batchMerge(dirp,Cadj)

%% Compare 
figure()
% Original red channel
subplot(5,1,1) 
imshow(B)
title('Original red channel')
% Red channel after adjusted
subplot(5,1,2)
imshow(Cadj)
title('Red channel after adjusted')
% Green channel after normalized
subplot(5,1,3)
Anorm = imhistmatch(A,Cadj,65536); 
imshow(Anorm)
title('Green channel after normalized')
% Red channel after normalized
subplot(5,1,4)
Bnorm = imhistmatch(B,Cadj,65536);
imshow(Bnorm)
title('Red channel after normalized')
% Merged channel after normalized
subplot(5,1,5)
M1 = 0.5*(im2double(Anorm)+im2double(Bnorm));
M2 = im2int16(M1)
imshow(M2)
title('Merged normalized channels')