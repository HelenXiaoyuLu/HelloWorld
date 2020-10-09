% Load img
Igfp = io.nd2.read('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200910 Spectra 24 well plate\Spectra 10nm spacing 700-1080\20200910_165518_103_NDExp_Well0000_Point0001_Seq0041.nd2');
Ijedi = io.nd2.read('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200910 Spectra 24 well plate\Spectra 10nm spacing 700-1080\20200910_165518_103_NDExp_Well0001_Point0000_Seq0164.nd2');
Igfp2 = io.nd2.read('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200910 Spectra 24 well plate\Spectra 10nm spacing 700-1080\20200910_165518_103_NDExp_Well0000_Point0001_Seq0080.nd2');

%% Tranditional method
bggfp = img.bglevel(Igfp);
bgjedi = img.bglevel(Ijedi);
  bgjedi = img.bglevel(Ijedi);

Igfpbgcrt = img.bgcrt(Igfp);
Ijedibgcrt = img.bgcrt(Ijedi);
%% 
figure(),subplot(1,2,1),imshow(Igfp,[0, 4095]);
subplot(1,2,2),imshow(Igfpbgcrt,[0, 4095])
H = imgaussfilt(Igfp, 2);  % apply before adding Inf values
figure(),imshow(H,[])
% Rolling ball
I = Ijedi;
Icross = zeros(size(I,1),1);
Icross = I(256,:);
flatidx = stdfilt(Icross);

%%
figure(1),clf,
subplot(2,1,1)
plot(Icross)
subplot(2,1,2)
plot(diff(Icross))
%%Rolling ball

%%
s = stdfilt(Ijedi, true(5));
figure(),imshow(s,[])
BW2 = imfill(s ,'holes');
figure(),imshow(BW2,[])
figure(),histogram(BW2,'BinWidth',1)