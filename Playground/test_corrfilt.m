%% Load data
clc
clear
input = struct();
input.path = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20201008 double blind\Benchmarking repeat\2pFstim data';
input.tgtpath = 'video';
input.refpath = 'image';
% Load reference
tgtdir = dir(fullfile(input.path, input.tgtpath, '*.nd2'));
refdir = dir(fullfile(input.path, input.refpath, '*.nd2'));
imgLib = table();
for i = 1:numel(tgtdir)
    imgLib.refCh{i} = io.nd2.read(fullfile(refdir(i).folder, refdir(i).name), 'ch', 2, 't', 1:20);
    imgLib.tgtCh{i} = io.nd2.read(fullfile(tgtdir(i).folder, tgtdir(i).name), 'ch', 1, 't', 0);
end

%%
i = 3;
binSize = 2;
detrendOrder = 2;
figshow = true;
bgthreshold = 40;
bgthresholdR = 40;
groupSize = 200;
I = squeeze(imgLib.tgtCh{i});
Iref = mean(squeeze(imgLib.refCh{i}),3);

% Get saturation mask
reader = bfGetReader(fullfile(refdir(i).folder, refdir(i).name));
nd2info = io.nd2.getOptics(reader);
reader.close();
SATVAR = nd2info.satValue;
satMaskChref = (Iref == SATVAR);
satMaskChtgt = any(I == SATVAR, 3);
satMask = satMaskChref | satMaskChtgt;

% Bin image
I = img.binsum(I.*(~satMask), binSize, binSize)./binSize./binSize;    
Iref = img.binsum(Iref.*(~satMask), binSize, binSize)./binSize./binSize;    

% background correction
bg = img.bglevel(I, nd2info.bitDepth);
Ibgcrt = img.bgcrt(I, bg);
bgref = img.bglevel(Iref, nd2info.bitDepth);
Irefbgcrt = img.bgcrt(Iref, bgref);

% Detrend 
baselinetf = math.isbaseline(reshape(mean(Ibgcrt, 1:2), [], 1), 'mode');
Idetbgcrt = img.detrend(Ibgcrt, 3, 'idleIdx', baselinetf', ...
     'figshow', figshow, 'order', detrendOrder);
 
 %% masking
I1 = Ibgcrt(:,:,1);
mask = any(I1 > bgthreshold, 3) & (Irefbgcrt > bgthresholdR);
IdetLbg = reshape(Idetbgcrt.*mask, [], size(Ibgcrt, 3));

template = mean(IdetLbg(mask(:), 1:end),1);
corrmu = corr2dv1d(IdetLbg, template, 2);
corrthreshold = 0;
corrmask = reshape(corrmu, size(mask)) > corrthreshold;
Irefmsked = Irefbgcrt.*corrmask;
Itgtmsked = I1.*corrmask;
Igr = I1./Irefbgcrt;
Igrmsked = Igr.*corrmask;
figure()
subplot(5,1,1),title('brightness mask'),
imshow(mask,[]);
subplot(5,1,2),title('correlation mask'),
imshow(corrmask,[]);
subplot(5,1,3)
imshow(Irefmsked,[]);
subplot(5,1,4)
imshow(Itgtmsked,[]);
subplot(5,1,5)
imshow(Irefbgcrt,[]);
% [maskt, resp] = img.corrFilt(...
%     IdetLbg(mask(:), 10:end), 2, 'groupSize', 20, ...    % remove first 10 frames 
%     'selectRatio', [0.005, 0.5], 'figshow', figshow);
tr = img.groupfun(Idetbgcrt, corrmask);
tr = tr./max(tr(:));
figure(),plot(tr(2,:))
%% brightness correlation
Irefmatch = imhistmatch(Irefbgcrt, I1);
pcc = peasons2d(I1.*mask, Irefmatch.*mask);
figure(4),subplot(3,1,1), imshow(pcc>0,[])
pccmask = pcc > 0;
[mcc1,mcc2] = manders2d(I1, Irefbgcrt);
figure(3),subplot(3,1,2), imshow(mcc1,[])
figure(3),subplot(3,1,3), imshow(mcc2,[])

%%
figure(), hold on,
scatter(reshape(I1,1,[]), reshape(Irefbgcrt,1,[]),5,'filled')
scatter(reshape(I1.*mask,1,[]), reshape(Irefbgcrt.*mask,1,[]),10)
% scatter(reshape(I1.*corrmask,1,[]), reshape(Igr.*corrmask,1,[]),15)
% scatter(reshape(I1.*pccmask,1,[]), reshape(Igr.*pccmask,1,[]),20)
xlabel("Green channel")
ylabel("Red channel")
legend("all","bgmask","corrmask","pccmask")

%%
function r = corr2dv1d(a, b, dim)
    % this is equivalent to r = 1 - pdist2(a, b, 'correlation') but faster. 
    % r = corr2dv1d(a, b, dim)
    % INPUT
    %   a: n by x matrix or x by n. 
    %   b: n by 1 vector or 1 by n. 
    %   dim: dimension of n. 
    % OUTPUT
    %   r: 1 by x or x by 1 vector. 
    % 
    % Zhuohe Liu, St-Pierre Lab, Nov. 2018
    % harry.liu@rice.edu
    
    if nargin < 3
        if size(a, 1) == size(b, 1)
            dim = 1;
        elseif size(a, 2) == size(b, 2)
            dim = 2;
        else
            error('Incompatible a, b dimension. ');
        end
    end
    n = size(a, dim);
    a = a - sum(a, dim)./n;
    b = b - sum(b)./n;
    r = sum(a.*b, dim)./sqrt(sum(a.*a, dim).*sum(b.*b));
end 

  

function [M1, M2] = manders2d(A, B)
    M1 = A./sum(A(:));
    M2 = B./sum(B(:));
end
