%% Getting responding map
clear
%dirp = 'D:\OneDrive - rice.edu\Francois\Platfrom\SingleCell\20200304 Fstim single cell mock co transfection\Fstim 1000Hz\Sel2\20200304_115900_387_NDExp_Well0000_Point0000_Seq0001.nd2';
dirp = 'D:\OneDrive - rice.edu\Francois\Platfrom\SingleCell\SC-X2-ASAP2s-60pct-2-1-3\X2_ASAP2s_60pct_2_1-3_trace.nd2';
% FastAnalysis2P perform fast analysis of field stimulation 2P
    % recordings. 
    % 
    % Zhuohe Liu, St-Pierre Lab, Jan. 2020
    % harry.liu@rice.edu
path = dirp;
figshow = false;
detrendOrder = 2;
ch = 1;
binSize = 1;
denoise = false;
    
    [~, fName] = fileparts(path);
    reader = bfGetReader(path);
    I = nd2.read(reader, 'ch', ch, 't', 0);
    t = nd2.getTimePoints(reader, 't', 0, 'dimOrder', 'TCSZ');
    nd2info = nd2.getOptics(reader);
    reader.close();
    SATVAR = nd2info.satValue;
    I = squeeze(I);
    I = img.binsum(I, binSize, binSize)./binSize./binSize;    

    Id = I;
    
    baselinetf = math.isbaseline(reshape(mean(Id, 1:2), [], 1), 'mode');
    Idet = img.detrend(Id, 3, t, 'idleIdx', baselinetf, ...
         'figshow', figshow, 'order', detrendOrder);
    
    % background correction
    bg = img.bglevel(img.binsum(mean(Idet, 3), 2, 2)/4, nd2info.bitDepth);
    Idetbgcrt = img.bgcrt(Idet, bg);
    IdetLbg = reshape(Idetbgcrt, [], size(Idet, 3));
    
    % masking
    [~,~,noisestd] = msp.ephys.estNoiseModel(I/SATVAR, figshow);
    mask = sum(Idet >= SATVAR, 3) == 0 & ...        % not saturated
           mean(Idet, 3) > bg & ...                 % not background (by brightness)
           stdfilt(mean(Idet, 3), true(15)) > noisestd*SATVAR*0.5;    % not background by std
    mask0 = mask;
    [maskt, resp] = img.corrFilt(...
        IdetLbg(mask(:), 10:end), 2, 'groupSize', 20, ...    % remove first 10 frames %  double(idleIdx(10:end))'
        'selectRatio', [0.005, 0.5], 'figshow', figshow);
    mask(mask(:)) = maskt;
    snrl = math.SNR(resp, 1)';
    
    if sum(mask(:)) > 0
        evalout = struct('Fmean', []);
        evalout.Fmean = cat(2, squeeze(mean(I, 1:2)), ...
                               squeeze(mean(Idet, 1:2)), ...
                               squeeze(mean(Idetbgcrt, 1:2)), ...
                               nanmean(IdetLbg(mask(:),:), 1)');  % original, detrended, bgcrt, masked
        legends = {'Full Frame Mean', 'Detrended', 'Background Corrected', 'Selected'};
        w = 600;
        h = 600;
        hs = figure('name', path);
        set(gcf,'units','points','position',[0 0 h w]);
        subplot(3,2,1);
        hold on
        for i = 1:2
            plot(t, evalout.Fmean(:,i), 'DisplayName', legends{i});
        end
        axis tight
        grid on
        ylabel('Intensity');
        legend('show');
        xlabel('Time [s]');
        title('Step 1: Detrend');
        subplot(3,2,2);
        plot(t, evalout.Fmean(:, 3), 'DisplayName', legends{3});
        hold on
        plot(t, evalout.Fmean(:, 4), 'DisplayName', legends{4});
        axis tight
        grid on
        ylabel('Intensity');
        xlabel('Time [s]');
        legend('show', 'Location', 'southoutside');
        title('Step 2: Masking');
        subplot(3,2,4);
        normidx = find(~isoutlier(evalout.Fmean(:, 4)) == 1);
        plot(t, math.normbase(evalout.Fmean(:, 4), normidx) - 1, 'DisplayName', 'dF/F0'); %#ok<FNDSB>
        xlabel('Time [s]');
        ylabel('dF/F0');
        title('Step 3: Infer dF/F0');
        axis tight
        grid on
        subplot(3,2,3)
        plot(snrl, 'DisplayName', 'SNR');
        xlabel(sprintf('Pixel Number (x %d)', 50));
        ylabel('SNR');
        axis tight
        grid on
        title('Step 2: Masking');
        subplot(3,2,5)
        imagesc(mean(Idet, 3));
        axis image
        title('Brightness');
        set(gca, 'XTickLabel', [], 'XTick',[], ...
                 'YTickLabel', [], 'YTick',[]);
        subplot(3,2,6);
        imagesc(mask + mask0);
        axis image
        title('Selected Pixels (Yellow)');
        set(gca, 'XTickLabel', [], 'XTick',[], ...
                 'YTickLabel', [], 'YTick',[]);
    else
        warning('No pixel selected in file: %s ', fName);
    end
    if ~usejava('desktop')  % close MATLAB if in no-desktop mode
        waitfor(hs);
        exit;
    end
[FILEPATH,NAME,EXT] = fileparts(dirp);
% save(strcat(FILEPATH,NAME,'rMap'))

%%
tl = length(t);
stimFile = [1,1,1,2,2,2,3,4,5];

nStimType = max(stimFile);
figure(2), imagesc(mask), axis equal
[Rx, Ry] = find(mask);
% Respn: number of Responding pixels; tl: time length 
Respn = length(Rx);
odd = [];
RespNor = zeros(Respn, nStimType);
figure(7),clf
for i = 1:Respn
    [RespNor(i,:), odd] = Response(Idetbgcrt(Rx(i), Ry(i),:),i,odd,tl,stimFile);
end
[RespNSort, idx] = sortrows(RespNor);
RespNornz = RespNSort(1:Respn-length(odd),:);
RxSort = Rx(idx);
RySort = Ry(idx);
RxNz = RxSort(1:Respn-length(odd),:);
RyNz = RySort(1:Respn-length(odd),:);

%% CLustering
RespNornz1 = 0.001+0.99*normalize(-RespNornz,1,'range');
MaxMax = max( max(abs(RespNornz1(:,:))));
Min = min(min(RespNornz1(:,:)));
Max = max(max(RespNornz1(:,:)));
Mid = mean([Min, Max]);

figure(8),clf, hold on

for stcount = 1:nStimType
    subplot(nStimType,1,stcount), hold on, axis equal, axis([0 512 -50 0]),
    for i = 1:length(RespNornz1)
        plot(RyNz(i),-RxNz(i),'.', 'Color', [RespNornz1(i,stcount)/MaxMax,1-RespNornz1(i,stcount)/MaxMax,RespNornz1(i,stcount)/MaxMax])
    end
    plot(10, -10, '.', 'Color', [Min/MaxMax,1-Min/MaxMax,Min/MaxMax])
    plot(10, -20, '.', 'Color', [Mid/MaxMax,1-Mid/MaxMax,Mid/MaxMax])
    plot(10, -30, '.', 'Color', [Max/MaxMax,1-Max/MaxMax,Max/MaxMax])
    title(strcat('Response ',num2str(stcount))) 
end

RespMap = zeros(size(I,1),size(I,2),nStimType);
for i = 1:length(RespNornz1) 
    for j = 1:nStimType
        RespMap(RxNz(i),RyNz(i),j) = RespNornz1(i,j);
    end
end

save(fullfile(FILEPATH,'RespMap.mat'),'RespMap')

%% Get single cell mask (manual method)
colormap colorcube	
clearvars -except FILEPATH RespMap nStimType stimFile
dirp = fullfile(FILEPATH,'20200304_115900_387_Well0000_Point0000_Seq0000.nd2');
I = nd2.read(dirp,'ch',2);
[mask, ROIs, ~] = ui.cellSel(I);
Savep_tif = fullfile(FILEPATH,'mskTiff.tiff');
imwrite(mask,Savep_tif)
Savep_mat = fullfile(FILEPATH,'mask.mat');
save(Savep_mat,'mask');
figure()
imshow(mask,[])

%% Apply convolution filter
selRespMap = RespMap(:,:,1);
figure(3)
clf
sgtitle('Convolution Filter')
subplot(4,1,1)
imshow(selRespMap,[])
title('short dF/F0 mask')
subplot(4,1,2)
Iblur = imgaussfilt(selRespMap,3,'FilterSize',21);
imshow(Iblur,[])
title('Convolve with Gaussian Filter')
subplot(4,1,3)
maskbw = logical(mask);
maskResp = maskbw.*Iblur;
imshow(maskResp,[])
title('Overlay with mask')
subplot(4,1,4)
nROI = max(mask(:));
ROIResp = zeros(1, nROI);
ROIRespMap = zeros(size(mask));
for i = 1:nROI
    idx = find(mask == i);
    ROIResp(i) = mean(Iblur(idx));
    ROIRespMap(idx) = ROIResp(i);
end
imshow(ROIRespMap,[])
title('ROI ranking')
colormap default
[rpConv,rkConv] = sort(ROIResp,'descend','MissingPlacement','last');

%% Nearest Neighbor Voting
% Create the gallery for the center of each ROI
Centroids = zeros(nROI,2);
for i = 1:nROI
     [xlist,ylist] = find(mask == i);
     [xCenter, yCenter]= findCenter(xlist, ylist);
     Centroids(i,:) = [yCenter, xCenter];
end
figure(4)
clf
sgtitle("Nearest Neighbor Voting")
subplot(4,1,1)
imshow(selRespMap,[],'InitialMagnification','fit')
title('short dF/F0 mask')

subplot(4,1,2) 
imshow(mask,[],'InitialMagnification','fit')
hold on 
scatter(Centroids(:,1),Centroids(:,2),'*r')
title('nucleaus centers')

% each Responding pixel vote for its nearest center
subplot(4,1,2) % show voting points gallery
RespPixelPos = zeros(nnz(selRespMap),2);
[RespPixelPos(:,2),RespPixelPos(:,1),RespPixelVal] = find(selRespMap);
hold on
scatter(RespPixelPos(:,1),RespPixelPos(:,2),1,'.y')
[k,dist] = dsearchn(Centroids,RespPixelPos);
subplot(4,1,3) % show voting pixels
Colors = 0.001+0.99*normalize(k*[1,1,1],'range');
imshow(mask,[],'InitialMagnification','fit')
hold on
title('Assign responding pixels to nucleaus centers')
scatter(RespPixelPos(:,1),RespPixelPos(:,2),1,Colors,'.')
subplot(4,1,4) % show voting result (per roi value)
ROIRespVotMap = double(mask);
ROIVoltResp = zeros(1, nROI);
imshow(mask,[],'InitialMagnification','fit')
hold on
for i = 1:nROI
    idx = find(k == i);
    xpos = RespPixelPos(idx,2);
    ypos = RespPixelPos(idx,1);
    ind = sub2ind(size(selRespMap),xpos,ypos);
    ROIVoltResp(i) = mean(selRespMap(ind),'all');
    nnz(ROIRespVotMap);
    if ~isnan(ROIVoltResp(i))
        ROIRespVotMap(find(ROIRespVotMap == i)) = ROIVoltResp(i);
    else
        ROIRespVotMap(find(ROIRespVotMap == i)) = 0;
    end
end
imshow(ROIRespVotMap,[],'InitialMagnification','fit')
title('rank ROI by mean dF/F0 of assigned pixels')
colormap default
[rpNNV,rkNNV] = sort(ROIVoltResp,'descend','MissingPlacement','last');

%% KNN weighted voting
figure(5)
clf
sgtitle('K-nearest neighbor voting')
subplot(3,1,1)
imshow(selRespMap,[],'InitialMagnification','fit')
title('short dF/F0 mask')

subplot(3,1,2) 
imshow(mask,[],'InitialMagnification','fit')
hold on 
scatter(Centroids(:,1),Centroids(:,2),'*r')
scatter(RespPixelPos(:,1),RespPixelPos(:,2),1,'.y')
title('nucleaus centers')

% KNN search
k = 2; % this number can be derive from the sparsity of the image, 2-5
[kIdx,Distance] = knnsearch(Centroids,RespPixelPos,'K',k);
voltScore = (1./sum(1./Distance,2)).*RespPixelVal*ones(1,k).*(1./Distance);
ROIRespkNNMap = double(mask);
ROIkNNResp = zeros(1, nROI);
for i = 1:nROI
    [kidn,kidk] = find(kIdx == i);
    ind = sub2ind(size(voltScore),kidn,kidk);
    ROIkNNResp(i) = mean(voltScore(ind),'all');
    if ~isnan(ROIkNNResp(i))
        ROIRespkNNMap(find(ROIRespkNNMap == i)) = ROIkNNResp(i);
    else
        ROIRespkNNMap(find(ROIRespkNNMap == i)) = 0;
    end
end
subplot(3,1,3) 

imshow(mask,[],'InitialMagnification','fit')
hold on
title('rank ROI by mean dF/F0 of assigned pixels')
imshow(ROIRespkNNMap,[],'InitialMagnification','fit')
colormap default
[rpKNN,rkKNN] = sort(ROIkNNResp,'descend','MissingPlacement','last');

%% Import manual ROI selection
manualROIraw = readmatrix("D:\OneDrive - rice.edu\Francois\Platfrom\SingleCell\SC-X2-ASAP2s-60pct-2-1-3\SC_manual0316_roi2");
manualROIresp = manualROIraw(:,find(manualROIraw(1,6:end-1)~=1)+5);
manualROIresp = manualROIresp(:,1:end-1); % exclude background roi
nMROI = size(manualROIresp,2);
MROIResp = zeros(nMROI, nStimType);
odd = [];
% figure(8)
% clf
% hold on
% for i = 1:nMROI
%     subplot(3,8,i)
%     a = wdenoise(manualROIresp(:,i),5);
%     arev = a/mean(a(1:300));
%     anorm = -arev+1;
%     plot(arev)
%     [~,locs] = findpeaks(anorm,'MinPeakHeight',0.1,...
%                                         'MinPeakDistance',400,'NPeaks',9);
%     hold on
%     scatter(locs,arev(locs),'r^')
% end

for i = 1:nMROI
    [MROIResp(i,:),odd] = Response(manualROIresp(:,i),i,odd,length(manualROIresp),stimFile);
end
[rpMROI,rkMROI] = sort(MROIResp(:,1),'ascend','MissingPlacement','last');


%% Comparison
I  = nd2.read('D:\OneDrive - rice.edu\Francois\Platfrom\SingleCell\SC-X2-ASAP2s-60pct-2-1-3\X2_ASAP2s_60pct_2_1-3_trace_roi12.nd2');
figure()
subplot(2,1,1)
imshow(I,[])
subplot(2,1,2)
imshow(mask,[])
m = 1:21;
seg = [7,9,8,11,18,19,12,6,15,14,1,2,3,5,4,13,16,10,17,21,20];
rpMROIad = -rpMROI';
rpMROIseq = zeros(1,nROI);
for i = 1:nROI
    rkMROIad(i) = seg(rkMROI(i));
end
figure()
subplot(1,2,1)
hold on
plot(rpMROIad,rpConv,'.','MarkerSize',20);
plot(rpMROIad,rpNNV,'.','MarkerSize',20);
plot(rpMROIad,rpKNN,'.','MarkerSize',20);
legend('Gussian Conv','Nearest Neighbor Voting','k-Nearest Neighbor')

subplot(1,2,2)
hold on
plot(rkMROIad,rkConv,'.','MarkerSize',20);
plot(rkMROIad,rkNNV,'.','MarkerSize',20);
plot(rkMROIad,rkKNN,'.','MarkerSize',20);
plot(m,m,'-y')
legend('Gussian Conv','Nearest Neighbor Voting','k-Nearest Neighbor')

% dF/F0 comparison ranked by manual mask
rpConvad = zeros(size(m));
rpNNVad = zeros(size(m));
rpkNNad = zeros(size(m));
for i = 1:nROI
    rpConvad(i) = rpConv(find(rkConv == rkMROIad(i)));
    rpNNVad(i) = rpNNV(find(rkNNV == rkMROIad(i)));
    rpKNNad(i) = rpKNN(find(rkKNN == rkMROIad(i)));   
end

rkMROIcrt = zeros(size(m));
rkNNVcrt = zeros(size(m));
rkKNNcrt = zeros(size(m));
rkConvcrt = zeros(size(m));
for i = 1:nROI
    rkMROIcrt(i) = i;
    rkConvcrt(i) = find(rkMROI == find(seg == rkConv(i)));
    rkNNVcrt(i) = find(rkMROI == find(seg == rkNNV(i)));
    rkKNNcrt(i) = find(rkMROI == find(seg == rkKNN(i)));
end
DivConv = sqrt(sum((rkConvcrt-rkMROIcrt).^2));
DivNNV = sqrt(sum((rkNNVcrt-rkMROIcrt).^2));
DivKNN = sqrt(sum((rkKNNcrt-rkMROIcrt).^2));
rpDivConv = sqrt(sum((normalize(-rpMROI','range')-normalize(rpConvad,'range')).^2));
rpDivNNV = sqrt(sum((normalize(-rpMROI','range')-normalize(rpNNVad,'range')).^2));
rpDivKNN = sqrt(sum((normalize(-rpMROI','range')-normalize(rpKNNad,'range')).^2));

figure()
hold on
title('Compare normalized dF/F0 per ROI')
rpConvad(isnan(rpConvad))=0;
rpNNVad(isnan(rpNNVad))=0;
rpKNNad(isnan(rpKNNad))=0;
plot(normalize(-rpMROI,'range'),normalize(rpConvad,'range'))
plot(normalize(-rpMROI,'range'),normalize(rpNNVad,'range'))
plot(normalize(-rpMROI,'range'),normalize(rpKNNad,'range'))
plot(0:0.01:1,0:0.01:1,'--k')
xlabel('Normalized -dF/F0 from manual segmentation')
ylabel('Normalized -dF/F0 from nucleaus segmentation')
legend(strcat('Gaussian Convolution Filter, Deviation = ',num2str(rpDivConv)),...
    strcat('Nearest Neighbor Voting, Deviation = ',num2str(rpDivNNV)),...
     strcat('k-Nearest Neighbor, Deviation = ',num2str(rpDivKNN)),'Location','SouthOutside')

figure()
sgtitle('Compare ROI ranking (by s1 dF/F0)')
subplot(1,3,1)
hold on
title('Gussian Conv')
scatter(m,rkMROIcrt,'o')
scatter(m,rkConvcrt,20,'filled')
xlabel("manual segmentation")
ylabel("nucleaus segmentation")

subplot(1,3,2)
hold on
title('Nearest Neighbor Voting')
scatter(m,rkMROIcrt,'o')
scatter(m,rkNNVcrt,20,'filled')
xlabel("manual segmentation")
ylabel("nucleaus segmentation")

subplot(1,3,3)
hold on
title('k-Nearest Neighbor')
scatter(m,rkMROIcrt,'o')
scatter(m,rkKNNcrt,20,'filled')
xlabel("manual segmentation")
ylabel("nucleaus segmentation")
%%
function [Rs,odd1] = Response(Idetbgcrt,numodd,odd,tLen,protocol)
    a4 = reshape(Idetbgcrt,[tLen,1]);
    a4wde = wdenoise(a4,5);
    a4wdenorm = -a4wde/mean(a4wde(1:300))+1;
    if tLen == 7650
        st = 300;
        md = 3000;
        ed = 6000;        
        nPeak = [6,3];
    elseif tLen == 4600      
        st = 800;
        md = 2500;
        ed = 4000;
        nPeak = [4,1];
    end
    % Find peaks correspond to short stimulations
    [~,locs_p1] = findpeaks(a4wdenorm(st+1:md),'MinPeakHeight',0.05,...
                                            'MinPeakDistance',400,'NPeaks',nPeak(1));
    % Find peaks correspond to short stimulations
    locs_p1 = locs_p1+st; 
    [~,locs_p2] = findpeaks(a4wdenorm(md+1:ed),'MinPeakHeight',0.1,...
                                            'MinPeakDistance',900,'NPeaks',nPeak(2));
    locs_p2 = locs_p2+md;     
    locs_p = cat(1, locs_p1, locs_p2);
    if length(locs_p) == sum(nPeak)
        figure(6), hold on, plot(-a4wdenorm); plot(locs_p,-a4wdenorm(locs_p),'r^','MarkerFaceColor','r')
        for k = 1: max(protocol)
            Rs(k) = - mean([a4wdenorm(locs_p(find( protocol == k)))])
        end
        odd1 = odd;
    else
        figure(7), hold on, plot(-a4wdenorm); plot(locs_p,-a4wdenorm(locs_p),'r^','MarkerFaceColor','r')
        disp(strcat('odd = ',num2str(numodd)));
        odd1 = [odd, numodd];
        Rs = zeros(1,max(protocol));
    end
end

function [xcoor,ycoor] = findCenter(xlist, ylist)
    xcoor = mean(xlist);
    ycoor = mean(ylist);
end
