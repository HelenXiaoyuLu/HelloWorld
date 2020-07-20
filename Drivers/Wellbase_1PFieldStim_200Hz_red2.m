close all
dirp = 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20190310-0318 Jun cp\20190318 Jun cp 1000Hz\1000Hz\FastAnalysis2P';
files = 'X2-cpmTagRFP1-170-169_1*';
flist = dir(fullfile(dirp, files));
llist = struct();
for i = 1:length(flist)
    F = openfig(fullfile(flist(i).folder, flist(i).name));
    handles=findobj(F,'Type','line');
    llist(i).XData = handles(2).XData;
    llist(i).YData = handles(2).YData;
    h = subplot(3,2,6);3
    I = getimage(h);
    llist(i).npixel = length(find(I(:) == 2));
end

%% 
i = 4;
trWell.Name{i} = "X2-cpmTagRFP1-170-169_2";
trWell.Raw{i} = llist;
Ylist = zeros(length(llist),length(llist(1).YData));
Nlist = zeros(length(llist),1);
for k = 1:length(llist)
    Nlist(k) = llist(i).npixel;
    Ylist(k,:) = Nlist(k)*llist(k).YData;  
end
trWell.weightedTrace{i} = [llist(1).XData;sum(Ylist,1)/sum(Nlist(:))];
trWell.detrendedTrace{i} = [llist(1).XData;detrend(sum(Ylist,1)/sum(Nlist(:)),10)];
figure()
plot(llist(1).XData,detrend(sum(Ylist,1)/sum(Nlist(:)),10))

%%
clear                     
dirp = 'D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20190310-0318 Jun cp\20190308 Jun cp 200Hz';
path = fullfile(dirp,'200Hz\X2-cpmTagRFP1-170-169_2_1-1.nd2');
load(fullfile(dirp,'Stim_Dim2Bright Long'),'stim');
import io.*
path = char(path);
figshow = false;
detrendOrder = 2;
ch = 1;
binSize = 2;    
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
[~,~,noisestd] = msp.ephys.estNoiseModel(I/SATVAR, figshow);
mask = sum(Idet >= SATVAR, 3) == 0 & ...        % not saturated
       mean(Idet, 3) > bg & ...                 % not background (by brightness)
       stdfilt(mean(Idet, 3), true(15)) > noisestd*SATVAR*0.5;
s = [];
s.y = mean(IdetLbg(mask(:), 10:end),1);
s.x = t(10:end);
delay = img.temporalAlignment(s.y, s.x, 'method', 'template', ...
'template', stim, 'figshow', false);
s.lag=delay;
tsin = timeseries(stim.y,stim.x-delay);
tsout = resample(tsin,t);
RSstim = tsout.Data;
RSstim(isnan(RSstim))=0;
mask1 = mask;
[maskt, resp] = img.corrFilt(...
    IdetLbg(mask(:), 10:end), 2,-RSstim(10:700), 'groupSize', 200, ...    % remove first 10 frames %  double(idleIdx(10:end))'
    'selectRatio', [0.05, 0.8], 'figshow', figshow);
evalout = struct('Fmean', []);
mask1(mask1(:)) = maskt;
evalout.Fmean = cat(2, squeeze(mean(I, 1:2)), ...
                               squeeze(mean(Idet, 1:2)), ...
                               squeeze(mean(Idetbgcrt, 1:2)), ...
                               nanmean(IdetLbg(mask1(:),:), 1)');  % original, detrended, bgcrt, masked
normidx = find(~isoutlier(evalout.Fmean(:, 4)) == 1);
mask2 = mask;
[maskt2, resp2] = img.corrFilt(...
    IdetLbg(mask(:), 10:end), 2, -RSstim(10:700), 'groupSize', 20, ...    % remove first 10 frames %  double(idleIdx(10:end))'
    'selectRatio', [0.005, 0.5], 'figshow', figshow);
mask2(mask2(:)) = maskt2;
evalout.Fmean2 = cat(2, squeeze(mean(I, 1:2)), ...
                               squeeze(mean(Idet, 1:2)), ...
                               squeeze(mean(Idetbgcrt, 1:2)), ...
                               nanmean(IdetLbg(mask2(:),:), 1)');  % original, detrended, bgcrt, masked
normidx2 = find(~isoutlier(evalout.Fmean2(:, 4)) == 1);
figure()
subplot(2,1,1)
plot(t, math.normbase(evalout.Fmean(:, 4), normidx) - 1, 'DisplayName', 'dF/F0'); 
subplot(2,1,2)
plot(t, math.normbase(evalout.Fmean2(:, 4), normidx2) - 1, 'DisplayName', 'dF/F0'); 
figure()
plot(t, evalout.Fmean2(:, 4)- 1, 'DisplayName', 'dF/F0'); 
