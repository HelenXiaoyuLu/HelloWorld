%% Load images
function [T_group, T_well, T_fov] = Wellbase_2PFieldStim(cfg)
    %% Wellbase_2PFieldStim is a pipeline main function. 
    % [T_group, T_well, T_fov] = Wellbase_2PFieldStim(cfg)
    % INPUT
    %   cfg (optional): job configuration object. Default, test dataset is used. 
    %       See the source code for default values. 
    % OUTPUT
    %   T_group, T_well, T_fov: result tables. 
    % 
    % Zhuohe Liu, St-Pierre Lab, June 2020
    % harry.liu@rice.edu
    
    arguments
        cfg (1,1) db.jobConfig = db.jobConfig()
    end
    
    % Default Inputs
    cfg.default = db.jobConfig(...
        'path', 'D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\BenchmarkingPlots\20190227_benchmarking\20190227_benchmarking_plate1', ...
        'tgtCh\path', '2p photobleaching video', ...
        'refCh\path', '2p photobleaching image',...
        'tgtCh\sampleFrame',10);
    
%% Phase 1: Import Files and Settings
    lib = sp.libData('2P photobleaching video');
    lib.protocol.path = cfg.default.getConfig('path');
    lib.loadFrame(...
        'path', cfg.default.getConfig('tgtCh\path'), ...
        'ch', 1, ...
        'ncpu', Inf);
    lib.loadFOV();
    lib.regroup('method', 'perwell');      % regroup
    
% 2nd library: load reference images (not the same size)
    lib2 = sp.libData('2P photobleaching 1P ref image');
    lib2.protocol.path = cfg.default.getConfig('path');
    lib2.loadFrame(...
        'path', cfg.default.getConfig('refCh\path'), ...
        'ch', 1, ...
        'ncpu', Inf);    
    lib2.loadFrame(...
        'path', cfg.default.getConfig('refCh\path'), ...
        'ch', 2, ...
        'asch',2,...
        'ncpu', Inf); 
    lib2.loadFOV()
    lib.regroup('method', 'perwell');  
    % temporal alignment
    lib2.getChildren('f').syncChannel('ch', 1);
    T = cfg.default.getConfig('tgtCh\sampleFrame');

%% Phase 2: Image Processing
    if T == 0 % Default: do background correction and masking on last frame
        T = lib.getChildren('f',1).nt;
    end 
    T_fov = table();
    for i = 1:numel(lib.getChildren('f'))
        FOV2P = lib.getChildren('f',i);
        FOV1P = lib2.getChildren('f',i);
        mSet2P = FOV2P.getImageSet("name",'rawFrames');
        mSet1P = FOV1P.getImageSet("name",'rawFrames');
        trace2P = mSet2P.getPixel();
        img1P = mSet1P.getPixel();
        sample2P = squeeze(trace2P(:,:,T));
        gfp1P = squeeze(img1P(100:400,:,1));
        rfp1P = squeeze(img1P(100:400,:,2));
        [tform, overlapratio, ~] = img.imalign(mat2gray(squeeze(trace2P(:,:,T))), ...
            imhistmatch(mat2gray(gfp1P),mat2gray(squeeze(trace2P(:,:,T)))),'modality','2p'); % base, mask
        imggfpreg = imwarp(gfp1P, tform,'OutputView', imref2d(size(sample2P)));
        imgrfpreg = imwarp(rfp1P, tform,'OutputView', imref2d(size(sample2P)));
        maskSat = mSet2P.getTag("C",1).getSaturationMask(squeeze(trace2P(:,:,1)));
        gfp2Pbgcrt = img.bgcrt(squeeze(trace2P(:,:,1)), img.bglevel(img.nanmasking(squeeze(trace2P(:,:,T)), ~maskSat), 12));
        trace2Pbgcrt = img.bgcrt(trace2P,img.bglevel(img.nanmasking(squeeze(trace2P(:,:,T)), ~maskSat), 12));
        gfp1Pbgcrt = img.bgcrt(imggfpreg,img.bglevel(imggfpreg, 16));
        rfp1Pbgcrt = img.bgcrt(imgrfpreg,img.bglevel(imgrfpreg, 16));
        finalMask = ~maskSat & gfp2Pbgcrt > 20;
        meangfp2P = mean(gfp2Pbgcrt(finalMask));
        meangfp1P = mean(gfp1Pbgcrt(finalMask));
        meanrfp1P = mean(rfp1Pbgcrt(finalMask));
        trace = img.groupfun(trace2Pbgcrt, finalMask, @nanmean, 'toCell', true, 'ignoreGroup', 0);
        T_fov.("Name")(i) = string(FOV2P.Name);
        T_fov.("Final Mask"){i} = finalMask;
        T_fov.("GFP 2P sample"){i} = gfp2Pbgcrt;
        T_fov.("GFP 1P aligned"){i} = gfp1Pbgcrt;
        T_fov.("RFP 1P aligned"){i} = rfp1Pbgcrt;
        T_fov.("Alignment Overlap Ratio")(i) = overlapratio;
        T_fov.("GFP 2P Mean Brightness")(i) = meangfp2P;
        T_fov.("GFP 1P Mean Brightness")(i) = meangfp1P;
        T_fov.("RFP 1P Mean Brightness")(i) = meanrfp1P;
        T_fov.("Trace bgcrt"){i} = trace;
    end
save(fullfile(cfg.default.getConfig('path'),'T_fov.mat'),'T_fov');  
save(fullfile(cfg.default.getConfig('path'),'lib.mat'),'lib','lib2');

%% Phase3: Merge and interpret results
% T_fov = [T_fov_plate1;T_fov_plate2];
T_fov_plate1 = sortrows(T_fov_plate1,'Name','ascend');
T_fov_plate1.("GFP2P/RFP1P")= T_fov_plate1.("GFP 2P Mean Brightness")./T_fov_plate1.("RFP 1P Mean Brightness");
T_fov_plate1.("GFP1P/RFP1P")= T_fov_plate1.("GFP 1P Mean Brightness")./T_fov_plate1.("RFP 1P Mean Brightness");
T_fov_plate1.("GFP2P/GFP1P")= T_fov_plate1.("GFP 2P Mean Brightness")./T_fov_plate1.("GFP 1P Mean Brightness");

vargroups = findgroups(T_fov_plate1.Name);
T_group_plate1 = table();
T_group_plate1.Name = splitapply(@unique, T_fov_plate1.Name, vargroups);
for i = 1:height(T_group_plate1)
    T_group_plate1.Name(i) = libplot.fzmatch(T_group_plate1.Name{i});
end
T_group_plate1.("n FOV") = splitapply(@numel, T_fov_plate1.Name, vargroups);
T_group_plate1.("GFP 2P Mean Brightness") = splitapply(@mean, T_fov_plate1.("GFP 2P Mean Brightness"), vargroups);
T_group_plate1.("GFP 2P Mean Brightness STD") = splitapply(@std, T_fov_plate1.("GFP 2P Mean Brightness"), vargroups);
T_group_plate1.("GFP 1P Mean Brightness") = splitapply(@mean, T_fov_plate1.("GFP 1P Mean Brightness"), vargroups);
T_group_plate1.("GFP 1P Mean Brightness STD") = splitapply(@std, T_fov_plate1.("GFP 1P Mean Brightness"), vargroups);
T_group_plate1.("RFP 1P Mean Brightness") = splitapply(@mean, T_fov_plate1.("RFP 1P Mean Brightness"), vargroups);
T_group_plate1.("RFP 1P Mean Brightness STD") = splitapply(@std, T_fov_plate1.("RFP 1P Mean Brightness"), vargroups);
T_group_plate1.("GFP2P/RFP1P") = splitapply(@mean, T_fov_plate1.("GFP2P/RFP1P"), vargroups);
T_group_plate1.("GFP2P/RFP1P STD") = splitapply(@std, T_fov_plate1.("GFP2P/RFP1P"), vargroups);
T_group_plate1.("GFP1P/RFP1P") = splitapply(@mean, T_fov_plate1.("GFP1P/RFP1P"), vargroups);
T_group_plate1.("GFP1P/RFP1P STD") = splitapply(@std, T_fov_plate1.("GFP1P/RFP1P"), vargroups);
T_group_plate1.("GFP2P/GFP1P") = splitapply(@mean, T_fov_plate1.("GFP2P/GFP1P"), vargroups);
T_group_plate1.("GFP2P/GFP1P STD") = splitapply(@std, T_fov_plate1.("GFP2P/GFP1P"), vargroups);

% Plate 2
T_fov_plate2 = sortrows(T_fov_plate2,'Name','ascend');
T_fov_plate2.("GFP2P/RFP1P")= T_fov_plate2.("GFP 2P Mean Brightness")./T_fov_plate2.("RFP 1P Mean Brightness");
T_fov_plate2.("GFP1P/RFP1P")= T_fov_plate2.("GFP 1P Mean Brightness")./T_fov_plate2.("RFP 1P Mean Brightness");
T_fov_plate2.("GFP2P/GFP1P")= T_fov_plate2.("GFP 2P Mean Brightness")./T_fov_plate2.("GFP 1P Mean Brightness");

vargroups = findgroups(T_fov_plate2.Name);
T_group_plate2 = table();
T_group_plate2.Name = splitapply(@unique, T_fov_plate2.Name, vargroups);
for i = 1:height(T_group_plate2)
    T_group_plate2.Name(i) = libplot.fzmatch(T_group_plate2.Name{i});
end
T_group_plate2.("n FOV") = splitapply(@numel, T_fov_plate2.Name, vargroups);
T_group_plate2.("GFP 2P Mean Brightness") = splitapply(@mean, T_fov_plate2.("GFP 2P Mean Brightness"), vargroups);
T_group_plate2.("GFP 2P Mean Brightness STD") = splitapply(@std, T_fov_plate2.("GFP 2P Mean Brightness"), vargroups);
T_group_plate2.("GFP 1P Mean Brightness") = splitapply(@mean, T_fov_plate2.("GFP 1P Mean Brightness"), vargroups);
T_group_plate2.("GFP 1P Mean Brightness STD") = splitapply(@std, T_fov_plate2.("GFP 1P Mean Brightness"), vargroups);
T_group_plate2.("RFP 1P Mean Brightness") = splitapply(@mean, T_fov_plate2.("RFP 1P Mean Brightness"), vargroups);
T_group_plate2.("RFP 1P Mean Brightness STD") = splitapply(@std, T_fov_plate2.("RFP 1P Mean Brightness"), vargroups);
T_group_plate2.("GFP2P/RFP1P") = splitapply(@mean, T_fov_plate2.("GFP2P/RFP1P"), vargroups);
T_group_plate2.("GFP2P/RFP1P STD") = splitapply(@std, T_fov_plate2.("GFP2P/RFP1P"), vargroups);
T_group_plate2.("GFP1P/RFP1P") = splitapply(@mean, T_fov_plate2.("GFP1P/RFP1P"), vargroups);
T_group_plate2.("GFP1P/RFP1P STD") = splitapply(@std, T_fov_plate2.("GFP1P/RFP1P"), vargroups);
T_group_plate2.("GFP2P/GFP1P") = splitapply(@mean, T_fov_plate2.("GFP2P/GFP1P"), vargroups);
T_group_plate2.("GFP2P/GFP1P STD") = splitapply(@std, T_fov_plate2.("GFP2P/GFP1P"), vargroups);

%% Plot results plate 1
sel = {'ASAP2s','JEDI-1P'};
sel = cellstr(T_group_plate1.Name);
libplot.formattedBar(T_group_plate1, 'GFP2P/RFP1P', ...
  'GFP1P/RFP1P', 'pp1std','GFP2P/RFP1P STD','pp2std','GFP1P/RFP1P STD', ...
  'sort_on', 'GFP2P/RFP1P', 'sort_direction','descend','Name','Name','ylabel', ...
  'Relative Brightness','selectNames',sel);
title('Plate 1')
libplot.formattedBar(T_group_plate1, 'GFP2P/RFP1P', ...
  'GFP1P/RFP1P', 'pp1std','GFP2P/RFP1P STD','pp2std','GFP1P/RFP1P STD', ...
  'sort_on', 'GFP2P/RFP1P', 'sort_direction','descend','Name','Name','ylabel', ...
  'Relative Brightness Normalized','norm2ctrl','ASAP2s','selectNames',sel);
title('Plate 1')

c = colormap(lines);
libplot.formattedBar(T_group_plate1, 'GFP 2P Mean Brightness', ...
  'RFP 1P Mean Brightness', 'pp1std','GFP 2P Mean Brightness STD','pp2std',...
  'RFP 1P Mean Brightness STD','ylabel','Mean Brightness [a.u.]', ...
  'sort_on', 'GFP 2P Mean Brightness', 'sort_direction','descend','Name',...
  'Name','selectNames',sel,'colormap', c([5,2],:));
title('Plate 1')
libplot.formattedBar(T_group_plate1, 'GFP 2P Mean Brightness', ...
  'RFP 1P Mean Brightness', 'pp1std','GFP 2P Mean Brightness STD','pp2std',...
  'RFP 1P Mean Brightness STD','ylabel','Mean Brightness Normalized [a.u.]', ...
  'sort_on', 'GFP 2P Mean Brightness', 'sort_direction','descend','Name',...
  'Name','selectNames',sel,'colormap', c([5,2],:),'norm2ctrl','ASAP2s');
title('Plate 1')

%% Plot results plate 2
T_group_plate2 = sortrows(T_group_plate2,'Name','ascend');
sel = cellstr(T_group_plate2.Name);
libplot.formattedBar(T_group_plate2, 'GFP2P/RFP1P', ...
  'GFP1P/RFP1P', 'pp1std','GFP2P/RFP1P STD','pp2std','GFP1P/RFP1P STD', ...
  'sort_on', 'GFP2P/RFP1P', 'sort_direction','descend','Name','Name','ylabel', ...
  'Relative Brightness','selectNames',sel);
title('Plate 2')
libplot.formattedBar(T_group_plate2, 'GFP2P/RFP1P', ...
  'GFP1P/RFP1P', 'pp1std','GFP2P/RFP1P STD','pp2std','GFP1P/RFP1P STD', ...
  'sort_on', 'GFP2P/RFP1P', 'sort_direction','descend','Name','Name','ylabel', ...
  'Relative Brightness Normalized','norm2ctrl','ASAP2s','selectNames',sel);
title('Plate 2')

c = colormap(lines);
libplot.formattedBar(T_group_plate2, 'GFP 2P Mean Brightness', ...
  'RFP 1P Mean Brightness', 'pp1std','GFP 2P Mean Brightness STD','pp2std',...
  'RFP 1P Mean Brightness STD','ylabel','Mean Brightness [a.u.]', ...
  'sort_on', 'GFP 2P Mean Brightness', 'sort_direction','descend','Name',...
  'Name','selectNames',sel,'colormap', c([5,2],:));
title('Plate 2')
libplot.formattedBar(T_group_plate2, 'GFP 2P Mean Brightness', ...
  'RFP 1P Mean Brightness', 'pp1std','GFP 2P Mean Brightness STD','pp2std',...
  'RFP 1P Mean Brightness STD','ylabel','Mean Brightness Normalized [a.u.]', ...
  'sort_on', 'GFP 2P Mean Brightness', 'sort_direction','descend','Name',...
  'Name','selectNames',sel,'colormap', c([5,2],:),'norm2ctrl','ASAP2s');
title('Plate 2')

i = 3;
figure(),imshow(T_fov_plate1.("Final Mask"){i})
figure(),subplot(2,1,1),imshow(T_fov_plate1.("GFP 2P sample"){i},[]),subplot(2,1,2),imshow(imhistmatch(T_fov_plate1.("GFP 1P aligned"){i},T_fov_plate1.("GFP 2P sample"){i}),[])
figure(), imshowpair(T_fov_plate1.("GFP 2P sample"){i},imhistmatch(T_fov_plate1.("RFP 1P aligned"){i},T_fov_plate1.("GFP 2P sample"){i}))
title('GFP 2P vs RFP 1P')
figure(),imshow(T_fov_plate1.("Final Mask"){i})

%% Correlate plate1 vs plate2
figure(6)
hold on
cmap = colorcube(60);
for i = 1:height(T_group_plate2)
    scatter(T_group_plate1.("GFP1P/RFP1P")(i),T_group_plate2.("GFP1P/RFP1P")(i),[],cmap(i,:),'filled')
end
xlabel('1P GFP/ 1P RFP (plate 1)')
ylabel('1P GFP/ 1P RFP (plate 2)')
libplot.lnrFitting(T_group_plate1.("GFP1P/RFP1P"),T_group_plate2.("GFP1P/RFP1P"),false,true,6);
legend(T_group_plate2.Name,'Location','eastoutside')
axis([0 1.5 0 1.8])
axis square

figure(7)
hold on
cmap = colorcube(60);
for i = 1:height(T_group_plate2)
    scatter(T_group_plate1.("GFP2P/RFP1P")(i),T_group_plate2.("GFP2P/RFP1P")(i),[],cmap(i,:),'filled')
end
xlabel('2P GFP/ 1P RFP (plate 1)')
ylabel('2P GFP/ 1P RFP (plate 2)')
libplot.lnrFitting(T_group_plate1.("GFP2P/RFP1P"),T_group_plate2.("GFP2P/RFP1P"),false,true,7);
legend(T_group_plate2.Name,'Location','eastoutside')
axis square

%%
T_well_stats = sortrows(T_well_stats,'Name','ascend');
T_well_stats.Properties.RowNames = T_well_stats.Name;
T_group_plate1.Properties.RowNames = T_group_plate1.Name;
figure()
hold on
cmap = colorcube(60);
for i = 1:height(T_well_stats)
    scatter(T_well_stats.GRRatioMean(i),T_group_plate1.("GFP1P/RFP1P")(i),[],cmap(i,:),'filled')
end
xlabel('1P G/R from 1P Fstim')
ylabel('1P G/F from 1P ref image')
libplot.lnrFitting(T_well_stats.GRRatioMean,T_group_plate1.("GFP1P/RFP1P"),false,true,1);
legend(T_group_plate1.Name,'Location','eastoutside')
axis([0 0.35 0 1.4])
axis square

%% Correlate plate1 vs 2P cyOFP
figure(9)
clf
hold on
cmap = colorcube(60);
nameList = T_group_plate1.Name;
for i = 1:length(nameList)
    scatter(T_group_plate1.("GFP1P/RFP1P")(nameList(i)),T_well_stats_1P.GRratioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_plate1.("GFP1P/RFP1P")(nameList);
y = T_well_stats_1P.GRratioMean_Mean(nameList);
xlabel('1P GFP/1P RFP mCherry')
ylabel('1P GFP/1P RFP cyOFP')
axis square
[k, R2, pVal]=libplot.lnrFitting(x,y,true,true,9,'polyfit') 
legend(nameList,'Location','eastoutside')
title('All benchmarking constructs')

figure(10)
clf
hold on
cmap = colorcube(60);
nameList = ["ASAP1";"ASAP1-EGFP";"ASAP1-dpOPT";"ASAP2s";"ASAP2s-H152E";"ASAP3";"JEDI-1P";"JEDI-2P";"JEDI-beta"];
for i = 1:length(nameList)
    scatter(T_group_plate1.("GFP1P/RFP1P")(nameList(i)),T_well_stats_1P.GRratioMean_Mean(nameList(i)),[],cmap(i,:),'filled')
end
x = T_group_plate1.("GFP1P/RFP1P")(nameList);
y = T_well_stats_1P.GRratioMean_Mean(nameList);
xlabel('1P GFP/1P RFP mCherry')
ylabel('1P GFP/1P RFP cyOFP')
libplot.lnrFitting(x,y,true,true,10,'polyfit') 
legend(nameList,'Location','eastoutside')
axis square
title('ASAP-like constructs')

