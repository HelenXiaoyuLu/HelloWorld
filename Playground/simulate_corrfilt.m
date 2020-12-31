%% Load data 
clc
clear
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\20200126 Benchmarking\plate2';
j0.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask foreground no thresholding'));
j1.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask foreground'));
j2.output = load(fullfile(dirp, '2pFstimSingleCell 20200126 20201026190510 T_roi max with manual mask background'));
j3.output = load(fullfile(dirp, '2pFieldStimCorrFilt 20200126 20201029201318 T_out (P2)'));  % correlation filter
j4.output = load(fullfile(dirp, '2pFstim 20200126 T_out'));  % thresholding method
j5.output = load(fullfile(dirp, '2PFieldStimCorrFilt 20200126 20201101005346 (P2)'));  % correlation filter ver_2
T_fovcorr = j3.output.T_fov;
T_fovcorr_v2 = j5.output.T_fov;
T_fovcorr.Name = splitapply(@parsename, T_fovcorr.Path(:), (1:height(T_fovcorr))'); 
T_fovcorr_v2.Name = splitapply(@parsename, T_fovcorr_v2.Path(:), (1:height(T_fovcorr_v2))'); 
T_thres = j4.output.T_fov;
T_thres.Name = splitapply(@parsename, T_thres.Path(:), (1:height(T_thres))'); 
selvarnames = ["ASAP1-EGFP-cyOFP"; "ASAP1-cyOFP"; "ASAP1-dpOPT-cyOFP"; ...
    "ASAP2s-H152E-cyOFP"; "ASAP2s-T207H-cyOFP"; "ASAP2s-cyOFP"; ...
    "ASAP3-cyOFP"; "JEDI-1P-cyOFP"; "JEDI-2P-cyOFP"];
isvarselected = @(s) any(contains(selvarnames, s));
T_fovcorr.("selected") = splitapply(isvarselected, T_fovcorr.Name(:), (1:height(T_fovcorr))'); 
T_fovcorr = T_fovcorr(T_fovcorr.("selected") == 1, :);
T_fovcorr_v2.("selected") = splitapply(isvarselected, T_fovcorr_v2.Name(:), (1:height(T_fovcorr_v2))'); 
T_fovcorr_v2  = T_fovcorr_v2(T_fovcorr_v2.("selected") == 1, :);
T_thres.("selected") = splitapply(isvarselected, T_thres.Name(:), (1:height(T_thres))'); 
T_thres = T_thres(T_thres.("selected") == 1, :);

%%
T = table();
T.Name = j1.output.T_fov.FOV_Name;
grouping = findgroups(T.Name);
cmap = lines;
T.("Foreground Mask") = j1.output.T_fov.ManualMask;
T.("Background Mask") = j2.output.T_fov.ManualMask;
desaturate = @(M1, M2) {M1{:} & ~M2{:}};
T.("Responding Mask") = splitapply(desaturate, T_fovcorr.finalMask(:), T_fovcorr.satMask(:), (1:height(T_fovcorr))');
T.("Responding Mask 2") = splitapply(desaturate, T_fovcorr_v2.finalMask(:), T_fovcorr_v2.satMask(:), (1:height(T_fovcorr_v2))');
T.("final Mask") = T_thres.finalMask;
T.("Bin G") = T_fovcorr.GMap;
T.("Bin R") = T_fovcorr.RMap;
T.("Bin final G") = T_thres.GMap;
T.("Bin final R") = T_thres.RMap;
T.G = j0.output.T_fov.GMap;
T.R = j0.output.T_fov.RMap;
for i = 1:height(T)
    % number of pixels in each map
    T.("Foreground n_pixel")(i) = sum(T.("Foreground Mask"){i}, 'all');
    T.("Background n_pixel")(i) = sum(T.("Background Mask"){i}, 'all');
    T.("Responding n_pixel")(i) = sum(T.("Responding Mask"){i}, 'all');
    T.("Responding n_pixel 2")(i) = sum(T.("Responding Mask 2"){i}, 'all');
    T.("final n_pixel")(i) = sum(T.("final Mask"){i}, 'all');
    % G/R after masking
    T.("Foreground G"){i} = T.("Foreground Mask"){i}.*T.G{i};
    T.("Foreground R"){i} = T.("Foreground Mask"){i}.*T.R{i};    
    T.("Background G"){i} = T.("Background Mask"){i}.*T.G{i};
    T.("Background R"){i} = T.("Background Mask"){i}.*T.R{i};
    T.("Responding G"){i} = T.("Responding Mask"){i}.*T.("Bin G"){i};
    T.("Responding R"){i} = T.("Responding Mask"){i}.*T.("Bin R"){i};
    T.("Responding G 2"){i} = T.("Responding Mask 2"){i}.*T.("Bin G"){i};
    T.("Responding R 2"){i} = T.("Responding Mask 2"){i}.*T.("Bin R"){i};
    T.("final G"){i} = T.("final Mask"){i}.*T.("Bin final G"){i};
    T.("final R"){i} = T.("final Mask"){i}.*T.("Bin final R"){i};
    forpcc = pearsons2d(T.("Foreground G"){i}, T.("Foreground R"){i});
    backpcc = pearsons2d(T.("Background G"){i}, T.("Background R"){i});
    resppcc = pearsons2d(T.("Responding G"){i}, T.("Responding R"){i});
    resppcc2 = pearsons2d(T.("Responding G 2"){i}, T.("Responding R 2"){i});
    finalpcc = pearsons2d(T.("final G"){i}, T.("final R"){i});
    allpcc  = pearsons2d(T.G{i},T.R{i});
    T.("Foreground PCC"){i} = forpcc;
    T.("Background PCC"){i} = backpcc;
    T.("Responding PCC"){i} = resppcc;
    T.("Responding PCC 2"){i} = resppcc2;
    T.("final PCC"){i} = finalpcc;
    T.("PCC"){i} = allpcc;
    T.("Foreground PCC Sum")(i) = sum(forpcc(:));
    T.("Background PCC Sum")(i) = sum(backpcc(:));
    T.("Responding PCC Sum")(i) = sum(resppcc(:));
    T.("Responding PCC Sum 2")(i) = sum(resppcc2(:));
    T.("final PCC Sum")(i) = sum(finalpcc(:));
    T.("PCC Mean")(i) = sum(allpcc(:));
    T.("Foreground GMean")(i) = mean(nonzeros(T.("Foreground G"){i}),'all');
    T.("Background GMean")(i) = mean(nonzeros(T.("Background G"){i}),'all');
    T.("Responding GMean")(i) = mean(nonzeros(T.("Responding G"){i}),'all');
    T.("Responding GMean 2")(i) = mean(nonzeros(T.("Responding G 2"){i}),'all');
    T.("final GMean")(i) = mean(nonzeros(T.("final G"){i}),'all');
    T.("all GMean")(i) = mean(nonzeros(T.G{i}),'all');
    T.("Foreground RMean")(i) = mean(nonzeros(T.("Foreground R"){i}),'all');
    T.("Background RMean")(i) = mean(nonzeros(T.("Background R"){i}),'all');
    T.("Responding RMean")(i) = mean(nonzeros(T.("Responding R"){i}),'all');
    T.("Responding RMean 2")(i) = mean(nonzeros(T.("Responding R 2"){i}),'all');
    T.("final RMean")(i) = mean(nonzeros(T.("final R"){i}),'all');
    T.("all RMean")(i) = mean(nonzeros(T.R{i}),'all');   
end
T.("Foreground GRratio") = T.("Foreground GMean")./T.("Foreground RMean");
T.("Background GRratio") = T.("Background GMean")./T.("Background RMean");
T.("Responding GRratio") = T.("Responding GMean")./T.("Responding RMean");
T.("Responding GRratio 2") = T.("Responding GMean 2")./T.("Responding RMean 2");
T.("final GRratio") = T.("final GMean")./T.("final RMean");
T.("all GRratio") = T.("all GMean")./T.("all RMean"); 

%%
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
t.Title.String = 'Channel intensity distribution';
for groupn = 1:max(grouping)
    ax = nexttile(t);
    hold(ax,'on')
    groupidx = find(grouping == groupn);
%     forR = [];
%     forG = [];
%     backR = [];
    Gpool = [];
    Rpool = [];
    for j = min(groupidx):max(groupidx)
        s = gobjects(1,1);
        forR = reshape(T.("Foreground R"){j},1,[]);
        forG = reshape(T.("Foreground G"){j},1,[]);
        backR = reshape(T.("Background R"){j},1,[]);
        backG = reshape(T.("Background G"){j},1,[]);
        respR = reshape(T.("Responding R"){j},1,[]);
        respG = reshape(T.("Responding G"){j},1,[]);
        respR2 = reshape(T.("Responding R 2"){j},1,[]);
        respG2 = reshape(T.("Responding G 2"){j},1,[]);
        finalR = reshape(T.("final R"){j},1,[]);
        finalG = reshape(T.("final G"){j},1,[]);
        allR = reshape(T.R{j},1,[]);
        allG = reshape(T.G{j},1,[]);
%         Gpool = [Gpool, respG];
%         Rpool = [Rpool, respR];
%         s(1) = scatter(allR, allG, 1, cmap(3,:),'filled', ...
% %              'MarkerFaceAlpha', 0.5, 'DisplayName','All');   
%         s(2) = scatter(forR, forG, 1, cmap(1,:),'filled',...
%             'DisplayName','Foreground');
%         s(3) = scatter(backR, backG, 1, cmap(2,:),...
%             'DisplayName','Background');    
%         s(4) = scatter(finalR, finalG, 1, cmap(6,:),'filled',...
%             'DisplayName','Foreground gating');    
%         s(1) = scatter(respR(respR>0), respG(respR>0), 1, cmap(5,:),'filled',...
%             'DisplayName','corrFilt gating');  
        s(1) = scatter(respR2(respR2>0), respG2(respR2>0), 1, cmap(4,:),'filled',...
            'DisplayName','corrFilt gating v2');  
        %         values = hist3([respR', respG'], [20, 20]);
%         imshow(values.', []);
%         view(2)
    end
%     cvalues = heatscatter(Rpool', Gpool', 40);
%     s(1) = scatter(ax, Rpool, Gpool, 5, cvalues, 'filled',...
%             'DisplayName','Foreground gating');  
%     colormap default
%     colorbar
    ax.Title.String = T.Name(j);
    ax.XAxis.Label.String = 'Red Intensity';
    ax.YAxis.Label.String = 'Green Intensity';
    ax.XAxis.Limits = [0,4096];
    ax.YAxis.Limits = [0,4096];
end

%%
% f = figure('Color', [1, 1, 1]);
% t = tiledlayout(f, 'flow', 'Padding', 'none');
% t.Title.String = 'Channel intensity distribution';
% ax = nexttile(t);
% hold(ax,'on')
% s = gobjects(1,5);
% for groupn = 1:max(grouping)   
%     groupidx = find(grouping == groupn);
%     s(1) = scatter(T.("all RMean"), T.("all GMean"), 10, cmap(3,:),'filled', ...
%              'MarkerFaceAlpha', 0.5, 'DisplayName','All');        
%     s(2) = scatter(T.("Foreground RMean"),T.("Foreground GMean"), 10, cmap(1,:),'filled',...
%         'DisplayName','Foreground');
%     s(3) = scatter(T.("Background RMean"), T.("Background GMean"), 10, cmap(2,:),...
%         'DisplayName','Background');    
%     s(4) = scatter(T.("final RMean"), T.("final GMean"), 10, cmap(6,:),'filled',...
%         'DisplayName','final');    
%     s(5) = scatter(T.("Responding RMean"), T.("Responding GMean"), 10, cmap(5,:),'filled',...
%         'DisplayName','Responding');       
% end
% ax.Title.String = T.Name(j);
% ax.XAxis.Label.String = 'Red Intensity';
% ax.YAxis.Label.String = 'Green Intensity';
% ax.XAxis.Limits = [0,4096];
% ax.YAxis.Limits = [0,4096];
%  legend(s)

%% Pull stats from fov to group
grouping = findgroups(T.Name);
weightedMean = @(x, w) sum(x.*w, 'all')./sum(w, 'all');
T_fov2group = table();
T_fov2group.Name = splitapply(@unique, T.Name, grouping);
T_fov2group.("Foreground n_pixel") = splitapply(@sum, T.("Foreground n_pixel"), grouping);
T_fov2group.("Background n_pixel") = splitapply(@sum, T.("Background n_pixel"), grouping);
T_fov2group.("Responding n_pixel") = splitapply(@sum, T.("Responding n_pixel"), grouping);
T_fov2group.("final n_pixel") = splitapply(@sum, T.("final n_pixel"), grouping);
T_fov2group.("Foreground PCC mean") = splitapply(weightedMean, T.("Foreground PCC Sum"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background PCC mean") = splitapply(weightedMean, T.("Background PCC Sum"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding PCC mean") = splitapply(weightedMean, T.("Responding PCC Sum"), T.("Responding n_pixel"), grouping);
T_fov2group.("final PCC mean") = splitapply(weightedMean, T.("final PCC Sum"), T.("final n_pixel"), grouping);
T_fov2group.("Foreground GMean") = splitapply(weightedMean, T.("Foreground GMean"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background GMean") = splitapply(weightedMean, T.("Background GMean"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding GMean") = splitapply(weightedMean, T.("Responding GMean"), T.("Responding n_pixel"), grouping);
T_fov2group.("final GMean") = splitapply(weightedMean, T.("final GMean"), T.("final n_pixel"), grouping);
T_fov2group.("Foreground RMean") = splitapply(weightedMean, T.("Foreground RMean"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background RMean") = splitapply(weightedMean, T.("Background RMean"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding RMean") = splitapply(weightedMean, T.("Responding RMean"), T.("Responding n_pixel"), grouping);
T_fov2group.("final RMean") = splitapply(weightedMean, T.("final RMean"), T.("final n_pixel"), grouping);

T_fov2group.("Foreground PCC std") = splitapply(@std, T.("Foreground PCC Sum"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background PCC std") = splitapply(@std, T.("Background PCC Sum"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding PCC std") = splitapply(@std, T.("Responding PCC Sum"), T.("Responding n_pixel"), grouping);
T_fov2group.("final PCC std") = splitapply(@std, T.("final PCC Sum"), T.("final n_pixel"), grouping);
T_fov2group.("Foreground Gstd") = splitapply(@std, T.("Foreground GMean"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background Gstd") = splitapply(@std, T.("Background GMean"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding Gstd") = splitapply(@std, T.("Responding GMean"), T.("Responding n_pixel"), grouping);
T_fov2group.("final Gstd") = splitapply(@std, T.("final GMean"), T.("final n_pixel"), grouping);
T_fov2group.("Foreground Rstd") = splitapply(@std, T.("Foreground RMean"), T.("Foreground n_pixel"), grouping);
T_fov2group.("Background Rstd") = splitapply(@std, T.("Background RMean"), T.("Background n_pixel"), grouping);
T_fov2group.("Responding Rstd") = splitapply(@std, T.("Responding RMean"), T.("Responding n_pixel"), grouping);
T_fov2group.("final Rstd") = splitapply(@std, T.("final RMean"), T.("final n_pixel"), grouping);

T_fov2group.("Foreground GRratio") = T_fov2group.("Foreground GMean")./T_fov2group.("Foreground RMean");
T_fov2group.("Background GRratio") = T_fov2group.("Background GMean")./T_fov2group.("Background RMean");
T_fov2group.("Responding GRratio") = T_fov2group.("Responding GMean")./T_fov2group.("Responding RMean");
T_fov2group.("final GRratio") = T_fov2group.("final GMean")./T_fov2group.("final RMean");

%% bar plot compare masks
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
t.Title.String = 'Channel intensity distribution (Group)';


for groupn = 1:height(T_fov2group)
    ax = nexttile(t);
    hold(ax,'on')
    e = gobjects(1,4);        
    e(1) = errorbar(T_fov2group.("Foreground RMean")(groupn),T_fov2group.("Foreground GMean")(groupn),...
        T_fov2group.("Foreground Gstd")(groupn), T_fov2group.("Foreground Gstd")(groupn), ...
        T_fov2group.("Foreground Rstd")(groupn), T_fov2group.("Foreground Rstd")(groupn),  ...
        'Color', cmap(1,:), 'LineStyle', 'none', 'DisplayName','Foreground');
    e(2) = errorbar(T_fov2group.("Background RMean")(groupn),T_fov2group.("Background GMean")(groupn),...
        T_fov2group.("Background Gstd")(groupn), T_fov2group.("Background Gstd")(groupn), ...
        T_fov2group.("Background Rstd")(groupn), T_fov2group.("Background Rstd")(groupn),  ...
        'Color', cmap(2,:), 'LineStyle', 'none', 'DisplayName','Background');   
    e(3) = errorbar(T_fov2group.("final RMean")(groupn),T_fov2group.("final GMean")(groupn),...
        T_fov2group.("final Gstd")(groupn), T_fov2group.("final Gstd")(groupn), ...
        T_fov2group.("final Rstd")(groupn), T_fov2group.("final Rstd")(groupn),  ...
        'Color', cmap(6,:), 'LineStyle', 'none', 'DisplayName','Foreground gating');     
    e(4) = errorbar(T_fov2group.("Responding RMean")(groupn),T_fov2group.("Responding GMean")(groupn),...
        T_fov2group.("Responding Gstd")(groupn), T_fov2group.("Responding Gstd")(groupn), ...
        T_fov2group.("Responding Rstd")(groupn), T_fov2group.("Responding Rstd")(groupn),  ...
        'Color', cmap(5,:), 'LineStyle', 'none', 'DisplayName','Correlation filter'); 
    ax.Title.String = T_fov2group.Name(groupn);
    ax.XAxis.Label.String = 'Red Intensity';
    ax.YAxis.Label.String = 'Green Intensity';
    ax.XAxis.Limits = [0,2000];
    ax.YAxis.Limits = [0,2000]; 
    legend(e);
end

 

%% Compare mask
n = 6;
idx = 1;
figure('Name', T.Name(idx))
sgtitle(T.Name(idx))
subplot(n,1,1)
imshow(T.("Bin final G"){idx}, [20 300]);
title("GMap")

subplot(n,1,2)
imshow(T.("Background Mask"){idx}, []);
title("Background Mask")

subplot(n,1,3)
imshow(T.("final Mask"){idx}, []);
title("Foreground gating Mask")

subplot(n,1,4)
imshow(T.("Foreground Mask"){idx}, []);
title("Manual Mask")

subplot(n,1,5)
imshow(T.("Responding Mask"){idx}, []);
title("Responding Mask")

subplot(n,1,6)
imshow(T.("Responding Mask 2"){idx}, []);
title("Responding Mask 2")

%% Simulation
sz = [32,512];
Ig = zeros(sz);
Ir = zeros(sz);
Igroup = zeros(sz);

for i = 1:numel(Ig)
    rnd = randi(100);
    if rnd < 60
        n_group = 1;
    elseif (60 <= rnd) && (rnd < 95)
        n_group = 2;
    else 
        n_group = 3;
    end
    [Ig(i), Ir(i)] = intensity(n_group);
    Igroup(i) = n_group;
end
%
[Ig,idx] = sort(Ig(:));
Ir = Ir(idx);
Igroup = Igroup(idx);

figure()
subplot(4,1,1)
imshow(reshape(Ig, sz),[])
subplot(4,1,2)
imshow(reshape(Ir, sz),[])
subplot(4,1,3)
imshow(reshape(Igroup, sz),[])
subplot(4,1,4), hold on
s(1) = scatter(Ig(Igroup == 1),Ir(Igroup == 1),5,'filled','DisplayName','background');
s(2) = scatter(Ig(Igroup == 2),Ir(Igroup == 2),5,'filled','DisplayName','membrane');
s(3) = scatter(Ig(Igroup == 3),Ir(Igroup == 3),5,'filled','DisplayName','clump');
legend(s,'location','northwest')
xlim([0,4096])
ylim([0,4096])
function [g, r] = intensity(N)
    g2r = 0.001;
    r2g = 0.001;
    switch N
        case 1 % background
            mu = 20;
            sigma = 10;
            g0 = normrnd(mu, sigma);
            r0 = normrnd(mu, sigma);
            g = g0 + r2g * r0;
            r = r0 + g2r * g0;
        case 2 % membrane
            mu = 180;
            sigma = 200;
            alpha = 2;
            g0 = normrnd(mu, sigma);
            r0 = g0*alpha + abs(normrnd(0, g0*alpha));
            g = g0 + r2g * r0;
            r = r0 + g2r * g0; 
        case 3 % clump
            mu = 100;
            sigma = 200;
            alpha = 10;
            g0 = normrnd(mu, sigma);
            r0 = g0*alpha + abs(normrnd(0, g0*alpha));
            g = g0 + r2g * r0;
            r = r0 + g2r * g0;  
    end
    r = abs(r);
    g = abs(g);
end

function R = pearsons2d(A, B)
    meanA = mean(A(:));
    meanB = mean(B(:));
    R = (A - meanA ).*(B - meanB)./sqrt(sum((A - meanA).^2,'all') .* sum((B - meanB).^2,'all'));
end

function name = parsename(s)
    parsedname = split(s,{'\2p Fstim video\','_P2'});
    name = parsedname(2);
end

function scatter_COL = heatscatter(X, Y, numbins)
    [values, centers] = hist3([X Y], [numbins, numbins]);
    centers_X = centers{1,1};
    centers_Y = centers{1,2};
    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);
    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end
    scatter_COL = zeros(length(X), 1);
    onepercent = round(length(X) / 100);
      for i = 1:length(X)
        if (mod(i,onepercent) == 0)
            fprintf('.');
        end            
        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;
        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1));
        if (~isempty(last_lower_X))
            last_lower_X = last_lower_X(end);
        else
            last_higher_X = find(c_X <= bins_X(:,2));
            if (~isempty(last_higher_X))
                last_higher_X = last_higher_X(1);
            end
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end
        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;
        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1));
        if (~isempty(last_lower_Y))
            last_lower_Y = last_lower_Y(end);
        else
            last_higher_Y = find(c_Y <= bins_Y(:,2));
            if (~isempty(last_higher_Y))
                last_higher_Y = last_higher_Y(1);
            end
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end
        scatter_COL(i) = values(id_X, id_Y);
    
      end
    
end