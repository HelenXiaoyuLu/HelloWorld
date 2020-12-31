%% Load job output 
clear;
clc;
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI';
j.output = load(fullfile(dirp, '20200206 bg roi mask 197-300'));
 
%% Pixel distribution
T_fov = j.output.T_fov;
% T_fov.Name = splitapply(@parsename, T_fov.Path(:), (1:height(T_fov))'); 
for i = 1:height(T_fov)
    T_fov.fgG{i} = mean(nonzeros(T_fov.GMap{i}.*T_fov.finalMask{i}),'all');
    T_fov.fgR{i} = mean(nonzeros(T_fov.RMap{i}.*T_fov.finalMask{i}),'all');
%     T_fov.bgG{i} = mean(nonzeros(T_fov.GMap{i}.*double(T_fov.bgMask{i})),'all');
%     T_fov.bgR{i} = mean(nonzeros(T_fov.RMap{i}.*double(T_fov.bgMask{i})),'all');
%     T_fov.clpG{i} = mean(nonzeros(T_fov.GMap{i}.*double(T_fov.clumpMask{i})),'all');
%     T_fov.clpR{i} = mean(nonzeros(T_fov.RMap{i}.*double(T_fov.clumpMask{i})),'all');
end
grouping = findgroups(T_fov.FOV_Name);

%% Plot
cmap = lines;
f = figure('Color', [1, 1, 1]);
t = tiledlayout(f, 'flow', 'Padding', 'none');
t.Title.String = 'Channel intensity distribution';
for groupn = 1:max(grouping)
    ax = nexttile(t);
    hold(ax,'on')
    groupidx = find(grouping == groupn);
    for j = min(groupidx):max(groupidx)
        s = gobjects(1,1);
        forR = reshape(T_fov.RMap{j}.*T_fov.finalMask{j},1,[]);
        forG = reshape(T_fov.GMap{j}.*T_fov.finalMask{j},1,[]);
%         backR = reshape(T_fov.RMap{j}.*double(T_fov.bgMask{j}),1,[]);
%         backG = reshape(T_fov.GMap{j}.*double(T_fov.bgMask{j}),1,[]);
%         clpR = reshape(T_fov.RMap{j}.*double(T_fov.clumpMask{j}),1,[]);
%         clpG = reshape(T_fov.GMap{j}.*double(T_fov.clumpMask{j}),1,[]);
        s(1) = scatter(forR, forG, 5, cmap(5,:),'filled', ...
              'MarkerFaceAlpha', 0.5, 'DisplayName','Foreground');   
%         s(2) = scatter(backR, backG, 5, cmap(1,:),'filled',...
%             'DisplayName','Background');
%         s(3) = scatter(clpR, clpG, 5, cmap(2,:),'filled',...
%             'DisplayName','Clump');    
    end
    legend(s)
    ax.Title.String = T_fov.FOV_Name(j);
    ax.XAxis.Label.String = 'Red Intensity';
    ax.YAxis.Label.String = 'Green Intensity';
    ax.XAxis.Limits = [0,1000];
    ax.YAxis.Limits = [0,200];
end

 