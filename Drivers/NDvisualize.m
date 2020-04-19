% Visualize the N-D screening results
clear
close all
CBheaders = readcell('D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20200102 Summary of Libray 001-017','FileType','spreadsheet','Sheet','mAppleGRAB_combine_rank');
CBscores = readmatrix('D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\20200102 Summary of Libray 001-017','FileType','spreadsheet','Sheet','mAppleGRAB_combine_rank');
varNames = CBheaders(2:end,2);
CBheaders = CBheaders(1,1:end);

% Group data by construct name
Controls = {'X2-cpmAppleGRAB ';'X2-cpmAppleGRAB-A148D ';'X2-cpmAppleGRAB-A148G ';...
    'X2-cpmAppleGRAB-A148S '; 'X2-cpmAppleGRAB-A148T ';'X2-cpmAppleGRAB-P149G '; 'X2-cpmAppleGRAB-L392G ';...
    'X2-cpmAppleGRAB-L392N '; 'X2-cpmAppleGRAB-L392Q '; 'X2-cpmAppleGRAB-L392S '; 'X2-cpmOrange2-174-173-N-1-C-1 '};
Libs = {'A148A-';'A148D-';'A148G-';'A148S-';'A148T-';...
    'P149P-';'P149G-';'L392L_';'L392G_';'L392N_'; 'L392Q_'; 'L392S_'};
for i = 1:length(Controls)
    idx = find(cellfun(@(s) ~isempty(strfind(s,Controls{i})), varNames));
    Controls{i,2} = CBscores(idx,4); % Ch2/Ch3 brightness
    Controls{i,3} = CBscores(idx,6); % 1ms dF/F0
    Controls{i,4} = CBscores(idx,7); % 20ms dF/F0
end

for i = 1:length(Libs)
    idx = find(cellfun(@(s) ~isempty(strfind(s,Libs{i})), varNames));
    Libs{i,2} = CBscores(idx,4); % Ch2/Ch3 brightness
    Libs{i,3} = CBscores(idx,6); % 1ms dF/F0
    Libs{i,4} = CBscores(idx,7); % 20ms dF/F0
end

%% 3D plot
figure()
cmap = lines;
subplot(1,3,1)
hold on
for l = 1:5
    scatter3(Libs{l,3},Libs{l,4},Libs{l,2},40,cmap(l,:),'.');
end

subplot(1,3,2)
hold on
for l = 6:7
    scatter3(Libs{l,3},Libs{l,4},Libs{l,2},40,cmap(l,:),'.');
end

subplot(1,3,3)
hold on
for l = 8:12
    scatter3(Libs{l,3},Libs{l,4},Libs{l,2},40,cmap(l,:),'.');
end

%% 3D surface plot
figure()
hold on
for l = 8:12
  %  plot3(Libs{l,3},Libs{l,4},Libs{l,2},'.-');
    tri = delaunay(Libs{l,3},Libs{l,4});
  %  plot(Libs{l,3},Libs{l,4},'.-')
    [r,c] = size(tri);
    trisurf(tri, Libs{l,3}, Libs{l,4}, Libs{l,2},'FaceColor',cmap(l,:));
    axis vis3d
end
shading interp
colorbar EastOutside

%% 2D plot F-dF20ms
figure()
cmap = lines;
subplot(1,3,1)
hold on
for l = 1:5
    scatter(Libs{l,2},Libs{l,4},40,cmap(l,:),'.');
end

subplot(1,3,2)
hold on
for l = 6:7
    scatter(Libs{l,2},Libs{l,4},40,cmap(l,:),'.');
end

subplot(1,3,3)
hold on
for l = 8:12
    scatter(Libs{l,2},Libs{l,4},40,cmap(l,:),'.');
end

%% 2D plot F-dF1ms
figure()
cmap = lines;
subplot(1,3,1)
hold on
for l = 1:5
    scatter(Libs{l,2},Libs{l,3},40,cmap(l,:),'.');
end

subplot(1,3,2)
hold on
for l = 6:7
    scatter(Libs{l,2},Libs{l,3},40,cmap(l,:),'.');
end

subplot(1,3,3)
hold on
for l = 8:12
    scatter(Libs{l,2},Libs{l,3},40,cmap(l,:),'.');
end

%% 2D plot DI1ms DI20ms
figure()
cmap = lines;
subplot(1,3,1)
hold on
for l = 1:5
    scatter(Libs{l,3}.*sqrt(Libs{l,2}),Libs{l,4}.*sqrt(Libs{l,2}),40,cmap(l,:),'.');
end

subplot(1,3,2)
hold on
for l = 6:7
    scatter(Libs{l,3}.*sqrt(Libs{l,2}),Libs{l,4}.*sqrt(Libs{l,2}),40,cmap(l,:),'.');
end

subplot(1,3,3)
hold on
for l = 8:12
    scatter(Libs{l,3}.*sqrt(Libs{l,2}),Libs{l,4}.*sqrt(Libs{l,2}),40,cmap(l,:),'.');
end