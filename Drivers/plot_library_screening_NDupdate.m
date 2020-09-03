%% Load data
clear
close all
clc
% Load lookup table for libraries
dirp = 'D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Fstim\202002-03 repeated screening towards JEDI';
LookUp = readtable(fullfile(dirp,'Lookup.xlsx'));

% Load lookup table for sequencing result of library screening
LookUpSeq = readtable(fullfile(dirp,'Lookup_seq.xlsx'));

% Load benchmarking result
j = load(fullfile(dirp,'20200126 2pFstim plate1'));
T_group_benchmarking = j.T_out.T_group;
T_well_benchmarking = j.T_out.T_well;
T_fov_benchmarking = j.T_out.T_fov;
CtrlFP = 'cyOFP';

% Load library screening result
files = dir(strcat(dirp,'\2pFstim*.mat'));
T_screening_2pFstim = table();
% for each screening plate
for i=1:length(files)
% Fill in plate information
    j.output = load(fullfile(files(i).folder,files(i).name),'T_lib','T_fov','T_well','T_group');
    T_fov = j.output.T_fov;
    T_well = j.output.T_well;
    T_group = j.output.T_group;
    numplate = split(j.output.T_lib.path,{'plate'});
    nwell = height(T_well);
% Calculate Photostability
    TraceRawNormSync = T_fov.TraceRaw.synchronize('dt', 0.001);
    T_fov.TraceRawNorm = TraceRawNormSync ./ TraceRawNormSync.integral(0, 0.01) / 100;   
    grouping = findgroups(T_fov.Well);
    T_fov2Well = table();
    T_fov2Well.Well = splitapply(@unique, T_fov.Well, grouping);
    T_fov2Well.TraceRawNormMean = splitapply(@mean, T_fov.TraceRawNorm, T_fov.Area, grouping);
    T_well = innerjoin(T_fov2Well, T_well, 'Keys', 'Well');    
    tmax = min(arrayfun(@(s) max(s.x), T_well.TraceRawNormMean));
    T_well.Photostability = T_well.TraceRawNormMean.integral(0, tmax) / tmax;    
    
% Calculate DI
    T_well.("Detectability Index Short Stim") = abs(T_well.("dFF0 Short Stim").*sqrt(T_well.GRRatioMean));
    T_well.("Detectability Index Long Stim") =  abs(T_well.("dFF0 Long Stim").*sqrt(T_well.GRRatioMean));
    plateid = str2double(numplate{end});
    if isnan(plateid)
        numplate = split(files(i).name,{'Plate','.mat'});
        plateid = str2double(numplate{end-1});
    end
    T_well.PlateID(:) = repmat(plateid,nwell,1);
    
% Fill in the sequence if known
    issequenced = find(LookUpSeq.Plate == plateid);
    isonpath = intersect(find(LookUpSeq.Plate == plateid),find(LookUpSeq.IsInEvolPath == 1));
    T_well.Sequence = findseq(LookUpSeq(issequenced,:), ...
        T_well.wellName);
% Check if the protein is on evolution path
    T_well.IsOnEvolPath = findevol(LookUpSeq(isonpath,:),...
        T_well.wellName);
% Merge tables
    T = innerjoin(T_group,T_well,'key','Group');
    T.ConstructName = strcat(T.Name,'_P',num2str(T.PlateID),T.wellName);
    L(i).T_fov = T_fov;
    L(i).T_well = T;
    L(i).T_group = T_group;
    T_screening_2pFstim = [T_screening_2pFstim; T];    
end

%% Standardize names and regroup data
% sort rows by construct name
T_screening = sortrows(T_screening_2pFstim,'ConstructName','ascend');
T_screening.Properties.RowNames = T_screening.ConstructName;

% Find parent for each library
idx.GFPLib = find(contains(T_screening.ConstructName,'-CT'));
idx.VSDLib = find(contains(T_screening.ConstructName,'-N124-'));
idx.ctrl = setdiff(1:height(T_screening),[idx.GFPLib;idx.VSDLib]);
T_screening.Parent(idx.ctrl) = splitapply(@libplot.fzmatch,...
    T_screening.ConstructName(idx.ctrl)',1:numel(idx.ctrl));
splitnames = split(T_screening.ConstructName(idx.GFPLib),'-CT');
T_screening.Parent(idx.GFPLib) = splitapply(@libplot.fzmatch,...
    splitnames(:,1)', 1:numel(idx.GFPLib));
splitnames = split(T_screening.ConstructName(idx.VSDLib),'-N124');
T_screening.Parent(idx.VSDLib) = splitapply(@libplot.fzmatch,...
    splitnames(:,1)',1:numel(idx.VSDLib));

% Standardized the name for each library 
findLegalName = @(x) string(LookUp.StandardName{strcmp(x,LookUp.ConstructName)});
T_screening.Library(idx.ctrl) = T_screening.Parent(idx.ctrl);
T_screening.Library([idx.GFPLib;idx.VSDLib]) = splitapply(findLegalName,...
    T_screening.ConstructName([idx.GFPLib;idx.VSDLib])', ...
    1:length([idx.GFPLib;idx.VSDLib]));

%% Pull stats for benchmarking 
TraceRawNormSync = T_fov_benchmarking.TraceRaw.synchronize('dt', 0.001);
T_fov_benchmarking.TraceRawNorm = TraceRawNormSync ./ TraceRawNormSync.integral(0, 0.01) / 100;

% Pull from FOV to Well
grouping = findgroups(T_fov_benchmarking.Well);
T_fov2Well = table(); 
T_fov2Well.Well = splitapply(@unique, T_fov_benchmarking.Well, grouping);
T_fov2Well.TraceRawNormMean = splitapply(@mean, T_fov_benchmarking.TraceRawNorm, T_fov_benchmarking.Area, grouping);
T_well_benchmarking = innerjoin(T_fov2Well, T_well_benchmarking, 'Keys', 'Well');
tmax = min(arrayfun(@(s) max(s.x), T_well_benchmarking.TraceRawNormMean));
T_well_benchmarking.Photostability = T_well_benchmarking.TraceRawNormMean.integral(0, tmax) / tmax;
T_well_benchmarking.("Detectability Index Long Stim") = abs(T_well_benchmarking.("dFF0 Long Stim").*sqrt(T_well_benchmarking.GRRatioMean));
T_well_benchmarking.("Detectability Index Short Stim") = abs(T_well_benchmarking.("dFF0 Short Stim").*sqrt(T_well_benchmarking.GRRatioMean));

% Pull from Well to Group
grouping = findgroups(T_well_benchmarking.Group);
weightedMean = @(x, w) sum(x.*w, 'all')./sum(w, 'all');
T_well2Group = table();
T_well2Group.Group = splitapply(@unique, T_well_benchmarking.Group, grouping);
T_well2Group.PhotostabilityWellMean = splitapply(@mean, T_well_benchmarking.Photostability, grouping);
T_well2Group.PhotostabilityWellStd = splitapply(@std, T_well_benchmarking.Photostability, grouping);
T_well2Group.Group = splitapply(@unique, T_well_benchmarking.Group, grouping);
T_well2Group.("dFF0 Short Stim") = splitapply(weightedMean, T_well_benchmarking.("dFF0 Short Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("dFF0 Short Stim STD") = splitapply(@std, T_well_benchmarking.("dFF0 Short Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("dFF0 Long Stim") = splitapply(weightedMean, T_well_benchmarking.("dFF0 Long Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("dFF0 Long Stim STD") = splitapply(@std, T_well_benchmarking.("dFF0 Long Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("Detectability Index Short Stim") = splitapply(weightedMean, T_well_benchmarking.("Detectability Index Short Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("Detectability Index Short Stim STD") = splitapply(@std, T_well_benchmarking.("Detectability Index Short Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("Detectability Index Long Stim") = splitapply(weightedMean, T_well_benchmarking.("Detectability Index Long Stim"), T_well_benchmarking.Area, grouping);
T_well2Group.("Detectability Index Long Stim STD") = splitapply(@std, T_well_benchmarking.("Detectability Index Long Stim"), T_well_benchmarking.Area, grouping);

T_group_benchmarking = innerjoin(T_well2Group, T_group_benchmarking, 'Keys', 'Group');
% Label each row
for i = 1:height(T_group_benchmarking)
    groupname = split(T_group_benchmarking.Name(i),['-',CtrlFP]);
    T_group_benchmarking.("Legal Name")(i) = groupname(1);
    T_group_benchmarking.Properties.RowNames(i) = groupname(1);
end

%% Plot evolutionary path from benchmarking
% clearvars e l 
sel = ["ASAP1","ASAP2s","ASAP2s-H152E","ASAP2s-T207H","JEDI-1P","JEDI-2P"];
f = figure();
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
pp1 = "Detectability Index Short Stim";
pp1std = "Detectability Index Short Stim STD";
pp2 = "PhotostabilityWellMean";
pp2std = "PhotostabilityWellStd";
norm2ctrl = "ASAP1";
DIAUC = [1, 2, 4];
% AUCoperation = "DI*$$\sqrt{AUC}$$=";
AUCoperation = "DI*AUC=";
startnodes = [1,2,2,3,3,4];
endnodes = [2,3,4,5,6,6];
note = ["VSD","Linker","GFP","VSD+Linker","VSD+Linker+GFP","VSD+Linker"];
xval = zeros(1,numel(sel));
xstd = zeros(1,numel(sel));
yval = zeros(1,numel(sel));
ystd = zeros(1,numel(sel));
for i = 1:numel(sel)  
    xval(i) = T_group_benchmarking.(pp1)(sel(i))./T_group_benchmarking.(pp1)(norm2ctrl);
    xstd(i) = T_group_benchmarking.(pp1std)(sel(i))./T_group_benchmarking.(pp1)(norm2ctrl);
    yval(i) = T_group_benchmarking.(pp2)(sel(i))./T_group_benchmarking.(pp2)(norm2ctrl);
    ystd(i) = T_group_benchmarking.(pp2std)(sel(i))./T_group_benchmarking.(pp2)(norm2ctrl);
    e(i) = errorbar(ax,xval(i),yval(i),ystd(i),...
        ystd(i),xstd(i),xstd(i),'.','MarkerSize',20);
end

% Plot countour line
for i = 1:numel(DIAUC)
    x = min(xval)-0.5:0.1:max(xval+0.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

% Plot evolution path
xevol = xval(endnodes)-xval(startnodes);
yevol = yval(endnodes)-yval(startnodes);
quiver(ax,xval(startnodes)+0.2*xevol,yval(startnodes)+0.2*yevol,...
    xevol,yevol,'r','MaxHeadSize',0.1,...
    'AutoScaleFactor',0.9)
for i = 1:length(startnodes)
    text(ax, mean([xval(startnodes(i)),xval(endnodes(i))]),...
        mean([yval(startnodes(i)),yval(endnodes(i))]),note(i),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square
axis(ax, [0 5 0 2]);
legend(e,sel,'location','eastoutside')

%% Scatter libraries & pull evolution path from benchmarking
clearvars e l s
sel = ["ASAP1","ASAP2s","ASAP2s-H152E","ASAP2s-T207H","JEDI-1P","JEDI-2P"];
f = figure();
c = lines;
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');

% Scatter library screening by parent
selLib = ["ASAP2s","ASAP2s-H152E","ASAP2s-T207H","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
ctrlidx = strcmp(T_screening.Parent,"ASAP1");
for i = 1:length(selLib)
    libidx = strcmp(T_screening.Parent,selLib(i));
%     xlibval = T_screening.(pp1lib)(libidx)./T_group.(pp1)(norm2ctrl);
%     ylibval = T_screening.(pp2lib)(libidx)./T_group.(pp2)(norm2ctrl);
    xlibval = T_screening.(pp1lib)(libidx)./nanmean(T_screening.(pp1lib)(ctrlidx));
    ylibval = T_screening.(pp2lib)(libidx)./nanmean(T_screening.(pp2lib)(ctrlidx));   
    s(i) = scatter(ax, xlibval, ylibval,'filled','MarkerFaceAlpha',0.5);
    s(i).CData = c(i+1,:);
    s(i).DataTipTemplate.DataTipRows(1) = dataTipTextRow('Parent',cellstr(T_screening.Parent(libidx)));
    s(i).DataTipTemplate.DataTipRows(2) = dataTipTextRow('Library Name',cellstr(T_screening.Library(libidx)));
    s(i).DataTipTemplate.DataTipRows(3) = dataTipTextRow('Well',cellstr(T_screening.wellName(libidx)));
    s(i).DataTipTemplate.DataTipRows(4) = dataTipTextRow('dF/F0 Short',T_screening.("dFF0 Short Stim")(libidx));
    s(i).DataTipTemplate.DataTipRows(5) = dataTipTextRow('Sequence',T_screening.Sequence(libidx));
end

% Pull data for evolutaionary path
pp1 = "Detectability Index Short Stim";
pp1std = "Detectability Index Short Stim STD";
pp2 = "PhotostabilityWellMean";
pp2std = "PhotostabilityWellStd";
norm2ctrl = "ASAP1";
DIAUC = [1, 2, 4, 8, 16];
% AUCoperation = "DI*$$\sqrt{AUC}$$=";
AUCoperation = "DI*AUC=";
startnodes = [1,2,2,3,3,4];
endnodes = [2,3,4,5,6,6];
note = ["VSD","Linker","GFP","VSD+Linker","VSD+Linker+GFP","VSD+Linker"];
xval = zeros(1,numel(sel));
xstd = zeros(1,numel(sel));
yval = zeros(1,numel(sel));
ystd = zeros(1,numel(sel));
for i = 1:numel(sel)  
    xval(i) = T_group_benchmarking.(pp1)(sel(i))./T_group_benchmarking.(pp1)(norm2ctrl);
    xstd(i) = T_group_benchmarking.(pp1std)(sel(i))./T_group_benchmarking.(pp1)(norm2ctrl);
    yval(i) = T_group_benchmarking.(pp2)(sel(i))./T_group_benchmarking.(pp2)(norm2ctrl);
    ystd(i) = T_group_benchmarking.(pp2std)(sel(i))./T_group_benchmarking.(pp2)(norm2ctrl);
    e(i) = errorbar(ax,xval(i),yval(i),ystd(i),...
        ystd(i),xstd(i),xstd(i),'.','MarkerSize',40,'color',c(i,:));
end
% Plot countour line
for i = 1:numel(DIAUC)
    x = min(xval)-0.5:0.1:max(xval+2.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

% Plot evolution path
xevol = xval(endnodes)-xval(startnodes);
yevol = yval(endnodes)-yval(startnodes);
quiver(ax,xval(startnodes)+0.2*xevol,yval(startnodes)+0.2*yevol,...
    xevol,yevol,'r','MaxHeadSize',0.1,...
    'AutoScaleFactor',0.9)
for i = 1:length(startnodes)
    text(ax, mean([xval(startnodes(i)),xval(endnodes(i))]),...
        mean([yval(startnodes(i)),yval(endnodes(i))]),note(i),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 9 0 2]);
title({'Evolution path pulled from benchmarking'; ...
    'Unnormalized for plate-to-plate variation'})
legend([e,s],[sel,selLib],'location','eastoutside')

%% Plot controls from different library
T_onboardctrl = table();
onboardctrls = ["ASAP1","ASAP2s","ASAP2s-T207H","JEDI-1P"];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
for i = 1:numel(L)
    for j = 1:length(onboardctrls)
        boolctrl = strcmp(L(i).T_well.Sequence,onboardctrls(j));
        T_ctrl_temp = L(i).T_well(boolctrl,:);
        if any(boolctrl)
            T_onboardctrl = [T_onboardctrl;T_ctrl_temp];
        end
    end
end
T_onboardctrl = rmmissing(T_onboardctrl,'DataVariables',{char(pp1lib), char(pp2lib)});
T_onboardctrl(find(T_onboardctrl.("dFF0 Short Stim")>-0.03),:) = [];
   
f = figure();
clearvars l s
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
groupbyplate = findgroups(T_onboardctrl.PlateID);
c = lines;
for i = 1:max(groupbyplate)
    idx = find(groupbyplate==i);
    l(i) = plot(ax,T_onboardctrl.(pp1lib)(idx),...
        T_onboardctrl.(pp2lib)(idx),'--','color',c(i,:));
end
groupbyctrl = findgroups(T_onboardctrl.Sequence);
for i = 1:numel(onboardctrls)
    idx = find(groupbyctrl==i);
    s(i) = scatter(ax,T_onboardctrl.(pp1lib)(idx),...
        T_onboardctrl.(pp2lib)(idx),'filled','MarkerFaceAlpha',0.5);
end
legend(ax,[s,l],[splitapply(@unique,T_onboardctrl.Sequence,groupbyctrl);...
    splitapply(@unique,T_onboardctrl.PlateID,groupbyplate)],...
    'location','eastoutside');
xlabel(ax,{pp1});
ylabel(ax,{pp2});
axis(ax, [0 inf 0 inf]);
title(ax,{'On-board Control Performance';'NaN discarded'; ...
    'Unnormalized for plate-to-plate variation'});

%% Plot controls from different library, normalize to a ref GEVI
ref = "ASAP2s-T207H";
f = figure();
clearvars l s 
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
c = lines;
for i = 1:numel(L)
    idx = find(T_onboardctrl.PlateID==i);
    refidx = intersect(find(T_onboardctrl.Sequence==ref),idx);
    if ~isempty(refidx)
        T_onboardctrl.(strcat("Normalized ",pp1lib))(idx) = T_onboardctrl.(pp1lib)(idx)./T_onboardctrl.(pp1lib)(refidx);
        T_onboardctrl.(strcat("Normalized ",pp2lib))(idx) = T_onboardctrl.(pp2lib)(idx)./T_onboardctrl.(pp2lib)(refidx);
        l(i) = plot(ax,T_onboardctrl.(strcat("Normalized ",pp1lib))(idx),...
        T_onboardctrl.(strcat("Normalized ",pp2lib))(idx),'--','color',c(i,:));
    elseif any(idx)
%         if i-1 > 0
%             samedayrefidx = intersect(find(T_onboardctrl.Sequence==ref),find(T_onboardctrl.PlateID==i-1));
%         else
            samedayrefidx = intersect(find(T_onboardctrl.Sequence==ref),find(T_onboardctrl.PlateID==i+1));
%         end
        T_onboardctrl.(strcat("Normalized ",pp1lib))(idx) = T_onboardctrl.(pp1lib)(idx)./T_onboardctrl.(pp1lib)(samedayrefidx);
        T_onboardctrl.(strcat("Normalized ",pp2lib))(idx) = T_onboardctrl.(pp2lib)(idx)./T_onboardctrl.(pp2lib)(samedayrefidx);
        l(i) = plot(ax,T_onboardctrl.(strcat("Normalized ",pp1lib))(idx),...
        T_onboardctrl.(strcat("Normalized ",pp2lib))(idx),'--','color',c(i,:));        
    end
end
l = l(find(isgraphics(l(:))));
groupbyctrl = findgroups(T_onboardctrl.Sequence);
for i = 1:numel(onboardctrls)
    idx = find(groupbyctrl==i);
    s(i) = scatter(ax,T_onboardctrl.(strcat("Normalized ",pp1lib))(idx),...
        T_onboardctrl.(strcat("Normalized ",pp2lib))(idx),...
        'filled','MarkerFaceAlpha',0.5);
end
xlabel(ax,{pp1;strcat("Normalized to ",ref)})
ylabel(ax,{pp2;strcat("Normalized to ",ref)})
axis(ax,'equal')
title(ax,{'On-board Control Performance';'NaN discarded ';strcat("Normalized to ",ref)})
legend(ax,[s,l],[splitapply(@unique,T_onboardctrl.Sequence,groupbyctrl);...
    splitapply(@unique,T_onboardctrl.PlateID,groupbyplate)],...
    'location','eastoutside');

%% Scatter libraries & pull evolution path from library screening 
% Pull controls on evolutionary path to a table
T_ctrlinlib = T_screening(find(T_screening.IsOnEvolPath),:);
T_ctrlinlib(find(T_ctrlinlib.("dFF0 Short Stim")>-0.03),:) = [];

% T_ctrlinlib = T_onboardctrl;
selpp = {'Area_T_well', 'GMean', 'RMean', 'GRRatioMean',...
    'dFF0 Short Stim', 'dFF0 Long Stim', 'Photostability',...
    'Detectability Index Short Stim','Detectability Index Long Stim', 'Sequence'};
T_ctrlinlib_groupall = libplot.groupWells(T_ctrlinlib,'groupOn','Sequence',...
    'selectProperties',selpp);
T_ctrlinlib_groupall.Properties.RowNames = T_ctrlinlib_groupall.Sequence;
selectedctrls = find(T_ctrlinlib_groupall.n_well>1);
T_ctrlinlib_group = T_ctrlinlib_groupall(selectedctrls,:);

% Plot evol path
f = figure();
clearvars e l s
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
pp1 = "Detectability Index Short Stim_Mean";
pp1std = "Detectability Index Short Stim_STD";
pp2 = "Photostability_Mean";
pp2std = "Photostability_STD";
norm2ctrl = "ASAP1";
DIAUC = [1, 2, 4, 8, 16];
% AUCoperation = "DI*$$\sqrt{AUC}$$=";
AUCoperation = "DI*AUC=";
xval = zeros(1,numel(sel));
xstd = zeros(1,numel(sel));
yval = zeros(1,numel(sel));
ystd = zeros(1,numel(sel));
nevol = height(T_ctrlinlib_group);
for i = 1:nevol 
    xval(i) = T_ctrlinlib_group.(pp1)(i)./T_ctrlinlib_group.(pp1)(norm2ctrl);
    xstd(i) = T_ctrlinlib_group.(pp1std)(i)./T_ctrlinlib_group.(pp1)(norm2ctrl);
    yval(i) = T_ctrlinlib_group.(pp2)(i)./T_ctrlinlib_group.(pp2)(norm2ctrl);
    ystd(i) = T_ctrlinlib_group.(pp2std)(i)./T_ctrlinlib_group.(pp2)(norm2ctrl);
    e(i) = errorbar(ax,xval(i),yval(i),ystd(i),...
        ystd(i),xstd(i),xstd(i),'.','MarkerSize',20);
end

% Plot countour line
for i = 1:numel(DIAUC)
    x = min(xval)-0.5:0.1:max(xval+2.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

%% Scatter library screening by parent
selLib = ["ASAP2s","ASAP2s-T207H","ASAP2s-H152E","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
ctrlidx = strcmp(T_screening.Parent,"ASAP1");
for i = 1:length(selLib)
    libidx = strcmp(T_screening.Parent,selLib(i));
%     xlibval = T_screening.(pp1lib)(libidx)./T_group.(pp1)(norm2ctrl);
%     ylibval = T_screening.(pp2lib)(libidx)./T_group.(pp2)(norm2ctrl);
    xlibval = T_screening.(pp1lib)(libidx)./nanmean(T_screening.(pp1lib)(ctrlidx));
    ylibval = T_screening.(pp2lib)(libidx)./nanmean(T_screening.(pp2lib)(ctrlidx));   
    s(i) = scatter(ax, xlibval, ylibval,'filled','MarkerFaceAlpha',0.5);
    s(i).DataTipTemplate.DataTipRows(1) = dataTipTextRow('Parent',cellstr(T_screening.Parent(libidx)));
    s(i).DataTipTemplate.DataTipRows(2) = dataTipTextRow('Library Name',cellstr(T_screening.Library(libidx)));
    s(i).DataTipTemplate.DataTipRows(3) = dataTipTextRow('Well',cellstr(T_screening.wellName(libidx)));
    s(i).DataTipTemplate.DataTipRows(4) = dataTipTextRow('dFF0 Short Stim',T_screening.("dFF0 Short Stim")(libidx));
    s(i).DataTipTemplate.DataTipRows(5) = dataTipTextRow('Sequence',T_screening.Sequence(libidx));
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 9 0 2]);
title({'Evolution path pulled from libraries'; ...
    'Unnormalized for plate-to-plate variation'})
legend([e,s],[sel,selLib],'location','eastoutside')

%% Scatter libraries & pull evolution path from library screening,...
%  normalize to selected on-board control 
ref = "ASAP2s-T207H";
refOrigin = "ASAP1";
f = figure();
clearvars l s 
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');
pp1lib = "Detectability Index Short Stim_Mean";
pp1libstd = "Detectability Index Short Stim_STD";
pp2lib = "TgtAUC_sig_NormaltoRandT_Mean";
pp2libstd = "TgtAUC_sig_NormaltoRandT_STD";
norm2ctrl = ref;
DIAUC = [1 2 4 8];
% AUCoperation = "DI*$$\sqrt{AUC}$$=";
AUCoperation = "DI*AUC=";
xval = zeros(1,numel(sel));
xstd = zeros(1,numel(sel));
yval = zeros(1,numel(sel));
ystd = zeros(1,numel(sel));
nevol = height(T_ctrlinlib_group);
xamplify = T_ctrlinlib_group.(pp1lib)(ref)/T_ctrlinlib_group.(pp1lib)(refOrigin);
yamplify = T_ctrlinlib_group.(pp2lib)(ref)/T_ctrlinlib_group.(pp2lib)(refOrigin);

for i = 1:nevol 
    
    xval(i) = T_ctrlinlib_group.(pp1lib)(i)./T_ctrlinlib_group.(pp1lib)(norm2ctrl)*xamplify;
    xstd(i) = T_ctrlinlib_group.(pp1libstd)(i)./T_ctrlinlib_group.(pp1lib)(norm2ctrl)*xamplify;
    yval(i) = T_ctrlinlib_group.(pp2lib)(i)./T_ctrlinlib_group.(pp2lib)(norm2ctrl)*yamplify;
    ystd(i) = T_ctrlinlib_group.(pp2libstd)(i)./T_ctrlinlib_group.(pp2lib)(norm2ctrl)*yamplify;
    e(i) = errorbar(ax,xval(i),yval(i),ystd(i),...
        ystd(i),xstd(i),xstd(i),'.','MarkerSize',20);
end

% Plot countour line
for i = 1:numel(DIAUC)
    x = max(min(xval)-0.5,0.1):0.1:max(xval+0.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

pp1lib = "Detectability Index Short Stim";
pp2lib = "TgtAUC_sig_NormaltoRandT";
 
% scatter libraries (normalized)
for i = 1:numel(L)
    idx = find(T_screening.Plate==i);
    refidx = intersect(find(T_screening.Sequence==ref),idx);
    if ~isempty(refidx)
        T_screening.(strcat("Normalized ",pp1lib))(idx) = T_screening.(pp1lib)(idx)./T_screening.(pp1lib)(refidx)*xamplify;
        T_screening.(strcat("Normalized ",pp2lib))(idx) = T_screening.(pp2lib)(idx)./T_screening.(pp2lib)(refidx)*yamplify;
    elseif any(idx)
        samedayrefidx = intersect(find(T_screening.Sequence==ref),find(T_screening.Plate==i-1));
        T_screening.(strcat("Normalized ",pp1lib))(idx) = T_screening.(pp1lib)(idx)./T_screening.(pp1lib)(samedayrefidx)*xamplify;
        T_screening.(strcat("Normalized ",pp2lib))(idx) = T_screening.(pp2lib)(idx)./T_screening.(pp2lib)(samedayrefidx)*yamplify;
    end
end

% Scatter library screening by parent
selLib = ["ASAP2s","ASAP2s-T207H","ASAP2s-H152E","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = strcat("Normalized ",pp1lib);
pp2lib = strcat("Normalized ",pp2lib);
for i = 1:length(selLib)
    libidx = strcmp(T_screening.Parent,selLib(i));
    xlibval = T_screening.(pp1lib)(libidx);
    ylibval = T_screening.(pp2lib)(libidx);   
    s(i) = scatter(ax, xlibval, ylibval,'filled','MarkerFaceAlpha',0.5);
    s(i).DataTipTemplate.DataTipRows(1) = dataTipTextRow('Parent',cellstr(T_screening.Parent(libidx)));
    s(i).DataTipTemplate.DataTipRows(2) = dataTipTextRow('Library Name',cellstr(T_screening.Library(libidx)));
    s(i).DataTipTemplate.DataTipRows(3) = dataTipTextRow('Well',cellstr(T_screening.wellName(libidx)));
    s(i).DataTipTemplate.DataTipRows(4) = dataTipTextRow('dF/F0 Short',T_screening.dFF0_shortStim(libidx));
    s(i).DataTipTemplate.DataTipRows(5) = dataTipTextRow('Sequence',T_screening.Sequence(libidx));
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 8 0 5]);
title({'Evolution path pulled from libraries'; ...
    'Normalized for plate-to-plate variation'})
legend([e,s],[sel,selLib],'location','eastoutside')

%% Fun collections
function seq = findseq(T_lookup,welllist)
    seq = strings(numel(welllist),1);
    for i = 1:numel(welllist)
        wellname = welllist(i);
        if iscell(wellname)
            wellname = string(wellname);
        end
        findwell = strcmp(wellname,T_lookup.Well);
        if any(findwell)
            seq(i) = string(T_lookup.SequencingResult{findwell});
        else
            seq(i) = "";
        end   
    end
end

function pname = findevol(T_lookup,welllist)
    pname = zeros(numel(welllist),1);
    for i = 1:numel(welllist)
        wellname = welllist(i);
        if iscell(wellname)
            wellname = string(wellname);
        end
        findwell = strcmp(wellname,T_lookup.Well);
        if any(findwell)
            pname(i) = 1;
        end    
    end
end