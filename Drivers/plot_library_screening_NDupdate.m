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

%% Load library screening result
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
pp1 = "Detectability Index Short Stim";
pp1std = "Detectability Index Short Stim STD";
pp2 = "PhotostabilityWellMean";
pp2std = "PhotostabilityWellStd";
norm2ctrl = "ASAP1";
startnodes = [1,2,2,3,3,4];
endnodes = [2,3,4,5,6,6];
note = ["VSD","Linker","GFP","VSD+Linker","VSD+Linker+GFP","VSD+Linker"];
T_plot = T_group_benchmarking(sel,:);
[ax,xval,yval,xstd,ystd,e] = plotevolpath(T_plot,pp1,pp2,'pp1std',pp1std,...
    'pp2std',pp2std,'norm2ctrl',norm2ctrl,'startnodes',startnodes,...
    'endnodes',endnodes,'note',note);
legend(e,sel,'location','eastoutside')

DIAUC = [1, 3, 5];
operation = "X.*Y";
x = min(xval)-0.5:0.1:max(xval+2.5);
y = min(yval)-1:0.1:max(yval+2.5);
[ax,C,h] = plotcontour(x,y,'ax',ax,'operation',operation,'LevelList',DIAUC,'note',"DI*AUC");
xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square
axis(ax, [0 5 0 2]);
title('Evolution path pulled from benchmarking');
    
%% Scatter libraries & pull evolution path from *benchmarking*
clearvars e l s
sel = ["ASAP1","ASAP2s","ASAP2s-H152E","ASAP2s-T207H","JEDI-1P","JEDI-2P"];
f = figure();
c = lines;
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');

% Scatter library screening by parent
T_plot = T_screening(find(T_screening.("dFF0 Short Stim")<-0.03),:);
selLib = ["ASAP2s","ASAP2s-H152E","ASAP2s-T207H","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
ctrlidx = find(strcmp(T_plot.Sequence,"ASAP1"));
stip = ["ConstructName","dFF0 Short Stim","GRRatioMean","Photostability","Sequence"];
[ax,s] = scatterbyparent(T_plot,selLib,pp1lib,pp2lib,'ax',ax,'ctrlidx',ctrlidx,...
    'tip',stip,'cmap',c(2:end,:));

% Plot evol path from benchmarking
T_plot = T_group_benchmarking(sel,:);
[ax,xval,yval,xstd,ystd,e] = plotevolpath(T_plot,pp1,pp2,'pp1std',pp1std,...
    'pp2std',pp2std,'ax',ax,'norm2ctrl',norm2ctrl,'startnodes',startnodes,...
    'endnodes',endnodes,'note',note);
legend([e,s],[sel,selLib],'location','eastoutside')

% Plot contour
DIAUC = [1, 3, 5];
operation = "X.*Y";
x = min(xval)-0.5:0.1:max(xval+2.5);
y = min(yval)-1:0.1:max(yval+2.5);
[ax,C,h] = plotcontour(x,y,'ax',ax,'operation',operation,'LevelList',DIAUC,'note',"DI*AUC");

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 9 0 2]);
title({'Evolution path pulled from benchmarking'; ...
    'Unnormalized for plate-to-plate variation'})

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

%% pull evolution path from library screening 
ref = "ASAP2s-T207H";
refOrigin = "ASAP1";
T_ctrlinlib = T_screening(find(T_screening.IsOnEvolPath),:);
T_ctrlinlib(find(T_ctrlinlib.("dFF0 Short Stim")>-0.03),:) = [];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
for i = 1:numel(L)
    idx = find(T_ctrlinlib.PlateID==i);
    refidx = intersect(find(T_ctrlinlib.Sequence==ref),idx);
    if ~isempty(refidx)
        T_ctrlinlib.(strcat("Normalized ",pp1lib))(idx) = T_ctrlinlib.(pp1lib)(idx)./T_ctrlinlib.(pp1lib)(refidx);
        T_ctrlinlib.(strcat("Normalized ",pp2lib))(idx) = T_ctrlinlib.(pp2lib)(idx)./T_ctrlinlib.(pp2lib)(refidx);
    elseif any(idx)
%         if i-1 > 0
%             samedayrefidx = intersect(find(T_onboardctrl.Sequence==ref),find(T_onboardctrl.PlateID==i-1));
%         else
            samedayrefidx = intersect(find(T_ctrlinlib.Sequence==ref),find(T_ctrlinlib.PlateID==i+1));
%         end
        T_ctrlinlib.(strcat("Normalized ",pp1lib))(idx) = T_ctrlinlib.(pp1lib)(idx)./T_ctrlinlib.(pp1lib)(samedayrefidx);
        T_ctrlinlib.(strcat("Normalized ",pp2lib))(idx) = T_ctrlinlib.(pp2lib)(idx)./T_ctrlinlib.(pp2lib)(samedayrefidx);       
    end
end
selpp = {'Area_T_well', 'GMean', 'RMean', 'GRRatioMean',...
    'dFF0 Short Stim', 'dFF0 Long Stim', 'Photostability',...
    'Detectability Index Short Stim','Detectability Index Long Stim', ...
    strcat("Normalized ",pp1lib),'Normalized Photostability','Sequence'};
T_ctrlinlib_groupall = libplot.groupWells(T_ctrlinlib,'groupOn','Sequence',...
    'selectProperties',selpp);
T_ctrlinlib_groupall.Properties.RowNames = T_ctrlinlib_groupall.Sequence;
selectedctrls = find(T_ctrlinlib_groupall.n_well>1);
T_plot = T_ctrlinlib_groupall(selectedctrls,:);

f = figure();
clearvars e l s
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');

% Method 1: normalized to ASAP2s-T207H
pp1 = "Normalized Detectability Index Short Stim_Mean";
pp1std = "Normalized Detectability Index Short Stim_STD";
pp2 = "Normalized Photostability_Mean";
pp2std = "Normalized Photostability_STD";
[ax,nxval,nyval,nxstd,nystd,e] = plotevolpath(T_plot,pp1,pp2,'pp1std',pp1std,...
    'pp2std',pp2std,'norm2ctrl',norm2ctrl,'startnodes',startnodes,...
    'endnodes',endnodes,'note',note,'ax',ax);
legend(e,T_plot.Sequence,'location','eastoutside')

% Method 1 contour
x = min(nxval)-0.5:0.1:max(nxval+1);
y = min(nyval)-1:0.1:max(nyval+1);
[ax,C,h] = plotcontour(x,y,'ax',ax,'operation',operation,'LevelList',DIAUC,'note',"DI*AUC");

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 8 0 2]);
title(ax,{'Evolution path pulled from libraries'; ...
    strcat("normalized to ",ref, " for plate-to-plate variation")})

% Method 2: non-normalized
f2 = figure();
t2 = tiledlayout(f2, 'flow', 'Padding', 'none');
ax2 = nexttile(t2);
hold(ax2, 'on');
pp1 = "Detectability Index Short Stim_Mean";
pp1std = "Detectability Index Short Stim_STD";
pp2 = "Photostability_Mean";
pp2std = "Photostability_STD";
[ax2,nnxval,nnyval,nnxstd,nnystd,e2] = plotevolpath(T_plot,pp1,pp2,'pp1std',pp1std,...
    'pp2std',pp2std,'norm2ctrl',norm2ctrl,'startnodes',startnodes,...
    'endnodes',endnodes,'note',note,'ax',ax2);
legend(e2,T_plot.Sequence,'location','eastoutside')

% Method 2 contour
x = min(nxval)-0.5:0.1:max(nxval+1);
y = min(nyval)-1:0.1:max(nyval+1);
[ax2,C2,h2] = plotcontour(x,y,'ax',ax2,'operation',operation,'LevelList',DIAUC,'note',"DI*AUC");

xlabel(ax2,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax2,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax2, [0 7 0 2]);
title(ax2,{'Evolution path pulled from libraries'; ...
    'Unnormalized for plate-to-plate variation'})

%% Scatter libraries & pull evolution path from library screening *unnormalized*
% Plot evol path
f = figure();
c = lines;
clearvars e l s
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');

% Scatter library screening by parent
refOrigin = "ASAP1";
selLib = ["ASAP2s","ASAP2s-H152E","ASAP2s-T207H","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
ctrlidx = strcmp(T_screening.Parent,refOrigin);
for i = 1:length(selLib)
    libidx = strcmp(T_screening.Parent,selLib(i));
%     xlibval = T_screening.(pp1lib)(libidx)./T_group.(pp1)(norm2ctrl);
%     ylibval = T_screening.(pp2lib)(libidx)./T_group.(pp2)(norm2ctrl);
    xlibval = T_screening.(pp1lib)(libidx)./nanmean(T_screening.(pp1lib)(ctrlidx));
    ylibval = T_screening.(pp2lib)(libidx)./nanmean(T_screening.(pp2lib)(ctrlidx));   
    s(i) = scatter(ax, xlibval, ylibval,'filled','MarkerFaceAlpha',0.5);
    s(i).CData = c(i+1,:);
    s(i).DataTipTemplate.DataTipRows(1) = dataTipTextRow('Parent',cellstr(T_screening.Parent(libidx)));
    s(i).DataTipTemplate.DataTipRows(2) = dataTipTextRow('Construct Name',cellstr(T_screening.ConstructName(libidx)));
    s(i).DataTipTemplate.DataTipRows(3) = dataTipTextRow('dF/F0 Short Stim',T_screening.("dFF0 Short Stim")(libidx));
    s(i).DataTipTemplate.DataTipRows(4) = dataTipTextRow('G/R',T_screening.("GRRatioMean")(libidx));
    s(i).DataTipTemplate.DataTipRows(5) = dataTipTextRow('Photostability (AUC)',T_screening.("Photostability")(libidx));
    s(i).DataTipTemplate.DataTipRows(6) = dataTipTextRow('Sequence',T_screening.Sequence(libidx));
end

% Plot evolutionary path (from previous section)
for i = 1:nevol 
    e(i) = errorbar(ax,nnxval(i),nnyval(i),nnystd(i),...
        nnystd(i),nnxstd(i),nnxstd(i),'.','MarkerSize',30,'color',c(i,:));
end
% Plot evolution path
nnxevol = nnxval(endnodes)-nnxval(startnodes);
nnyevol = nnyval(endnodes)-nnyval(startnodes);
quiver(ax,nnxval(startnodes)+0.2*nnxevol,nnyval(startnodes)+0.2*nnyevol,...
    nnxevol,nnyevol,'r','MaxHeadSize',0.1,...
    'AutoScaleFactor',0.98)

% Plot countour line
for i = 1:numel(DIAUC)
    x = min(xval)-0.2:0.1:max(xval+2.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 9 0 2]);
title({'Evolution path pulled from libraries'; ...
    'Unnormalized for plate-to-plate variation'})
legend([e,s],[sel,selLib],'location','eastoutside')

%% Scatter libraries & pull evolution path from library screening, *normalize* to selected on-board control 
ref = "ASAP2s-T207H";
refOrigin = "ASAP1";
f = figure();
clearvars e l s
t = tiledlayout(f, 'flow', 'Padding', 'none');
ax = nexttile(t);
hold(ax, 'on');

pp1lib = "Detectability Index Short Stim";
pp2lib = "Photostability";
 
% Normalize libraries to on-board control
for i = 1:numel(L)
    idx = find(T_screening.PlateID==i);
    refidx = intersect(find(T_screening.Sequence==ref),idx);
    if ~isempty(refidx) && T_screening.("dFF0 Short Stim")(refidx)<=-0.03
        T_screening.(strcat("Normalized ",pp1lib))(idx) = T_screening.(pp1lib)(idx)./T_screening.(pp1lib)(refidx);
        T_screening.(strcat("Normalized ",pp2lib))(idx) = T_screening.(pp2lib)(idx)./T_screening.(pp2lib)(refidx);
    elseif any(idx)
        samedayrefidx = intersect(find(T_screening.Sequence==ref),find(T_screening.PlateID==i+1));
        T_screening.(strcat("Normalized ",pp1lib))(idx) = T_screening.(pp1lib)(idx)./T_screening.(pp1lib)(samedayrefidx);
        T_screening.(strcat("Normalized ",pp2lib))(idx) = T_screening.(pp2lib)(idx)./T_screening.(pp2lib)(samedayrefidx);
    end
end

% Scatter library screening by parent
selLib = ["ASAP2s","ASAP2s-H152E","ASAP2s-T207H","ASAP2s-H152E-T207H","ASAP2s-H152E-Q397H"];
pp1lib = strcat("Normalized ",pp1lib);
pp2lib = strcat("Normalized ",pp2lib);
stip = ["ConstructName","dFF0 Short Stim","GRRatioMean","Photostability","Sequence"];
T_plot = T_screening(find(T_screening.("dFF0 Short Stim")<-0.03),:);
ctrlidx = find(strcmp(T_plot.Sequence,refOrigin));
[ax,s] = scatterbyparent(T_plot,selLib,pp1lib,pp2lib,'ax',ax,'ctrlidx',ctrlidx,...
    'tip',stip,'cmap',c(2:end,:));


% Plot evolutionary path (from previous section)
for i = 1:nevol 
    e(i) = errorbar(ax,nxval(i),nyval(i),nystd(i),...
        nystd(i),nxstd(i),nxstd(i),'.','MarkerSize',30,'color',c(i,:));
end
% Plot evolution path
nxevol = nxval(endnodes)-nxval(startnodes);
nyevol = nyval(endnodes)-nyval(startnodes);
quiver(ax,nxval(startnodes)+0.2*nxevol,nyval(startnodes)+0.2*nyevol,...
    nxevol,nyevol,'r','MaxHeadSize',0.1,...
    'AutoScaleFactor',0.98)
% Plot countour line
for i = 1:numel(DIAUC)
    x = min(xval)-0.2:0.1:max(xval+2.5);
    y = (DIAUC(i)./x).^1;
    l(i) = plot(ax,x,y,'--k');
    text(ax,x(end),y(end),strcat(AUCoperation,num2str(DIAUC(i))),...
        'Interpreter','latex');
end

xlabel(ax,{pp1,strcat("Normalized to ",norm2ctrl)})
ylabel(ax,{pp2,strcat("Normalized to ",norm2ctrl)})
axis square               
axis(ax, [0 9 0 2]);
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

function [ax,xval,yval,xstd,ystd,e] = plotevolpath(T_plot,pp1,pp2,p)
    arguments 
        T_plot (:,:) table
        pp1 (1,1) string {mustBeVarofTable(pp1, T_plot)} % x-value
        pp2 (1,1) string {mustBeVarofTable(pp2, T_plot)} % y-value
        p.pp1std (1,1) string {mustBeVarofTable(p.pp1std, T_plot)} = ""
        p.pp2std (1,1) string {mustBeVarofTable(p.pp2std, T_plot)} = ""
        p.ax (1,1) handle = nexttile(tiledlayout(figure(), 'flow', 'Padding', 'none'))
        p.norm2ctrl (1,1) string {mustBeRowofTable(p.norm2ctrl, T_plot)} = ""        
        p.startnodes (1,:) double = []
        p.endnodes (1,:) double = []
        p.notetxt (1,:) string = ""
        p.cmap (:,3) double = lines
    end
    pp1std = p.pp1std;
    pp2std = p.pp2std;
    ax = p.ax;
    norm2ctrl = p.norm2ctrl;
    startnodes = p.startnodes;
    endnodes = p.endnodes;
    note = p.notetxt;
    cmap = p.cmap;
    
    % Scatter plot of variants on evol path
    hold(ax, 'on');
    nvar = height(T_plot);
    xval = zeros(1,nvar);
    xstd = zeros(1,nvar);
    yval = zeros(1,nvar);
    ystd = zeros(1,nvar);
    for i = 1:nvar
        xval(i) = T_plot.(pp1)(i)./T_plot.(pp1)(norm2ctrl);
        yval(i) = T_plot.(pp2)(i)./T_plot.(pp2)(norm2ctrl);
        if ~isempty(pp1std)
            xstd(i) = T_plot.(pp1std)(i)./T_plot.(pp1)(norm2ctrl);
        else
            xstd(i) = 0;
        end
        if ~isempty(pp2std)
            ystd(i) = T_plot.(pp2std)(i)./T_plot.(pp2)(norm2ctrl);
        end
        e(i) = errorbar(ax,xval(i),yval(i),ystd(i),...
            ystd(i),xstd(i),xstd(i),'.','MarkerSize',20,'color',cmap(i,:));
    end
    
    % Plot evolution path
    if ~isempty(startnodes) && ~isempty(endnodes)       
        xevol = xval(endnodes)-xval(startnodes);
        yevol = yval(endnodes)-yval(startnodes);
        quiver(ax,xval(startnodes)+0.2*xevol,yval(startnodes)+0.2*yevol,...
            xevol,yevol,'r','MaxHeadSize',0.1,...
            'AutoScaleFactor',0.98)
    end
    
    % Add notes to evolution arrows
    if ~isempty(note)
        for i = 1:length(startnodes)
        text(ax, mean([xval(startnodes(i)),xval(endnodes(i))]),...
            mean([yval(startnodes(i)),yval(endnodes(i))]),note(i),...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
end

function [ax,C,h] = plotcontour(x,y,p)
    arguments 
        x (1,:) double
        y (1,:) double
        p.ax (1,1) handle = nexttile(tiledlayout(figure(), 'flow', 'Padding', 'none'))
        p.operation (1,1) string = "X.*Y"
        p.LevelList (1,:) double = []
        p.note (1,:) string = ""
    end
    ax = p.ax;
    operation = p.operation;
    note = p.note;
    LevelList = p.LevelList;
    hold(ax, 'on');
    [X,Y] = meshgrid(x,y);
    Z = eval(operation);
    switch note
        case "" % default: text numbers)
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z,'ShowText','on');
            else
                [C,h] = contour(X,Y,Z,'ShowText','on','LevelList',LevelList);
            end
        case "off"
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z);
            else
                [C,h] = contour(X,Y,Z,'LevelList',LevelList);
            end
        otherwise
            if isempty(LevelList)
                [C,h] = contour(X,Y,Z,'ShowText','on');     
                h.DisplayName = note;
            else
                [C,h] = contour(X,Y,Z,'ShowText','on','LevelList',LevelList);
                h.DisplayName = note;
            end
    end        
end

function [ax,s] = scatterbyparent(T_plot,parent,pp1,pp2,p)
    arguments
        T_plot (:,:) table
        parent (1,:) string
        pp1 (1,1) string {mustBeVarofTable(pp1, T_plot)} % x-value
        pp2 (1,1) string {mustBeVarofTable(pp2, T_plot)} % y-value       
        p.ax (1,1) handle = nexttile(tiledlayout(figure(), 'flow', 'Padding', 'none'))
        p.norm2ctrl (1,1) string = ""    
        p.ctrlidx (1,:) double = []   
        p.tip (:,:) string
        p.cmap (:,3) double = lines
    end
    ax = p.ax;
    norm2ctrl = p.norm2ctrl;
    ctrlidx = p.ctrlidx;
    tip = p.tip;
    cmap = p.cmap;
    if ~isempty(ctrlidx)
        normx = nanmean(T_plot.(pp1)(ctrlidx));
        normy = nanmean(T_plot.(pp2)(ctrlidx));
    elseif ~isempty(norm2ctrl)
        ctrlidx = find(strcmp(T_plot.Sequence,norm2ctrl));
        normx = nanmean(T_plot.(pp1)(ctrlidx));
        normy = nanmean(T_plot.(pp2)(ctrlidx));
    else
        normx = 1;
        normy = 1;
    end
    for i = 1:length(parent)
        libidx = strcmp(T_plot.Parent,parent(i));
        xlibval = T_plot.(pp1)(libidx)./normx;
        ylibval = T_plot.(pp2)(libidx)./normy;   
        s(i) = scatter(ax, xlibval, ylibval,'filled','MarkerFaceAlpha',0.5);
        s(i).CData = cmap(i,:);
        for j = 1:length(tip)
            roweg = T_plot.(tip(i))(libidx);
            if iscell(roweg)
            s(i).DataTipTemplate.DataTipRows(j) = dataTipTextRow(tip(j),...
                cellstr(T_plot.(tip(j))(libidx)));
            elseif isnumeric(roweg) | isstring(roweg)
            s(i).DataTipTemplate.DataTipRows(j) = dataTipTextRow(tip(j),...
                T_plot.(tip(j))(libidx));  
            else
                disp("cannot show " + tip(j) + " in datatip");
            end
        end
    end
end
    
function mustBeVarofTable(s,T)
% mustBeVarofTable validate string is a member of the table variables
    arguments
        s (:,1) string
        T (:,:) table
    end
    if ~ismember(s,T.Properties.VariableNames)
        throwAsCaller(MException('validator:mustBeVarofTable', ...
            'Input property name must be variables of the table'));
    end
end

function mustBeRowofTable(s,T)
% mustBeVarofTable validate string is a member of the table variables
    arguments
        s (:,1) string
        T (:,:) table
    end
    if ~ismember(s,T.Properties.RowNames)
        throwAsCaller(MException('validator:mustBeRowofTable', ...
            'Input name must be ROW NAMES of the table'));
    end
end