%% Extract data from Ephys 3-test figureand export to target excel 
clear
clc
close all
path = "D:\OneDrive - rice.edu\Paper\Data\Ephys";
figName = "201901310403_ASAP2s-40x-HEK293A.fig";
tarTable = "Ephys_201901310403_Three-test_One-photon_ASAP2s.xlsx";
GEVIName = "ASAP2s"; 
f1 = openfig(fullfile(path,figName));
%% Get voltage clamp trace
h1 = findobj(subplot(3,3,2), 'Type', 'line');
Var1 = cell(2*numel(h1),1);
Var1(:) = {GEVIName};
Var2 = zeros(2*numel(h1),1);
Var3 = cell(2*numel(h1),1);
Var4 = cell(2*numel(h1),1);
Var5 = NaN(2*numel(h1),1300);
i = 1;
j = 1;
while j <= numel(h1)    
    spName = split(h1(j).DisplayName,'(');
    Var2(i) = str2num(strcat('(',spName{2}));
    Var2(i+1) = str2num(strcat('(',spName{2}));
    Var3{i} = spName{1};
    Var3{i+1} = spName{1};
    Var4{i} = 't';
    Var4{i+1} = 'f';
    tLen = length(h1(j).XData);
    Var5(i,1:tLen) = h1(j).XData;
    Var5(i+1,1:tLen) = h1(j).YData;
    j = j+1;
    i = i+2;
end
VClamp = table(Var1,Var2,Var3,Var4,Var5);
writetable(VClamp,fullfile(path,tarTable),'FileType','spreadsheet','Sheet','VoltageClamp');

%% Get dF/F - V curve
h2 = findobj(subplot(3,3,3), 'Type', 'line');
Var1 = cell(numel(h2)+1,1);
Var1(:) = {GEVIName};
Var2 = cell(numel(h2)+1,1);
Var2{1} = 'Voltage';
Var3 = zeros(numel(h2)+1,4);
Var3(1,:) = [-100, -70, -40, 30];
for i = 1:numel(h2)
    Var2{i+1} = h2(i).DisplayName;
    Var3(i+1,:) = h2(i).YData;
end
VFcurve = table(Var1,Var2,Var3);
writetable(VFcurve,fullfile(path,tarTable),'FileType','spreadsheet','Sheet','TransitionCurve');

%% Get 2Hz APs
h3 = findobj(subplot(3,3,5), 'Type', 'line');
Var1 = cell(2*numel(h3),1);
Var1(:) = {GEVIName};
Var2 = cell(2*numel(h3),1);
Var3 = cell(2*numel(h3),1);
Var4 = NaN(2*numel(h3),150);
i = 1;
j = 1;
while j <= numel(h3)    
    Var2{i} = h3(j).DisplayName;
    Var2{i+1} = h3(j).DisplayName;
    Var3{i} = 't';
    Var3{i+1} = 'f';
    tLen = length(h3(j).XData);
    Var4(i,1:tLen) = h3(j).XData;
    Var4(i+1,1:tLen) = h3(j).YData;
    j = j+1;
    i = i+2;
end
shortAPs = table(Var1,Var2,Var3,Var4);
writetable(shortAPs,fullfile(path,tarTable),'FileType','spreadsheet','Sheet','2HzAPs');

%% Get 100Hz APs
h4 = findobj(subplot(3,3,6), 'Type', 'line');
Var4 = NaN(2*numel(h4),600);
i = 1;
j = 1;
while j <= numel(h3)    
    tLen = length(h4(j).XData);
    Var4(i,1:tLen) = h4(j).XData;
    Var4(i+1,1:tLen) = h4(j).YData;
    j = j+1;
    i = i+2;
end
longAPs = table(Var1,Var2,Var3,Var4);
writetable(longAPs,fullfile(path,tarTable),'FileType','spreadsheet','Sheet','100HzAPs');
