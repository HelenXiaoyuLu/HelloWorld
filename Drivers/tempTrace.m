% Access data from old fstimanalysis
close all

cmap = colormap(lines);
wellsel = lib.getChildren('g',63).getChildren('w',1);
traceT = table();
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    traceT.Name{i} = strcat(wellsel.Name,'fov',num2str(i));
    for j = 1:fov.Children.n        
        time = wellsel.getChildren('event',j).trel;
        traceT.("Stim_" + j){i}  = [time'; mean(wellsel.getChildren('event',j).normData,1)];
    end
end

figure()
hold on
croptime = 0.1;
cropframe = croptime/0.001;
StimShort = zeros(16,cropframe*2+1);
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    for j = 1:1       
        tr = traceT.("Stim_" + j){i};
        idx0 = find(tr(1,:) == 0);
        StimShort((i-1)*4+j,:) = tr(2,idx0-cropframe:idx0+cropframe);
    end
end
plot(tr(1,idx0-cropframe:idx0+cropframe),mean(StimShort,1),'Color',cmap(1,:))
% plot(tr(1,idx0-cropframe:idx0+cropframe),1-mean(StimShort,1),'Color',cmap(1,:))
% axis([-croptime,croptime,0.4, 1.05])
% xlabel('Time [s]')      
% ylabel('Normalized dF/F0')
figure()
croptime = [-0.1,0.25];
cropframe = croptime/0.001;
StimLong = zeros(length(wellsel.Children),cropframe(2)-cropframe(1)+1);
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    tr = traceT.("Stim_" + 2){i};
    idx0 = find(tr(1,:) == 0);
    StimLong(i,:) = tr(2,idx0+cropframe(1):idx0+cropframe(2));
end
% deltaT = tend - tr(1,idx0+cropframe(1))+0.001;
% plot(deltaT+tr(1,idx0+cropframe(1):idx0+cropframe(2)),1-mean(StimLong,1),'Color',cmap(1,:))
plot(tr(1,length(StimLong)),mean(StimLong,1),'Color',cmap(1,:))
axis([-tend,deltaT + croptime(2),-inf, 0.7])
xlabel('Time [s]')      
ylabel('Normalized dF/F0')
%%
close all

cmap = colormap(lines);
wellsel = lib.getChildren('g',63).getChildren('w',1);

traceT = table();
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    traceT.Name{i} = strcat(wellsel.Name,'fov',num2str(i));
    for j = 1:fov.Children.n        
        time = wellsel.getChildren('event',j).trel;
        traceT.("Stim_" + j){i}  = [time'; mean(wellsel.getChildren('event',j).normData,1)];
    end
end
Stim1 = [];
figure()
hold on
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    tr = traceT.("Stim_" + 1){i};
    idx0 = find(tr(1,:) == 0);
    Stim1(i,:) = tr(2,:);
end
plot(tr(1,:),mean(Stim1,1),'Color',cmap(1,:))
title('Stim_1')

Stim2 = [];
figure()
hold on
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    tr = traceT.("Stim_" + 2){i};
    idx0 = find(tr(1,:) == 0);
    Stim2(i,:) = tr(2,:);
end
plot(tr(1,:),mean(Stim2,1),'Color',cmap(1,:))
title('Stim_2')
Stim3 = [];
figure()
hold on
for i = 1:length(wellsel.Children)
    fov = wellsel.Children(i);
    tr = traceT.("Stim_" + 3){i};
    idx0 = find(tr(1,:) == 0);
    Stim3(i,:) = tr(2,:);
end
plot(tr(1,:),mean(Stim3,1),'Color',cmap(1,:))
title('Stim_3')

%%
traces = readmatrix('D:\OneDrive - rice.edu\Francois\Grant\20200428 U01 Milestones\Fig\greenD2B\Trace_bettervariants.xlsx');
figure()
hold on
plot(traces(:,8),traces(:,9)-1)
plot(traces(:,15),traces(:,16)-1)
plot(traces(:,22),traces(:,23)-1)
plot(traces(:,36),traces(:,37)-1)
plot(traces(:,43),traces(:,44)-1)

axis([-0.17, 0.17, -0.05, 0.25])
xlabel('Time/s')
ylabel('Response (% \DeltaF/F0)')
set(gcf, 'Position',  [110, 100, 400, 400])
box on
legend('N150V','N150D','V151X','A148E','JEDI-alpha')
axis square

figure()
hold on
plot(traces(:,12),traces(:,13)-1)
plot(traces(:,19),traces(:,20)-1)
plot(traces(:,26),traces(:,27)-1)
plot(traces(:,40),traces(:,41)-1)
plot(traces(:,47),traces(:,48)-1)
axis([-0.2, 0.3, -0.05, 0.7])
axis square
xlabel('Time/s')
ylabel('Response (% \DeltaF/F0)')
set(gcf, 'Position',  [110, 100, 450, 450])
box on
legend('N150V','N150D','V151X','A148E','JEDI-alpha','location','best')

