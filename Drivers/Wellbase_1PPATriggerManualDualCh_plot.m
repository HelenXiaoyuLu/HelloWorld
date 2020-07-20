% Wellbase_1PPATriggerManualdualCh
groupgroups = findgroups(T_fov.Name);
T_group = table();
T_group.Name = splitapply(@unique, T_fov.Name, groupgroups);
n_group = max(groupgroups);

%% Scatch raw traces
figure(1)
clf
sgtitle('Raw Traces Red (up) vs Green (down)')
cmap = colormap(lines);
hold on
for i = 1:height(T_fov)
    subplot(2,n_group,groupgroups(i))
  %  hold on
    plot(T_fov.TraceRawRed(i).x', T_fov.TraceRawRed(i).y', 'color', cmap(2,:));
    xlabel('t [s]')
    ylabel({'Background corrected','fluorescence [a.u]'})
    title(T_fov.Name(i))
    xlim([0.5 15])
    axis square
    
    subplot(2,n_group,groupgroups(i)+n_group)
  %  hold on
    plot(T_fov.TraceRawGreen(i).x', T_fov.TraceRawGreen(i).y', 'color', cmap(5,:));
    xlabel('t [s]')
    ylabel({'Background corrected','fluorescence [a.u]'})
    title(T_fov.Name(i))  
    xlim([0.5 15])
    axis square
end

%% Scatch Final traces
figure(2)
clf
sgtitle('Raw Traces Red (up) vs Green (down)')
cmap = colormap(lines);
hold on
for i = 1:height(T_fov)
    subplot(2,n_group,groupgroups(i))
  % hold on
    plot(T_fov.TraceFinalRed(i).x', T_fov.TraceFinalRed(i).y', 'color', cmap(2,:));
    xlabel('t [s]')
    ylabel({'Background corrected','fluorescence [a.u]'})
    title(T_fov.Name(i))
    xlim([0.5 15])
    ylim([0.8 1.6])
    axis square
    
    subplot(2,n_group,groupgroups(i)+n_group)
  % hold on
    plot(T_fov.TraceFinalGreen(i).x', T_fov.TraceFinalGreen(i).y', 'color', cmap(5,:));
    xlabel('t [s]')
    ylabel({'Background corrected','fluorescence [a.u]'})
    title(T_fov.Name(i))  
    xlim([0.5 15])
    ylim([0.5 2.8])
    axis square
end



%% Plot stimulation
load('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\20200708 photoactivation train\stim.mat')
cmap = colormap(lines);
figure()
hold on
pw = 10;
plot([0; stim.x; 15],[0; stim.y/20*pw; 0],'color',cmap(5,:))
plot([0 15],[100 100],'color',cmap(2,:))
ylim([0 105])
xlabel('t [s]')
ylabel({'Excitaiton Power (%)'})
legend('470nm','555nm')

