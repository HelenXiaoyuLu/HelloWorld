% Load
load('D:\OneDrive - rice.edu\Francois\ASAPScreening\Wet\Data\20200710 JEDI photoactivation\20200710 TrainStim photoactivation result.mat');

%%
TrialNames = fieldnames(Tgroup);
cmap = colormap(lines);
figure(1)
clf
hold on
for i = 1:length(TrialNames)-1
    T_out = Tfov.(TrialNames{i});
    for j = 1:height(T_out)
        tr = T_out.TraceRaw(j);
        ynorm = normalize(tr.y,'range');
        subplot(2,6,j)
        title(T_out.Name(j),'Interpreter','none')
        hold on
        plot(tr.x,ynorm,'color',cmap(i,:))
        xlabel('t [s]')
        ylabel('Normalized Fluorescence')
        axis square
    end
end
legend(TrialNames(1:end-1),'location','best')

%%
figure(2)
clf
hold on
for i = 1:length(TrialNames)-1
    T_out = Tfov.(TrialNames{i});
    for j = 1:height(T_out)
        tr = T_out.TraceRaw(j);
        ynorm = normalize(tr.y,'range');
        subplot(2,6,j)
        title(T_out.Name(j),'Interpreter','none')
        hold on
        plot(tr.x,tr.y,'color',cmap(i,:))
        xlabel('t [s]')
        ylabel(' Fluorescence')
        axis square
    end
end
legend(TrialNames(1:end-1))
%%
figure(3)
clf
hold on
for i = 1:length(TrialNames)-1
    T_out = Tfov.(TrialNames{i});
    for j = 1:height(T_out)
        tr = T_out.TraceRaw(j);
        ynorm = normalize(tr.y,'range');
        subplot(2,6,j)
        title(T_out.Name(j),'Interpreter','none')
        hold on
        plot(tr.x,ynorm,'color',cmap(i,:))
        xlabel('t [s]')
        ylabel('Normalized Fluorescence')
        ymax = min(ynorm(end)+0.2,1);
        ymin = max(ynorm(end)-0.2,0);
        ylim([ymin ymax]);
        axis square
    end
end
legend(TrialNames(1:end-1),'location','best')

%%
cmap2 = colormap(jet(6));
figure(4)
clf
hold on
for i = 1:length(TrialNames)-1
    T_out = Tfov.(TrialNames{i});
    for j = 1:height(T_out)
        tr = T_out.TraceFinal(j);
        ynorm = normalize(tr.y,'range');
        subplot(2,6,j)
        title(T_out.Name(j),'Interpreter','none')
        hold on
        plot(tr.x,tr.y,'color',cmap2(i,:))
        xlabel('t [s]')
        ylabel('Normalized Fluorescence')
        ymax = min(ynorm(end)+0.2,1);
        ymin = max(ynorm(end)-0.2,0);
        axis square
    end
end
legend(TrialNames(1:end-1),'location','best')