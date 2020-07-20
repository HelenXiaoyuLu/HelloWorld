% Create the distribution
tpnum = 50;
fpnum = 47000;
tpcenter = 4; %10^4
fpcenter = 0; %10^0
tpdistribution = 10.^normrnd(tpcenter,0.75,[1,tpnum]);
fpdistribution = 10.^normrnd(fpcenter,0.75,[1,fpnum]);
figure()
hold on
plot(tpdistribution, normrnd(1,1,length(tpdistribution)),'g+')
plot(fpdistribution, ones(1,length(fpdistribution)),'b+')
set(gca,'Xscale','log') % yeah, that is indeed easy!
set(gcf, 'Position',  [50, 10, 500, 500])

%%
ttdistribution  = sort([tpdistribution,fpdistribution],2,'ascend');
precision = [];
recall = [];
for gate = 0.1:0.01:3.5
    selcells = sum(ttdistribution>10^gate);
    tp = sum(tpdistribution>10^gate);
    fp = sum(fpdistribution>10^gate);
    precision = [precision,tp/selcells];
    recall = [recall,tp/pAnum];
end
figure()
plot(recall, precision)
xlabel('recall')
ylabel('precision')
set(gcf, 'Position',  [50, 10, 500, 500])

%%
RecallM = [0.522959184	0.558418367	0.557142857	0.575510204	0.580102041	0.582397959	0.564540816	0.493112245	0.461989796	0.408418367 nan nan nan nan nan nan nan;
    nan 0.522101326 0.510880653	0.513090785	0.511334013	0.515527598	0.519607843	0.519607843	0.519607843	0.519607843	0.509803922	0.509803922	0.5	0.480392157	0.480392157	0.480392157	0.460784314;
    nan	0.561729412	0.551793419	0.516891187	0.531734791	0.502980798	0.486523281	0.480375043	0.449447207	0.422680412	0.412371134	0.402061856	0.381443299	0.360824742	0.340206186	0.340206186	0.329896907];
RecallY = [0.596577381	0.47577381	0.405535714	0.343392857	0.291666667	0.229166667	0.177083333 nan;
    0.409579792	0.3513875	0.325369861	0.277546181	0.222036944 nan nan nan;
0.298364878	0.227477073	0.19456	0.1608	0.136515122	0.107317073	0.07804878	0.058536585];
PrecisionM = [0.077622113	0.22453585	0.493670886	0.573170732	0.89036805	0.91576414	0.913330582	0.902006533	0.962785752	0.958108917 nan nan nan nan nan nan nan;
    nan	0.446561	0.559659796	0.635636057	0.732981316	0.882518432	1	1	1	1	1	1	1	1	1	1	1; ...
    nan	0.608885028	0.690418298	0.725189076	0.89579403	0.942072795	0.959343609	0.978989998	0.977576655	1	1	1	1	1	1	1	1];
PrecisionY = [0.065808697	0.317901603	0.628621517	0.785539216	0.933333333	0.956521739	0.944444444 nan;
    0.101692372	0.230418243	0.510088156	0.799866511	0.864767351 nan nan nan;
0.120879051	0.309579321	0.6053718	0.804706572	0.874943725	0.956521739	0.941176471	0.923076923];

% Assumptions: 
% In non-pA population,  is popGRate
% In pA population, the GFP ratio is pAGRate
trialName = 'Mammalian Trial 3'
pAnum = 102; % Total photoactivated cell number
pAsuccessrate = 0.55; % Successful rate of PA, also includes retrival rate
Ntotal = 47276; % Total number of cells 
npAnum = Ntotal-pAnum; % Number of non-photoactivated cells
g = 2.4:0.01:3.8; % Threshold of gates
sigma0 = 0.78; % standard deviation of mCherry signal in non-photoactivated population
sigma1 = 0.78; % standard deviation of mCherry signal in photoactivated population
popGRate = 0.02; % the GFP ratio in non-photoactivated population
pAGRate = 1; % the GFP ratio in photoactivated population
mu0 = 0; % mean mCherry signal of non-photoactivated population (log-scale)
mu1 = 4; % mean mCherry signal of photoactivated population (log-scale); 

c = colormap(lines);
c2 = colormap(flipud(jet));

% True positive = GREEN+/RED+
truep = (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*pAGRate*pAsuccessrate;

% False positive = GREEN-/RED+ 
% This is either from bleed-through from non-photoactivated cells,
% or photoactivated GREEN- cells by mistake
falsep = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*sigma0))))*npAnum*(1-popGRate)+...
    (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*(1-pAGRate);

precision = truep./(truep+falsep);
recall = truep./pAnum;
figure(4)
clf
sgtitle(trialName)
subplot(1,2,1)
hold on
pAval = normrnd(mu1,sigma1,[1,pAnum]);
yval = normrnd(1,0.5,[1,pAnum]);
npAval = normrnd(mu0,sigma0,[1,npAnum]);
npAyval = normrnd(1,0.5,[1,npAnum]);
plot(10.^npAval,npAyval,'.','color',[200 200 200]/256);
plot(10.^pAval,yval,'.','color',c(2,:));
set(gca,'Xscale','log')
xlabel('mCherry signal')
legend('non-photoactivated','photoactivated','location','northeast')
legend('boxon')

subplot(1,2,2)
hold on
plot(recall, precision,'color',c(2,:))
xlabel('recall')
ylabel('precision')

plot(RecallM(3,:),PrecisionM(3,:),'+','color',c(1,:))
% plot(RecallM(2,:),precisionM(2,:),'+','color',c(2,:))
% plot(RecallM(3,:),precisionM(3,:),'+','color',c(3,:))
annotation('textbox', [0.8, 0.8, 0.1, 0.1], 'String', {strcat('sigma0: ',num2str(sigma0));...
    strcat('sigma1: ',num2str(sigma1));strcat('mu1: ',num2str(mu1));...
    strcat('pAGRate: ',num2str(pAGRate));strcat('pAsuccessrate',num2str(pAsuccessrate))},'LineStyle','none')
set(gcf, 'Position',  [50, 10, 1000, 400])

% subplot(1,3,3)
% precision = [];
% recall = [];
% for gate = g
%     tp = sum(pAval>gate)*pAGRate;
%     fp = sum(npAval>gate)*(1-popGRate);
%     precision = [precision,tp/(tp+fp)];
%     recall = [recall,tp/pAnum];
% end
% hold on
% plot(recall, precision)
% plot(RecallM(1,:),PrecisionM(1,:),'+','color',c(1,:))
% xlabel('recall')
% ylabel('precision')
% set(gcf, 'Position',  [110, 100, 1000, 450])

figure(3)
clf
hold on
plot(g,recall)
plot(g,precision)
legend('recall','precision')