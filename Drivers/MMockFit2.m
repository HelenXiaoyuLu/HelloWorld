%%
% Assumptions: 
% In non-pA population, the GFP ratio is popGRate
% In pA population, the GFP ratio is pAGRate
trialName = 'Yeast Trial 3'
pAnum = 205;
Ntotal = 621987;
npAnum = Ntotal-pAnum;
g = 2:0.01:7;
sigma0 = 0.83;
sigma1 = 1.5;
popGRate = 1 - 1/(1+0.008571429);
pAGRate = 0.96;
mu0 = 0;
mu1 = 3;

c = colormap(lines);
c2 = colormap(flipud(jet));
truep = (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*pAGRate;
   % (1-0.5*(1+erf((g-mu0)/(sqrt(2)*sigma0))))*npAnum*popGRate;
falsep = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*sigma0))))*npAnum*(1-popGRate)+...
    (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*(1-pAGRate);
falsep2 = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*sigma0))))*npAnum*(1-popGRate)+...
    (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum - truep;
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
%set(gcf, 'Position',  [50, 10, 500, 500])
RecallM = [0.362631758 0.522101326 0.510880653	0.513090785	0.511334013	0.515527598	0.519607843	0.519607843	0.519607843	0.519607843	0.509803922	0.509803922	0.5	0.480392157	0.480392157	0.480392157	0.460784314;
    0.493659629	0.561729412	0.551793419	0.516891187	0.531734791	0.502980798	0.486523281	0.480375043	0.449447207	0.422680412	0.412371134	0.402061856	0.381443299	0.360824742	0.340206186	0.340206186	0.329896907;
    0.522959184	0.558418367	0.557142857	0.575510204	0.580102041	0.582397959	0.564540816	0.493112245	0.461989796	0.408418367 nan nan nan nan nan nan nan];
RecallY = [0.596577381	0.47577381	0.405535714	0.343392857	0.291666667	0.229166667	0.177083333 nan;
    0.409579792	0.3513875	0.325369861	0.277546181	0.222036944 nan nan nan;
0.298364878	0.227477073	0.19456	0.1608	0.136515122	0.107317073	0.07804878	0.058536585];
precisionM = [0.167377259	0.446561	0.559659796	0.635636057	0.732981316	0.882518432	1	1	1	1	1	1	1	1	1	1	1; ...
    0.31527135	0.608885028	0.690418298	0.725189076	0.89579403	0.942072795	0.959343609	0.978989998	0.977576655	1	1	1	1	1	1	1	1;...
    0.077622113	0.22453585	0.493670886	0.573170732	0.89036805	0.91576414	0.913330582	0.902006533	0.962785752	0.958108917 nan nan nan nan nan nan nan];
PrecisionY = [0.065808697	0.317901603	0.628621517	0.785539216	0.933333333	0.956521739	0.944444444 nan;
    0.101692372	0.230418243	0.510088156	0.799866511	0.864767351 nan nan nan;
0.120879051	0.309579321	0.6053718	0.804706572	0.874943725	0.956521739	0.941176471	0.923076923];

plot(RecallY(3,:),PrecisionY(3,:),'+','color',c(1,:))
% plot(RecallM(2,:),precisionM(2,:),'+','color',c(2,:))
% plot(RecallM(3,:),precisionM(3,:),'+','color',c(3,:))
set(gcf, 'Position',  [110, 100, 1000, 450])
annotation('textbox', [0.8, 0.8, 0.1, 0.1], 'String', {strcat('sigma0: ',num2str(sigma0));...
    strcat('sigma1: ',num2str(sigma1));strcat('mu1: ',num2str(mu1));...
    strcat('pAGRate: ',num2str(pAGRate))},'LineStyle','none')