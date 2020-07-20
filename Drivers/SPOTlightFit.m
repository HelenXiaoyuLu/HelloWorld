%% Data from FACS
%SPfit = table();
SPfit.Name = {'Mammalian 1';'Mammalian 2';'mammalian 3';'Yeast 1';'Yeast 2';'Yeast 3'};
SPfit.Properties.RowNames = {'Mammalian 1';'Mammalian 2';'Mammalian 3';'Yeast 1';'Yeast 2';'Yeast 3'};
SPfit.Recall = {[nan	0.558418367	0.557142857	0.575510204	0.580102041	0.582397959	0.564540816	0.493112245	0.461989796	0.408418367 nan nan nan nan nan nan nan];   
    [nan	0.561729412	0.551793419	0.516891187	0.531734791	0.502980798	0.486523281	0.480375043	0.449447207	0.422680412	0.412371134	0.402061856	0.381443299	0.360824742	0.340206186	0.340206186	0.329896907];
[nan 0.522101326 0.510880653	0.513090785	0.511334013	0.515527598	0.519607843	0.519607843	0.519607843	0.519607843	0.509803922	0.509803922	0.5	0.480392157	0.480392157	0.480392157	0.460784314];
[0.596577381	0.47577381	0.405535714	0.343392857	0.291666667	0.229166667	0.177083333 nan nan nan nan nan nan nan nan nan nan];
    [0.496997361	0.422350069 0.409579792	0.3513875	0.325369861	0.277546181	0.222036944 nan nan nan nan nan nan nan nan nan nan];
[0.298364878	0.227477073	0.19456	0.1608	0.136515122	0.107317073	0.07804878	0.058536585 nan nan nan nan nan nan nan nan nan]};

SPfit.Precision = {[nan	0.22453585	0.493670886	0.573170732	0.89036805	0.91576414	0.913330582	0.902006533	0.962785752	0.958108917 nan nan nan nan nan nan nan];
    [nan	0.608885028	0.690418298	0.725189076	0.89579403	0.942072795	0.959343609	0.978989998	0.977576655	1	1	1	1	1	1	1	1];
    [nan	0.446561	0.559659796	0.635636057	0.732981316	0.882518432	1	1	1	1	1	1	1	1	1	1	1];
[0.065808697	0.317901603	0.628621517	0.785539216	0.933333333	0.956521739	0.944444444 nan nan nan nan nan nan nan nan nan nan];
   [0.015940871	0.03801582 0.101692372	0.230418243	0.510088156	0.799866511	0.864767351 nan nan nan nan nan nan nan nan nan nan];
[0.120879051	0.309579321	0.6053718	0.804706572	0.874943725	0.956521739	0.941176471	0.923076923 nan nan nan nan nan nan nan nan nan]};
SPfit.('Total cell number') = {109298;86904;47276;621066;724460;621987};
SPfit.('PA cell number') = {56;97;102;96;144;205};
SPfit.('Population GFP ratio') = {0.016463415;0.039401185;0.02;0.005432596;0.004;0.0052};
SPfit.('mu0') = {0; 0; 0; 0; 0; 0};
SPfit.('mu1') = {3; 4; 4; 3; 3; 3};

%% Input
trialName = 'Mammalian 1'; % Mammalian or Yeast, trial 1-3
pAnum = SPfit.('PA cell number'){trialName}; % Total photoactivated cell number
Ntotal = SPfit.('Total cell number'){trialName}; % Total number of cells 
npAnum = Ntotal-pAnum; % Number of non-photoactivated cells
rec = SPfit.Recall{trialName}; 
prec = SPfit.Precision{trialName};
popGRate = SPfit.('Population GFP ratio'){trialName}; % the GFP ratio in non-photoactivated population
mu0 = SPfit.('mu0'){trialName}; % mean mCherry signal of non-photoactivated population (log-scale)
mu1 = SPfit.('mu1'){trialName}; % mean mCherry signal of photoactivated population (log-scale); 
pAsuccessrate = max(rec)+0.2; % Successful rate of PA, also includes retrival rate
pAGRate = max(prec); % the GFP ratio in photoactivated population

% Fitting 
currentR = inf;
currentSigma = [0 0];
g = 0:0.01:7; % Threshold of gates
for sigma0 = 0.3:0.01:1 % standard deviation of mCherry signal in non-photoactivated population
    for sigma1 = 0.3:0.01:2 % standard deviation of mCherry signal in photoactivated population
        % True positive = GREEN+/RED+
        truep = (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*pAGRate*pAsuccessrate;
        % False positive = GREEN-/RED+ 
        % This is either from bleed-through from non-photoactivated cells,
        % or photoactivated GREEN- cells by mistake
        falsep = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*sigma0))))*npAnum*(1-popGRate)+...
            (1-0.5*(1+erf((g-mu1)/(sqrt(2)*sigma1))))*pAnum*(1-pAGRate)*pAsuccessrate;
        precision = truep./(truep+falsep);
        recall = truep./pAnum;
        [sortedrec,idx] = sort(rec);
        sortedprec = prec(idx);
        [srecall1,id1] = sort(recall);
        [srecall,id] = unique(srecall1);
        sprecision = precision(id1);
        sprecision = sprecision(id);
        sortedrec = rmmissing(sortedrec);
        sortedprec = rmmissing(sortedprec);
        precinterp = interp1(srecall,sprecision,sortedrec);
        R2 = nansum((precinterp - sortedprec).^2);
        if R2 < currentR
            currentR = R2;
            currentSigma = [sigma0,sigma1];
        end
    end
end
SPfit.('Sigma0'){trialName} = currentSigma(1);
SPfit.('Sigma1'){trialName} = currentSigma(2);

%% Plot
c = colormap(lines);
c2 = colormap(flipud(jet));
g = 0.6:0.01:3.8; % Threshold of gates
truep = (1-0.5*(1+erf((g-mu1)/(sqrt(2)*currentSigma(2)))))*pAnum*pAGRate*pAsuccessrate;
falsep = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*currentSigma(1)))))*npAnum*(1-popGRate)+...
            (1-0.5*(1+erf((g-mu1)/(sqrt(2)*currentSigma(2)))))*pAnum*(1-pAGRate)*pAsuccessrate;
precision = truep./(truep+falsep);
recall = truep./pAnum;
figure(1)
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
%plot(recall, precision2,'color',c(2,:))
xlabel('recall')
ylabel('precision')
plot(rec,prec,'+','color',c(1,:));
SPfit.('Gate Threshold'){trialName} = g;

annotation('textbox', [0.8, 0.8, 0.1, 0.1], 'String', {strcat('sigma0: ',num2str(sigma0));...
    strcat('sigma1: ',num2str(sigma1));strcat('mu0: ',num2str(mu0));strcat('mu1: ',num2str(mu1));...
    strcat('pAGRate: ',num2str(pAGRate));strcat('pAsuccessrate: ',num2str(pAsuccessrate))},'LineStyle','none')
set(gcf, 'Position',  [50, 10, 1000, 400])

% save('D:\OneDrive - rice.edu\Francois\Paper\SPOTlight\CurveFitting\spfit2.mat','SPfit')

%%
figure(2)
hold on
for i = 1:6
    rec = SPfit.Recall{i};
    prec = SPfit.Precision{i};
    plot(rec,prec,'+','color',c(i,:));
    mu0 = SPfit.mu0{i};
    mu1 = SPfit.mu1{i};
    s0 = SPfit.Sigma0{i};
    s1 = SPfit.Sigma1{i};
    g = SPfit.('Gate Threshold'){i};
    pAnum = SPfit.('PA cell number'){i};
    npAnum = SPfit.('Total cell number'){i} - pAnum;
    popGRate = SPfit.('Population GFP ratio'){i}
    pAGRate = max(prec); % the GFP ratio in photoactivated population 
    pAsuccessrate = max(rec);
    truep = (1-0.5*(1+erf((g-mu1)/(sqrt(2)*s1))))*pAnum*pAGRate*pAsuccessrate;
    falsep = (1-0.5*(1+erf((g-mu0)/(sqrt(2)*s0))))*npAnum*(1-popGRate)+...
            (1-0.5*(1+erf((g-mu1)/(sqrt(2)*s1))))*pAnum*(1-pAGRate)*pAsuccessrate;
    precision = truep./(truep+falsep);
    recall = truep./pAnum;
    plot(recall,precision,'color',c(i,:))
end
xlabel('recall')
ylabel('precision')
legend(SPfit.Properties.RowNames)

