figure(1)
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

