% [k, R2, pVal] = libplot.lnrFitting(x,y,intercept,figureTrue,figNumber)
% Generate a linear fitting based on the give x and y vector and intercept,
% and plot on the specified figure (optional). The fitting method is FITLM function in MATLAB
% 
% INPUT(all required):
%   x: vector for x axis data for fitting
%   y: vector for y axis data for fitting
%   intercept: indicate whether the fitting come from origin, false=come from 0,true=not from 0
%   figureTrue: whether we need a figure, true/false
%   figNumber: need to assign figure a number
% OUTPUT:
%   k,R2,pVal: y = kx, R-square and pValue.
%   
%   Example usage:
%   [k, R2, pVal]= libplot.lnrFitting(x,y,false,true,6);
%   perform a liner fitting on x and y from 0 origin with a figure at figure(6)
function [k, R2, pVal] = lnrFitting(x,y,intercept,figureTrue,axh,fitfun)
    switch fitfun
        case 'polyfit'
            k = polyfit(x,y,1);
            SSError = sum(((k(1)*x+k(2))-y).^2);
            SST = sum((y-mean(y)).^2);
            R2 = 1-(SSError/(numel(x)-2))/(SST/(numel(x)-1)); % adjusted R2
            pVal = nan;
            RMSE = sqrt(SSError/(numel(x)));
        case 'fitlm'
            p = fitlm(x,y,'Intercept',intercept);   
            k = p.Coefficients.('Estimate'); % k(1) = intercept, k(2) = slope
            if numel(k)==1
                k = [k(1),0]; % k(1) = slope, k(2) = intercept = 0
            else
                k = [k(2),k(1)]; % k(1) = slope, k(2) = intercept
            end
            R2 = p.Rsquared.Adjusted;
            RMSE = p.RMSE;
            pVal = p.Coefficients.('pValue');
    end
    if figureTrue
        if min(x) <0
            linestart = min(x)-0.1*min(x);
        else
            linestart = 0;
        end
        lineend = max(x)+0.1*max(x);
        hold(axh, 'on');
        xplot = linestart:0.01:lineend;        
        plot(axh, xplot,xplot*k(1)+k(2) ,'--k')
        text(0.05,0.9,{strcat('k = ',num2str(k(1))),...
            strcat('R2 = ',num2str(R2)),strcat('RMSE = ',num2str(RMSE)),...
            strcat('Intercept = ',num2str(k(2)))},'Color',[0 0 0],'units','normalized');
    end
end