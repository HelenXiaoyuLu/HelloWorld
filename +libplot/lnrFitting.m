% A quick way to generate a linear fitting based on the give x and y value 
% and intercept, and plot on the specified figure (optional). The fitting
% method is FITLM function in MATLAB
% Example usage:
%   [k, R2] = lnrFitting(x,y,false,true,6);
% INPUT:
%   x/y: values to be fitted
%   intercept: specify a value for intecept, default = false
%   figureTrue/figNumber: if plot is needed/which figure to plot on
% OUTPUT:
%   k,R2,pVal: y = kx, R-square and pValue.

function [k, R2, pVal] = lnrFitting(x,y,intercept,figureTrue,figNumber)
    p = fitlm(x,y,'Intercept',intercept);   
    k = p.Coefficients.('Estimate');
    R2 = p.Rsquared.Adjusted;
    pVal = p.Coefficients.('pValue');
    if figureTrue
        figure(figNumber)
        if min(x) <0
            linestart = min(min(x),min(y));
        else
            linestart = 0;
        end
        xplot = linestart:0.01:max(max(x),max(y));        
        plot(xplot,xplot*p.Coefficients(1,1).Variables,'--k')
        text(linestart+0.01,max(y),{strcat('k = ',num2str(k)),strcat('R2 = ',num2str(R2)),...
            strcat('p-value = ',num2str(pVal))},'Color',[0 0 0])
    end
end