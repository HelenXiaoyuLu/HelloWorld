% A Helper for barplotall to plot horizontal lines and generate a text box

function pltHorizon(val,no,legend,varargin)
    p = inputParser;
    p.addRequired('val', @isnumeric);
    p.addRequired('no', @isnumeric);
    p.addRequired('legend', @iscell);
    p.addParameter('arrow',false, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));
    p.addParameter('LineSpec', '--k') % assign the control for scatter
    p.parse(val,no,legend,varargin{:});
    arrow = p.Results.arrow;  
    LineSpec = p.Results.LineSpec; 
   
    % Plot horizontal threshold line
    X = 1:no;
    plot(X,ones(1,no)*val,LineSpec)
    text(X(end)+1,val,legend','HorizontalAlignment','left','VerticalAlignment', 'baseline')        
end