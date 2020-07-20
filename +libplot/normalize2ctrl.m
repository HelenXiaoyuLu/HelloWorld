function libTable = normalize2ctrl(libTable,pp1,ctrl,varargin)
    p = inputParser;
    p.addRequired('libTable', @istable);
    p.addRequired('pp1', @(x) any(strcmp(x, libTable.Properties.VariableNames)));
    p.addRequired('ctrl', @ischar);
    p.addParameter('RefCol', 'plibName', @ischar);
    p.parse(libTable,pp1,ctrl,varargin{:})
    RefCol = p.Results.RefCol;    
    pp1col = find(strcmp(pp1,libTable.Properties.VariableNames));
    normctrlsel = strmatch(ctrl, libTable.(RefCol),'exact') ;
    Tnormctrl = libTable(normctrlsel,:);
    meanctrlppl = nanmean(Tnormctrl{:,pp1col});
    newVarname = strcat(pp1,'Norm2',ctrl);
    libTable.(newVarname)(:) = libTable{:,pp1col}/meanctrlppl;   
end