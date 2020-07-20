% parse library screening result and generate a standardized table for
% plotting function. 
% USAGE:
%   Screening = libplot.tableparse(Screening,'parse_on','Legal Name', ...
%   'mask_on','Area_OnePhoton','mskth', 50);
% INPUT:
%   REQUIRED: Table (tab)
%   PARAMETERS (Name-value pairs): 
%       mask_on: do masking on which property, default is 'masked'
%       mskth: (lower) threshold for mask_on
%       parseName: boolean value, if parsing names is needed or not.
%       parse_on: specify the column of names to parse on, this should be
%       in the formats (e.g.) "ASAP1-R415 (P0A7)" or "[JEDI-1P] (P0A7)"
%       (more works TBD for univeral name parsing
%       squareName: if the name of GEVI is round by []
% OUTPUT:
%    tab: a table with standardized name/platewell information, if
%    parseName = true;
%    emptywells: a list of the index of empty wells

function [tab, emptywells] = tableparse(tab,varargin)
    p = inputParser;
    p.addRequired('tab', @istable);   
    p.addParameter('mask_on', 'masked', @ischar);
    p.addParameter('mskth', 0, @isnumeric);
    p.addParameter('parse_on', 'Name', @(x) any(strcmp(x, tab.Properties.VariableNames)));
    p.addParameter('parseName', true, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));  
    p.addParameter('squareName', false, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));    
    p.parse(tab,varargin{:});
    mask_on = p.Results.mask_on;
    mskth = p.Results.mskth;
    parse_on = p.Results.parse_on;
    parseName = p.Results.parseName;
    squareName = p.Results.squareName;
    emptywells = [];
    switch squareName
        case 1
        switch parseName
            case 1
                for w = 1:height(tab)
                    nameparse = split(tab.(parse_on){w},["[",']']);
                    np2 = split(nameparse{end},["-",'(',')'])
                    tab.pPlateWell{w} = np2{end-1};
                    tab.plibName{w} = strcat(np2{1},'-', np2{2});
                    tab.pMutation{w} = np2{2};
                    if length(nameparse) > 1
                        tab.pParent{w} = nameparse{2};
                    else
                        tab.pParent{w} = nameparse{2};
                    end                    
                end
            otherwise
                for w = 1:height(tab) % loop wells
                    if ~any(find(emptywells == w)) && ~isempty(tab.(mask_on)(w)) && tab.(mask_on)(w)<=mskth
                        warning(['entry ' num2str(w) ' removed for less than ', num2str(mskth), ' pixels'])
                        emptywells = [emptywells,w];
                    end                
                end  
        end
        otherwise
        switch parseName
            case 1
                for w = 1:height(tab) % loop wells
                    nameparse = split(tab.(parse_on){w},["-",'(',')']);
                    if length(nameparse)>3
                        tab.pPlateWell{w} = nameparse{end-1};
                        tab.plibName{w} = strcat(nameparse{1},'-', nameparse{2});
                        tab.pParent{w} = nameparse{1};
                        tab.pMutation{w} = nameparse{2}; 
                    elseif length(nameparse)==3
                        tab.pPlateWell{w} = nameparse{end-1};
                        tab.pParent{w} = nameparse{1};
                        tab.plibName{w} = nameparse{1};
                        tab.pMutation{w} = ''; 
                    else
                        warning(['entry ' num2str(w) ' removed for empty or unparsable name'])
                        emptywells = [emptywells,w];
                    end
                    if ~any(find(emptywells == w)) && ~isempty(tab.(mask_on)(w)) && tab.(mask_on)(w)<=mskth 
                        warning(['entry ' num2str(w) ' removed for less than ', num2str(mskth), ' pixels'])
                        emptywells = [emptywells,w];
                    end
                end
            otherwise
                for w = 1:height(tab) % loop wells
                    if ~any(find(emptywells == w)) && ~isempty(tab.(mask_on)(w)) && tab.(mask_on)(w)<=mskth
                        warning(['entry ' num2str(w) ' removed for less than ', num2str(mskth), ' pixels'])
                        emptywells = [emptywells,w];
                    end                
                end
        end
    tab(emptywells,:) = [];
end