% Group wells based on names, and do average and std for selected
% properties
% Used for benchmarking experiments (with repeated wells)
% USAGE:
%   T_well = libplot.groupWells(T_well,'groupOn','Legal Name','selectProperties',selpp1P);
% INPUT:
%   REQUIRED: Table(T_well)
%   PARAMETERS (Name-value pairs): 
%   groupOn: specify a name column, this is the criteria for grouping 
%   selectProperties: which properties to be average/std. string values 
%   will be pooled as a cell array 
% Xiaoyu updated on May 20th, 2020

function T_group = groupWells(T_well,varargin)
    p = inputParser;
    p.addRequired('T_well', @istable);
    p.addParameter('groupOn', 'Name', @(x) any(strcmp(x, T_well.Properties.VariableNames)));
    p.addParameter('selectProperties', T_well.Properties.VariableNames, @iscell); % group specified properties
    p.parse(T_well,varargin{:});
    groupOn = p.Results.groupOn;
    selpp = p.Results.selectProperties;  
    T_group = table();
    Wellgroups = findgroups(T_well.(groupOn));
    T_group.n_well = splitapply(@numel, T_well.(groupOn), Wellgroups);
    mergeFcn = @(strs) {str2arr(strs)};
    for i = 1:length(selpp)
        pp = selpp{i};
        if strmatch(pp,groupOn,'exact'); % Assume this column is standardized name
            T_group.(pp)(:) = splitapply(@unique, T_well.(pp), Wellgroups);
        elseif isnumeric(T_well.(pp)(1))
            T_group.(strcat(pp,'_Mean'))(:) = splitapply(@nanmean, T_well.(pp), Wellgroups);
            T_group.(strcat(pp,'_STD'))(:) = splitapply(@nanstd, T_well.(pp), Wellgroups);
        elseif isstring(T_well.(pp)(1))
            T_group.(pp)(:) = splitapply(mergeFcn, T_well.(pp), Wellgroups);
%         else
%             T_group.(strcat(pp,'_Mean'))(:) = splitapply(@mean, T_well.(pp), Wellgroups);
%             T_group.(strcat(pp,'_STD'))(:) = splitapply(@std, T_well.(pp), Wellgroups);
        end
    end
    
end
function merged = str2arr(strs)
    merged = strings(size(strs));
    for i = 1:length(strs)
        merged(i) = strs(i);
    end
end