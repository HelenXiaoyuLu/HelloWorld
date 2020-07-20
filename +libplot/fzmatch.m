% fzmatch return the closest legal name of a query GEVI name
% USAGE
%   lgname = fzmatch('JEDI1p');
%   lgname = fzmatch('Ace2N','guess','Ace');
%   lgname = fzmatch('Voltron','new',{'Voltron'}); ???
%   lgname = fzmatch('JEDI-1P-H152 (P0E2)','exempt',true); 
% INPUT
%   
% guess: make an initial guess of one of the commonly used nicknames, to
% avoid goint through calculating distances

function [lgname,minD] = fzmatch(query,varargin)
    p = inputParser;
    addRequired(p,'query',@(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p,'guess',@(x) ischar(x) || isstring(x) || iscell(x));
    p.addParameter('new', {},@iscell);
    p.addParameter('returnopt','string', @ischar);
    p.addParameter('exempt',false, @(n) validateattributes(n, ...
        {'logical'},{'scalar'}));
    p.parse(query,varargin{:});
    query = char(query);
    guess = char(p.Results.guess);
    new = p.Results.new;
    returnopt = p.Results.returnopt;
    exempt = p.Results.exempt;
    persistent LegalName NickNames
    LegalName = {'ASAP1','ASAP2s','ASAP3','ASAP2s-H152E','ASAP2s-T207H',...
        'ASAP2s-H152E-Q397H','JEDI-1P','JEDI-2P','Bongwoori-Pos6',...
        'Bongwoori-R3','ArcLight-A242','Ace2N-mNeon','MacQ-mCitrine',...
        'JEDI-alpha','JEDI-beta','JEDI-alpha-N124I-R406K','Marina','Lyn-EGFP','ASAP1-dpOPT','ASAP1-EGFP'};
    NickNames = {{'ASAP1'},{'ASAP2s','ASAP2'},{'ASAP3'},{'ASAP2s-H152E','GV3'},...
        {'ASAP2s-T207H','GV53'},{'ASAP2s-H152E-Q397H','GV75'},{'JEDI-1P','CBGV75','X2'},...
        {'JEDI-2P','CBGV79','X3'},{'Bongwoori-Pos6','BWR-P6'},{'Bongwoori-R3','BWR-R3'},...
        {'ArcLight-A242','Arclight'},{'Ace2N-mNeon','Ace'},{'MacQ-mCitrine','MacQ'},...
        {'JEDI-alpha','JEDIa'},{'JEDI-beta','JEDIb'},{'JEDI-alpha-N124I-R406K','JEDI-alpha-X1'},...
        {'Marina'},{'Lyn-EGFP','lyn'},{'ASAP1-dpOPT'},{'ASAP1-EGFP'}};
    if ~isempty(new)
        LegalName{end+1} = new{1};
        NickNames{end+1} = new(end);
    end
    
    lgname = '';
    minD = 0;
    dict = containers.Map(LegalName,NickNames);
    if isKey(dict,query) || exempt 
        lgname = query; % return the query name if it's legal
    elseif ~isempty(guess) % try find 
        entryidx = find(cellfun(@(c) any(strcmp(c, guess)), NickNames));
        if ~isempty(entryidx)
            lgname = LegalName{entryidx};            
        end
    end
    
    % if query or guess name wasn't found, try matching the query to
    % a neighbor with minimum Levenshtein distance
    if isempty(lgname)
        minD = inf;
        entryidx = inf;
        for i = 1:length(NickNames)
            for j = 1:length(NickNames{i})
                currentD = LVdistance(NickNames{i}{j},query);
                if currentD < minD
                    minD = currentD;
                    entryidx = i;
                end
            end
        end
        lgname = LegalName{entryidx};
        disp(strcat('Replacing',{' '},query,{' '},'to',{' '},lgname));
    end 
    
    if strcmp(returnopt ,'string')
        lgname = string(lgname);
    end    
end

% Calculate LVdistance using Wagner and Fischer algorithm (see link)
% https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm
function [d,LVmat] = LVdistance(s1,s2)
    LVmat = zeros(length(s1),length(s2));
    LVmat(:,1) = 0:length(s1)-1;
    LVmat(1,:) = 0:length(s2)-1;
    for i = 2:length(s1)
        for j = 2:length(s2)
            if s1(i) == s2(j)
                subscost = 0;
            elseif isnumeric(s1(i)) 
                subscost = 2; % add more weight to numbers 
            else 
                subscost = 1;
            end
            LVmat(i,j) = min([LVmat(i-1,j)+1,LVmat(i,j-1)+1,LVmat(i-1,j-1)+subscost]);
        end
    end
    d = LVmat(end,end);
end