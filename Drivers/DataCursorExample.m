close all

% first plot
f1 = figure;
hold on
p1 = plot(1:10,ones(1,10));
row = dataTipTextRow('Name',1:10);
datacursormode on; % enable datatip mode
c1 = datacursormode(f1); % get the cursor mode
d1 = c1.createDatatip(p1); % create a new datatip
p1.DataTipTemplate.DataTipRows(end+1) = row;

% second plot
f2 = figure;
hold on
p2 = plot(1:10,sin(1:10));
row = dataTipTextRow('Name',1:10);
datacursormode on;
c2 = datacursormode(f2);
d2 = c2.createDatatip(p2);
p2.DataTipTemplate.DataTipRows(end+1) = row;

% register the function to execute when the datatip changes.    
set(d1,'UpdateFcn',@(cursorMode,eventData) onDataTipUpdate(cursorMode,eventData, d2))
set(d2,'UpdateFcn',@(cursorMode,eventData) onDataTipUpdate(cursorMode,eventData, d1))

% callback function when the datatip changes
function displayText = onDataTipUpdate(cursorMode,eventData, d)
   pos = get(eventData,'Position'); % the new position of the datatip
   displayText = {['X: ',num2str(pos(1))], ...
                      ['Y: ',num2str(pos(2))]}; % construct the datatip text
   d.Position(1) = pos(1); % update the location of the other datatip.
%    plot(pos(1),pos(2),'-s','MarkerSize',50)
end