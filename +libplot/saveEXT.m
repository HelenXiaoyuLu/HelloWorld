% Save figures (formatted)

function saveEXT(savePath)
    if ~isempty(savePath)
        [~,~,EXT] = fileparts(savePath);
        switch EXT
            case {'.png','.jpg'}
                exportgraphics(gcf,savePath,'Resolution','600');
            case '.pdf'
                exportgraphics(gcf,savePath,'ContentType',...
                    'vector','BackgroundColor','none');
            otherwise
                saveas(gcf,savePath);
        end
    end
end