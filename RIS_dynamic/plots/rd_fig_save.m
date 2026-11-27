function files = rd_fig_save(fig,name,o)
%RD_FIG_SAVE Export a figure in the requested formats, then close it.
files = strings(0,1);
drawnow;
for j = 1:numel(o.formats)
    format = lower(o.formats{j});
    path = fullfile(o.outputFolder,[name '.' format]);
    switch format
        case 'png'
            exportgraphics(fig,path,'Resolution',o.resolution);
        case 'pdf'
            exportgraphics(fig,path,'ContentType','vector');
        case 'fig'
            set(fig,'CreateFcn','set(gcbo,''Visible'',''on'')');
            savefig(fig,path,'compact');
        otherwise
            error('rd:Format','Unsupported figure format "%s".',format);
    end
    files(end+1,1) = string(path); %#ok<AGROW>
end
if strcmp(o.visible,'off'), close(fig); end
end
