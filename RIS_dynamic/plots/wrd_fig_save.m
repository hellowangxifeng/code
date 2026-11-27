function paths = wrd_fig_save(fig,name,opts)
%WRD_FIG_SAVE Write <name>.png and an editable <name>.fig, then close.
paths = [string(fullfile(opts.outputFolder,[name '.png'])); ...
    string(fullfile(opts.outputFolder,[name '.fig']))];
drawnow;
exportgraphics(fig,paths(1),'Resolution',opts.resolution);
% A figure saved while invisible would reopen invisible; force it visible.
set(fig,'CreateFcn','set(gcbo,''Visible'',''on'')');
savefig(fig,paths(2),'compact');
if strcmp(opts.visible,'off'), close(fig); end
end
