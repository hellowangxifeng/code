function file=rp_s2_finish(f,c,stem)
%RP_S2_FINISH Compact physical figures, 10-point ticks, editable FIG plus PNG.
axesHandles=findall(f,'Type','axes');
set(axesHandles,'FontUnits','points','FontSize',10,'LineWidth',.75,'Box','on');
set(findall(f,'Type','colorbar'),'FontSize',10);
% Re-evaluate corners AFTER final title/layout/font changes, so saved legends
% are tested at the actual physical size used for PNG/FIG export.
legends=findall(f,'Type','legend');
for j=1:numel(legends)
    lg=legends(j); info=lg.UserData;
    if ~isstruct(info) || ~isfield(info,'axisTag'), continue; end
    ax=findobj(f,'Type','axes','Tag',info.axisTag); handles=gobjects(1,numel(info.sourceTags));
    for k=1:numel(handles), handles(k)=findobj(ax,'Type','line','Tag',info.sourceTags{k}); end
    appearance=struct('Color',lg.Color,'TextColor',lg.TextColor,'Box',lg.Box,'FontSize',lg.FontSize);
    rp_s2_legend(ax,handles,lg.String,lg.NumColumns,lg.FontSize,lg.ItemTokenSize,appearance);
end
file=fullfile(c.output,[stem '.png']);
savefig(f,fullfile(c.output,[stem '.fig']));
exportgraphics(f,file,'Resolution',c.plot.resolution,'BackgroundColor','white');
close(f);
end
