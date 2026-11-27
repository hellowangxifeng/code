function rd_fig_axes(ax,xLimits,o,floorValue,label,showX,showY)
%RD_FIG_AXES Log risk axes with the exact pointing (FOV) floor as a dotted line.
set(ax,'YScale','log','XLim',xLimits,'YLim',o.yLimits,'FontName',o.fontName, ...
    'FontSize',o.fontSize,'XMinorTick','on','Box','on','Layer','top');
grid(ax,'on');
if isfinite(floorValue) && floorValue>o.yLimits(1)
    yline(ax,floorValue,':k','LineWidth',0.8,'HandleVisibility','off');
end
if showX, xlabel(ax,'$P_t$ (dBm)','Interpreter','latex','FontSize',o.fontSize+1); end
if showY, ylabel(ax,label,'FontName',o.fontName,'FontSize',o.fontSize+1); end
end
