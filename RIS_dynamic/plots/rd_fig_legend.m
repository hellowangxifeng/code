function lg = rd_fig_legend(layout,handles,names,o,columns)
%RD_FIG_LEGEND One shared legend above the tiles, never on top of a curve.
valid = isgraphics(handles);
if nargin<5, columns = min(4,nnz(valid)); end
lg = legend(handles(valid),names(valid),'Orientation','horizontal', ...
    'NumColumns',max(1,columns),'Box','off','FontName',o.fontName,'FontSize',o.fontSize);
lg.Layout.Tile = 'north';
end
